#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 the mmcif-tools authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
Copies the label identifiers of an mmCIF file onto the author ones.

By default copies both chains and residue numbers; an option copies only
one of them:

    -chain   auth_asym_id <- label_asym_id. Every molecule gets a chain of
             its own: in PDB entries, ligands and waters share the author
             chain of their protein, but each has its own label_asym_id.
    -res     auth_seq_id <- label_seq_id, for polymer residues. They are
             then numbered 1, 2, ... along the deposited sequence, so
             missing residues keep their place (unlike mmcif_reres), and
             insertion codes are removed. Ligands and waters have no
             label_seq_id and keep their numbers.

This is the mmCIF counterpart of pdb_chainxseg / pdb_segxchain and, with
-chain, of pdb_chainbows. Rows of _atom_site_anisotrop and _struct_conn
follow. Other categories that mention author chains or residue numbers
(_pdbx_poly_seq_scheme, _struct_conf, ...) are not updated. A file without
auth_asym_id or auth_seq_id columns (e.g. an integrative model) gets them;
residues without a label_seq_id are then numbered '?'. Nothing is
written if two residues would end up with the same chain and number, e.g.
with -res when a ligand of the chain already has one of the new numbers.

Usage:
    mmcif_label2auth [-chain | -res] <cif file>

Example:
    mmcif_label2auth 1abc.cif          # chains and residue numbers
    mmcif_label2auth -chain 1abc.cif   # one chain per molecule

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cif, _cli

NULLS = ('.', '?')
OPTIONS = ('chain', 'res')


def check_input(args):
    """Return (file handle, copy chains, copy residue numbers)."""
    fhandle, option, _ = _cli.parse_args(args, __doc__)
    if option is not None and option not in OPTIONS:
        _cli.fail(__doc__, 'Invalid option: {!r} (use -chain or -res)'.format('-' + option))
    return fhandle, option != 'res', option != 'chain'


def _new_ids(loop, chains, residues):
    """A function: atom values -> (old, new) (chain, number, insertion code)."""
    def pos(name, needed):
        if name in loop:
            return loop.index(name)
        if needed:
            raise _cif.CIFError('{} has no {} column (needed here)'.format(loop.category, name),
                                loop.lineno)
        return None

    auth = _atoms.getter(pos('auth_asym_id', False), pos('auth_seq_id', False),
                         pos('pdbx_PDB_ins_code', False))
    label_chain = pos('label_asym_id', chains)
    label_seq = pos('label_seq_id', residues)

    def new_ids(values):
        old = auth(values)
        chain, number, icode = old
        if chains:
            chain = values[label_chain]
            if chain in NULLS:
                raise _cif.CIFError('An atom has no label_asym_id: {}'.format(values))
        if residues and _atoms.as_int(values[label_seq]) is not None:
            number = values[label_seq]
            if icode is not None:
                icode = '?'
        return old, (chain, number, icode)
    return new_ids


def residue_map(source, chains, residues):
    """{old (chain, number, insertion code): new one}, after checking that no
    two residues get the same new identifiers. Old (chain, number) pairs are
    keys too, for rows without an insertion code column (None when two
    residues share them)."""
    mapping, owner = {}, {}
    for ev in _cif.parse(source.lines(), _atoms.ATOM_SITE):
        if type(ev) is _cif.Row:
            values = ev.values
            old, new = new_ids(values)
            if old not in mapping:
                mapping[old] = new
                short = old[:2]
                mapping[short] = new if short not in mapping else None
            here = (model(values), new)
            ident = (old, values[label_chain]) if label_chain is not None else old
            first = owner.setdefault(here, ident)
            if first != ident:
                raise _cif.CIFError(
                    'Residues {} and {} would both become {} {}{}'.format(
                        _name(first[0] if label_chain is not None else first),
                        _name(old), new[0], new[1], _icode(new[2])), ev.lineno)
        elif type(ev) is _cif.Loop:
            new_ids = _new_ids(ev, chains, residues)
            model = _atoms.getter(ev.index('pdbx_PDB_model_num')
                                  if 'pdbx_PDB_model_num' in ev else None)
            label_chain = ev.index('label_asym_id') if 'label_asym_id' in ev else None
    return mapping


def _icode(icode):
    return '' if icode is None or icode in NULLS else icode


def _name(ids):
    return '{} {}{}'.format(ids[0], ids[1], _icode(ids[2]))


def run(fhandle, chains=True, residues=True, scheme=None):
    """Yield the events of the file with label identifiers copied."""
    source = _atoms.as_source(fhandle)
    mapping = residue_map(source, chains, residues)
    kind = None
    for ev in _cif.parse(source.lines(), (_atoms.ATOM_SITE, _atoms.ANISOTROP, _atoms.STRUCT_CONN)):
        etype = type(ev)
        if etype is _cif.Row:
            if kind == _atoms.ATOM_SITE:
                old, new = new_ids(ev.values)
                if target is not None:  # with the added author columns ('?' at first)
                    ev = ev.moved(target)
                for column, before, after in zip(names, old, new):
                    if column is not None and after is not None and after != before:
                        ev[column] = after
            elif kind is not None:
                for positions in refs:
                    values = ev.values
                    old = tuple(values[p] if p is not None else None for p in positions)
                    new = mapping.get(old if positions[2] is not None else old[:2])
                    if new is None:
                        continue
                    for p, before, after in zip(positions, old, new):
                        if p is not None and after != before:
                            ev[ev.loop.columns[p]] = after
            yield ev
            continue
        if etype is _cif.Loop:
            kind = ev.category.lower()
            if kind == _atoms.ATOM_SITE:
                new_ids = _new_ids(ev, chains, residues)
                missing = [c for c, wanted in (('auth_asym_id', chains), ('auth_seq_id', residues))
                           if wanted and c not in ev]
                target = ev.with_columns(ev.columns + missing) if missing else None
                names = [c if c in ev or c in missing else None
                         for c in ('auth_asym_id', 'auth_seq_id', 'pdbx_PDB_ins_code')]
                if target is not None:
                    ev = target
            else:
                refs = [r[:3] for r in _atoms._reference_columns(ev, 'auth')]
        elif etype is _cif.LoopEnd:
            kind = None
        yield ev


def main():
    fhandle, chains, residues = check_input(sys.argv[1:])
    _cli.output(run(fhandle, chains, residues), fhandle)


if __name__ == '__main__':
    main()
