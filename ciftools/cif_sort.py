#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 the cif-tools authors
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
Sorts the atoms (_atom_site rows) of an mmCIF file.

As in pdb_sort, ATOM records come before HETATM records, and each group is
sorted by:
    chain, residue number, insertion code, altloc, atom id   (default, -CR)
    chain, altloc, atom id                                    (-C)
    chain order of appearance, residue number, insertion
    code, altloc, atom id                                     (-R)
so atoms keep their order within a residue. Chains sort as text. Unlike
pdb_sort, multi-model files are supported (each model is sorted on its
own, models stay in order), and nothing is removed. The rows of
_atom_site_anisotrop are put in the new order of their atoms.

Chains and residue numbers are auth identifiers (label ones with --label).

Usage:
    cif_sort [-<option>] <cif file>

Example:
    cif_sort 1abc.cif      # sorts by chain and residues
    cif_sort -C 1abc.cif   # sorts by chain only
    cif_sort -R 1abc.cif   # sorts by residue number/icode only

This program is part of the cif-tools suite, a port of pdb-tools to mmCIF.
"""

import sys

from ciftools import _atoms, _cif, _cli

NULLS = ('.', '?')


def check_input(args):
    """Return (file handle, sort keys, identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__)
    keys = set((option or 'CR').upper())
    for key in keys:
        if key not in 'CR':
            _cli.fail(__doc__, 'Sorting key is not valid: {!r}'.format(key))
    return fhandle, ''.join(sorted(keys)), scheme


def _number(value):
    num = _atoms.as_int(value)
    return (0, num, '') if num is not None else (1, 0, value)


def _key_function(fields, keys, chain_order):
    """A function giving each row its sort key."""
    model, group, chain, resnum = fields.model, fields.group, fields.chain, fields.resnum
    icode, altloc, id_pos = fields.icode, fields.altloc, fields.id
    by_chain, by_residue = 'C' in keys, 'R' in keys

    def key(values, position):
        chain_id = values[chain] if chain is not None else ''
        if chain_id not in chain_order:
            chain_order[chain_id] = len(chain_order)
        out = [chain_id if by_chain else chain_order[chain_id]]
        if by_residue:
            code = values[icode] if icode is not None else '?'
            out.append(_number(values[resnum]) if resnum is not None else (1, 0, ''))
            out.append('' if code in NULLS else code)
        loc = values[altloc] if altloc is not None else '.'
        out.append('' if loc in NULLS else loc)
        out.append(_number(values[id_pos]) if id_pos is not None else (0, position, ''))
        hetatm = group is not None and values[group] == 'HETATM'
        return (values[model] if model is not None else '',), hetatm, tuple(out)
    return key


def run(fhandle, keys='CR', scheme=None):
    """Yield the events of the file with the atoms sorted."""
    kind = None
    rank = {}  # atom id -> new position, for _atom_site_anisotrop
    buffered = []
    for ev in _cif.parse(fhandle, (_atoms.ATOM_SITE, _atoms.ANISOTROP)):
        etype = type(ev)
        if etype is _cif.Row and kind is not None:
            if kind == _atoms.ATOM_SITE:
                model, hetatm, order = key(ev.values, len(buffered))
                if model not in model_order:
                    model_order[model] = len(model_order)
                atom_id = ev.values[id_pos] if id_pos is not None else None
                buffered.append(((model_order[model], hetatm, order), atom_id, ev.lead, ev.text))
            else:
                atom_id = ev.values[aniso_id] if aniso_id is not None else None
                buffered.append((rank.get(atom_id, len(rank)), atom_id, ev.lead, ev.text))
            continue
        if etype is _cif.Loop:
            kind = ev.category.lower()
            loop, buffered = ev, []
            if kind == _atoms.ATOM_SITE:
                fields = _atoms.Fields(ev, scheme)
                key = _key_function(fields, keys, {})
                id_pos, model_order = fields.id, {}
            else:
                aniso_id = ev.index('id') if 'id' in ev else None
        elif etype is _cif.LoopEnd and kind is not None:
            buffered.sort(key=lambda item: item[0])
            if kind == _atoms.ATOM_SITE:
                rank = {atom_id: i for i, (_, atom_id, _, _) in enumerate(buffered)}
            for _, _, lead, text in buffered:
                yield _cif.Row(loop, (), _lead(lead), text)
            buffered, kind = [], None
        yield ev


def _lead(lead):
    """Keep comment lines before a moved row, but not mid-line spacing."""
    return lead if lead.endswith('\n') else ''


def main():
    fhandle, keys, scheme = check_input(sys.argv[1:])
    _cli.output(run(fhandle, keys, scheme), fhandle)


if __name__ == '__main__':
    main()
