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
Deletes insertion codes (_atom_site.pdbx_PDB_ins_code) and shifts the
residue numbering (auth_seq_id) to make room.

As in pdb_fixinsert, the first residue carrying an insertion code just
loses it, and every later residue whose number was already seen moves the
numbering of the rest of the chain up by one: 51, 52A, 52B, 53 becomes
51, 52, 53, 54. With an option, only the listed residues lose their codes
(the numbering still shifts for the rest of the chain).

The shift restarts at each new chain instance (a change of chain or of
label_asym_id, e.g. between a protein and its ligands or waters: where
pdb_fixinsert meets a TER record) and in each model (pdb_fixinsert never
restarts between models). Chains may have several
characters and numbers may be negative: -A9,B12,AB-3. Insertion codes only
exist with author identifiers, so --label is not accepted. The matching
columns of _atom_site_anisotrop and _struct_conn are updated too. A file
without a pdbx_PDB_ins_code column (e.g. an integrative model) has no
insertion codes and is written unchanged.

Usage:
    mmcif_fixinsert [-<option>] <cif file>

Example:
    mmcif_fixinsert 1abc.cif            # deletes ALL insertion codes
    mmcif_fixinsert -A9,B12 1abc.cif    # deletes insertion codes of residue 9
                                        # of chain A and 12 of chain B

This program is part of the mmcif-tools suite.
"""

import re
import sys

from mmciftools import _atoms, _cif, _cli

NULLS = ('.', '?')
RESIDUE = re.compile(r'(\S+?)(-?[0-9]+)\Z')


def parse_residues(option):
    """'A9,B12' -> {('A', 9), ('B', 12)} (ValueError if invalid)."""
    residues = set()
    for item in option.split(','):
        match = RESIDUE.match(item.strip())
        if not match:
            raise ValueError('Option invalid: {!r}'.format(item))
        residues.add((match.group(1), int(match.group(2))))
    return residues


def check_input(args, doc=__doc__):
    """Return (file handle, residues or empty set, identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, doc)
    if scheme == 'label':
        _cli.fail(doc, 'Insertion codes only exist with author identifiers (--auth).')
    try:
        residues = parse_residues(option) if option else set()
    except ValueError as err:
        _cli.fail(doc, err)
    return fhandle, residues, scheme


def _fixer(fields, residues):
    """A function giving each atom row its (new number, new insertion code
    or None to keep it), or None to leave the row alone."""
    fields.need('chain', 'resnum', 'icode')
    model, chain, resnum, icode = fields.model, fields.chain, fields.resnum, fields.icode
    uid_of = _atoms.getter(fields.resname, fields.chain, fields.resnum, fields.icode)
    # The shift restarts where pdb_fixinsert would meet a TER: at each new
    # chain instance (label_asym_id), e.g. between a protein and its ligands.
    loop = fields.loop
    segment_of = _atoms.getter(chain, loop.index('label_asym_id') if 'label_asym_id' in loop else None)
    current_model, current_segment = object(), object()
    seen, prev, offset, clean = set(), None, 0, False

    def fix(values):
        nonlocal current_model, current_segment, seen, prev, offset, clean
        if model is not None and values[model] != current_model:
            current_model, current_segment = values[model], segment_of(values)
            seen, prev, offset = set(), None, 0
        segment = segment_of(values)
        if segment != current_segment:
            current_segment, offset = segment, 0
        num = _atoms.as_int(values[resnum])
        if num is None:
            return None
        has_code = values[icode] not in NULLS
        residue = (values[chain], num)
        uid = uid_of(values)
        if uid != prev:
            clean = (has_code or residue in seen) and (not residues or residue in residues)
            prev = uid
            if residue in seen:
                offset += 1
        seen.add(residue)
        return str(num + offset), ('?' if clean and has_code else None)
    return fix


def _key(chain, resnum, icode, resname):
    return chain, resnum, None if icode in NULLS else icode, resname


def _has_codes(source):
    """True if the _atom_site loop has an insertion code column."""
    for ev in _cif.parse(source.lines(), _atoms.ATOM_SITE):
        if type(ev) is _cif.Loop:
            return 'pdbx_PDB_ins_code' in ev
    return True  # no atoms: let the usual checks report it


def run(fhandle, residues=frozenset(), scheme=None):
    """Yield the events of the file with insertion codes removed."""
    source = _atoms.as_source(fhandle)
    if not _has_codes(source):  # no insertion code column: nothing to remove
        for ev in _cif.parse(source.lines(), ()):
            yield ev
        return

    # First look: new (number, code) of the residues of the first model, for
    # _struct_conn (before _atom_site) and _atom_site_anisotrop.
    changes, first_model = {}, None
    for ev in _cif.parse(source.lines(), _atoms.ATOM_SITE):
        if type(ev) is _cif.Row:
            values = ev.values
            if fields.model is not None:
                if first_model is None:
                    first_model = values[fields.model]
                elif values[fields.model] != first_model:
                    break
            new = fix(values)
            if new is not None:
                key = _key(*where(values))
                changes[key] = new
                # also without the residue name, for rows that have none (e.g. written by gemmi)
                changes.setdefault(key[:3] + (None,), new)
        elif type(ev) is _cif.Loop:
            fields = _atoms.Fields(ev, 'auth')
            fix = _fixer(fields, residues)
            where = _atoms.getter(fields.chain, fields.resnum, fields.icode, fields.resname)

    # Anisotropic rows follow their atom (by id): models may hold different residues.
    has_aniso = any(type(ev) is _cif.Loop for ev in _cif.parse(source.lines(), _atoms.ANISOTROP))
    by_id = {}
    kind = None
    for ev in _cif.parse(source.lines(), (_atoms.ATOM_SITE, _atoms.ANISOTROP, _atoms.STRUCT_CONN)):
        etype = type(ev)
        if etype is _cif.Row:
            if kind == _atoms.ATOM_SITE:
                new = fix(ev.values)
                if new is not None:
                    if new[0] != ev.values[resnum_pos]:
                        ev[resnum_col] = new[0]
                    if new[1] is not None:
                        ev[icode_col] = new[1]
                    if has_aniso and id_pos is not None:
                        by_id[ev.values[id_pos]] = new
            elif kind == _atoms.ANISOTROP and by_id and aniso_id is not None and ev.values[aniso_id] in by_id:
                new = by_id[ev.values[aniso_id]]
                c_pos, r_pos, i_pos, n_pos = refs[0]
                if r_pos is not None and new[0] != ev.values[r_pos]:
                    ev[ev.loop.columns[r_pos]] = new[0]
                if i_pos is not None and new[1] is not None and new[1] != ev.values[i_pos]:
                    ev[ev.loop.columns[i_pos]] = new[1]
            elif kind is not None:
                for c_pos, r_pos, i_pos, n_pos in refs:
                    if r_pos is None:
                        continue
                    values = ev.values
                    key = _key(*(values[p] if p is not None else None
                                 for p in (c_pos, r_pos, i_pos, n_pos)))
                    new = changes.get(key)
                    if new is None and n_pos is not None:  # partner without a residue name match
                        new = changes.get(key[:3] + (None,))
                    if new is not None:
                        if new[0] != values[r_pos]:
                            ev[ev.loop.columns[r_pos]] = new[0]
                        if new[1] is not None and i_pos is not None:
                            ev[ev.loop.columns[i_pos]] = new[1]
        elif etype is _cif.Loop:
            kind = ev.category.lower()
            if kind == _atoms.ATOM_SITE:
                fields = _atoms.Fields(ev, 'auth')
                fix = _fixer(fields, residues)
                resnum_pos, id_pos = fields.resnum, fields.id
                resnum_col, icode_col = ev.columns[fields.resnum], ev.columns[fields.icode]
            else:
                refs = _atoms._reference_columns(ev, 'auth')
                aniso_id = ev.index('id') if kind == _atoms.ANISOTROP and 'id' in ev else None
        elif etype is _cif.LoopEnd:
            kind = None
        yield ev


def main(doc=__doc__):
    fhandle, residues, scheme = check_input(sys.argv[1:], doc)
    _cli.output(run(fhandle, residues, scheme), fhandle)


if __name__ == '__main__':
    main()
