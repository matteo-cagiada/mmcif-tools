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
Renames atoms sequentially (C1, C2, O1, ...) for each HETATM residue.

As in pdb_uniqname, only HETATM atoms are renamed (waters included), the
name is the element (_atom_site.type_symbol, see mmcif_element) followed by a
counter per element that restarts at each residue, and a missing element is
an error.

Atom names are the same chemical identity in both identifier schemes, so
both auth_atom_id and label_atom_id are changed, together with the atom
names of _atom_site_anisotrop and _struct_conn (pdb_uniqname leaves ANISOU
records alone). Residues are told apart by model, name, chain, number and
insertion code (auth identifiers, or label ones with --label).

Usage:
    mmcif_uniqname <cif file>

Example:
    mmcif_uniqname 1abc.cif

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cif, _cli

NULLS = ('.', '?')
NAME_COLUMNS = ('auth_atom_id', 'label_atom_id')


def check_input(args):
    """Return (file handle, identifier scheme)."""
    fhandle, _, scheme = _cli.parse_args(args, __doc__, takes_option=False)
    return fhandle, scheme


def _namer(fields):
    """A function giving each atom row its new name (None: not renamed)."""
    fields.need('group', 'element')
    group, element, id_pos = fields.group, fields.element, fields.id
    residue = _atoms.getter(fields.model, fields.resname, fields.chain, fields.resnum, fields.icode)
    prev, counts = object(), {}

    def name(row):
        nonlocal prev, counts
        values = row.values
        if values[group] != 'HETATM':
            return None
        symbol = values[element]
        if symbol in NULLS:
            raise _cif.CIFError('No element found for atom {} (see mmcif_element)'.format(
                values[id_pos] if id_pos is not None else '?'), row.lineno)
        key = residue(values)
        if key != prev:
            prev, counts = key, {}
        counts[symbol] = counts.get(symbol, 0) + 1
        return symbol + str(counts[symbol])
    return name


def run(fhandle, scheme=None):
    """Yield the events of the file with HETATM atoms renamed."""
    source = _atoms.as_source(fhandle)

    # First look: the new names, by atom id (for _atom_site_anisotrop) and by
    # link identity (for _struct_conn, which comes before _atom_site).
    by_id, renamed = {}, []
    conn_loop = atom_loop = None
    for ev in _cif.parse(source.lines(), (_atoms.STRUCT_CONN, _atoms.ATOM_SITE)):
        if type(ev) is _cif.Row:
            if ev.loop is atom_loop:
                new = name(ev)
                if new is not None:
                    renamed.append((ev.values, new))
                    if id_pos is not None:
                        by_id[ev.values[id_pos]] = new
        elif type(ev) is _cif.Loop:
            if ev.category.lower() == _atoms.STRUCT_CONN:
                conn_loop = ev
            elif atom_loop is None:
                atom_loop = ev
                fields = _atoms.Fields(ev, scheme)
                name, id_pos = _namer(fields), fields.id
    links, columns = {}, None
    if conn_loop is not None and atom_loop is not None:
        columns = _atoms._link_columns(atom_loop, conn_loop)
        if columns:
            atom_key = _atoms.getter(*[atom_loop.index(a) for a, _ in columns])
            links = {atom_key(values): new for values, new in renamed}

    kind = None
    for ev in _cif.parse(source.lines(), (_atoms.ATOM_SITE, _atoms.ANISOTROP, _atoms.STRUCT_CONN)):
        etype = type(ev)
        if etype is _cif.Row:
            if kind == _atoms.ATOM_SITE:
                new = name(ev)
                if new is not None:
                    for col in name_cols:
                        if ev[col] != new:
                            ev[col] = new
            elif kind == _atoms.ANISOTROP:
                new = by_id.get(ev.get('id'))
                if new is not None:
                    for col in aniso_cols:
                        ev[col] = new
            elif kind == _atoms.STRUCT_CONN and links:
                for partner, cols in partners:
                    new = links.get(partner(ev.values))
                    if new is not None:
                        for col in cols:
                            ev[col] = new
        elif etype is _cif.Loop:
            kind = ev.category.lower()
            if kind == _atoms.ATOM_SITE:
                name = _namer(_atoms.Fields(ev, scheme))
                name_cols = [c for c in NAME_COLUMNS if c in ev]
            elif kind == _atoms.ANISOTROP:
                aniso_cols = [c for c in ('pdbx_auth_atom_id', 'pdbx_label_atom_id') if c in ev]
            elif kind == _atoms.STRUCT_CONN and links:
                first, second = _atoms._partner_getters(ev, columns)
                partners = [(getter, [c.format(n) for c in ('ptnr{}_label_atom_id', 'ptnr{}_auth_atom_id')
                                      if c.format(n) in ev])
                            for n, getter in ((1, first), (2, second))]
        elif etype is _cif.LoopEnd:
            kind = None
        yield ev


def main():
    fhandle, scheme = check_input(sys.argv[1:])
    _cli.output(run(fhandle, scheme), fhandle)


if __name__ == '__main__':
    main()
