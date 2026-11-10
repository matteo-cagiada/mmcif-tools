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
Assigns the elements (_atom_site.type_symbol) of an mmCIF file from the
atom names.

By default only missing elements ('?', '.', or no type_symbol column) are
assigned, because deposited elements are authoritative; -all re-assigns
every atom (pdb_element always overwrites).

mmCIF atom names carry no column alignment (pdb_element tells ' CA ', an
alpha carbon, from 'CA  ', calcium, by it), so the rules are:
    1. an atom named like its residue is that element (ZN in ZN, CA in CA)
    2. SE in MSE / SEC is selenium
    3. in ligands (HETATM, not a standard residue) a name starting with a
       common two-letter element is that element (CL1, BR2, FE, ZN1, MG ...;
       not CA, NA, HG, CD, which in ligands are mostly carbon, nitrogen or
       hydrogen names, e.g. the heme nitrogens NA-ND)
    4. otherwise the first letter of the name, after leading digits
Atom names are auth_atom_id (label_atom_id with --label). Atoms for which no
element can be found get '?', and a warning is printed. An element new to the
file gets a row in _atom_type, the table that lists the elements.

Usage:
    mmcif_element [-all] <cif file>

Example:
    mmcif_element 1abc.cif        # fills in missing elements
    mmcif_element -all 1abc.cif   # re-assigns every element

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cif, _cli

NULLS = ('.', '?')
ATOM_TYPE = '_atom_type'


def check_input(args):
    """Return (file handle, reassign all, identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__)
    if option is not None and option != 'all':
        _cli.fail(__doc__, 'You provided an invalid option: {!r}'.format(option))
    return fhandle, option == 'all', scheme


def run(fhandle, reassign=False, scheme=None, unknown=None):
    """Yield the events of the file with elements assigned.

    `unknown`, if given, is a list that receives the ids of atoms for which
    no element was found.
    """
    source = _atoms.as_source(fhandle)
    missing = _missing_types(source, reassign, scheme)
    for ev in _assign(_cif.parse(source.lines(), (_atoms.ATOM_SITE, ATOM_TYPE) if missing else _atoms.ATOM_SITE),
                      reassign, scheme, unknown):
        if missing and type(ev) is _cif.LoopEnd and ev.loop.category.lower() == ATOM_TYPE:
            for symbol in missing:  # a new element is listed in _atom_type
                values = ['?'] * len(ev.loop.columns)
                values[ev.loop.index('symbol')] = symbol
                yield _cif.Row(ev.loop, values)
            missing = []
        yield ev


def _missing_types(source, reassign, scheme):
    """The elements new to the file that _atom_type does not list (none when
    the file has no _atom_type)."""
    listed = None
    for ev in _cif.parse(source.lines(), ATOM_TYPE):
        if type(ev) is _cif.Loop and listed is None and 'symbol' in ev:
            listed, symbol = set(), ev.index('symbol')
        elif type(ev) is _cif.Row and listed is not None:
            listed.add(ev.values[symbol])
    if listed is None:
        return []
    known = set(listed) | set(NULLS)  # and the elements the atoms already have

    def atoms():
        for ev in _cif.parse(source.lines(), _atoms.ATOM_SITE):
            if type(ev) is _cif.Row:
                known.add(ev.get('type_symbol'))
            yield ev
    new = {}
    for ev in _assign(atoms(), reassign, scheme):
        if type(ev) is _cif.Row:
            new[ev['type_symbol']] = True
    return [element for element in new if element not in known]


def _assign(events, reassign, scheme, unknown=None):
    """The events with the elements of the atoms assigned."""
    kind = None
    for ev in events:
        etype = type(ev)
        if etype is _cif.Row and kind == _atoms.ATOM_SITE:
            values = ev.values
            current = values[element_pos] if element_pos is not None else '?'
            if reassign or current in NULLS:
                hetatm = group is not None and values[group] == 'HETATM'
                new = _atoms.guess_element(values[name], values[resname] if resname is not None
                                           else '', hetatm) or '?'
                if new == '?' and unknown is not None:
                    unknown.append(values[id_pos] if id_pos is not None else '?')
                if element_pos is None:
                    ev = ev.moved(loop, {'type_symbol': new})
                elif new != current:
                    ev[element_col] = new
        elif etype is _cif.Loop:
            kind = ev.category.lower()
            if kind != _atoms.ATOM_SITE:
                yield ev
                continue
            fields = _atoms.Fields(ev, scheme)
            fields.need('atom')
            name, resname, group, id_pos = fields.atom, fields.resname, fields.group, fields.id
            element_pos = fields.element
            if element_pos is None:  # add the column
                loop = ev.with_columns(ev.columns + ['type_symbol'])
                ev = loop
            else:
                element_col = ev.columns[element_pos]
        yield ev


def main():
    fhandle, reassign, scheme = check_input(sys.argv[1:])
    unknown = []

    def events():
        for ev in run(fhandle, reassign, scheme, unknown):
            yield ev
        if unknown:
            sys.stderr.write('WARNING: no element found for {} atom(s), e.g. id {}\n'.format(
                len(unknown), ', '.join(unknown[:5])))
    _cli.output(events(), fhandle)


if __name__ == '__main__':
    main()
