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
element can be found get '?', and a warning is printed.

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
    for ev in _cif.parse(fhandle):
        etype = type(ev)
        if etype is _cif.Row:
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
