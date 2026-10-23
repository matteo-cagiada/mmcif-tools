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
Renumbers the atom ids (_atom_site.id) of an mmCIF file starting from a
given number (default 1).

The ids of _atom_site_anisotrop are changed to match. In mmCIF the atom id
must be unique in the whole file, so numbering continues across models
(pdb_reatom restarts at every MODEL), and there is no upper limit.

Usage:
    cif_reatom -<number> <cif file>

Example:
    cif_reatom -10 1abc.cif   # renumbers from 10
    cif_reatom --1 1abc.cif   # renumbers from -1

This program is part of the cif-tools suite, a port of pdb-tools to mmCIF.
"""

import sys

from ciftools import _atoms, _cif, _cli


def check_input(args):
    """Return (file handle, first number, identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__)
    try:
        start = 1 if option is None else int(option)
    except ValueError:
        _cli.fail(__doc__, 'You provided an invalid atom serial number: {!r}'.format(option))
    return fhandle, start, scheme


def run(fhandle, start, scheme=None):
    """Yield the events of the file with atoms renumbered from start."""
    new_ids = {}  # old id -> new id, for _atom_site_anisotrop
    number = start
    kind = None
    for ev in _cif.parse(fhandle, (_atoms.ATOM_SITE, _atoms.ANISOTROP)):
        etype = type(ev)
        if etype is _cif.Row and id_col is not None:
            old = ev.values[id_pos]
            if kind == _atoms.ATOM_SITE:
                new = str(number)
                number += 1
                new_ids[old] = new
                if new != old:
                    ev[id_col] = new
            else:
                new = new_ids.get(old)
                if new is not None and new != old:
                    ev[id_col] = new
        elif etype is _cif.Loop:
            kind = ev.category.lower()
            id_pos = ev.index('id') if 'id' in ev else None
            id_col = ev.columns[id_pos] if id_pos is not None else None
        yield ev


def main():
    fhandle, start, scheme = check_input(sys.argv[1:])
    _cli.output(run(fhandle, start, scheme), fhandle)


if __name__ == '__main__':
    main()
