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
Sets the temperature factor (_atom_site.B_iso_or_equiv) of all atoms
(default 10.0).

Values are written with two decimals; a file without the column gets it.
Anisotropic factors (_atom_site_anisotrop) are not changed, as in pdb_b.

Usage:
    mmcif_b -<bfactor> <cif file>

Example:
    mmcif_b -10.0 1abc.cif

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cif, _cli


def check_input(args):
    """Return (file handle, b-factor, identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__)
    try:
        bfactor = 10.0 if option is None else float(option)
    except ValueError:
        _cli.fail(__doc__, 'You provided an invalid b-factor value: {!r}'.format(option))
    return fhandle, bfactor, scheme


def run(fhandle, bfactor, scheme=None):
    """Yield the events of the file with every B-factor set to bfactor."""
    value = '{:.2f}'.format(bfactor)
    col, target = None, None
    for ev in _cif.parse(fhandle):
        if type(ev) is _cif.Row:
            if target is not None:  # the column is added
                ev = ev.moved(target, {'B_iso_or_equiv': value})
            elif col is not None:
                ev[col] = value
        elif type(ev) is _cif.Loop:
            fields = _atoms.Fields(ev, scheme)
            target = None
            if fields.bfactor is None:
                target = ev.with_columns(ev.columns + ['B_iso_or_equiv'])
                ev = target
            else:
                col = ev.columns[fields.bfactor]
        yield ev


def main():
    fhandle, bfactor, scheme = check_input(sys.argv[1:])
    _cli.output(run(fhandle, bfactor, scheme), fhandle)


if __name__ == '__main__':
    main()
