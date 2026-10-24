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
Sets the occupancy (_atom_site.occupancy) of all atoms (default 1.0).

Values are written with two decimals.

Usage:
    mmcif_occ -<occupancy> <cif file>

Example:
    mmcif_occ -1.0 1abc.cif

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cif, _cli


def check_input(args):
    """Return (file handle, occupancy, identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__)
    try:
        occupancy = 1.0 if option is None else float(option)
    except ValueError:
        _cli.fail(__doc__, 'You provided an invalid occupancy value: {!r}'.format(option))
    return fhandle, occupancy, scheme


def run(fhandle, occupancy, scheme=None):
    """Yield the events of the file with every occupancy set to occupancy."""
    value = '{:.2f}'.format(occupancy)
    col = None
    for ev in _cif.parse(fhandle):
        if type(ev) is _cif.Row:
            if col is not None:
                ev[col] = value
        elif type(ev) is _cif.Loop:
            fields = _atoms.Fields(ev, scheme)
            fields.need('occupancy')
            col = ev.columns[fields.occupancy]
        yield ev


def main():
    fhandle, occupancy, scheme = check_input(sys.argv[1:])
    _cli.output(run(fhandle, occupancy, scheme), fhandle)


if __name__ == '__main__':
    main()
