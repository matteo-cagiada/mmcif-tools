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
Keeps the first N atoms (_atom_site rows) of an mmCIF file.

Atoms are counted over the whole file, across models. Unlike pdb_head, which
outputs nothing but the first N ATOM/HETATM lines, the rest of the file
(metadata categories) is kept, so the output is still a complete mmCIF file;
rows of _atom_site_anisotrop and _struct_conn that refer to removed atoms
are removed.

Usage:
    cif_head -<num> <cif file>

Example:
    cif_head -100 1abc.cif   # first 100 atoms of the file

This program is part of the cif-tools suite, a port of pdb-tools to mmCIF.
"""

import sys

from ciftools import _atoms, _cli


def check_input(args):
    """Return (file handle, number of atoms, identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__, needs_option=True)
    try:
        count = int(option)
        if count < 1:
            raise ValueError
    except ValueError:
        _cli.fail(__doc__, 'Option must be a positive numerical value: {!r}'.format(option))
    return fhandle, count, scheme


def run(fhandle, count, scheme=None):
    """Yield the events of the file, keeping only the first `count` atoms."""
    def make_keep(fields):
        seen = 0

        def keep(values):
            nonlocal seen
            seen += 1
            return seen <= count
        return keep

    return _atoms.filter_atoms(fhandle, make_keep, scheme)


def main():
    fhandle, count, scheme = check_input(sys.argv[1:])
    _cli.output(run(fhandle, count, scheme), fhandle)


if __name__ == '__main__':
    main()
