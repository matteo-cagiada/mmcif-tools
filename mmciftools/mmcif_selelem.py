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
Selects all atoms of the given element(s) in an mmCIF file.

Elements are read from _atom_site.type_symbol and compared ignoring case
(pdb_selelem only matches upper-case elements). Rows of
_atom_site_anisotrop and _struct_conn that refer to removed atoms are
removed too.

Usage:
    mmcif_selelem -<element>[,<element>...] <cif file>

Example:
    mmcif_selelem -H 1abc.cif     # selects all hydrogens
    mmcif_selelem -N 1abc.cif     # selects all nitrogens
    mmcif_selelem -H,N 1abc.cif   # selects all hydrogens and nitrogens

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cli


def check_input(args):
    """Return (file handle, elements, identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__, needs_option=True)
    elements = [e.upper() for e in _atoms.parse_ids(option)]
    if not elements:
        _cli.fail(__doc__, 'Element set cannot be empty')
    return fhandle, elements, scheme


def run(fhandle, elements, scheme=None):
    """Yield the events of the file, keeping only atoms of the given elements."""
    elements = set(e.upper() for e in elements)

    def make_keep(fields):
        fields.need('element')
        element = fields.element
        return lambda values: values[element].upper() in elements

    return _atoms.filter_atoms(fhandle, make_keep, scheme)


def main():
    fhandle, elements, scheme = check_input(sys.argv[1:])
    _cli.output(run(fhandle, elements, scheme), fhandle)


if __name__ == '__main__':
    main()
