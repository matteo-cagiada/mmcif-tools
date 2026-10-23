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
Selects all residues matching the given name(s) in an mmCIF file.

Residue names are matched on auth_comp_id (label_comp_id with --label);
the names given are upper-cased, as in pdb_selresname. Rows of
_atom_site_anisotrop and _struct_conn that refer to removed atoms are
removed too.

Usage:
    cif_selresname -<name>[,<name>...] <cif file>

Example:
    cif_selresname -ALA 1abc.cif      # keeps only alanines
    cif_selresname -ASP,GLU 1abc.cif  # keeps (-) charged residues

This program is part of the cif-tools suite, a port of pdb-tools to mmCIF.
"""

import sys

from ciftools import _atoms, _cli


def check_input(args):
    """Return (file handle, residue names, identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__, needs_option=True)
    names = [n.upper() for n in _atoms.parse_ids(option)]
    if not names:
        _cli.fail(__doc__, 'Residue name set cannot be empty')
    return fhandle, names, scheme


def run(fhandle, names, scheme=None):
    """Yield the events of the file, keeping only residues with the given names."""
    names = set(names)

    def make_keep(fields):
        fields.need('resname')
        resname = fields.resname
        return lambda values: values[resname] in names

    return _atoms.filter_atoms(fhandle, make_keep, scheme)


def main():
    fhandle, names, scheme = check_input(sys.argv[1:])
    _cli.output(run(fhandle, names, scheme), fhandle)


if __name__ == '__main__':
    main()
