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
Selects all HETATM records (ligands, ions, waters, ...) in an mmCIF file.

Atoms are selected on _atom_site.group_PDB. Unlike pdb_selhetatm, which
drops every other record, the rest of the file (metadata categories) is
kept; rows of _atom_site_anisotrop and _struct_conn that refer to removed
atoms are removed.

Usage:
    mmcif_selhetatm <cif file>

Example:
    mmcif_selhetatm 1abc.cif

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cli


def check_input(args):
    """Return (file handle, identifier scheme)."""
    fhandle, _, scheme = _cli.parse_args(args, __doc__, takes_option=False)
    return fhandle, scheme


def run(fhandle, scheme=None):
    """Yield the events of the file, keeping only HETATM atoms."""
    def make_keep(fields):
        fields.need('group')
        group = fields.group
        return lambda values: values[group] == 'HETATM'

    return _atoms.filter_atoms(fhandle, make_keep, scheme)


def main():
    fhandle, scheme = check_input(sys.argv[1:])
    _cli.output(run(fhandle, scheme), fhandle)


if __name__ == '__main__':
    main()
