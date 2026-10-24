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
Deletes atoms by entity: by entity number (label_entity_id) or by type.

Entity types are polymer, non-polymer (ligands and ions), branched
(glycans), macrolide and water; see mmcif_selentity. Unlike
mmcif_delhetatm, this never cuts modified residues (HETATM records inside
a polymer) out of their chain: -non-polymer,water removes the ligands, ions
and waters and nothing else. Rows of _atom_site_anisotrop and _struct_conn
of removed atoms are removed too.

Usage:
    mmcif_delentity -<entity id or type>[,...] <cif file>

Example:
    mmcif_delentity -water 1abc.cif
    mmcif_delentity -non-polymer,water 1abc.cif   # keeps chains and glycans
    mmcif_delentity -2 1abc.cif                   # removes entity 2

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cli
from mmciftools.mmcif_selentity import matcher, parse_entities


def check_input(args):
    """Return (file handle, entity ids, entity types)."""
    fhandle, option, _ = _cli.parse_args(args, __doc__, needs_option=True)
    ids, types = parse_entities(option, __doc__)
    return fhandle, ids, types


def run(fhandle, ids, types, scheme=None):
    """Yield the events of the file without the atoms of the entities."""
    source = _atoms.as_source(fhandle)
    make_match = matcher(source, ids, types)

    def make_keep(fields):
        match = make_match(fields)
        return lambda values: not match(values)
    return _atoms.filter_atoms(source, make_keep, scheme)


def main():
    fhandle, ids, types = check_input(sys.argv[1:])
    _cli.output(_cli.lazy(run, fhandle, ids, types), fhandle)


if __name__ == '__main__':
    main()
