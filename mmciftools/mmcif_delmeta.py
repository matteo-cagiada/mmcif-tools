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
Deletes the metadata that edits make stale, keeping the coordinates and the
metadata that stays true.

PDB entries carry dozens of categories that refer to residues by chain and
number: validation reports, secondary structure, sequence schemes, sites,
assemblies, ... The other tools do not update them, so after an edit they no
longer match the atoms, and most programs that read mmCIF ignore them anyway.
This tool keeps:

    _atom_site                 the atoms
    _atom_site_anisotrop       their anisotropic B-factors
    _struct_conn               the links (disulfides, metal sites, ...)
    _entry                     the entry id
    _cell, _symmetry           the unit cell and space group
    _entity, _entity_poly,     the molecules and their sequences
    _entity_poly_seq
    _chem_comp, _atom_type     the chemical components and elements

and deletes every other category. With -all, only the atoms, anisotropic
B-factors and links are kept (mmcif_keepcoord also drops the anisotropic
B-factors, as pdb_keepcoord drops ANISOU records).

Run it first in a pipeline, before the edits. Two kept items can still
disagree after edits: _entity keeps describing molecules whose atoms a later
tool removes, and _entity_poly.pdbx_strand_id names author chains, which
chain renames do not update.

Usage:
    mmcif_delmeta [-all] <cif file>

Example:
    mmcif_delmeta 1abc.cif > 1abc_clean.cif
    mmcif_delmeta 1abc.cif | mmcif_delentity -non-polymer,water > protein.cif
    mmcif_delmeta -all 1abc.cif > atoms_only.cif

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cif, _cli

COORDINATES = (_atoms.ATOM_SITE, _atoms.ANISOTROP, _atoms.STRUCT_CONN)
KEPT = ('_entry', '_cell', '_symmetry', '_entity', '_entity_poly', '_entity_poly_seq', '_chem_comp',
        '_atom_type')


def check_input(args):
    """Return (file handle, delete all metadata, identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__)
    if option not in (None, 'all'):
        _cli.fail(__doc__, 'Invalid option: {!r} (the only option is -all)'.format('-' + option))
    return fhandle, option == 'all', scheme


def run(fhandle, everything=False, scheme=None):
    """Yield the events of the file without its metadata (all of it if
    `everything`, else all but KEPT)."""
    categories = COORDINATES + (() if everything else KEPT)
    for ev in _cif.parse(fhandle, categories):
        if type(ev) is not _cif.Verbatim:  # every other category is Verbatim here
            yield ev


def main():
    fhandle, everything, scheme = check_input(sys.argv[1:])
    _cli.output(run(fhandle, everything, scheme), fhandle)


if __name__ == '__main__':
    main()
