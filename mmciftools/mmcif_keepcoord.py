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
Removes all non-coordinate categories from an mmCIF file.

Keeps the data_ line, _atom_site and _struct_conn (the mmCIF counterparts of
the ATOM/HETATM and CONECT records that pdb_keepcoord keeps); every other
category is removed, _atom_site_anisotrop included (pdb_keepcoord drops
ANISOU records too).

With -basic, it also keeps the categories that programs reading mmCIF often
use and that edits do not make wrong, as they do not point at atoms or
residue numbers:

    _entry                     the entry id
    _cell, _symmetry           the unit cell and space group
    _entity, _entity_poly,     the molecules and their sequences
    _entity_poly_seq
    _chem_comp, _atom_type     the chemical components and elements

Two things can still disagree after edits: _entity keeps describing
molecules whose atoms a later tool removes, and _entity_poly.pdbx_strand_id
names author chains, which chain renames do not update. Everything else
(validation reports, secondary structure, sequence schemes, sites,
assemblies, ...) refers to residues by number, the other tools do not update
it, and after edits it no longer matches the atoms.

Usage:
    mmcif_keepcoord [-basic] <cif file>

Example:
    mmcif_keepcoord 1abc.cif          # coordinates only
    mmcif_keepcoord -basic 1abc.cif   # also entry, cell, entities, components

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cif, _cli

BASIC = ('_entry', '_cell', '_symmetry', '_entity', '_entity_poly', '_entity_poly_seq', '_chem_comp',
         '_atom_type')


def check_input(args):
    """Return (file handle, keep the basic categories too, identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__)
    if option not in (None, 'basic'):
        _cli.fail(__doc__, 'Invalid option: {!r} (the only option is -basic)'.format('-' + option))
    return fhandle, option == 'basic', scheme


def run(fhandle, basic=False, scheme=None):
    """Yield the events of the file without the non-coordinate categories
    (keeping the BASIC ones too if `basic`)."""
    categories = (_atoms.ATOM_SITE, _atoms.STRUCT_CONN) + (BASIC if basic else ())
    for ev in _cif.parse(fhandle, categories):
        if type(ev) is not _cif.Verbatim:  # every other category is Verbatim here
            yield ev


def main():
    fhandle, basic, scheme = check_input(sys.argv[1:])
    _cli.output(run(fhandle, basic, scheme), fhandle)


if __name__ == '__main__':
    main()
