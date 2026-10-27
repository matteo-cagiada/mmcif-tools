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
ANISOU records too). To keep the metadata that edits do not make stale
(entry, unit cell, entities, sequences, chemical components), use
mmcif_delmeta instead.

Usage:
    mmcif_keepcoord <cif file>

Example:
    mmcif_keepcoord 1abc.cif

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cif, _cli


def check_input(args):
    """Return (file handle, identifier scheme)."""
    fhandle, _, scheme = _cli.parse_args(args, __doc__, takes_option=False)
    return fhandle, scheme


def run(fhandle, scheme=None):
    """Yield the events of the file without the non-coordinate categories."""
    for ev in _cif.parse(fhandle, (_atoms.ATOM_SITE, _atoms.STRUCT_CONN)):
        if type(ev) is not _cif.Verbatim:  # every other category is Verbatim here
            yield ev


def main():
    fhandle, scheme = check_input(sys.argv[1:])
    _cli.output(run(fhandle, scheme), fhandle)


if __name__ == '__main__':
    main()
