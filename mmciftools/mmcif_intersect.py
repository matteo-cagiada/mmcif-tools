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
Keeps only the atoms of the first mmCIF file that are present in all the
others.

As in pdb_intersect, atoms are judged equal if their name, altloc, residue
name, residue number, insertion code and chain are the same (auth
identifiers, or label ones with --label; models are not compared), and
coordinates are taken from the first file. Unlike pdb_intersect, the output
is the first file with the other atoms removed (metadata kept; rows of
_atom_site_anisotrop and _struct_conn of removed atoms removed).

Usage:
    mmcif_intersect <cif file> <cif file> [...]

Example:
    mmcif_intersect 1xyz.cif 1abc.cif

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cif, _cli


def check_input(args):
    """Return (paths, identifier scheme)."""
    paths, _, scheme = _cli.parse_files(args, __doc__)
    return paths, scheme


def _identity(fields):
    fields.need('chain', 'resnum', 'resname', 'atom')
    return _atoms.getter(fields.atom, fields.altloc, fields.resname, fields.chain,
                         fields.resnum, fields.icode)


def atoms_of(path, scheme=None):
    """Set of atom identities of a file."""
    found = set()
    with _cif.open_input(path) as fh:
        for ev in _cif.parse(fh):
            if type(ev) is _cif.Row:
                found.add(identity(ev.values))
            elif type(ev) is _cif.Loop:
                identity = _identity(_atoms.Fields(ev, scheme))
    return found


def run(paths, scheme=None):
    """Yield the events of the first file, keeping the common atoms."""
    common = None
    for path in paths[1:]:
        found = atoms_of(path, scheme)
        common = found if common is None else common & found

    def make_keep(fields):
        identity = _identity(fields)
        if common is None:  # a single file
            return lambda values: True
        return lambda values: identity(values) in common

    with _cif.open_input(paths[0]) as fh:
        for ev in _atoms.filter_atoms(fh, make_keep, scheme):
            yield ev


def main():
    paths, scheme = check_input(sys.argv[1:])
    _cli.output(run(paths, scheme))


if __name__ == '__main__':
    main()
