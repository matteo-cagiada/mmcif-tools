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
Merges several mmCIF files into one multi-model (ensemble) file.

As pdb_mkensemble, the n-th file becomes model n (all its atoms, whatever
their model numbers), the provenance of each model is written at the top
(as comment lines), the other categories are dropped except the
_struct_conn of the first file, and the models are not checked against
each other (see mmcif_chkensemble).

The _atom_site loop has the columns of all files (missing values are
'?'), pdbx_PDB_model_num is added if needed, and atom ids are renumbered
from 1 over the whole ensemble, since mmCIF atom ids must be unique.

Usage:
    mmcif_mkensemble <cif file> <cif file> [...]

Example:
    mmcif_mkensemble 1abc.cif 1xyz.cif > ensemble.cif

This program is part of the mmcif-tools suite.
"""

import os
import sys

from mmciftools import _atoms, _cif, _cli

MODEL = 'pdbx_PDB_model_num'


def check_input(args):
    """Return (paths, identifier scheme)."""
    paths, _, scheme = _cli.parse_files(args, __doc__)
    return paths, scheme


def union_columns(paths):
    """_atom_site columns of all files, in order of first appearance."""
    columns, seen = [], set()
    for path in paths:
        found = _atoms.atom_columns(path)
        if found is None:
            raise _cif.CIFError('{} has no _atom_site category'.format(path))
        for column in found:
            if column.lower() not in seen:
                seen.add(column.lower())
                columns.append(column)
    return columns


def run(paths, scheme=None):
    """Yield the events of the ensemble."""
    columns = union_columns(paths)
    if MODEL.lower() not in {c.lower() for c in columns}:
        columns.append(MODEL)
    for n, path in enumerate(paths, 1):
        yield '# MODEL {} FROM {}\n'.format(n, os.path.basename(path))

    loop, number = None, 0
    links = []  # _struct_conn of the first file, written after the atoms
    for n, path in enumerate(paths, 1):
        source_loop = None
        with _cif.open_input(path) as fh:
            for ev in _cif.parse(fh, (_atoms.ATOM_SITE, _atoms.STRUCT_CONN)):
                etype = type(ev)
                if etype is _cif.DataBlock:
                    if n == 1 and loop is None:
                        yield ev
                    continue
                if etype is _cif.Verbatim:
                    continue
                kind = ev.loop.category.lower() if etype is not _cif.Loop else ev.category.lower()
                if kind == _atoms.STRUCT_CONN:
                    if n == 1:
                        links.append(ev)
                elif etype is _cif.Loop:
                    if source_loop is None:
                        source_loop = ev
                        if loop is None:
                            loop = ev.with_columns(columns)
                            yield '#\n'
                            yield loop
                elif etype is _cif.Row and ev.loop is source_loop:
                    number += 1
                    yield ev.moved(loop, {'id': str(number), MODEL: str(n)})
    if loop is not None:
        yield _cif.LoopEnd(loop, '#\n')
    for ev in links:
        yield ev


def main():
    paths, scheme = check_input(sys.argv[1:])
    _cli.output(run(paths, scheme))


if __name__ == '__main__':
    main()
