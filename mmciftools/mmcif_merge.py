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
Merges the atoms of several mmCIF files into one.

pdb_merge concatenates files; concatenated mmCIF files would just be
several data blocks, so instead the output is the first file with the atoms
of all the others added to its _atom_site loop (which gets the columns of
all files; missing values are '?'). The atoms of the first file keep their
ids, those of the next files are renumbered after them. Only _atom_site is
taken from the other files (from their first data block, added to the
first data block of the first file): their metadata, _atom_site_anisotrop and
_struct_conn rows are not merged, and _entity / _struct_asym are not
updated.

Chains are not renamed: a warning is printed when two files use the same
chain (auth_asym_id, or label_asym_id with --label); see mmcif_rplchain.
Pipe the result through mmcif_tidy to align the columns.

Usage:
    mmcif_merge <cif file> <cif file> [...]

Example:
    mmcif_merge 1abc.cif 1xyz.cif > merged.cif

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cif, _cli
from mmciftools.mmcif_mkensemble import union_columns


def check_input(args):
    """Return (paths, identifier scheme)."""
    paths, _, scheme = _cli.parse_files(args, __doc__)
    return paths, scheme


def run(paths, scheme=None, warnings=None):
    """Yield the events of the merged file.

    `warnings`, if given, is a list that receives the chain clash messages.
    """
    columns = union_columns(paths)
    first_atoms = None
    max_id = 0
    chains_by_file = []

    def note(values, chain, chains):
        nonlocal max_id
        chains.add(values[chain])
        if id_pos is not None:
            num = _atoms.as_int(values[id_pos])
            if num is not None and num > max_id:
                max_id = num

    with _cif.open_input(paths[0]) as fh:
        for ev in _cif.parse(fh):
            etype = type(ev)
            if etype is _cif.Loop and first_atoms is None:
                first_atoms = ev
                loop = ev if [c.lower() for c in columns] == [c.lower() for c in ev.columns] \
                    else ev.with_columns(columns)
                fields = _atoms.Fields(ev, scheme)
                fields.need('chain')
                chain, id_pos = fields.chain, fields.id
                chains = set()
                yield loop
                continue
            if etype is _cif.Row and ev.loop is first_atoms:
                note(ev.values, chain, chains)
                yield ev if loop is first_atoms else ev.moved(loop)
                continue
            if etype is _cif.LoopEnd and ev.loop is first_atoms:
                chains_by_file.append(chains)
                for ev2 in _other_atoms(paths[1:], loop, scheme, chains_by_file, lambda: max_id):
                    yield ev2
                yield ev
                continue
            yield ev

    seen = {}
    for n, chains in enumerate(chains_by_file):
        for chain_id in sorted(chains):
            seen.setdefault(chain_id, []).append(paths[n])
    for chain_id, files in sorted(seen.items()):
        if len(files) > 1 and warnings is not None:
            warnings.append('WARNING: chain {} is in several files: {}\n'.format(
                chain_id, ', '.join(files)))


def _other_atoms(paths, loop, scheme, chains_by_file, last_id):
    """Rows of the atoms of the other files, moved into `loop`, renumbered."""
    number = last_id()
    for path in paths:
        chains = set()
        with _cif.open_input(path) as fh:
            source_loop = None
            for ev in _cif.parse(fh):
                if type(ev) is _cif.Loop and source_loop is None:
                    source_loop = ev
                    fields = _atoms.Fields(ev, scheme)
                    fields.need('chain')
                    chain = fields.chain
                elif type(ev) is _cif.Row and ev.loop is source_loop:
                    chains.add(ev.values[chain])
                    number += 1
                    yield ev.moved(loop, {'id': str(number)})
        chains_by_file.append(chains)


def main():
    paths, scheme = check_input(sys.argv[1:])
    warnings = []

    def events():
        for ev in run(paths, scheme, warnings):
            yield ev
        for message in warnings:
            sys.stderr.write(message)
    _cli.output(events())


if __name__ == '__main__':
    main()
