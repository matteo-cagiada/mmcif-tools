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
Splits an mmCIF file into several, each containing one chain.

As pdb_splitchain, writes <name>_<chain>.cif files in the current directory
(<name> is the input file name without its extension, or 'splitchains' when
reading stdin), overwriting existing ones, and prints nothing. Unlike
pdb_splitchain, each file is a complete mmCIF file: the metadata
categories are kept, the data block is renamed <block>_<chain>, and the
rows of other chains are removed (_atom_site_anisotrop, _struct_conn,
_struct_conf, _pdbx_poly_seq_scheme, the molecules with no atoms left, ...),
as mmcif_selchain does.

Chains are auth_asym_id (label_asym_id with --label).

File names are unique ignoring case: on macOS and Windows, 1abc_A.cif and
1abc_a.cif would be the same file, so a chain whose file name differs from
an earlier one only in case gets a number (1abc_a~2.cif), with a note on
stderr. pdb_splitchain overwrites the file instead.

Usage:
    mmcif_splitchain <cif file>

Example:
    mmcif_splitchain 1abc.cif   # writes 1abc_A.cif, 1abc_B.cif, ...

This program is part of the mmcif-tools suite.
"""

import os
import sys

from mmciftools import _atoms, _cif, _cli, _refs


def check_input(args):
    """Return (file handle, output name base, identifier scheme)."""
    fhandle, _, scheme = _cli.parse_args(args, __doc__, takes_option=False)
    path = next((a for a in args if not a.startswith('-')), None)
    return fhandle, _atoms.base_name(path, 'splitchains'), scheme


BATCH = 100  # output files written in one reading of the input


def chains_of(source, scheme=None):
    """Chain ids in order of first appearance."""
    return survey(source, scheme)[0]


def survey(source, scheme=None):
    """(chain ids in order, _refs.Removal with the chains as owners, whether
    the file has anisotropic rows)."""
    has_aniso = any(type(ev) is _cif.Loop for ev in _cif.parse(source.lines(), _atoms.ANISOTROP))
    chains = {}

    def make_owner(loop):
        fields = _atoms.Fields(loop, scheme)
        fields.need('chain')
        chain = fields.chain

        def owner(values):
            chains.setdefault(values[chain], None)
            return values[chain]
        return owner
    removal = _refs.Removal(source, make_owner)
    return list(chains), removal, has_aniso


def _write_batch(source, batch, scheme, removal, has_aniso):
    """Write the files of `batch` ({chain: path}) in one reading of the input."""
    handles = {c: open(p, 'w', encoding=_cif.ENCODING, errors=_cif.ERRORS, newline='')
               for c, p in batch.items()}
    writers = {c: _cif.Writer(fh) for c, fh in handles.items()}
    everyone = list(writers.values())
    views = dict(zip(batch, removal.views(list(batch))))
    following = set(_refs.FOLLOWING)
    try:
        kind, id_chain = None, {}
        events = _cif.parse(source.lines(), (_atoms.ATOM_SITE, _atoms.ANISOTROP) + _refs.FOLLOWING)
        for ev in events:
            etype = type(ev)
            if etype is _cif.Row:
                if kind == _atoms.ATOM_SITE:
                    values = ev.values
                    if has_aniso and id_pos is not None:
                        id_chain[values[id_pos]] = values[chain_pos]
                    writer = writers.get(values[chain_pos])
                    if writer is not None:
                        writer.write(ev)
                    continue
                if kind == _atoms.ANISOTROP and aniso_id is not None:
                    owner = id_chain.get(ev.values[aniso_id])
                    if owner is not None:  # the row of one atom goes with its chain
                        if owner in writers:
                            writers[owner].write(ev)
                        continue
                elif kind in following:
                    candidates = removal.candidates(ev)
                    for c, writer in writers.items():
                        if candidates is None or c in candidates:
                            row = _cif.Row(ev.loop, ev.values, ev.lead, ev.text, ev.lineno)  # each file edits its own
                            if views[c].row(row):
                                writer.write(row)
                    continue
                for writer in everyone:
                    writer.write(ev)
                continue
            if etype is _cif.DataBlock:
                for c, writer in writers.items():
                    name = ev.name + '_' + c
                    writer.write(_cif.DataBlock(name, 'data_{}\n'.format(name)))
                continue
            if etype is _cif.Loop:
                kind = ev.category.lower()
                if kind == _atoms.ATOM_SITE:
                    fields = _atoms.Fields(ev, scheme)
                    chain_pos, id_pos = fields.chain, fields.id
                elif kind == _atoms.ANISOTROP:
                    aniso_id = ev.index('id') if 'id' in ev else None
            elif etype is _cif.LoopEnd:
                kind = None
            for writer in everyone:
                writer.write(ev)
        for writer in everyone:
            writer.close()
    finally:
        for fh in handles.values():
            fh.close()


def file_names(base, chains):
    """{chain: file name}, unique ignoring case."""
    names, taken = {}, set()
    for chain_id in chains:
        name, n = '{}_{}'.format(base, chain_id), 1
        while name.lower() in taken:
            n += 1
            name = '{}_{}~{}'.format(base, chain_id, n)
        taken.add(name.lower())
        names[chain_id] = name + '.cif'
    return names


def run(fhandle, base, scheme=None, outdir='.', notes=None):
    """Write one file per chain; return their paths. `notes`, if given, is a
    list that receives a message for each file name changed to stay unique."""
    source = _atoms.as_source(fhandle)
    chains, removal, has_aniso = survey(source, scheme)
    names = file_names(base, chains)
    if notes is not None:
        for chain_id, name in names.items():
            if '~' in name[len(base) + len(chain_id) + 1:]:
                notes.append('chain {} written to {}: its name differs from another chain only '
                             'in case'.format(chain_id, name))
    paths = {c: os.path.join(outdir, n) for c, n in names.items()}
    for start in range(0, len(chains), BATCH):
        batch = {c: paths[c] for c in chains[start:start + BATCH]}
        _write_batch(source, batch, scheme, removal, has_aniso)
    return [paths[c] for c in chains]


def main():
    fhandle, base, scheme = check_input(sys.argv[1:])
    notes = []
    try:
        run(fhandle, base, scheme, notes=notes)
        for note in notes:
            sys.stderr.write('Note: {}\n'.format(note))
    except _cif.CIFError as err:
        sys.stderr.write('ERROR!! {}\n'.format(err))
        sys.exit(1)
    finally:
        fhandle.close()
    sys.exit(0)


if __name__ == '__main__':
    main()
