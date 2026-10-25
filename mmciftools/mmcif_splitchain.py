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
categories are kept, the data block is renamed <block>_<chain>, and rows of
_atom_site_anisotrop and _struct_conn of other chains are removed.

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

from mmciftools import _atoms, _cif, _cli


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
    """(chain ids in order, link routing, whether the file has anisotropic rows).

    Link routing is (columns, {partner identity: chains of the atoms that
    have it}) for the partners of _struct_conn rows, as filter_atoms judges
    them; None when there are no links to route.
    """
    conn_loop, conns, has_aniso = None, [], False
    for ev in _cif.parse(source.lines(), (_atoms.STRUCT_CONN, _atoms.ANISOTROP)):
        if type(ev) is _cif.Row and ev.loop.category.lower() == _atoms.STRUCT_CONN:
            conns.append(ev)
        elif type(ev) is _cif.Loop:
            if ev.category.lower() == _atoms.STRUCT_CONN:
                conn_loop = ev
            else:
                has_aniso = True
    chains, owners, columns = {}, {}, None
    for ev in _cif.parse(source.lines()):
        if type(ev) is _cif.Row:
            values = ev.values
            chains.setdefault(values[chain], None)
            if atom_key is not None:
                key = atom_key(values)
                if key in wanted:
                    owners.setdefault(key, set()).add(values[chain])
        elif type(ev) is _cif.Loop:
            fields = _atoms.Fields(ev, scheme)
            fields.need('chain')
            chain = fields.chain
            atom_key = None
            if conns and columns is None:
                columns = _atoms._link_columns(ev, conn_loop)
                if columns:
                    atom_key = _atoms.getter(*[ev.index(a) for a, _ in columns])
                    first, second = _atoms._partner_getters(conn_loop, columns)
                    wanted = {g(r.values) for r in conns for g in (first, second)}
    routing = (columns, owners) if conns and columns else None
    return list(chains), routing, has_aniso


def _write_batch(source, batch, scheme, routing, has_aniso):
    """Write the files of `batch` ({chain: path}) in one reading of the input."""
    handles = {c: open(p, 'w', encoding=_cif.ENCODING, errors=_cif.ERRORS, newline='')
               for c, p in batch.items()}
    writers = {c: _cif.Writer(fh) for c, fh in handles.items()}
    everyone = list(writers.values())
    if routing is not None:
        columns, owners = routing
    try:
        kind, id_chain = None, {}
        events = _cif.parse(source.lines(), (_atoms.ATOM_SITE, _atoms.ANISOTROP, _atoms.STRUCT_CONN))
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
                elif kind == _atoms.STRUCT_CONN and routing is not None:
                    allowed = None  # chains whose files keep this link (None: all)
                    for partner in (first(ev.values), second(ev.values)):
                        found = owners.get(partner)
                        if found is not None:  # a partner that exists in the input
                            allowed = found if allowed is None else allowed & found
                    for c, writer in writers.items():
                        if allowed is None or c in allowed:
                            writer.write(ev)
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
                elif kind == _atoms.STRUCT_CONN and routing is not None:
                    first, second = _atoms._partner_getters(ev, columns)
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
    chains, routing, has_aniso = survey(source, scheme)
    names = file_names(base, chains)
    if notes is not None:
        for chain_id, name in names.items():
            if '~' in name[len(base) + len(chain_id) + 1:]:
                notes.append('chain {} written to {}: its name differs from another chain only '
                             'in case'.format(chain_id, name))
    paths = {c: os.path.join(outdir, n) for c, n in names.items()}
    for start in range(0, len(chains), BATCH):
        batch = {c: paths[c] for c in chains[start:start + BATCH]}
        _write_batch(source, batch, scheme, routing, has_aniso)
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
