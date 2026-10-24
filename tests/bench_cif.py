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
Benchmark the mmCIF streaming layer on a synthetic structure.

The structure is generated in memory, line by line (nothing is downloaded
or written to disk), so the memory numbers reflect the parser alone.

Usage:
    python tests/bench_cif.py [n_atoms]     (default: 1,000,000)
"""

import os
import sys
import time
import tracemalloc

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mmciftools import _cif  # noqa: E402

COLUMNS = [
    'group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id', 'label_comp_id',
    'label_asym_id', 'label_entity_id', 'label_seq_id', 'pdbx_PDB_ins_code', 'Cartn_x',
    'Cartn_y', 'Cartn_z', 'occupancy', 'B_iso_or_equiv', 'pdbx_formal_charge',
    'auth_seq_id', 'auth_comp_id', 'auth_asym_id', 'auth_atom_id', 'pdbx_PDB_model_num',
]
ATOMS = [('N', 'N'), ('C', 'CA'), ('C', 'C'), ('O', 'O'), ('C', 'CB')]


def structure(n_atoms):
    """Yield the lines of a wwPDB-style file with n_atoms atoms."""
    for line in ('data_BENCH\n', '#\n', '_entry.id BENCH\n', '#\n', 'loop_\n'):
        yield line
    for col in COLUMNS:
        yield '_atom_site.{} \n'.format(col)
    for i in range(n_atoms):
        res = i // len(ATOMS) + 1
        chain = 'ABCDEFGH'[(res // 500) % 8]
        elem, name = ATOMS[i % len(ATOMS)]
        x, y, z = (i % 977) * 0.137, (i % 613) * 0.211, (i % 389) * 0.307
        yield ('ATOM   {:<7d} {} {:<3s} . ALA {} 1 {:<5d} ? {:<8.3f} {:<8.3f} {:<8.3f} 1.00 '
               '{:<6.2f} ? {:<5d} ALA {} {:<3s} 1 \n').format(
                   i + 1, elem, name, chain, res, x, y, z, 20 + i % 30, res, chain, name)
    yield '#\n'


class NullSink:
    """A text stream that only counts what is written."""

    def __init__(self):
        self.chars = 0

    def write(self, text):
        self.chars += len(text)


def edit_bfactors(events):
    """Like mmcif_b: one value changed in every row."""
    for ev in events:
        if isinstance(ev, _cif.Row):
            ev['B_iso_or_equiv'] = '10.00'
        yield ev


def edit_ids(events):
    """Like mmcif_rplchain + mmcif_reres: three values changed in every row."""
    for ev in events:
        if isinstance(ev, _cif.Row):
            ev['auth_asym_id'] = 'Z'
            ev['label_asym_id'] = 'Z'
            ev['auth_seq_id'] = str(int(ev['auth_seq_id']) + 1000)
        yield ev


MODES = {'copy': None, 'edit 1 column': edit_bfactors, 'edit 3 columns': edit_ids}


def run(n_atoms, mode):
    events = _cif.parse(structure(n_atoms))
    if MODES[mode] is not None:
        events = MODES[mode](events)
    sink = NullSink()
    _cif.write(events, sink)
    return sink.chars


def main():
    n_atoms = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000
    n_chars = sum(len(line) for line in structure(n_atoms))
    print('Synthetic structure: {:,} atoms, {:.0f} MB of CIF text'.format(n_atoms, n_chars / 1e6))

    t0 = time.perf_counter()
    for _ in structure(n_atoms):
        pass
    t_gen = time.perf_counter() - t0
    print('  {:<16s} {:6.2f} s   (generating the lines; subtracted below)'.format('baseline', t_gen))

    for mode in MODES:
        t0 = time.perf_counter()
        chars = run(n_atoms, mode)
        elapsed = time.perf_counter() - t0 - t_gen
        if mode == 'copy':
            assert chars == n_chars, 'copy round-trip changed the size'
        print('  {:<16s} {:6.2f} s   {:5.2f} us/atom'.format(mode, elapsed, 1e6 * elapsed / n_atoms))

    # Peak memory must not grow with the number of atoms (streaming).
    for n in (n_atoms // 10, n_atoms):
        tracemalloc.start()
        run(n, 'copy')
        peak = tracemalloc.get_traced_memory()[1]
        tracemalloc.stop()
        print('  peak memory, {:>9,} atoms: {:6.2f} MB'.format(n, peak / 1e6))


if __name__ == '__main__':
    main()
