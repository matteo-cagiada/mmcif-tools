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

"""Helpers shared by the test modules (not a test module itself)."""

import io
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(HERE, 'data')
REAL = os.path.join(HERE, 'real')
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from mmciftools import _cif  # noqa: E402

FIXTURES = sorted(f for f in os.listdir(DATA) if f.endswith('.cif'))
INVALID_CIF = {'empty_atom_site.cif'}  # a loop with no rows (seen in the wild)

MODES = [('atom_site', _cif.DEFAULT_CATEGORIES), ('all', _cif.ALL)]

# Replacement values for the edit-everything tests: shorter, longer, quoted,
# null, empty, keyword-like and multi-line values.
CHOICES = ['X', 'two words', "C1'", 'multi\nline', '?', '.', 'a_much_longer_value_than_before',
           '', '"q"', 'data_z', 'mixed \' and " quotes', '#hash', '-1.5']


def read(name):
    """Text of a fixture in tests/data."""
    with open(os.path.join(DATA, name), encoding='utf-8', newline='') as fh:
        return fh.read()


def real_entries():
    """Paths of the real PDB entries fetched by tests/fetch_real.py."""
    if not os.path.isdir(REAL):
        return []
    return sorted(os.path.join(REAL, f) for f in os.listdir(REAL) if f.endswith('.cif.gz'))


def read_path(path):
    with _cif.open_input(path) as fh:
        return fh.read()


def roundtrip(text, categories=_cif.DEFAULT_CATEGORIES, edit=None):
    """Parse `text`, optionally pass the events through `edit`, write back."""
    events = _cif.parse(io.StringIO(text, newline=''), categories)
    if edit is not None:
        events = edit(events)
    out = io.StringIO()
    _cif.write(events, out)
    return out.getvalue()


def all_rows(text):
    return [(ev.loop.category, list(ev.values))
            for ev in _cif.parse(io.StringIO(text, newline=''), _cif.ALL)
            if isinstance(ev, _cif.Row)]


def edit_everything(text):
    """Replace every value of every row; return (output, expected rows)."""
    expected = []

    def edit(events):
        for ev in events:
            if isinstance(ev, _cif.Row):
                n = len(expected)
                for j, col in enumerate(ev.loop.columns):
                    ev[col] = CHOICES[(n + j) % len(CHOICES)]
                expected.append((ev.loop.category, list(ev.values)))
            yield ev
    return roundtrip(text, _cif.ALL, edit), expected


def rows(text, category='_atom_site'):
    return [ev for ev in _cif.parse(io.StringIO(text, newline=''), category)
            if isinstance(ev, _cif.Row) and ev.loop.category.lower() == category.lower()]


def edit_rows(func, category='_atom_site'):
    """Build an `edit` callable applying func(row) -> row or None to rows."""
    def edit(events):
        for ev in events:
            if isinstance(ev, _cif.Row) and ev.loop.category == category:
                ev = func(ev)
                if ev is None:
                    continue
            yield ev
    return edit


def changed_lines(before, after):
    a, b = before.splitlines(True), after.splitlines(True)
    if len(a) != len(b):
        return None
    return [(x, y) for x, y in zip(a, b) if x != y]


def run_tool(module, text, *args, scheme=None):
    """Run a cif tool's run() on text; return its output text."""
    result = module.run(io.StringIO(text, newline=''), *args, scheme=scheme)
    if module.__name__.endswith(('mmcif_wc', 'mmcif_tofasta')):
        return ''.join(result)
    out = io.StringIO()
    _cif.write(result, out)
    return out.getvalue()


def reference_problems(text):
    """Set of messages about rows of _atom_site_anisotrop and _struct_conn
    that match no atom (empty when consistent)."""
    from mmciftools import _atoms
    return set(_atoms.reference_problems(io.StringIO(text, newline='')))


def our_tables(text):
    """{(block, category): [row values]} for every category, nulls as ''."""
    tables, block = {}, None
    for ev in _cif.parse(io.StringIO(text, newline=''), _cif.ALL):
        if isinstance(ev, _cif.DataBlock):
            block = ev.name
        elif isinstance(ev, _cif.Row):
            vals = ['' if v in ('.', '?') else v for v in ev.values]
            tables.setdefault((block, ev.loop.category.lower()), []).append(vals)
    return tables


def split_by_filtering(text, scheme=None):
    """{chain: file text} as mmcif_splitchain wrote it before phase 5: one
    filter_atoms pass per chain (the reference for the batched writer)."""
    from mmciftools import _atoms, mmcif_splitchain
    source = _atoms.Source(io.StringIO(text, newline=''))
    out = {}
    for chain_id in mmcif_splitchain.chains_of(source, scheme):
        def make_keep(fields, chain_id=chain_id):
            chain = fields.chain
            return lambda values: values[chain] == chain_id
        buf = io.StringIO()
        _cif.write(_atoms.rename_block(_atoms.filter_atoms(source, make_keep, scheme), '_' + chain_id), buf)
        out[chain_id] = buf.getvalue()
    return out


def split_files(text, scheme=None):
    """{chain: file text} written by mmcif_splitchain.run()."""
    import shutil
    import tempfile
    from mmciftools import _atoms, mmcif_splitchain
    tmp = tempfile.mkdtemp()
    try:
        paths = mmcif_splitchain.run(io.StringIO(text, newline=''), 'x', scheme, outdir=tmp)
        chains = mmcif_splitchain.chains_of(_atoms.Source(io.StringIO(text, newline='')), scheme)
        out = {}
        for chain_id, path in zip(chains, paths):
            with open(path, encoding='utf-8', errors='surrogateescape', newline='') as fh:
                out[chain_id] = fh.read()
        return out
    finally:
        shutil.rmtree(tmp)
