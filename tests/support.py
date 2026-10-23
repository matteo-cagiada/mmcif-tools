# -*- coding: utf-8 -*-
#
# Copyright 2026 the cif-tools authors
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

from ciftools import _cif  # noqa: E402

try:
    import gemmi
except ImportError:  # optional cross-check only
    gemmi = None

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
            if isinstance(ev, _cif.Row) and ev.loop.category == category]


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
    if module.__name__.endswith(('cif_wc', 'cif_tofasta')):
        return ''.join(result)
    out = io.StringIO()
    _cif.write(result, out)
    return out.getvalue()


def reference_problems(text):
    """Rows of _atom_site_anisotrop and _struct_conn that do not match an atom.

    Returns a set of strings describing them (empty when consistent).
    """
    from ciftools import _atoms
    events = list(_cif.parse(io.StringIO(text, newline=''),
                             (_atoms.ATOM_SITE, _atoms.ANISOTROP, _atoms.STRUCT_CONN)))
    atoms = [ev for ev in events if isinstance(ev, _cif.Row) and ev.loop.category == '_atom_site']
    problems = set()
    by_id = {a['id']: a for a in atoms if 'id' in a}
    for ev in events:
        if not isinstance(ev, _cif.Row):
            continue
        if ev.loop.category == '_atom_site_anisotrop':
            atom = by_id.get(ev.get('id'))
            if atom is None:
                problems.add('anisotrop {} has no atom'.format(ev.get('id')))
                continue
            for ref, col in (('pdbx_auth_asym_id', 'auth_asym_id'), ('pdbx_auth_seq_id', 'auth_seq_id'),
                             ('pdbx_label_asym_id', 'label_asym_id')):
                if ref in ev and col in atom and ev[ref] != atom[col]:
                    problems.add('anisotrop {} {}={} but atom has {}'.format(
                        ev['id'], ref, ev[ref], atom[col]))
    conns = [ev for ev in events if isinstance(ev, _cif.Row) and ev.loop.category == '_struct_conn']
    if conns and atoms:
        columns = _atoms._link_columns(atoms[0].loop, conns[0].loop)
        if columns:
            atom_key = _atoms.getter(*[atoms[0].loop.index(a) for a, _ in columns])
            present = {atom_key(a.values) for a in atoms}
            first, second = _atoms._partner_getters(conns[0].loop, columns)
            for conn in conns:
                for n, partner in ((1, first), (2, second)):
                    if partner(conn.values) not in present:
                        problems.add('struct_conn {} partner {} missing'.format(conn.get('id'), n))
    return problems


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


def gemmi_tables(text):
    """The same as our_tables, read by gemmi."""
    tables = {}
    for block in gemmi.cif.read_string(text):
        for item in block:
            if item.loop is not None:
                tags, loop = item.loop.tags, item.loop
                cat = tags[0].split('.')[0].lower()
                width = loop.width()
                for r in range(loop.length()):
                    vals = [gemmi.cif.as_string(loop[r, c]) for c in range(width)]
                    tables.setdefault((block.name, cat), []).append(vals)
            elif item.pair is not None:
                tag, value = item.pair
                cat = tag.split('.')[0].lower()
                rows_ = tables.setdefault((block.name, cat), [[]])
                rows_[0].append(gemmi.cif.as_string(value))
    return tables
