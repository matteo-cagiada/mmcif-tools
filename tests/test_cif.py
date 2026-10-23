#!/usr/bin/env python
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

"""
Unit tests for the shared mmCIF layer (ciftools/_cif.py).

Only the standard library is required; the gemmi cross-checks are skipped
when gemmi is not installed.
"""

import gzip
import io
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(HERE, 'data')
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
    with open(os.path.join(DATA, name), encoding='utf-8', newline='') as fh:
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


class TestRoundTrip(unittest.TestCase):
    """The phase-0 criterion: unmodified files come back byte-for-byte."""

    def assertRoundTrips(self, text):
        for mode, categories in MODES:
            with self.subTest(mode=mode):
                self.assertEqual(roundtrip(text, categories), text)

    def test_fixtures(self):
        for name in FIXTURES:
            with self.subTest(fixture=name):
                self.assertRoundTrips(read(name))

    def test_crlf(self):
        for name in ('wwpdb_like.cif', 'text_fields.cif', 'pairs_atom_site.cif'):
            with self.subTest(fixture=name):
                self.assertRoundTrips(read(name).replace('\n', '\r\n'))

    def test_no_final_newline(self):
        for name in FIXTURES:
            with self.subTest(fixture=name):
                self.assertRoundTrips(read(name).rstrip('\n'))

    def test_gzip_input(self):
        text = read('wwpdb_like.cif')
        tmp = tempfile.mkdtemp()
        try:
            path = os.path.join(tmp, 'x.cif.gz')
            with gzip.open(path, 'wt', encoding='utf-8', newline='') as fh:
                fh.write(text)
            with _cif.open_input(path) as fh:
                out = io.StringIO()
                _cif.write(_cif.parse(fh), out)
            self.assertEqual(out.getvalue(), text)
        finally:
            shutil.rmtree(tmp)

    def test_cli_stdin_stdout_bytes(self):
        for name in ('wwpdb_like.cif', 'text_fields.cif'):
            data = read(name).replace('\n', '\r\n').encode('utf-8')
            for args in ([], ['--all']):
                with self.subTest(fixture=name, args=args):
                    proc = subprocess.run(
                        [sys.executable, '-m', 'ciftools._cif'] + args,
                        input=data, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                        cwd=ROOT, check=True)
                    self.assertEqual(proc.stdout, data)
                    self.assertEqual(proc.stderr, b'')

    def test_cli_reports_errors(self):
        bad = b"data_x\nloop_\n_atom_site.id\n_atom_site.type_symbol\n1 'C\n"
        proc = subprocess.run([sys.executable, '-m', 'ciftools._cif'], input=bad,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=ROOT)
        self.assertEqual(proc.returncode, 1)
        self.assertIn(b'line 5', proc.stderr)

    def test_library_does_not_import_gemmi(self):
        code = 'import sys; import ciftools._cif; assert "gemmi" not in sys.modules'
        subprocess.run([sys.executable, '-c', code], cwd=ROOT, check=True)


class TestReader(unittest.TestCase):

    def test_tokens(self):
        line = "a 'b c' \"d e\" 'N,N'-x' a#b O5' # a comment 'x\n"
        toks = [(v, k) for _, _, v, k in _cif._tokens(line)]
        self.assertEqual([v for v, _ in toks], ['a', 'b c', 'd e', "N,N'-x", 'a#b', "O5'"])
        self.assertTrue(all(k == _cif.VALUE for _, k in toks))

    def test_token_kinds(self):
        kinds = [k for _, _, _, k in _cif._tokens('_a.b loop_ data_x save_y stop_ x\n')]
        self.assertEqual(kinds, [_cif.TAG, _cif.KEYWORD, _cif.KEYWORD, _cif.KEYWORD,
                                 _cif.KEYWORD, _cif.VALUE])

    def test_unterminated_quote_has_line_number(self):
        text = "data_x\nloop_\n_atom_site.id\n_atom_site.label_atom_id\n1 N\n2 'CA\n"
        with self.assertRaises(_cif.CIFError) as ctx:
            rows(text)
        self.assertEqual(ctx.exception.lineno, 6)

    def test_incomplete_row(self):
        text = "data_x\nloop_\n_atom_site.id\n_atom_site.label_atom_id\n1 N\n2\n#\n_cell.a 1\n"
        with self.assertRaises(_cif.CIFError) as ctx:
            rows(text)
        self.assertEqual(ctx.exception.lineno, 6)

    def test_unterminated_text_field(self):
        with self.assertRaises(_cif.CIFError):
            roundtrip("data_x\n_struct.title\n;never closed\n", _cif.ALL)

    def test_values(self):
        atoms = rows(read('wwpdb_like.cif'))
        self.assertEqual(len(atoms), 19)
        first = atoms[0]
        self.assertEqual(first['auth_seq_id'], '10')
        self.assertEqual(first['label_seq_id'], '1')
        self.assertEqual(first['Cartn_x'], '21.598')
        self.assertEqual(first['CARTN_X'], '21.598')  # tags are case-insensitive
        self.assertEqual([a['pdbx_PDB_ins_code'] for a in atoms if a['label_seq_id'] == '3'],
                         ['A'] * 4)
        self.assertEqual([a['label_alt_id'] for a in atoms if a['auth_atom_id'] == 'OG'],
                         ['A', 'B'])
        self.assertEqual(atoms[-1]['label_seq_id'], '.')
        self.assertEqual(first.lineno, 68)

    def test_quoted_atom_names(self):
        names = [a['auth_atom_id'] for a in rows(read('nucleic.cif'))]
        self.assertEqual(names[2:5], ["O5'", "C5'", "C4'"])

    def test_text_field_layouts(self):
        atoms = rows(read('text_fields.cif'))
        self.assertEqual([a['id'] for a in atoms], [str(i) for i in range(1, 9)])
        details = [a['details'] for a in atoms]
        self.assertEqual(details, [
            '.', 'alpha carbon', 'row split over two lines', '?', '?',
            'a text field\ninside the atom_site loop', "double 'quoted'", 'a#b'])
        self.assertEqual(atoms[4].lead, ' ')  # second row on a shared line

    def test_keywords_inside_text_fields_are_ignored(self):
        text = read('text_fields.cif')
        struct = rows(text, '_struct')
        self.assertTrue(struct[0]['title'].startswith('Crystal structure'))
        self.assertIn('loop_ and data_', struct[0]['title'])
        citation = rows(text, '_citation')
        self.assertEqual([c['title'] for c in citation],
                         ['A title that spans\ntwo lines', "Second, with a quote: it''s here"])

    def test_pairs_category(self):
        atoms = rows(read('pairs_atom_site.cif'))
        self.assertEqual(len(atoms), 1)
        self.assertTrue(atoms[0].loop.pairs)
        self.assertEqual(atoms[0]['auth_seq_id'], '501')

    def test_events(self):
        events = list(_cif.parse(io.StringIO(read('multi_block.cif'))))
        blocks = [ev.name for ev in events if isinstance(ev, _cif.DataBlock)]
        self.assertEqual(blocks, ['first', 'second'])
        loops = [ev for ev in events if isinstance(ev, _cif.Loop)]
        self.assertEqual(len(loops), 2)
        ends = [ev for ev in events if isinstance(ev, _cif.LoopEnd)]
        self.assertEqual([e.nrows for e in ends], [2, 1])

    def test_only_requested_categories_are_parsed(self):
        events = list(_cif.parse(io.StringIO(read('wwpdb_like.cif'))))
        cats = {ev.category for ev in events if isinstance(ev, _cif.Loop)}
        self.assertEqual(cats, {'_atom_site'})
        events = list(_cif.parse(io.StringIO(read('wwpdb_like.cif')), 'atom_site_anisotrop'))
        cats = {ev.category for ev in events if isinstance(ev, _cif.Loop)}
        self.assertEqual(cats, {'_atom_site_anisotrop'})

    def test_column_order_does_not_matter(self):
        atoms = rows(read('column_order.cif'))
        self.assertEqual(atoms[0]['id'], '1')
        self.assertEqual(atoms[0]['Cartn_x'], '11.480922')
        self.assertEqual(atoms[3]['auth_comp_id'], 'GLN')


class TestQuote(unittest.TestCase):

    def test_cases(self):
        cases = {
            'CA': 'CA', '.': '.', '?': '?', '': "''", "O5'": '"O5\'"',
            'a b': "'a b'", 'data_x': "'data_x'", 'LOOP_': "'LOOP_'", '_x': "'_x'",
            '#x': "'#x'", '[O-]': "'[O-]'", ';x': "';x'", '$x': "'$x'",
            "it's here": '"it\'s here"', 'a\tb': "'a\tb'", 'x;y': 'x;y', 'a#b': 'a#b',
        }
        for value, expected in cases.items():
            with self.subTest(value=value):
                self.assertEqual(_cif.quote(value), expected)

    def test_both_quote_kinds_use_a_text_field(self):
        value = 'it\' s "both" kinds" here'
        self.assertEqual(_cif.quote(value), ';' + value + '\n;')

    def test_unrepresentable(self):
        with self.assertRaises(_cif.CIFError):
            _cif.quote('first\n;second')

    def test_fuzz(self):
        """Random strings survive quote() -> write -> parse."""
        rng = random.Random(1234)
        alphabet = 'aZ09._-+ \t\'";#$[]\n'
        for _ in range(5000):
            value = ''.join(rng.choice(alphabet) for _ in range(rng.randint(0, 8)))
            if '\n;' in value:
                with self.assertRaises(_cif.CIFError):
                    _cif.quote(value)
                continue
            text = 'data_t\nloop_\n_t.a\n_t.b\n{} end\n'.format(_cif.quote(value))
            with self.subTest(value=value):
                parsed = rows(text, '_t')
                self.assertEqual(len(parsed), 1)
                self.assertEqual(parsed[0]['a'], value)
                self.assertEqual(parsed[0]['b'], 'end')


class TestEdits(unittest.TestCase):

    def test_same_length_edit_changes_one_value(self):
        text = read('wwpdb_like.cif')

        def bump(row):
            if row['id'] == '5':
                row['B_iso_or_equiv'] = '99.99'
            return row
        diff = changed_lines(text, roundtrip(text, edit=edit_rows(bump)))
        self.assertEqual(len(diff), 1)
        old, new = diff[0]
        self.assertEqual(new, old.replace('26.47', '99.99'))

    def test_columns_stay_aligned(self):
        text = read('wwpdb_like.cif')

        def rename(row):
            row['auth_comp_id'] = 'X'          # shorter
            row['label_atom_id'] = 'CA2'       # longer, padding absorbs it
            return row
        out = roundtrip(text, edit=edit_rows(rename))
        for old, new in changed_lines(text, out):
            # model number (the last value) keeps its column
            self.assertEqual(old.rstrip().rindex(' '), new.rstrip().rindex(' '))
        self.assertEqual({a['auth_comp_id'] for a in rows(out)}, {'X'})

    def test_values_needing_quotes(self):
        text = read('wwpdb_like.cif')
        values = ["C1'", 'two words', 'data_x', 'multi\nline', '']

        def assign(row):
            i = int(row['id'])
            if i <= len(values):
                row['auth_atom_id'] = values[i - 1]
            return row
        out = roundtrip(text, edit=edit_rows(assign))
        got = [a['auth_atom_id'] for a in rows(out)][:len(values)]
        self.assertEqual(got, values)
        self.assertEqual(len(rows(out)), 19)

    def test_non_string_values_are_converted(self):
        row = rows(read('wwpdb_like.cif'))[0]
        row['occupancy'] = 0.5
        self.assertEqual(row['occupancy'], '0.5')

    def test_unmodified_assignment_keeps_raw_text(self):
        text = read('wwpdb_like.cif')

        def same(row):
            row['Cartn_x'] = row['Cartn_x']
            return row
        self.assertEqual(roundtrip(text, edit=edit_rows(same)), text)

    def test_deleting_all_rows_drops_the_loop(self):
        text = read('wwpdb_like.cif')
        out = roundtrip(text, edit=edit_rows(lambda row: None))
        self.assertNotIn('_atom_site.id', out)
        self.assertIn('_atom_site_anisotrop.id', out)
        self.assertIn('_cell.length_a', out)
        self.assertEqual(rows(out), [])

    def test_empty_input_loop_is_kept(self):
        text = read('empty_atom_site.cif')
        self.assertEqual(roundtrip(text, edit=edit_rows(lambda row: None)), text)

    def test_deleting_rows_that_share_a_line(self):
        text = read('text_fields.cif')
        for drop in ('4', '5', '6', '3'):
            with self.subTest(drop=drop):
                out = roundtrip(text, edit=edit_rows(lambda r: None if r['id'] == drop else r))
                ids = [a['id'] for a in rows(out)]
                self.assertEqual(ids, [str(i) for i in range(1, 9) if str(i) != drop])

    def test_text_field_in_edited_rows(self):
        text = read('text_fields.cif')

        def assign(row):
            row['details'] = 'now\nmultiline' if row['id'] in ('4', '7') else 'flat'
            return row
        out = roundtrip(text, edit=edit_rows(assign))
        details = {a['id']: a['details'] for a in rows(out)}
        self.assertEqual(details['4'], 'now\nmultiline')
        self.assertEqual(details['5'], 'flat')
        self.assertEqual(details['6'], 'flat')
        self.assertEqual(details['7'], 'now\nmultiline')

    def test_new_rows(self):
        text = read('nmr_models.cif')

        def add(events):
            for ev in events:
                if isinstance(ev, _cif.Loop):
                    loop = ev
                if isinstance(ev, _cif.LoopEnd):
                    values = ['ATOM', '7', 'C', 'CB', '.', 'ALA', 'A', '1', '1', '?', '1', '2',
                              '3', '1', '0', '?', '1', 'ALA', 'A', 'CB', '2']
                    yield _cif.Row(loop, values)
                yield ev
        out = roundtrip(text, edit=add)
        atoms = rows(out)
        self.assertEqual(len(atoms), 7)
        self.assertEqual(atoms[-1]['auth_atom_id'], 'CB')
        self.assertTrue(out.endswith('ALA A CB 2\n'))

    def test_rows_from_another_loop_are_mapped_by_column(self):
        text = read('multi_block.cif')
        small = _cif.Loop('_atom_site', ['id', 'auth_asym_id', 'group_PDB'])

        def add(events):
            for ev in events:
                if isinstance(ev, _cif.LoopEnd) and ev.nrows == 1:
                    yield _cif.Row(small, ['9', 'Z', 'HETATM'])
                yield ev
        atoms = rows(roundtrip(text, edit=add))
        self.assertEqual(atoms[-1]['id'], '9')
        self.assertEqual(atoms[-1]['auth_asym_id'], 'Z')
        self.assertEqual(atoms[-1]['Cartn_x'], '?')

    def test_every_value_edited(self):
        for name in FIXTURES:
            with self.subTest(fixture=name):
                out, expected = edit_everything(read(name))
                self.assertEqual(all_rows(out), expected)
                self.assertEqual(roundtrip(out, _cif.ALL), out)  # output is stable

    def test_pairs_edit_changes_one_line(self):
        text = read('pairs_atom_site.cif')

        def move(row):
            row['Cartn_x'] = '-10.5'
            return row
        diff = changed_lines(text, roundtrip(text, edit=edit_rows(move)))
        self.assertEqual(diff, [('_atom_site.Cartn_x           1.000\n',
                                 '_atom_site.Cartn_x           -10.5\n')])

    def test_pairs_with_second_row_become_a_loop(self):
        text = read('pairs_atom_site.cif')

        def duplicate(events):
            for ev in events:
                if isinstance(ev, _cif.Row):
                    yield ev
                    copy = _cif.Row(ev.loop, list(ev.values))
                    copy['id'] = '2'
                    copy['auth_seq_id'] = '502'
                    yield copy
                    continue
                yield ev
        out = roundtrip(text, edit=duplicate)
        self.assertIn('loop_\n_atom_site.group_PDB\n', out)
        atoms = rows(out)
        self.assertEqual([a['auth_seq_id'] for a in atoms], ['501', '502'])
        self.assertFalse(atoms[0].loop.pairs)


class TestScheme(unittest.TestCase):

    def setUp(self):
        patcher = mock.patch.dict(os.environ)
        patcher.start()
        self.addCleanup(patcher.stop)
        os.environ.pop('CIFTOOLS_SCHEME', None)

    def test_default_is_auth(self):
        self.assertEqual(_cif.get_scheme(), 'auth')

    def test_environment(self):
        os.environ['CIFTOOLS_SCHEME'] = 'LABEL'
        self.assertEqual(_cif.get_scheme(), 'label')
        self.assertEqual(_cif.get_scheme('auth'), 'auth')  # explicit wins

    def test_bad_scheme(self):
        with self.assertRaises(ValueError):
            _cif.get_scheme('segid')

    def test_column(self):
        loop = rows(read('wwpdb_like.cif'))[0].loop
        self.assertEqual(_cif.column(loop, 'chain'), 'auth_asym_id')
        self.assertEqual(_cif.column(loop, 'chain', 'label'), 'label_asym_id')
        self.assertEqual(_cif.column(loop, 'resnum', 'label'), 'label_seq_id')

    def test_column_falls_back(self):
        loop = rows(read('af3_like.cif'))[0].loop  # no auth_comp_id / auth_atom_id
        self.assertEqual(_cif.column(loop, 'resname'), 'label_comp_id')
        self.assertEqual(_cif.column(loop, 'atom'), 'label_atom_id')
        self.assertEqual(_cif.column(loop, 'chain'), 'auth_asym_id')

    def test_column_missing(self):
        loop = _cif.Loop('_atom_site', ['id'])
        with self.assertRaises(_cif.CIFError):
            _cif.column(loop, 'chain')


@unittest.skipUnless(gemmi, 'gemmi not installed (optional cross-check)')
class TestAgainstGemmi(unittest.TestCase):
    """Compare with an independent parser."""

    @staticmethod
    def ours(text):
        """{(block, category): [row values]} for every category, nulls as ''."""
        tables, block = {}, None
        for ev in _cif.parse(io.StringIO(text, newline=''), _cif.ALL):
            if isinstance(ev, _cif.DataBlock):
                block = ev.name
            elif isinstance(ev, _cif.Row):
                vals = ['' if v in ('.', '?') else v for v in ev.values]
                tables.setdefault((block, ev.loop.category.lower()), []).append(vals)
        return tables

    @staticmethod
    def theirs(text):
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

    def test_same_values_as_gemmi(self):
        for name in FIXTURES:
            if name in INVALID_CIF:
                continue
            with self.subTest(fixture=name):
                text = read(name)
                self.assertEqual(self.ours(text), self.theirs(text))

    def test_every_value_edited_is_valid(self):
        for name in FIXTURES:
            if name in INVALID_CIF:
                continue
            with self.subTest(fixture=name):
                out, _ = edit_everything(read(name))
                self.assertEqual(self.ours(out), self.theirs(out))

    def test_edited_output_is_valid(self):
        def scramble(row):
            if 'B_iso_or_equiv' in row:
                row['B_iso_or_equiv'] = '1.5'
            if 'auth_atom_id' in row:
                row['auth_atom_id'] = row['auth_atom_id'] + "'"
            if 'details' in row:
                row['details'] = 'two\nlines'
            return row
        for name in FIXTURES:
            if name in INVALID_CIF:
                continue
            with self.subTest(fixture=name):
                out = roundtrip(read(name), edit=edit_rows(scramble))
                self.assertEqual(self.ours(out), self.theirs(out))


if __name__ == '__main__':
    unittest.main()
