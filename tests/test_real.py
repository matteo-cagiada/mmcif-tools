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
Tests on real PDB entries (tests/real_set.txt).

Skipped unless the entries have been fetched with `python tests/fetch_real.py`.
"""

import io
import os
import time
import unittest

from support import (  # noqa: E402
    _cif, all_rows, edit_everything, gemmi, gemmi_tables, our_tables, read_path,
    real_entries, reference_problems, roundtrip, rows, run_tool)

from ciftools import (  # noqa: E402
    _atoms, cif_b, cif_chain, cif_chkensemble, cif_delchain, cif_delelem, cif_delhetatm,
    cif_delres, cif_delresname, cif_element, cif_fixinsert, cif_frompdb, cif_gap, cif_head,
    cif_keepcoord, cif_mkensemble, cif_occ, cif_reatom, cif_reres, cif_rplchain, cif_rplresname,
    cif_selaltloc, cif_selatom, cif_selchain, cif_selelem, cif_selhetatm, cif_selmodel, cif_selres,
    cif_selresname, cif_shiftres, cif_sort, cif_splitmodel, cif_tidy, cif_tofasta, cif_topdb,
    cif_uniqname, cif_validate, cif_wc)

ENTRIES = real_entries()
SKIP = 'real entries not fetched (run: python tests/fetch_real.py)'
GEMMI_MAX_ATOMS = 30000  # gemmi re-reads of larger outputs take too long for a test run
INF = float('inf')

# Pairs of tools that must split the atoms between them.
COMPLEMENTS = [('cif_selchain', 'cif_delchain'), ('cif_selresname', 'cif_delresname'),
               ('cif_selhetatm', 'cif_delhetatm'), ('cif_selelem', 'cif_delelem')]
# Tools that keep every atom
KEEP_ALL = ('cif_chain', 'cif_rplchain', 'cif_b', 'cif_occ', 'cif_reres', 'cif_shiftres',
            'cif_keepcoord', 'cif_element -all', 'cif_rplresname', 'cif_fixinsert', 'cif_reatom',
            'cif_uniqname', 'cif_sort')
# Entries whose chains have no missing residues: PDB -> mmCIF gives back the label ids
COMPLETE = ('1d3z', '1ehz', '1ejg', '1ppb', '4hhb')


def name(path):
    return os.path.basename(path).split('.')[0]


def calls_for(atoms):
    """Tool invocations suited to a structure: [(label, module, args)]."""
    chains = list(dict.fromkeys(a['auth_asym_id'] for a in atoms))
    numbers = [n for n in (_atoms.as_int(a['auth_seq_id']) for a in atoms) if n is not None]
    first = min(numbers)
    return [
        ('cif_selchain', cif_selchain, ([chains[0]],)),
        ('cif_delchain', cif_delchain, ([chains[0]],)),
        ('cif_selres range', cif_selres, (set(), [(first, first + 10, 1)])),
        ('cif_selres step', cif_selres, ({first}, [(-INF, INF, 3)])),
        ('cif_delres range', cif_delres, ((first, first + 10, 1),)),
        ('cif_delres step', cif_delres, ((-INF, INF, 4),)),
        ('cif_selresname', cif_selresname, (['HOH'],)),
        ('cif_delresname', cif_delresname, (['HOH'],)),
        ('cif_selatom', cif_selatom, (['CA', 'P'],)),
        ('cif_selelem', cif_selelem, (['C'],)),
        ('cif_selhetatm', cif_selhetatm, ()),
        ('cif_delhetatm', cif_delhetatm, ()),
        ('cif_selaltloc', cif_selaltloc, (None,)),
        ('cif_selaltloc A', cif_selaltloc, ('A',)),
        ('cif_selmodel', cif_selmodel, ([1],)),
        ('cif_chain', cif_chain, ('Z',)),
        ('cif_rplchain', cif_rplchain, ((chains[0], 'ZZ'),)),
        ('cif_b', cif_b, (1.5,)),
        ('cif_occ', cif_occ, (0.5,)),
        ('cif_reres', cif_reres, (1,)),
        ('cif_shiftres', cif_shiftres, (-1000,)),
        ('cif_delelem', cif_delelem, (['C'],)),
        ('cif_keepcoord', cif_keepcoord, ()),
        ('cif_head', cif_head, (100,)),
        ('cif_element -all', cif_element, (True,)),
        ('cif_rplresname', cif_rplresname, (('HOH', 'WAT'),)),
        ('cif_fixinsert', cif_fixinsert, (set(),)),
        ('cif_reatom', cif_reatom, (1000,)),
        ('cif_uniqname', cif_uniqname, ()),
        ('cif_sort', cif_sort, ('CR',)),
    ]


@unittest.skipUnless(ENTRIES, SKIP)
class TestRealLayer(unittest.TestCase):
    """The shared mmCIF layer on real files."""

    @classmethod
    def setUpClass(cls):
        cls.texts = {name(p): read_path(p) for p in ENTRIES}

    def test_round_trip(self):
        for entry, text in self.texts.items():
            for categories in (_cif.DEFAULT_CATEGORIES, _cif.ALL):
                with self.subTest(entry=entry, categories=categories):
                    self.assertEqual(roundtrip(text, categories), text)

    def test_every_value_edited(self):
        for entry, text in self.texts.items():
            with self.subTest(entry=entry):
                out, expected = edit_everything(text)
                self.assertEqual(all_rows(out), expected)

    @unittest.skipUnless(gemmi, 'gemmi not installed (optional cross-check)')
    def test_same_values_as_gemmi(self):
        for entry, text in self.texts.items():
            with self.subTest(entry=entry):
                self.assertEqual(our_tables(text), gemmi_tables(text))

    @unittest.skipUnless(gemmi, 'gemmi not installed (optional cross-check)')
    def test_edited_output_is_valid(self):
        for entry, text in self.texts.items():
            with self.subTest(entry=entry):
                out, _ = edit_everything(text)
                self.assertEqual(our_tables(out), gemmi_tables(out))


@unittest.skipUnless(ENTRIES, SKIP)
class TestRealTools(unittest.TestCase):
    """Every tool on every real entry, checked for consistency."""

    @classmethod
    def setUpClass(cls):
        cls.texts = {name(p): read_path(p) for p in ENTRIES}
        cls.timings = []

    @classmethod
    def tearDownClass(cls):
        if os.environ.get('CIFTOOLS_TIMINGS'):
            for entry, label, seconds in sorted(cls.timings, key=lambda t: -t[2])[:15]:
                print('{:6.2f}s  {:5s} {}'.format(seconds, entry, label))

    def test_tools(self):
        for entry, text in self.texts.items():
            atoms = rows(text)
            before = reference_problems(text)
            outputs = {}
            for label, module, args in calls_for(atoms):
                with self.subTest(entry=entry, tool=label):
                    start = time.perf_counter()
                    out = run_tool(module, text, *args)
                    self.timings.append((entry, label, time.perf_counter() - start))
                    outputs[label] = out
                    self.assertEqual(roundtrip(out, _cif.ALL), out)  # output re-reads
                    self.assertLessEqual(reference_problems(out), before)
                    if label in KEEP_ALL:
                        self.assertEqual(len(rows(out)), len(atoms))
                    if gemmi and len(atoms) <= GEMMI_MAX_ATOMS:
                        self.assertEqual(our_tables(out), gemmi_tables(out))
            for sel, dele in COMPLEMENTS:
                with self.subTest(entry=entry, pair=sel):
                    self.assertEqual(len(rows(outputs[sel])) + len(rows(outputs[dele])), len(atoms))

    def test_reports(self):
        for entry, text in self.texts.items():
            with self.subTest(entry=entry):
                summary = run_tool(cif_wc, text, '')
                self.assertTrue(summary.startswith('No. models:\t'))
                fasta = run_tool(cif_tofasta, text, False, False)
                self.assertTrue(fasta.startswith('>' + entry.upper() + '|'))
                seqres = run_tool(cif_tofasta, text, True, True)
                self.assertTrue(seqres.startswith('>' + entry.upper() + '|'))

    def test_element_matches_deposited(self):
        for entry, text in self.texts.items():
            with self.subTest(entry=entry):
                self.assertEqual(run_tool(cif_element, text, True), text)

    def test_gap(self):
        for entry, text in self.texts.items():
            with self.subTest(entry=entry):
                report = run_tool(cif_gap, text)
                self.assertTrue(report.endswith(' gap(s) in the structure\n'))
        self.assertTrue(run_tool(cif_gap, self.texts['6vxx']).startswith('A:'))

    def test_pdb_round_trip(self):
        auth = ['auth_asym_id', 'auth_seq_id', 'auth_comp_id', 'auth_atom_id', 'pdbx_PDB_ins_code',
                'label_alt_id', 'pdbx_PDB_model_num', 'type_symbol']
        label = ['label_asym_id', 'label_entity_id', 'label_seq_id']
        for entry, text in self.texts.items():
            with self.subTest(entry=entry):
                try:
                    pdb = run_tool(cif_topdb, text)
                except _cif.CIFError:
                    self.assertEqual(entry, '4v6x')  # multi-character chains
                    continue
                back = ''.join(cif_frompdb.run(io.StringIO(pdb), entry))
                before, after = rows(text), rows(back)
                self.assertEqual(len(after), len(before))
                for col in auth + (label if entry in COMPLETE else []):
                    self.assertEqual([r[col] for r in after], [r[col] for r in before], col)

    def test_split_and_rebuild_ensemble(self):
        import shutil
        import tempfile
        tmp = tempfile.mkdtemp()
        try:
            text = self.texts['1d3z']
            paths = cif_splitmodel.run(io.StringIO(text), '1d3z', outdir=tmp)
            self.assertEqual(len(paths), 10)
            out = io.StringIO()
            _cif.write(cif_mkensemble.run(paths), out)
            self.assertEqual(len(rows(out.getvalue())), len(rows(text)))
            status = cif_chkensemble.run(io.StringIO(out.getvalue()))[2]
            self.assertEqual(status, 0)
        finally:
            shutil.rmtree(tmp)

    def test_validate(self):
        for entry, text in self.texts.items():
            with self.subTest(entry=entry):
                lines, errors = cif_validate.run(io.StringIO(text, newline=''))
                self.assertEqual(errors, 0, ''.join(lines))

    def test_tidy(self):
        for entry, text in self.texts.items():
            with self.subTest(entry=entry):
                out = run_tool(cif_tidy, text)
                self.assertEqual(our_tables(out), our_tables(text))  # same data
                self.assertEqual(run_tool(cif_tidy, out), out)       # already tidy


if __name__ == '__main__':
    unittest.main()
