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
Tests on real PDB entries (tests/real_set.txt).

Skipped unless the entries have been fetched with `python tests/fetch_real.py`.
"""

import io
import os
import time
import unittest

from support import (  # noqa: E402
    _cif, all_rows, edit_everything, our_tables, read_path,
    real_entries, reference_problems, roundtrip, rows, run_tool)

from mmciftools import (  # noqa: E402
    _atoms, mmcif_assembly, mmcif_auth2label, mmcif_b, mmcif_blocks, mmcif_cat, mmcif_chain, mmcif_chkensemble,
    mmcif_delcat, mmcif_delchain, mmcif_delentity, mmcif_delelem, mmcif_delhetatm, mmcif_delres, mmcif_delresname,
    mmcif_element, mmcif_fixinsert, mmcif_info, mmcif_frompdb, mmcif_gap, mmcif_head, mmcif_keepcat,
    mmcif_keepcoord, mmcif_label2auth, mmcif_mkensemble, mmcif_occ, mmcif_reatom, mmcif_reres,
    mmcif_rplchain, mmcif_rplresname, mmcif_selaltloc, mmcif_selatom, mmcif_selchain, mmcif_selelem,
    mmcif_selentity, mmcif_selhetatm, mmcif_selmodel, mmcif_selres, mmcif_selresname, mmcif_shiftres, mmcif_sort,
    mmcif_splitmodel, mmcif_tidy, mmcif_tofasta, mmcif_topdb, mmcif_uniqname, mmcif_validate,
    mmcif_wc)

ENTRIES = real_entries()
SKIP = 'real entries not fetched (run: python tests/fetch_real.py)'
INF = float('inf')

# Pairs of tools that must split the atoms between them.
COMPLEMENTS = [('mmcif_selentity', 'mmcif_delentity'), ('mmcif_selchain', 'mmcif_delchain'), ('mmcif_selresname', 'mmcif_delresname'),
               ('mmcif_selhetatm', 'mmcif_delhetatm'), ('mmcif_selelem', 'mmcif_delelem')]
# Tools that keep every atom
KEEP_ALL = ('mmcif_chain', 'mmcif_rplchain', 'mmcif_b', 'mmcif_occ', 'mmcif_reres', 'mmcif_shiftres',
            'mmcif_keepcoord', 'mmcif_element -all', 'mmcif_rplresname', 'mmcif_fixinsert', 'mmcif_reatom',
            'mmcif_uniqname', 'mmcif_sort', 'mmcif_label2auth', 'mmcif_label2auth -chain', 'mmcif_delcat')
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
        ('mmcif_selchain', mmcif_selchain, ([chains[0]],)),
        ('mmcif_delchain', mmcif_delchain, ([chains[0]],)),
        ('mmcif_selres range', mmcif_selres, (set(), [(first, first + 10, 1)])),
        ('mmcif_selres step', mmcif_selres, ({first}, [(-INF, INF, 3)])),
        ('mmcif_delres range', mmcif_delres, ((first, first + 10, 1),)),
        ('mmcif_delres step', mmcif_delres, ((-INF, INF, 4),)),
        ('mmcif_selresname', mmcif_selresname, (['HOH'],)),
        ('mmcif_delresname', mmcif_delresname, (['HOH'],)),
        ('mmcif_selatom', mmcif_selatom, (['CA', 'P'],)),
        ('mmcif_selelem', mmcif_selelem, (['C'],)),
        ('mmcif_selhetatm', mmcif_selhetatm, ()),
        ('mmcif_delhetatm', mmcif_delhetatm, ()),
        ('mmcif_selaltloc', mmcif_selaltloc, (None,)),
        ('mmcif_selaltloc A', mmcif_selaltloc, ('A',)),
        ('mmcif_selmodel', mmcif_selmodel, ([1],)),
        ('mmcif_chain', mmcif_chain, ('Z',)),
        ('mmcif_rplchain', mmcif_rplchain, ((chains[0], 'ZZ'),)),
        ('mmcif_b', mmcif_b, (1.5,)),
        ('mmcif_occ', mmcif_occ, (0.5,)),
        ('mmcif_reres', mmcif_reres, (1,)),
        ('mmcif_shiftres', mmcif_shiftres, (-1000,)),
        ('mmcif_delelem', mmcif_delelem, (['C'],)),
        ('mmcif_keepcoord', mmcif_keepcoord, ()),
        ('mmcif_head', mmcif_head, (100,)),
        ('mmcif_element -all', mmcif_element, (True,)),
        ('mmcif_rplresname', mmcif_rplresname, (('HOH', 'WAT'),)),
        ('mmcif_fixinsert', mmcif_fixinsert, (set(),)),
        ('mmcif_reatom', mmcif_reatom, (1000,)),
        ('mmcif_uniqname', mmcif_uniqname, ()),
        ('mmcif_sort', mmcif_sort, ('CR',)),
        ('mmcif_label2auth', mmcif_label2auth, (True, True)),
        ('mmcif_label2auth -chain', mmcif_label2auth, (True, False)),
        ('mmcif_delcat', mmcif_delcat, (['_pdbx_validate*', '_struct_conf'],)),
        ('mmcif_selentity', mmcif_selentity, (set(), {'polymer', 'branched'})),
        ('mmcif_delentity', mmcif_delentity, (set(), {'polymer', 'branched'})),
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


@unittest.skipUnless(ENTRIES, SKIP)
class TestRealTools(unittest.TestCase):
    """Every tool on every real entry, checked for consistency."""

    @classmethod
    def setUpClass(cls):
        cls.texts = {name(p): read_path(p) for p in ENTRIES}
        cls.timings = []

    @classmethod
    def tearDownClass(cls):
        if os.environ.get('MMCIFTOOLS_TIMINGS'):
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
            for sel, dele in COMPLEMENTS:
                with self.subTest(entry=entry, pair=sel):
                    self.assertEqual(len(rows(outputs[sel])) + len(rows(outputs[dele])), len(atoms))

    def test_reports(self):
        for entry, text in self.texts.items():
            with self.subTest(entry=entry):
                summary = run_tool(mmcif_wc, text, '')
                self.assertTrue(summary.startswith('No. models:\t'))
                fasta = run_tool(mmcif_tofasta, text, False, False)
                self.assertTrue(fasta.startswith('>' + entry.upper() + '|'))
                seqres = run_tool(mmcif_tofasta, text, True, True)
                self.assertTrue(seqres.startswith('>' + entry.upper() + '|'))

    def test_element_matches_deposited(self):
        for entry, text in self.texts.items():
            with self.subTest(entry=entry):
                self.assertEqual(run_tool(mmcif_element, text, True), text)

    def test_gap(self):
        for entry, text in self.texts.items():
            with self.subTest(entry=entry):
                report = run_tool(mmcif_gap, text)
                self.assertTrue(report.endswith(' gap(s) in the structure\n'))
        self.assertTrue(run_tool(mmcif_gap, self.texts['6vxx']).startswith('A:'))

    def test_pdb_round_trip(self):
        auth = ['auth_asym_id', 'auth_seq_id', 'auth_comp_id', 'auth_atom_id', 'pdbx_PDB_ins_code',
                'label_alt_id', 'pdbx_PDB_model_num', 'type_symbol']
        label = ['label_asym_id', 'label_entity_id', 'label_seq_id']
        for entry, text in self.texts.items():
            with self.subTest(entry=entry):
                try:
                    pdb = run_tool(mmcif_topdb, text)
                except _cif.CIFError:
                    self.assertEqual(entry, '4v6x')  # multi-character chains
                    continue
                back = ''.join(mmcif_frompdb.run(io.StringIO(pdb), entry))
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
            paths = mmcif_splitmodel.run(io.StringIO(text), '1d3z', outdir=tmp)
            self.assertEqual(len(paths), 10)
            out = io.StringIO()
            _cif.write(mmcif_mkensemble.run(paths), out)
            self.assertEqual(len(rows(out.getvalue())), len(rows(text)))
            status = mmcif_chkensemble.run(io.StringIO(out.getvalue()))[2]
            self.assertEqual(status, 0)
        finally:
            shutil.rmtree(tmp)

    def test_label_and_author_copies(self):
        for entry, text in self.texts.items():
            before, atoms = reference_problems(text), rows(text)
            with self.subTest(entry=entry, tool='mmcif_label2auth'):
                after = rows(run_tool(mmcif_label2auth, text, True, True))
                self.assertEqual([a['auth_asym_id'] for a in after], [a['label_asym_id'] for a in atoms])
                self.assertEqual([a['auth_seq_id'] for a in after],
                                 [a['label_seq_id'] if a['label_seq_id'] != '.' else a['auth_seq_id']
                                  for a in atoms])
            with self.subTest(entry=entry, tool='mmcif_label2auth -res'):
                try:
                    out = run_tool(mmcif_label2auth, text, False, True)
                except _cif.CIFError as err:
                    self.assertIn('would both become', str(err))
                else:
                    self.assertEqual(len(rows(out)), len(atoms))
                    self.assertLessEqual(reference_problems(out), before)
            with self.subTest(entry=entry, tool='mmcif_auth2label'):
                try:
                    out = run_tool(mmcif_auth2label, text)
                except _cif.CIFError as err:
                    self.assertIn('ligand or water', str(err))
                    continue
                after = rows(out)
                self.assertLessEqual(reference_problems(out), before)
                # one new label chain per old one: no molecules merged
                pairs = {(a['label_asym_id'], b['label_asym_id']) for a, b in zip(atoms, after)}
                self.assertEqual(len({p[0] for p in pairs}), len(pairs))
                self.assertEqual(len({p[1] for p in pairs}), len(pairs))
                polymer = [b for a, b in zip(atoms, after) if a['label_seq_id'] != '.']
                self.assertEqual([b['label_asym_id'] for b in polymer], [b['auth_asym_id'] for b in polymer])

    def test_category_tools(self):
        for entry, text in self.texts.items():
            with self.subTest(entry=entry):
                table = run_tool(mmcif_cat, text, '_atom_site').splitlines()
                self.assertEqual(len(table), len(rows(text)) + 1)
                self.assertEqual(run_tool(mmcif_blocks, text), entry.upper() + '\n')
                self.assertEqual(run_tool(mmcif_keepcat, text, ['_atom_site', '_struct_conn']),
                                 run_tool(mmcif_keepcoord, text))
                out = run_tool(mmcif_delcat, text, ['_pdbx_validate*'])
                self.assertFalse([line for line in out.splitlines()
                                  if line.lower().startswith('_pdbx_validate')])
                self.assertEqual(run_tool(mmcif_delcat, out, ['_pdbx_validate*']), out)

    def test_assembly(self):
        for entry, text in self.texts.items():
            with self.subTest(entry=entry):
                atoms = rows(text)
                if entry.startswith('af-'):  # AlphaFold DB files have no assemblies
                    with self.assertRaises(_cif.CIFError):
                        run_tool(mmcif_assembly, text, '1')
                    continue
                out = run_tool(mmcif_assembly, text, '1')
                after = rows(out)
                self.assertEqual(roundtrip(out, _cif.ALL), out)
                self.assertLessEqual(reference_problems(out), reference_problems(text))
                self.assertEqual(rows(out, '_atom_site_anisotrop'), [])
                self.assertEqual([a['id'] for a in after], [str(i) for i in range(1, len(after) + 1)])
                if entry == '1stp':  # one chain, four operations
                    self.assertEqual(len(after), 4 * len(atoms))
                    self.assertEqual(list(dict.fromkeys(a['auth_asym_id'] for a in after)),
                                     ['A', 'A-2', 'A-3', 'A-4'])
                else:  # identity only: the same atoms, chain by chain
                    key = lambda a: (a['pdbx_PDB_model_num'], a['label_asym_id'], a['label_atom_id'],  # noqa
                                     a['label_seq_id'], a['label_alt_id'], a['Cartn_x'], a['Cartn_z'])
                    self.assertEqual(sorted(map(key, after)), sorted(map(key, atoms)))

    def test_info(self):
        for entry, text in self.texts.items():
            with self.subTest(entry=entry):
                lines = run_tool(mmcif_info, text).splitlines()
                self.assertEqual(lines[0], 'Entry:\t' + entry.upper())
                items = {line.split(':', 1)[0] for line in lines}
                self.assertIn('Entity', items)
                self.assertIn('Model quality' if entry.startswith('af-') else 'Method', items)

    def test_validate(self):
        for entry, text in self.texts.items():
            with self.subTest(entry=entry):
                lines, errors = mmcif_validate.run(io.StringIO(text, newline=''))
                self.assertEqual(errors, 0, ''.join(lines))

    def test_tidy(self):
        for entry, text in self.texts.items():
            with self.subTest(entry=entry):
                out = run_tool(mmcif_tidy, text)
                self.assertEqual(our_tables(out), our_tables(text))  # same data
                self.assertEqual(run_tool(mmcif_tidy, out), out)       # already tidy


if __name__ == '__main__':
    unittest.main()
