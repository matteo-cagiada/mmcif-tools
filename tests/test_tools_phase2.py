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

"""Unit tests for the phase-2 tools, on the fixtures in tests/data."""

import io
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

from support import (  # noqa: E402
    DATA, ROOT, _cif, our_tables, read, reference_problems, rows, run_tool)

from ciftools import (  # noqa: E402
    _atoms, cif_chkensemble, cif_delelem, cif_delinsertion, cif_element, cif_fixinsert,
    cif_frompdb, cif_gap, cif_head, cif_intersect, cif_keepcoord, cif_merge, cif_mkensemble,
    cif_reatom, cif_rplresname, cif_sort, cif_splitchain, cif_splitmodel, cif_tidy, cif_topdb,
    cif_uniqname, cif_validate)

WWPDB = read('wwpdb_like.cif')
NMR = read('nmr_models.cif')


def ids(text, category='_atom_site'):
    return [r['id'] for r in rows(text, category)]


def column(text, name, category='_atom_site'):
    return [r[name] for r in rows(text, category)]


def chain_text(numbers, chain='A'):
    """A small protein: one CA per residue, 3.8 A apart; numbers like '52A'."""
    head = ('data_c\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n'
            '_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n'
            '_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.pdbx_PDB_ins_code\n'
            '_atom_site.auth_comp_id\n_atom_site.auth_atom_id\n_atom_site.Cartn_x\n'
            '_atom_site.Cartn_y\n_atom_site.Cartn_z\n')
    out = []
    for i, label in enumerate(numbers):
        num, code = (label[:-1], label[-1]) if label[-1].isalpha() else (label, '?')
        x = 3.8 * i if not label.startswith('!') else 3.8 * i + 20
        num = num.lstrip('!')
        out.append('ATOM {} C CA ALA A {} {} {} ALA CA {:.3f} 0.000 0.000\n'.format(
            i + 1, chain, num, code, x))
    return head + ''.join(out)


class TestSingleFileTools(unittest.TestCase):

    def test_delelem(self):
        out = run_tool(cif_delelem, WWPDB, ['o'])
        self.assertEqual(len(ids(out)), 11)
        self.assertEqual(rows(out, '_struct_conn'), [])  # SER OG is gone
        self.assertEqual(ids(out, '_atom_site_anisotrop'), ['1', '2', '16'])

    def test_keepcoord(self):
        out = run_tool(cif_keepcoord, WWPDB)
        cats = {k[1] for k in our_tables(out)}
        self.assertEqual(cats, {'_atom_site', '_struct_conn'})
        self.assertEqual(ids(out), ids(WWPDB))
        self.assertTrue(out.startswith('data_1XYZ\n'))

    def test_head(self):
        out = run_tool(cif_head, WWPDB, 3)
        self.assertEqual(ids(out), ['1', '2', '3'])
        self.assertEqual(ids(out, '_atom_site_anisotrop'), ['1', '2'])
        self.assertEqual(rows(out, '_struct_conn'), [])
        self.assertIn('_cell.length_a', out)

    def test_element_fills_missing(self):
        unset = WWPDB.replace(' 1  N N ', ' 1  ? N ').replace(' 16 S S ', ' 16 ? S ')
        self.assertNotEqual(unset, WWPDB)
        self.assertEqual(run_tool(cif_element, unset, False), WWPDB)

    def test_element_adds_the_column(self):
        text = read('multi_block.cif')
        no_element = text.replace('_atom_site.type_symbol\n', '').replace(' 1 N N ', ' 1 N ') \
            .replace(' 2 C CA ', ' 2 CA ').replace(' 1 N N ALA', ' 1 N ALA')
        out = run_tool(cif_element, no_element, False)
        self.assertEqual(column(out, 'type_symbol'), column(text, 'type_symbol'))

    def test_element_rules(self):
        cases = [('CA', 'ALA', False, 'C'), ('CA', 'CA', True, 'CA'), ('NA', 'HEM', True, 'N'),
                 ('FE', 'HEM', True, 'FE'), ('CL1', 'LIG', True, 'CL'), ('HG21', 'THR', False, 'H'),
                 ('1HB', 'ALA', False, 'H'), ("O5'", 'G', False, 'O'), ('SE', 'MSE', True, 'SE'),
                 ('ZN', 'ZN', True, 'ZN'), ('O', 'HOH', True, 'O'), ('1', 'X', True, None)]
        for name, resname, het, expected in cases:
            with self.subTest(name=name, resname=resname):
                self.assertEqual(_atoms.guess_element(name, resname, het), expected)

    def test_rplresname(self):
        out = run_tool(cif_rplresname, WWPDB, ('SER', 'SEP'))
        atoms = rows(out)
        self.assertEqual({a['auth_comp_id'] for a in atoms if a['label_seq_id'] == '2'}, {'SEP'})
        self.assertEqual({a['label_comp_id'] for a in atoms if a['label_seq_id'] == '2'}, {'SEP'})
        self.assertEqual(rows(out, '_struct_conn')[0]['ptnr1_label_comp_id'], 'SEP')
        self.assertEqual(run_tool(cif_rplresname, WWPDB, ('ser', 'X')), WWPDB)  # case-sensitive

    def test_fixinsert(self):
        out = run_tool(cif_fixinsert, chain_text(['51', '52A', '52B', '53']), set())
        self.assertEqual(column(out, 'auth_seq_id'), ['51', '52', '53', '54'])
        self.assertEqual(column(out, 'pdbx_PDB_ins_code'), ['?'] * 4)
        out = run_tool(cif_fixinsert, chain_text(['52', '52A', '53']), set())
        self.assertEqual(column(out, 'auth_seq_id'), ['52', '53', '54'])
        out = run_tool(cif_fixinsert, chain_text(['9', '9A', '10', '10A', '11']), {('A', 9)})
        self.assertEqual(list(zip(column(out, 'auth_seq_id'), column(out, 'pdbx_PDB_ins_code'))),
                         [('9', '?'), ('10', '?'), ('11', '?'), ('12', 'A'), ('13', '?')])

    def test_fixinsert_restarts_at_ligands_and_models(self):
        out = run_tool(cif_fixinsert, WWPDB, set())
        numbers = sorted(set(column(out, 'auth_seq_id')), key=int)
        self.assertEqual(numbers, ['10', '11', '12', '201', '301', '302'])
        self.assertEqual(run_tool(cif_delinsertion, WWPDB, set()), out)
        self.assertEqual(cif_fixinsert.parse_residues('A9,AB-3'), {('A', 9), ('AB', -3)})
        with self.assertRaises(ValueError):
            cif_fixinsert.parse_residues('9A')

    def test_reatom(self):
        self.assertEqual(ids(run_tool(cif_reatom, NMR, 5)), [str(i) for i in range(5, 11)])
        out = run_tool(cif_reatom, WWPDB, 100)
        self.assertEqual(ids(out, '_atom_site_anisotrop'), ['100', '101', '115'])
        self.assertEqual(reference_problems(out), set())

    def test_uniqname(self):
        out = run_tool(cif_uniqname, WWPDB)
        het = [(a['auth_atom_id'], a['label_atom_id']) for a in rows(out) if a['group_PDB'] == 'HETATM']
        self.assertEqual(het, [('S1', 'S1'), ('O1', 'O1'), ('O1', 'O1'), ('O1', 'O1')])
        self.assertEqual(rows(out, '_struct_conn')[0]['ptnr2_label_atom_id'], 'O1')
        self.assertEqual(column(out, 'auth_atom_id')[:15], column(WWPDB, 'auth_atom_id')[:15])
        with self.assertRaises(_cif.CIFError):
            run_tool(cif_uniqname, WWPDB.replace(' 16 S S ', ' 16 ? S '))

    def test_gap(self):
        out = run_tool(cif_gap, chain_text(['1', '2', '3', '5', '5A', '6', '!7']))
        self.assertEqual(out, 'A:ALA3 < Seq. Gap > A:ALA5\n'
                              'A:ALA6 <   23.80A > A:ALA7\n'
                              'Found 2 gap(s) in the structure\n')
        self.assertEqual(run_tool(cif_gap, WWPDB), 'Found 0 gap(s) in the structure\n')

    def test_sort(self):
        lines = WWPDB.splitlines(True)
        atom = [i for i, line in enumerate(lines) if line.startswith(('ATOM', 'HETATM'))]
        shuffled = lines[:atom[0]] + lines[atom[0]:atom[-1] + 1][::-1] + lines[atom[-1] + 1:]
        self.assertEqual(run_tool(cif_sort, ''.join(shuffled), 'CR'), WWPDB)
        aniso = [i for i, line in enumerate(lines) if line.startswith(('1  N', '2  C', '16 S'))]
        swapped = list(lines)
        swapped[aniso[0]], swapped[aniso[2]] = swapped[aniso[2]], swapped[aniso[0]]
        self.assertEqual(ids(run_tool(cif_sort, ''.join(swapped), 'CR'), '_atom_site_anisotrop'),
                         ['1', '2', '16'])

    def test_sort_keys(self):
        text = chain_text(['3', '1', '2'])
        self.assertEqual(column(run_tool(cif_sort, text, 'CR'), 'auth_seq_id'), ['1', '2', '3'])
        self.assertEqual(column(run_tool(cif_sort, text, 'C'), 'auth_seq_id'), ['3', '1', '2'])
        two = NMR.replace('ATOM 1 N', 'ATOM 9 N')
        self.assertEqual(column(run_tool(cif_sort, two, 'CR'), 'pdbx_PDB_model_num'),
                         ['1'] * 3 + ['2'] * 3)

    def test_tidy(self):
        for name in ('wwpdb_like.cif', 'text_fields.cif', 'ccd_like.cif', 'multi_block.cif'):
            with self.subTest(fixture=name):
                text = read(name)
                out = run_tool(cif_tidy, text)
                self.assertEqual(our_tables(out), our_tables(text))
                self.assertEqual(run_tool(cif_tidy, out), out)
        self.assertEqual(run_tool(cif_tidy, WWPDB.replace('\n', '\r\n')), run_tool(cif_tidy, WWPDB))

    def test_validate(self):
        lines, errors = cif_validate.run(io.StringIO(WWPDB))
        self.assertEqual(errors, 0)
        broken = WWPDB.replace('ATOM   2  C CA  . MET A 1 1 ? 22.767', 'ATOM   1  C CA  . MET A 1 1 ? x.767')
        lines, errors = cif_validate.run(io.StringIO(broken))
        self.assertEqual(errors, 3)  # duplicate id, coordinate, anisotrop row 2 without its atom
        self.assertTrue(any('anisotrop 2' in line for line in lines))
        lines, errors = cif_validate.run(io.StringIO("data_x\nloop_\n_a.b\n'1\n"))
        self.assertEqual(errors, 1)
        self.assertIn('syntax', lines[0])


class TestConversion(unittest.TestCase):

    def test_topdb(self):
        out = run_tool(cif_topdb, WWPDB)
        lines = out.splitlines()
        self.assertEqual(lines[0], 'ATOM      1  N   MET A  10      21.598   3.026  14.624  1.00 29.41           N  ')
        self.assertEqual({len(l) for l in lines if l.startswith(('ATOM', 'HETATM', 'END'))}, {80})
        self.assertEqual(lines[11][17:27], 'GLY A  11A')
        self.assertTrue(lines[-1].startswith('END'))
        nmr = run_tool(cif_topdb, NMR).splitlines()
        self.assertEqual([l for l in nmr if l.startswith(('MODEL', 'ENDMDL'))],
                         ['MODEL        1', 'ENDMDL', 'MODEL        2', 'ENDMDL'])
        with self.assertRaises(_cif.CIFError):
            run_tool(cif_topdb, read('multichar_chains.cif'))

    def test_round_trip(self):
        for name in ('wwpdb_like.cif', 'nmr_models.cif', 'nucleic.cif'):
            with self.subTest(fixture=name):
                text = read(name)
                back = run_tool_lines(cif_frompdb, run_tool(cif_topdb, text))
                cols = ['auth_asym_id', 'auth_seq_id', 'auth_comp_id', 'auth_atom_id',
                        'pdbx_PDB_model_num', 'label_asym_id', 'label_seq_id']
                for col in cols:
                    self.assertEqual(column(back, col), column(text, col), col)
                self.assertEqual(cif_validate.run(io.StringIO(back))[1], 0)

    def test_frompdb_labels_and_cell(self):
        pdb = ('HEADER    TEST' + ' ' * 48 + '9ZZZ\n'
               'CRYST1   50.000   60.000   70.000  90.00  90.00  90.00 P 21 21 21    4\n'
               'ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N\n'
               'HETATM    2  N   MSE A   2       2.000   2.000   3.000  1.00 10.00           N\n'
               'ATOM      3  N   GLY A   3       3.000   2.000   3.000  1.00 10.00           N\n'
               'TER       4      GLY A   3\n'
               'HETATM    5 ZN    ZN A 101       9.000   9.000   9.000  1.00 10.00          ZN2+\n'
               'HETATM    6  O   HOH A 201       8.000   9.000   9.000  1.00 10.00\n'
               'ATOM      7  N   ALA B   1       1.000   5.000   3.000  1.00 10.00           N\n'
               'END\n')
        out = ''.join(cif_frompdb.run(io.StringIO(pdb), 'x'))
        atoms = rows(out)
        self.assertEqual([(a['label_asym_id'], a['label_entity_id'], a['label_seq_id']) for a in atoms],
                         [('A', '1', '1'), ('A', '1', '2'), ('A', '1', '3'), ('B', '2', '.'),
                          ('C', '3', '.'), ('D', '4', '1')])
        self.assertEqual([a['type_symbol'] for a in atoms], ['N', 'N', 'N', 'ZN', 'O', 'N'])
        self.assertEqual(atoms[3]['pdbx_formal_charge'], '2')
        self.assertIn('data_9ZZZ\n', out)
        self.assertEqual(rows(out, '_cell')[0]['length_b'], '60.000')
        self.assertEqual(rows(out, '_symmetry')[0]['space_group_name_H-M'], 'P 21 21 21')


def run_tool_lines(module, text):
    return ''.join(module.run(io.StringIO(text, newline='')))


class TestFileTools(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp)

    def path(self, name, text):
        path = os.path.join(self.tmp, name)
        with open(path, 'w', newline='') as fh:
            fh.write(text)
        return path

    def test_splitchain(self):
        text = read('multichar_chains.cif')
        paths = cif_splitchain.run(io.StringIO(text), 'big', outdir=self.tmp)
        self.assertEqual([os.path.basename(p) for p in paths], ['big_Axp.cif', 'big_A1.cif', 'big_zz.cif'])
        with open(paths[0]) as fh:
            first = fh.read()
        self.assertTrue(first.startswith('data_4BIG_Axp\n'))
        self.assertEqual(ids(first), ['1', '2'])
        fhandle, base, _ = cif_splitchain.check_input(['--auth', os.path.join(DATA, 'nucleic.cif')])
        fhandle.close()
        self.assertEqual(base, 'nucleic')

    def test_splitmodel(self):
        paths = cif_splitmodel.run(io.StringIO(NMR), 'nmr', outdir=self.tmp)
        self.assertEqual([os.path.basename(p) for p in paths], ['nmr_1.cif', 'nmr_2.cif'])
        with open(paths[1]) as fh:
            self.assertEqual(ids(fh.read()), ['4', '5', '6'])
        single = cif_splitmodel.run(io.StringIO(WWPDB), 'one', outdir=self.tmp)
        self.assertEqual([os.path.basename(p) for p in single], ['one_1.cif'])

    def test_mkensemble_and_chkensemble(self):
        one = self.path('one.cif', WWPDB)
        two = self.path('two.cif', WWPDB)
        out = ''.join(_written(cif_mkensemble.run([one, two])))
        self.assertTrue(out.startswith('# MODEL 1 FROM one.cif\n# MODEL 2 FROM two.cif\n'))
        self.assertEqual(ids(out), [str(i) for i in range(1, 39)])
        self.assertEqual(column(out, 'pdbx_PDB_model_num'), ['1'] * 19 + ['2'] * 19)
        self.assertEqual(len(rows(out, '_struct_conn')), 1)
        self.assertNotIn('_cell', out)
        out_text, err, status = cif_chkensemble.run(io.StringIO(out))
        self.assertEqual((out_text, status), (['Ensemble of 2 models *seems* OK\n'], 0))
        mixed = ''.join(_written(cif_mkensemble.run([one, self.path('nmr.cif', NMR)])))
        out_text, err, status = cif_chkensemble.run(io.StringIO(mixed))
        self.assertEqual(status, 1)
        self.assertEqual(err[0], 'Models 1 and 2 differ:\n')
        self.assertIn('A:ALA1:N\n', err)

    def test_mkensemble_union_of_columns(self):
        out = ''.join(_written(cif_mkensemble.run([
            self.path('a.cif', read('af3_like.cif')), self.path('b.cif', NMR)])))
        atoms = rows(out)
        self.assertEqual(len(atoms), 9)
        self.assertEqual(atoms[0]['auth_comp_id'], '?')  # af3_like has no auth_comp_id
        self.assertEqual(atoms[3]['auth_comp_id'], 'ALA')

    def test_merge(self):
        warnings = []
        out = ''.join(_written(cif_merge.run([self.path('a.cif', WWPDB), self.path('b.cif', NMR)],
                                             warnings=warnings)))
        self.assertEqual(ids(out), [str(i) for i in range(1, 26)])
        self.assertEqual(ids(out, '_atom_site_anisotrop'), ['1', '2', '16'])
        self.assertIn('_cell.length_a', out)
        self.assertEqual(len(warnings), 1)
        self.assertIn('chain A', warnings[0])
        self.assertEqual(cif_validate.run(io.StringIO(out))[1], 0)

    def test_intersect(self):
        part = run_tool(cif_head, WWPDB, 5)
        out = ''.join(_written(cif_intersect.run([self.path('a.cif', WWPDB), self.path('b.cif', part)])))
        self.assertEqual(ids(out), ['1', '2', '3', '4', '5'])
        self.assertIn('_cell.length_a', out)

    def test_command_lines(self):
        nmr = os.path.join(DATA, 'nmr_models.cif')
        cases = [
            ('cif_merge', [], None),  # usage
            ('cif_merge', ['missing.cif'], b'File not found'),
            ('cif_head', ['-0', nmr], b'positive numerical value'),
            ('cif_fixinsert', ['--label', nmr], b'author identifiers'),
            ('cif_sort', ['-X', nmr], b'Sorting key is not valid'),
            ('cif_element', ['-some', nmr], b'invalid option'),
            ('cif_rplresname', ['-ALA', nmr], b'Invalid option value'),
        ]
        for name, args, message in cases:
            with self.subTest(tool=name, args=args):
                proc = subprocess.run([sys.executable, '-m', 'ciftools.' + name] + args,
                                      stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE, cwd=ROOT)
                self.assertEqual(proc.returncode, 1)
                if message:
                    self.assertIn(message, proc.stderr)
        proc = subprocess.run([sys.executable, '-m', 'ciftools.cif_splitchain',
                               os.path.join(DATA, 'multichar_chains.cif')],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.tmp,
                              env=dict(os.environ, PYTHONPATH=ROOT))
        self.assertEqual(proc.returncode, 0)
        self.assertEqual(sorted(os.listdir(self.tmp)),
                         ['multichar_chains_A1.cif', 'multichar_chains_Axp.cif', 'multichar_chains_zz.cif'])
        proc = subprocess.run([sys.executable, '-m', 'ciftools.cif_validate', nmr],
                              stdout=subprocess.PIPE, cwd=ROOT)
        self.assertEqual(proc.returncode, 0)


def _written(events):
    """Text of an event stream."""
    out = io.StringIO()
    _cif.write(events, out)
    return [out.getvalue()]


if __name__ == '__main__':
    unittest.main()
