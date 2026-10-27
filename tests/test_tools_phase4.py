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

"""Unit tests for the phase-4 tools (entities, assemblies, metadata)."""

import io
import os
import re
import subprocess
import sys
import unittest

from support import (  # noqa: E402
    DATA, FIXTURES, ROOT, _cif, read, read_path, reference_problems, roundtrip, rows, run_tool,
    split_by_filtering, split_files)

from mmciftools import (  # noqa: E402
    mmcif_assembly, mmcif_delcat, mmcif_delentity, mmcif_delhetatm, mmcif_delinsertion,
    mmcif_b, mmcif_chain, mmcif_fixinsert, mmcif_head, mmcif_info, mmcif_label2auth, mmcif_occ,
    mmcif_reres, mmcif_selaltloc, mmcif_selelem, mmcif_selentity, mmcif_splitchain, mmcif_tofasta)

WWPDB = read('wwpdb_like.cif')

# Two operations: the identity and a 2-fold axis along z through (5, 0, z)
ASSEMBLY = """data_asm
loop_
_pdbx_struct_assembly.id
_pdbx_struct_assembly.details
_pdbx_struct_assembly.oligomeric_count
1 'dimer' 2
2 'product' 4
3 'no identity' 1
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 1,2         A,B
2 '(1,2)(1,2)' A
3 2           A
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.type
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[3]
1 'identity operation'         1  0 0 0  0  1 0 0 0 0 1 0
2 'crystal symmetry operation' -1 0 0 10 0 -1 0 0 0 0 1 0
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr1_symmetry
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.ptnr2_symmetry
metalc1 metalc A HIS 1 NE2 A 1   1_555 B ZN . ZN A 101 1_555
metalc2 metalc A HIS 1 NE2 A 1   1_555 B ZN . ZN A 101 2_655
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM   1 N  NE2 HIS A 1 1 1.000 2.000 3.000 1   A 1
HETATM 2 ZN ZN  ZN  B 2 . 3.000 2.000 3.000 101 A 1
loop_
_atom_site_anisotrop.id
_atom_site_anisotrop.U[1][1]
1 0.1
"""


def drop_columns(text, names):
    """text without some _atom_site columns (as in integrative models)."""
    lines = text.split('\n')
    start = next(i for i, line in enumerate(lines) if line.startswith('_atom_site.'))
    end = start
    while lines[end].startswith('_atom_site.'):
        end += 1
    tags = [line.strip() for line in lines[start:end]]
    drop = {tags.index('_atom_site.' + n) for n in names}
    body = [' '.join(t for k, t in enumerate(line.split()) if k not in drop)
            if line.startswith(('ATOM', 'HETATM')) else line for line in lines[end:]]
    return '\n'.join(lines[:start] + [t for k, t in enumerate(tags) if k not in drop] + body)


def column_of(text, name, category='_atom_site'):
    return [r[name] for r in rows(text, category)]


def atoms(text):
    return [(a['label_asym_id'], a['auth_asym_id'], a['id'],
             tuple(round(float(a[c]), 3) for c in ('Cartn_x', 'Cartn_y', 'Cartn_z')))
            for a in rows(text)]


def cli(name, *args):
    return subprocess.run([sys.executable, '-m', 'mmciftools.' + name] + list(args),
                          stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, cwd=ROOT)


class TestEntities(unittest.TestCase):

    def count(self, module, *args, text=WWPDB):
        return len(rows(run_tool(module, text, *args)))

    def test_by_type_and_number(self):
        self.assertEqual(self.count(mmcif_selentity, set(), {'polymer'}), 15)
        self.assertEqual(self.count(mmcif_selentity, set(), {'water'}), 2)
        self.assertEqual(self.count(mmcif_selentity, {'2'}, {'water'}), 4)
        self.assertEqual(self.count(mmcif_delentity, set(), {'non-polymer', 'water'}), 15)
        for ids, types in ((set(), {'polymer'}), ({'2'}, set()), ({'3'}, {'non-polymer'})):
            with self.subTest(ids=ids, types=types):
                self.assertEqual(self.count(mmcif_selentity, ids, types)
                                 + self.count(mmcif_delentity, ids, types), len(rows(WWPDB)))

    def test_modified_residues_stay(self):
        text = re.sub(r'^ATOM(   [1-4] )', r'HETATM\1', WWPDB, flags=re.M)  # MET as HETATM
        self.assertEqual(self.count(mmcif_delentity, set(), {'non-polymer', 'water'}, text=text), 15)
        self.assertEqual(self.count(mmcif_delhetatm, text=text), 11)

    def test_without_entity_category(self):
        text = run_tool(mmcif_delcat, WWPDB, ['_entity'])
        self.assertEqual(self.count(mmcif_selentity, set(), {'water'}, text=text), 2)
        self.assertEqual(self.count(mmcif_selentity, set(), {'polymer'}, text=text), 15)
        self.assertEqual(self.count(mmcif_selentity, set(), {'non-polymer'}, text=text), 2)

    def test_links_follow(self):
        out = run_tool(mmcif_delentity, WWPDB, set(), {'non-polymer'})
        self.assertEqual(rows(out, '_struct_conn'), [])  # its partner SO4 is gone
        self.assertEqual(len(rows(out, '_atom_site_anisotrop')), 2)
        self.assertEqual(reference_problems(out), set())

    def test_command_line(self):
        path = os.path.join(DATA, 'wwpdb_like.cif')
        self.assertEqual(cli('mmcif_selentity', '-POLYMER,3', path).returncode, 0)
        proc = cli('mmcif_delentity', '-,', path)
        self.assertEqual(proc.returncode, 1)
        self.assertIn(b'at least one entity', proc.stderr)
        no_entity = WWPDB.replace('_atom_site.label_entity_id \n', '_atom_site.label_entity_idx \n')
        with self.assertRaises(_cif.CIFError):
            run_tool(mmcif_selentity, no_entity, {'1'}, set())


class TestAssembly(unittest.TestCase):

    def test_expressions(self):
        cases = [('1', [['1']]), ('1,2,5', [['1', '2', '5']]), ('(1-3)', [['1', '2', '3']]),
                 ('(1-2)(3,4)', [['1', '2'], ['3', '4']]), ('(X0)(1-2)', [['X0'], ['1', '2']]),
                 (' ( 1 - 2 ) ', [['1', '2']]), ('P', [['P']])]
        for text, groups in cases:
            with self.subTest(text=text):
                self.assertEqual(mmcif_assembly.parse_expression(text), groups)
        for bad in ('(1-2', '(a-b)', '', '()', '(1)x(2)'):
            with self.subTest(text=bad):
                with self.assertRaises(_cif.CIFError):
                    mmcif_assembly.parse_expression(bad)

    def test_dimer(self):
        out = run_tool(mmcif_assembly, ASSEMBLY, '1')
        self.assertEqual(atoms(out), [('A', 'A', '1', (1.0, 2.0, 3.0)), ('A-2', 'A-2', '2', (9.0, -2.0, 3.0)),
                                      ('B', 'A', '3', (3.0, 2.0, 3.0)), ('B-2', 'A-2', '4', (7.0, -2.0, 3.0))])
        links = rows(out, '_struct_conn')
        self.assertEqual([(r['id'], r['ptnr1_label_asym_id'], r['ptnr1_auth_asym_id'],
                           r['ptnr2_label_asym_id'], r['ptnr2_auth_asym_id']) for r in links],
                         [('metalc1', 'A', 'A', 'B', 'A'), ('metalc1-2', 'A-2', 'A-2', 'B-2', 'A-2')])
        self.assertEqual(rows(out, '_atom_site_anisotrop'), [])
        self.assertEqual(reference_problems(out), set())
        self.assertEqual(roundtrip(out, _cif.ALL), out)  # re-reads

    def test_product_and_naming(self):
        out = run_tool(mmcif_assembly, ASSEMBLY, '2')
        self.assertEqual([a[:2] + a[3:] for a in atoms(out)],
                         [('A', 'A', (1.0, 2.0, 3.0)), ('A-1x2', 'A-1x2', (9.0, -2.0, 3.0)),
                          ('A-2x1', 'A-2x1', (9.0, -2.0, 3.0)), ('A-2x2', 'A-2x2', (1.0, 2.0, 3.0))])
        self.assertEqual(rows(out, '_struct_conn'), [])  # B is not in the assembly
        out = run_tool(mmcif_assembly, ASSEMBLY, '3')
        self.assertEqual(atoms(out), [('A-2', 'A-2', '1', (9.0, -2.0, 3.0))])

    def test_errors(self):
        with self.assertRaises(_cif.CIFError) as ctx:
            run_tool(mmcif_assembly, ASSEMBLY, '9')
        self.assertIn('the file has 1, 2, 3', str(ctx.exception))
        with self.assertRaises(_cif.CIFError):
            run_tool(mmcif_assembly, WWPDB, '1')
        broken = ASSEMBLY.replace('1 1,2         A,B', '1 1,7         A,B')
        with self.assertRaises(_cif.CIFError) as ctx:
            run_tool(mmcif_assembly, broken, '1')
        self.assertIn('No operation 7', str(ctx.exception))

    def test_command_line(self):
        proc = cli('mmcif_assembly', os.path.join(DATA, 'wwpdb_like.cif'))
        self.assertEqual(proc.returncode, 1)
        self.assertIn(b'No assemblies', proc.stderr)


class TestInfo(unittest.TestCase):

    def test_summary(self):
        lines = run_tool(mmcif_info, WWPDB).splitlines()
        self.assertEqual(lines[0], 'Entry:\t1XYZ')
        entities = [line for line in lines if line.startswith('Entity:')]
        self.assertEqual(len(entities), 3)
        self.assertTrue(entities[0].startswith('Entity:\t1\tpolymer\t'), entities[0])
        self.assertIn('label chain A', entities[0])
        assembly = run_tool(mmcif_info, ASSEMBLY).splitlines()
        self.assertEqual([line for line in assembly if line.startswith('Assembly:')],
                         ['Assembly:\t1\t2 chains\tdimer', 'Assembly:\t2\t4 chains\tproduct',
                          'Assembly:\t3\t1 chain\tno identity'])

    def test_blocks_and_models(self):
        self.assertEqual(run_tool(mmcif_info, read('multi_block.cif')), 'Entry:\tfirst\n\nEntry:\tsecond\n')
        self.assertIn('Model quality:\t', run_tool(mmcif_info, read('af3_like.cif')))

    def test_command_line(self):
        proc = cli('mmcif_info', os.path.join(DATA, 'wwpdb_like.cif'))
        self.assertEqual(proc.returncode, 0)
        self.assertTrue(proc.stdout.startswith(b'Entry:\t1XYZ\n'))


class TestStressFixes(unittest.TestCase):
    """Problems found by the 200-system stress test."""

    def test_split_chains_differing_only_in_case(self):
        import io
        import shutil
        import tempfile
        text = WWPDB.replace('? 301 HOH A', '? 301 HOH a').replace('? 302 HOH A', '? 302 HOH a')
        tmp = tempfile.mkdtemp()
        try:
            notes = []
            paths = mmcif_splitchain.run(io.StringIO(text), '6ozp', outdir=tmp, notes=notes)
            self.assertEqual([os.path.basename(p) for p in paths], ['6ozp_A.cif', '6ozp_a~2.cif'])
            self.assertEqual(sorted(os.listdir(tmp)), ['6ozp_A.cif', '6ozp_a~2.cif'])
            self.assertEqual(sum(len(rows(read_path(p))) for p in paths), len(rows(WWPDB)))
            self.assertEqual(len(notes), 1)
            self.assertEqual(mmcif_splitchain.file_names('x', ['A', 'B', 'a', 'b', 'aa', 'AA', 'Aa']),
                             {'A': 'x_A.cif', 'B': 'x_B.cif', 'a': 'x_a~2.cif', 'b': 'x_b~2.cif',
                              'aa': 'x_aa.cif', 'AA': 'x_AA~2.cif', 'Aa': 'x_Aa~3.cif'})
        finally:
            shutil.rmtree(tmp)

    def test_batched_splitchain_writes_the_same_files(self):  # 7Y7A: 1792 chains
        texts = [read(n) for n in FIXTURES if n != 'empty_atom_site.cif']
        texts.append(WWPDB.replace('? 301 HOH A', '? 301 HOH a'))
        for text in texts:
            for scheme in (None, 'label'):
                self.assertEqual(split_files(text, scheme), split_by_filtering(text, scheme))
        batch = mmcif_splitchain.BATCH
        mmcif_splitchain.BATCH = 2  # several readings of the input
        try:
            for name in ('multichar_chains.cif', 'wwpdb_like.cif'):
                self.assertEqual(split_files(read(name), 'label'), split_by_filtering(read(name), 'label'))
        finally:
            mmcif_splitchain.BATCH = batch

    def test_files_without_author_residue_columns(self):  # 9A1M, an integrative model
        text = drop_columns(WWPDB, ['auth_seq_id', 'auth_asym_id', 'pdbx_PDB_ins_code'])
        self.assertEqual(run_tool(mmcif_fixinsert, text, set()), text)
        self.assertEqual(run_tool(mmcif_delinsertion, text, set()), text)
        out = rows(run_tool(mmcif_label2auth, text, True, True))
        self.assertEqual([(a['auth_asym_id'], a['auth_seq_id']) for a in out],
                         [(a['label_asym_id'], a['label_seq_id'] if a['label_seq_id'] != '.' else '?')
                          for a in out])
        chains_only = rows(run_tool(mmcif_label2auth, text, True, False))
        self.assertNotIn('auth_seq_id', chains_only[0].loop)
        self.assertEqual([a['auth_asym_id'] for a in chains_only], [a['label_asym_id'] for a in chains_only])

    def test_head_counts_across_data_blocks(self):  # round 2: two blocks
        two = WWPDB + WWPDB.replace('data_1XYZ', 'data_1XYZ_copy', 1)
        self.assertEqual(len(rows(run_tool(mmcif_head, two, 5))), 5)
        self.assertEqual(len(rows(run_tool(mmcif_head, two, 25))), 25)
        self.assertEqual(len(rows(run_tool(mmcif_head, WWPDB, 5))), 5)

    def test_fixinsert_links_without_residue_names(self):  # round 2: a file written by gemmi
        text = WWPDB.replace(
            "_struct_conn.pdbx_dist_value \nhydrog1 hydrog A SER 2 OG B SO4 . O1 2.95 \n",
            "_struct_conn.pdbx_dist_value \n_struct_conn.ptnr1_auth_asym_id \n_struct_conn.ptnr1_auth_seq_id \n"
            "_struct_conn.pdbx_ptnr1_PDB_ins_code \n_struct_conn.ptnr2_auth_asym_id \n"
            "_struct_conn.ptnr2_auth_seq_id \n_struct_conn.pdbx_ptnr2_PDB_ins_code \n"
            "hydrog1 hydrog A GLY 3 CA B SO4 . O1 2.95 A 11 A A 201 ? \n")
        out = run_tool(mmcif_fixinsert, text, set())
        link = rows(out, '_struct_conn')[0]
        self.assertEqual((link['ptnr1_auth_seq_id'], link['pdbx_ptnr1_PDB_ins_code']), ('12', '?'))
        self.assertEqual(reference_problems(out), set())

    def test_anisotropic_rows_follow_their_atom(self):  # round 2: 4PTJ, models that differ
        text = ('data_t\nloop_\n' + ''.join('_atom_site.{}\n'.format(c) for c in (
            'group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_comp_id', 'label_asym_id', 'label_seq_id',
            'auth_seq_id', 'auth_asym_id', 'Cartn_x', 'Cartn_y', 'Cartn_z', 'pdbx_PDB_model_num'))
                + 'ATOM 1 N N MET A 1 10 A 0 0 0 1\nATOM 2 N N SER A 2 11 A 0 0 0 1\n'
                'ATOM 3 N N SER A 2 11 A 0 0 0 2\n'
                'loop_\n_atom_site_anisotrop.id\n_atom_site_anisotrop.pdbx_auth_asym_id\n'
                '_atom_site_anisotrop.pdbx_auth_seq_id\n_atom_site_anisotrop.U[1][1]\n'
                '1 A 10 0.1\n2 A 11 0.1\n3 A 11 0.1\n')
        out = run_tool(mmcif_reres, text, 1)
        self.assertEqual([r['auth_seq_id'] for r in rows(out)], ['1', '2', '1'])  # model 2 has no MET
        self.assertEqual(column_of(out, 'pdbx_auth_seq_id', '_atom_site_anisotrop'), ['1', '2', '1'])
        self.assertEqual(reference_problems(out), set())
        coded = text.replace('_atom_site.pdbx_PDB_model_num\n', '_atom_site.pdbx_PDB_model_num\n'
                             '_atom_site.pdbx_PDB_ins_code\n').replace(' 0 0 0 1\n', ' 0 0 0 1 ?\n')
        coded = coded.replace('ATOM 3 N N SER A 2 11 A 0 0 0 2\n', 'ATOM 3 N N SER A 2 10 A 0 0 0 2 A\n')
        coded = coded.replace('3 A 11 0.1', '3 A 10 0.1')
        self.assertEqual(reference_problems(run_tool(mmcif_fixinsert, coded, set())), set())

    def test_edits_need_the_columns_of_their_scheme(self):  # round 2: files without label (or author) ids
        no_auth = drop_columns(WWPDB, ['auth_asym_id', 'auth_seq_id'])
        with self.assertRaises(_cif.CIFError) as ctx:
            run_tool(mmcif_chain, no_auth, 'Z')
        self.assertIn('no auth_asym_id column', str(ctx.exception))
        no_label = drop_columns(WWPDB, ['label_asym_id', 'label_seq_id'])
        with self.assertRaises(_cif.CIFError) as ctx:
            run_tool(mmcif_reres, no_label, 1, scheme='label')
        self.assertIn('no label_seq_id column', str(ctx.exception))

    def test_missing_optional_columns(self):  # round 2: modified copies
        no_alt = drop_columns(WWPDB, ['label_alt_id'])
        self.assertEqual(run_tool(mmcif_selaltloc, no_alt, None), no_alt)
        no_b = drop_columns(WWPDB, ['B_iso_or_equiv', 'occupancy'])
        self.assertEqual({r['B_iso_or_equiv'] for r in rows(run_tool(mmcif_b, no_b, 5.0))}, {'5.00'})
        self.assertEqual({r['occupancy'] for r in rows(run_tool(mmcif_occ, no_b, 0.5))}, {'0.50'})
        with self.assertRaises(_cif.CIFError) as ctx:
            run_tool(mmcif_selelem, drop_columns(WWPDB, ['type_symbol']), ['C'])
        self.assertIn('no type_symbol column', str(ctx.exception))

    def test_upper_case_tags(self):  # round 2: tags are case-insensitive
        import shutil
        import tempfile
        from mmciftools import mmcif_merge
        upper = re.sub(r'^(_[A-Za-z_]+\.\S+)', lambda m: m.group(1).upper(), WWPDB, flags=re.M)
        tmp = tempfile.mkdtemp()
        try:
            path = os.path.join(tmp, 'u.cif')
            with open(path, 'w', newline='') as fh:
                fh.write(upper)
            out = io.StringIO()
            _cif.write(mmcif_merge.run([path, path]), out)
            ids = [r['id'] for r in rows(out.getvalue())]
            self.assertEqual(len(ids), 2 * len(rows(WWPDB)))
            self.assertEqual(len(set(ids)), len(ids))  # renumbered, no duplicates
        finally:
            shutil.rmtree(tmp)
        no_b = re.sub(r'^(_[A-Za-z_]+\.\S+)', lambda m: m.group(1).upper(),
                      drop_columns(WWPDB, ['B_iso_or_equiv']), flags=re.M)
        self.assertEqual({r['B_iso_or_equiv'] for r in rows(run_tool(mmcif_b, no_b, 3.0))}, {'3.00'})

    def test_seqres_without_poly_seq_scheme(self):  # files not from the PDB; 3IRI (sugars only)
        text = ('data_x\nloop_\n_entity_poly.entity_id\n_entity_poly.pdbx_strand_id\n1 A,C\n2 B\n'
                'loop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n_entity_poly_seq.mon_id\n'
                '1 1 MET\n1 2 GLY\n1 2 ALA\n2 1 DA\n2 2 DC\n'
                'loop_\n_struct_asym.id\n_struct_asym.entity_id\nA 1\nB 2\nC 1\nD 3\n')
        self.assertEqual(run_tool(mmcif_tofasta, text, True, True), '>x|A\nMG\n>x|B\nAC\n>x|C\nMG\n')
        self.assertEqual(run_tool(mmcif_tofasta, text, False, True), '>x|ABC\nMGACMG\n')
        with self.assertRaises(_cif.CIFError) as ctx:
            run_tool(mmcif_tofasta, 'data_x\n_entry.id x\n', False, True)
        self.assertIn('no polymer sequence', str(ctx.exception))


if __name__ == '__main__':
    unittest.main()
