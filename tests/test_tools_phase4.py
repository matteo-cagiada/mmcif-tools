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

import os
import re
import subprocess
import sys
import unittest

from support import (  # noqa: E402
    DATA, ROOT, _cif, read, reference_problems, roundtrip, rows, run_tool)

from mmciftools import (  # noqa: E402
    mmcif_assembly, mmcif_delcat, mmcif_delentity, mmcif_delhetatm, mmcif_info, mmcif_selentity)

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


if __name__ == '__main__':
    unittest.main()
