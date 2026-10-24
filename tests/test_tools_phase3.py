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

"""Unit tests for the phase-3 tools and fixes, on the fixtures in tests/data."""

import io
import os
import subprocess
import sys
import unittest

from support import (  # noqa: E402
    DATA, ROOT, _cif, read, reference_problems, rows, run_tool)

from ciftools import (  # noqa: E402
    cif_auth2label, cif_blocks, cif_cat, cif_chain, cif_delcat, cif_frompdb, cif_keepcat,
    cif_keepcoord, cif_label2auth, cif_rplchain, cif_selaltloc, cif_topdb, cif_validate)

WWPDB = read('wwpdb_like.cif')
NMR = read('nmr_models.cif')
MULTICHAR = read('multichar_chains.cif')

# wwpdb_like.cif with author columns in _struct_conn (as in PDB entries)
CONN_AUTH = WWPDB.replace(
    "_struct_conn.pdbx_dist_value \nhydrog1 hydrog A SER 2 OG B SO4 . O1 2.95 \n",
    "_struct_conn.pdbx_dist_value \n_struct_conn.ptnr1_auth_asym_id \n"
    "_struct_conn.ptnr1_auth_seq_id \n_struct_conn.pdbx_ptnr1_PDB_ins_code \n"
    "_struct_conn.ptnr2_auth_asym_id \n_struct_conn.ptnr2_auth_seq_id \n"
    "_struct_conn.pdbx_ptnr2_PDB_ins_code \n"
    "hydrog1 hydrog A SER 2 OG B SO4 . O1 2.95 A 11 ? A 201 ? \n")


def column(text, name, category='_atom_site'):
    return [r[name] for r in rows(text, category)]


def values(text, category='_atom_site'):
    return [list(r.values) for r in rows(text, category)]


def pdb_atom(serial, name, resname, chain, resnum, x, group='ATOM', element=None):
    element = element or name[0]
    return '{:<6}{:>5} {:^4} {:>3} {:1}{:>4}    {:>8.3f}{:>8.3f}{:>8.3f}  1.00  0.00          {:>2}\n'.format(
        group, serial, name, resname, chain, resnum, x, 0.0, 0.0, element)


def convert(pdb):
    return ''.join(cif_frompdb.run(io.StringIO(pdb), 'x'))


def cli(name, *args):
    return subprocess.run([sys.executable, '-m', 'ciftools.' + name] + list(args),
                          stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, cwd=ROOT)


class TestFromPdbTer(unittest.TestCase):

    def test_ter_ends_a_molecule_in_a_blank_chain(self):
        pdb = (pdb_atom(1, 'CA', 'ALA', '', 1, 0.0) + pdb_atom(2, 'CA', 'ALA', '', 2, 3.8) + 'TER\n'
               + pdb_atom(3, 'CA', 'GLY', '', 1, 10.0) + pdb_atom(4, 'CA', 'GLY', '', 2, 13.8) + 'TER\n'
               + pdb_atom(5, 'NA', 'NA', '', 1, 20.0, 'HETATM', 'NA') + 'TER\n'
               + pdb_atom(6, 'O', 'HOH', '', 1, 30.0, 'HETATM') + 'END\n')
        atoms = rows(convert(pdb))
        self.assertEqual([(a['label_asym_id'], a['label_entity_id'], a['label_seq_id'], a['auth_asym_id'])
                          for a in atoms],
                         [('A', '1', '1', 'A'), ('A', '1', '2', 'A'), ('B', '2', '1', 'B'),
                          ('B', '2', '2', 'B'), ('C', '3', '.', 'C'), ('D', '4', '.', 'D')])

    def test_blank_chain_avoids_used_names(self):
        pdb = (pdb_atom(1, 'CA', 'ALA', 'A', 1, 0.0) + 'TER\n'
               + pdb_atom(2, 'CA', 'ALA', 'C', 1, 10.0) + 'TER\n'
               + pdb_atom(3, 'NA', 'NA', '', 1, 20.0, 'HETATM', 'NA') + 'END\n')
        atoms = rows(convert(pdb))
        self.assertEqual([(a['label_asym_id'], a['auth_asym_id']) for a in atoms],
                         [('A', 'A'), ('B', 'C'), ('C', 'B')])

    def test_later_models_without_ter(self):
        pdb = ('MODEL        1\n' + pdb_atom(1, 'CA', 'ALA', 'A', 1, 0.0) + 'TER\n'
               + pdb_atom(2, 'ZN', 'ZN', 'A', 101, 9.0, 'HETATM', 'ZN') + 'ENDMDL\n'
               'MODEL        2\n' + pdb_atom(1, 'CA', 'ALA', 'A', 1, 0.0)
               + pdb_atom(2, 'ZN', 'ZN', 'A', 101, 9.0, 'HETATM', 'ZN') + 'ENDMDL\nEND\n')
        atoms = rows(convert(pdb))
        self.assertEqual([(a['label_asym_id'], a['label_seq_id']) for a in atoms],
                         [('A', '1'), ('B', '.'), ('A', '1'), ('B', '.')])


class TestTopdbTer(unittest.TestCase):

    def test_ter_after_each_polymer(self):
        lines = run_tool(cif_topdb, WWPDB).splitlines()
        ter = [i for i, line in enumerate(lines) if line.startswith('TER')]
        self.assertEqual(len(ter), 1)
        self.assertEqual(lines[ter[0]].rstrip(), 'TER      16      GLY A  11A')
        self.assertEqual(len(lines[ter[0]]), 80)
        self.assertEqual(lines[ter[0] - 1][:6], 'ATOM  ')
        self.assertEqual(lines[ter[0] + 1][:11], 'HETATM   17')
        nmr = run_tool(cif_topdb, NMR).splitlines()
        self.assertEqual([line[:6] for line in nmr if not line.startswith('ATOM')],
                         ['MODEL ', 'TER   ', 'ENDMDL', 'MODEL ', 'TER   ', 'ENDMDL', 'END   '])

    def test_round_trip_keeps_label_ids(self):
        back = convert(run_tool(cif_topdb, WWPDB))
        for col in ('label_asym_id', 'label_seq_id', 'auth_asym_id', 'auth_seq_id'):
            self.assertEqual(column(back, col), column(WWPDB, col), col)


class TestQuotedNulls(unittest.TestCase):
    TEXT = ("data_x\nloop_\n_a.x\n_a.y\n'.' .\n\"?\" ?\n;\n.\n;\nz\n#\n"
            "_b.v '.'\n")

    def test_read_as_text(self):
        values = [ev.values for ev in _cif.parse(io.StringIO(self.TEXT), ('_a', '_b'))
                  if isinstance(ev, _cif.Row)]
        self.assertIsInstance(values[0][0], _cif.Text)
        self.assertNotIn(values[0][0], ('.', '?'))  # not a null
        self.assertIn(values[0][1], ('.', '?'))     # a null
        self.assertNotEqual(values[1][0], '?')
        self.assertEqual(values[1][0], _cif.Text('?'))
        self.assertEqual({_cif.Text('.'), '.'}, {'.', _cif.Text('.')})
        self.assertEqual(len({_cif.Text('.'), '.'}), 2)
        self.assertEqual(values[3], (_cif.Text('.'),))

    def test_written_in_quotes(self):
        self.assertEqual(_cif.quote(_cif.Text('.')), "'.'")
        self.assertEqual(_cif.quote(_cif.Text('?')), "'?'")
        self.assertEqual(_cif.quote('.'), '.')
        loop = _cif.Loop.new('_a', ['x', 'y'])
        out = io.StringIO()
        _cif.write([loop, _cif.Row(loop, [_cif.Text('.'), '.']), _cif.LoopEnd(loop)], out)
        self.assertEqual(out.getvalue(), "loop_\n_a.x\n_a.y\n'.' .\n")

    def test_edit_keeps_quotes(self):
        events = list(_cif.parse(io.StringIO(self.TEXT), '_a'))
        next(ev for ev in events if isinstance(ev, _cif.Row))['y'] = 'q'
        out = io.StringIO()
        _cif.write(events, out)
        self.assertEqual(out.getvalue(), self.TEXT.replace("'.' .", "'.' q"))


class TestSelaltlocLinks(unittest.TestCase):

    def test_links_follow_the_chosen_conformer(self):
        text = WWPDB.replace(
            "_struct_conn.pdbx_dist_value \nhydrog1 hydrog A SER 2 OG B SO4 . O1 2.95 \n",
            "_struct_conn.pdbx_dist_value \n_struct_conn.pdbx_ptnr1_label_alt_id \n"
            "_struct_conn.pdbx_ptnr2_label_alt_id \nhydrog1 hydrog A SER 2 OG B SO4 . O1 2.95 A . \n"
            "hydrog2 hydrog A SER 2 OG B SO4 . O1 3.10 B . \n")
        for altloc, kept in ((None, 'hydrog1'), ('B', 'hydrog2')):
            with self.subTest(altloc=altloc):
                out = run_tool(cif_selaltloc, text, altloc)
                links = rows(out, '_struct_conn')
                self.assertEqual([(r['id'], r['pdbx_ptnr1_label_alt_id']) for r in links],
                                 [(kept, '.')])

    def test_without_altloc_links(self):
        self.assertEqual(run_tool(cif_selaltloc, WWPDB, None).count('hydrog1'), 1)


class TestCat(unittest.TestCase):

    def test_table(self):
        lines = run_tool(cif_cat, WWPDB, '_atom_site', ['id', 'AUTH_ASYM_ID', 'Cartn_x']).splitlines()
        self.assertEqual(lines[0], 'id\tAUTH_ASYM_ID\tCartn_x')
        self.assertEqual(lines[1], '1\tA\t21.598')
        self.assertEqual(len(lines), len(rows(WWPDB)) + 1)
        full = run_tool(cif_cat, WWPDB, '_atom_site').splitlines()
        self.assertEqual(full[0].split('\t'), rows(WWPDB)[0].loop.columns)

    def test_pairs_and_escapes(self):
        self.assertEqual(run_tool(cif_cat, WWPDB, '_cell').splitlines()[1].split('\t')[:2],
                         ['1XYZ', '53.632'])
        text = "data_x\n_t.a\n;line 1\nline\t2\n;\n_t.b 'a\\b'\n_t.c '.'\n"
        self.assertEqual(run_tool(cif_cat, text, '_t'), 'a\tb\tc\nline 1\\nline\\t2\ta\\\\b\t.\n')

    def test_header_once_per_column_set(self):
        block = "data_{0}\n_c.id {0}\n_c.name x\n"
        text = block.format('a') + block.format('b') + "data_c\n_c.id c\n"
        self.assertEqual(run_tool(cif_cat, text, '_c'), 'id\tname\na\tx\nb\tx\nid\nc\n')

    def test_errors(self):
        with self.assertRaises(_cif.CIFError):
            run_tool(cif_cat, WWPDB, '_nothing')
        with self.assertRaises(_cif.CIFError):
            run_tool(cif_cat, WWPDB, '_atom_site', ['bogus'])


class TestDelKeepCat(unittest.TestCase):

    def test_delcat(self):
        out = run_tool(cif_delcat, WWPDB, ['_cell', '_atom_site_anisotrop'])
        self.assertEqual(rows(out, '_cell'), [])
        self.assertEqual(rows(out, '_atom_site_anisotrop'), [])
        self.assertEqual(len(rows(out)), len(rows(WWPDB)))
        self.assertNotIn('_cell.', out)
        self.assertNotIn('#\n#', out.replace(' ', ''))  # one separator between categories
        self.assertEqual(run_tool(cif_delcat, WWPDB, ['_nothing']), WWPDB)

    def test_prefix(self):
        out = run_tool(cif_delcat, WWPDB, ['_atom_site*'])
        self.assertNotIn('_atom_site', out)
        self.assertIn('_struct_conn.id', out)
        kept = run_tool(cif_keepcat, WWPDB, ['_STRUCT*'])
        self.assertEqual([line for line in kept.splitlines() if line.startswith('_')],
                         [line for line in WWPDB.splitlines() if line.startswith('_struct')])

    def test_keepcat(self):
        out = run_tool(cif_keepcat, WWPDB, ['_cell'])
        self.assertTrue(out.startswith('data_1XYZ\n'))
        self.assertEqual(values(out, '_cell'), values(WWPDB, '_cell'))
        self.assertEqual(rows(out), [])
        for text in (WWPDB, NMR, read('multi_block.cif')):
            self.assertEqual(run_tool(cif_keepcat, text, ['_atom_site', '_struct_conn']),
                             run_tool(cif_keepcoord, text))

    def test_bad_names(self):
        for option in ('-_', '-a*b', '-_cell.length_a', '-,'):
            with self.subTest(option=option):
                proc = cli('cif_delcat', option, os.path.join(DATA, 'wwpdb_like.cif'))
                self.assertEqual(proc.returncode, 1)
                self.assertIn(b'ERROR!!', proc.stderr)


class TestBlocks(unittest.TestCase):
    MULTI = read('multi_block.cif')

    def test_list_and_select(self):
        self.assertEqual(run_tool(cif_blocks, self.MULTI), 'first\nsecond\n')
        second = run_tool(cif_blocks, self.MULTI, ['SECOND'])
        self.assertTrue(second.startswith('data_second\n'))
        self.assertEqual(values(second), values(self.MULTI)[-len(values(second)):])
        both = run_tool(cif_blocks, self.MULTI, ['second', 'first'])
        self.assertEqual(run_tool(cif_blocks, both), 'first\nsecond\n')  # file order
        with self.assertRaises(_cif.CIFError):
            run_tool(cif_blocks, self.MULTI, ['third'])

    def test_stops_after_the_last_block(self):
        def lines():
            yield 'data_a\n_x.v 1\n'
            yield 'data_b\n'
            raise AssertionError('read past the requested blocks')
        out = io.StringIO()
        _cif.write(cif_blocks.run(lines(), ['a']), out)
        self.assertEqual(out.getvalue(), 'data_a\n_x.v 1\n')


class TestLabel2Auth(unittest.TestCase):

    def test_chains_and_numbers(self):
        out = run_tool(cif_label2auth, CONN_AUTH, True, True)
        atoms = rows(out)
        self.assertEqual([a['auth_asym_id'] for a in atoms], [a['label_asym_id'] for a in atoms])
        self.assertEqual([(a['auth_seq_id'], a['pdbx_PDB_ins_code']) for a in atoms
                          if a['label_seq_id'] != '.'],
                         [(a['label_seq_id'], '?') for a in atoms if a['label_seq_id'] != '.'])
        self.assertEqual([a['auth_seq_id'] for a in atoms if a['label_seq_id'] == '.'],
                         ['201', '201', '301', '302'])
        link = rows(out, '_struct_conn')[0]
        self.assertEqual([link[c] for c in ('ptnr1_auth_asym_id', 'ptnr1_auth_seq_id',
                                            'ptnr2_auth_asym_id', 'ptnr2_auth_seq_id')],
                         ['A', '2', 'B', '201'])
        self.assertEqual(reference_problems(out), set())

    def test_only_one(self):
        chains = rows(run_tool(cif_label2auth, CONN_AUTH, True, False))
        self.assertEqual([a['auth_seq_id'] for a in chains], column(CONN_AUTH, 'auth_seq_id'))
        self.assertEqual([a['auth_asym_id'] for a in chains], column(CONN_AUTH, 'label_asym_id'))
        numbers = rows(run_tool(cif_label2auth, CONN_AUTH, False, True))
        self.assertEqual([a['auth_asym_id'] for a in numbers], column(CONN_AUTH, 'auth_asym_id'))
        self.assertEqual(sorted({a['auth_seq_id'] for a in numbers}), ['1', '2', '201', '3', '301', '302'])

    def test_refuses_clashes(self):
        text = CONN_AUTH.replace('? 201 SO4 A', '? 2   SO4 A')
        with self.assertRaises(_cif.CIFError) as ctx:
            run_tool(cif_label2auth, text, False, True)
        self.assertIn('would both become A 2', str(ctx.exception))
        self.assertEqual(len(rows(run_tool(cif_label2auth, text, True, True))), len(rows(text)))

    def test_options(self):
        proc = cli('cif_label2auth', '-both', os.path.join(DATA, 'wwpdb_like.cif'))
        self.assertEqual(proc.returncode, 1)
        self.assertIn(b'Invalid option', proc.stderr)


class TestAuth2Label(unittest.TestCase):

    def test_polymer_chains_follow(self):
        renamed = run_tool(cif_rplchain, WWPDB, ('A', 'H'))
        out = run_tool(cif_auth2label, renamed)
        atoms = rows(out)
        self.assertEqual(sorted({(a['label_asym_id'], a['label_comp_id']) for a in atoms}),
                         [('B', 'SO4'), ('C', 'HOH'), ('H', 'GLY'), ('H', 'MET'), ('H', 'SER')])
        self.assertEqual(column(out, 'pdbx_label_asym_id', '_atom_site_anisotrop'), ['H', 'H', 'B'])
        self.assertEqual(rows(out, '_struct_conn')[0]['ptnr1_label_asym_id'], 'H')
        self.assertEqual(reference_problems(out), set())
        self.assertEqual(run_tool(cif_auth2label, WWPDB), WWPDB)  # already the same

    def test_multichar(self):
        out = run_tool(cif_auth2label, MULTICHAR)
        self.assertEqual(column(out, 'label_asym_id'), column(MULTICHAR, 'auth_asym_id'))

    def test_refuses_merges(self):
        with self.assertRaises(_cif.CIFError) as ctx:
            run_tool(cif_auth2label, run_tool(cif_rplchain, WWPDB, ('A', 'B')))
        self.assertIn('ligand or water', str(ctx.exception))
        with self.assertRaises(_cif.CIFError) as ctx:
            run_tool(cif_auth2label, run_tool(cif_chain, MULTICHAR, 'Z'))
        self.assertIn('holds two polymers', str(ctx.exception))

    def test_command_line(self):
        proc = cli('cif_auth2label', os.path.join(DATA, 'wwpdb_like.cif'))
        self.assertEqual(proc.returncode, 0)
        self.assertEqual(proc.stdout.decode(), WWPDB)


class TestCommandLine(unittest.TestCase):

    def test_help(self):
        for name in ('cif_cat', 'cif_delcat', 'cif_keepcat', 'cif_blocks', 'cif_label2auth',
                     'cif_auth2label'):
            with self.subTest(tool=name):
                proc = cli(name, '--help')
                self.assertEqual(proc.returncode, 0)
                self.assertIn(b'Usage:', proc.stdout)

    def test_outputs_validate(self):
        path = os.path.join(DATA, 'wwpdb_like.cif')
        for name, args in (('cif_delcat', ['-_cell']), ('cif_keepcat', ['-_atom_site']),
                           ('cif_label2auth', []), ('cif_blocks', ['-1XYZ'])):
            with self.subTest(tool=name):
                proc = cli(name, *(args + [path]))
                self.assertEqual(proc.returncode, 0, proc.stderr)
                lines, errors = cif_validate.run(io.StringIO(proc.stdout.decode(), newline=''))
                self.assertEqual(errors, 0, ''.join(lines))

    def test_missing_items(self):
        path = os.path.join(DATA, 'wwpdb_like.cif')
        for name, args, message in (('cif_cat', [path], b'option is required'),
                                    ('cif_cat', ['-_nothing', path], b'No _nothing category'),
                                    ('cif_blocks', ['-nope', path], b'No data block named nope')):
            with self.subTest(tool=name, args=args):
                proc = cli(name, *args)
                self.assertEqual(proc.returncode, 1)
                self.assertIn(message, proc.stderr)


if __name__ == '__main__':
    unittest.main()
