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

"""Unit tests for the cif_* tools, on the fixtures in tests/data."""

import gzip
import io
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
import urllib.error

from support import DATA, ROOT, _cif, gemmi, gemmi_tables, our_tables, read, rows  # noqa: E402

from ciftools import (  # noqa: E402
    _atoms,
    cif_b, cif_chain, cif_delchain, cif_delhetatm, cif_delres, cif_delresname, cif_fetch,
    cif_occ, cif_reres, cif_rplchain, cif_selaltloc, cif_selatom, cif_selchain, cif_selelem,
    cif_selhetatm, cif_selmodel, cif_selres, cif_selresname, cif_shiftres, cif_tofasta, cif_wc)

WWPDB = read('wwpdb_like.cif')


def tool(module, text, *args, scheme=None):
    """Run a tool's run() on text; return its output text."""
    result = module.run(io.StringIO(text, newline=''), *args, scheme=scheme)
    if module in (cif_wc, cif_tofasta):
        return ''.join(result)
    out = io.StringIO()
    _cif.write(result, out)
    return out.getvalue()


def ids(text, category='_atom_site'):
    return [r['id'] for r in rows(text, category)]


def column(text, name, category='_atom_site'):
    return [r[name] for r in rows(text, category)]


class TestSelection(unittest.TestCase):

    def test_selchain(self):
        out = tool(cif_selchain, WWPDB, ['A'])
        self.assertEqual(out, WWPDB)  # every atom is in auth chain A
        out = tool(cif_selchain, WWPDB, ['A'], scheme='label')
        self.assertEqual(ids(out), [str(i) for i in range(1, 16)])
        self.assertEqual(ids(out, '_atom_site_anisotrop'), ['1', '2'])
        self.assertEqual(rows(out, '_struct_conn'), [])  # its SO4 partner is gone
        self.assertIn('_cell.length_a', out)

    def test_selchain_multichar(self):
        out = tool(cif_selchain, read('multichar_chains.cif'), ['Axp', 'zz'])
        self.assertEqual(ids(out), ['1', '2', '5'])

    def test_delchain(self):
        out = tool(cif_delchain, WWPDB, ['C'], scheme='label')
        self.assertEqual(len(ids(out)), 17)
        self.assertEqual(len(rows(out, '_struct_conn')), 1)
        out = tool(cif_delchain, WWPDB, ['B'], scheme='label')
        self.assertEqual(ids(out, '_atom_site_anisotrop'), ['1', '2'])
        self.assertEqual(rows(out, '_struct_conn'), [])

    def test_selres(self):
        out = tool(cif_selres, WWPDB, {11}, [])
        self.assertEqual(ids(out), [str(i) for i in range(5, 16)])  # SER 11 and GLY 11A
        self.assertNotIn('_atom_site_anisotrop.id', out)  # all its rows are gone
        out = tool(cif_selres, WWPDB, set(), [(2, 2, 1)], scheme='label')
        self.assertEqual(column(out, 'label_comp_id'), ['SER'] * 7)

    def test_selres_step(self):
        # residue occurrences: 10, 11 (SER and 11A share chain + number), 201, 301, 302
        out = tool(cif_selres, WWPDB, set(), [(float('-inf'), float('inf'), 2)])
        self.assertEqual(sorted(set(column(out, 'auth_seq_id'))), ['10', '201', '302'])

    def test_selres_parsing(self):
        singles, ranges = cif_selres.parse_selection('1,4:6,10::2,-3:')
        self.assertEqual(singles, {1})
        self.assertEqual(ranges[0], (4, 6, 1))
        self.assertEqual(ranges[1], (10, float('inf'), 2))
        self.assertEqual(ranges[2], (-3, float('inf'), 1))
        for bad in ('a', '1:2:3:4', '5:1', '1:5:0', '1,,2'):
            with self.subTest(option=bad), self.assertRaises(ValueError):
                cif_selres.parse_selection(bad)

    def test_delres(self):
        out = tool(cif_delres, WWPDB, (11, 11, 1))
        self.assertEqual(sorted(set(column(out, 'auth_seq_id'))), ['10', '201', '301', '302'])
        out = tool(cif_delres, WWPDB, (float('-inf'), float('inf'), 2))
        self.assertEqual(sorted(set(column(out, 'auth_seq_id'))), ['11', '301'])

    def test_resname(self):
        self.assertEqual(ids(tool(cif_selresname, WWPDB, ['HOH'])), ['18', '19'])
        out = tool(cif_delresname, WWPDB, ['HOH', 'SO4'])
        self.assertEqual(len(ids(out)), 15)
        self.assertEqual(ids(out, '_atom_site_anisotrop'), ['1', '2'])

    def test_selatom(self):
        self.assertEqual(ids(tool(cif_selatom, WWPDB, ['CA'])), ['2', '6', '13'])
        nucleic = tool(cif_selatom, read('nucleic.cif'), ["O5'", 'P'])
        self.assertEqual(column(nucleic, 'auth_atom_id'), ['P', "O5'", 'P'])

    def test_selelem(self):
        out = tool(cif_selelem, WWPDB, ['o'])
        self.assertEqual(ids(out), ['4', '8', '10', '11', '15', '17', '18', '19'])

    def test_hetatm(self):
        self.assertEqual(ids(tool(cif_selhetatm, WWPDB)), ['16', '17', '18', '19'])
        out = tool(cif_delhetatm, WWPDB)
        self.assertEqual(len(ids(out)), 15)
        self.assertEqual(ids(out, '_atom_site_anisotrop'), ['1', '2'])
        self.assertEqual(rows(out, '_struct_conn'), [])

    def test_selmodel(self):
        nmr = read('nmr_models.cif')
        self.assertEqual(ids(tool(cif_selmodel, nmr, [2])), ['4', '5', '6'])
        self.assertEqual(ids(tool(cif_selmodel, nmr, [3])), [])
        no_models = read('multi_block.cif')
        self.assertEqual(tool(cif_selmodel, no_models, [5]), no_models)

    def test_links_with_auth_columns(self):
        """Links are matched on every identifier present in both categories."""
        text = WWPDB.replace(
            '_struct_conn.pdbx_dist_value \nhydrog1 hydrog A SER 2 OG B SO4 . O1 2.95 \n',
            '_struct_conn.pdbx_dist_value \n_struct_conn.ptnr1_auth_seq_id \n'
            '_struct_conn.ptnr2_auth_seq_id \nhydrog1 hydrog A SER 2 OG B SO4 . O1 2.95 11 201 \n'
            'hydrog2 hydrog A SER 2 OG C HOH . O 3.10 11 301 \n')
        out = tool(cif_delres, text, (301, 301, 1))
        self.assertEqual([r['id'] for r in rows(out, '_struct_conn')], ['hydrog1'])


class TestAltloc(unittest.TestCase):

    def test_highest_occupancy(self):
        out = tool(cif_selaltloc, WWPDB, None)
        atoms = rows(out)
        self.assertEqual(len(atoms), 18)
        og = [a for a in atoms if a['auth_atom_id'] == 'OG']
        self.assertEqual([(a['id'], a['label_alt_id'], a['occupancy']) for a in og],
                         [('10', '.', '0.60')])

    def test_given_altloc(self):
        og = [a for a in rows(tool(cif_selaltloc, WWPDB, 'B')) if a['auth_atom_id'] == 'OG']
        self.assertEqual([(a['id'], a['label_alt_id']) for a in og], [('11', '.')])
        self.assertEqual(tool(cif_selaltloc, WWPDB, 'Z'), WWPDB)

    def test_microheterogeneity_is_not_mixed(self):
        # residue 5 modelled as SER (A, 0.4) or THR (B, 0.6)
        rows_ = [
            ('1', 'N', '.', 'SER', '1.00'), ('2', 'CA', 'A', 'SER', '0.40'),
            ('3', 'OG', 'A', 'SER', '0.40'), ('4', 'CA', 'B', 'THR', '0.60'),
            ('5', 'OG1', 'B', 'THR', '0.60'), ('6', 'CG2', 'B', 'THR', '0.60'),
        ]
        text = ('data_x\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.label_atom_id\n'
                '_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.auth_asym_id\n'
                '_atom_site.auth_seq_id\n_atom_site.occupancy\n')
        text += ''.join('ATOM {} {} {} {} A 5 {}\n'.format(*r) for r in rows_)
        kept = [(a['id'], a['label_alt_id']) for a in rows(tool(cif_selaltloc, text, None))]
        self.assertEqual(kept, [('1', '.'), ('4', '.'), ('5', '.'), ('6', '.')])
        kept = [a['id'] for a in rows(tool(cif_selaltloc, text, 'A'))]
        self.assertEqual(kept, ['1', '2', '3', '5', '6'])  # OG1/CG2 have no A: kept as is


class TestEdits(unittest.TestCase):

    def test_chain(self):
        out = tool(cif_chain, WWPDB, 'Z')
        self.assertEqual(set(column(out, 'auth_asym_id')), {'Z'})
        self.assertEqual(column(out, 'label_asym_id'), column(WWPDB, 'label_asym_id'))
        out = tool(cif_chain, WWPDB, 'Q', scheme='label')
        self.assertEqual(set(column(out, 'label_asym_id')), {'Q'})
        self.assertEqual(column(out, 'pdbx_label_asym_id', '_atom_site_anisotrop'), ['Q'] * 3)
        link = rows(out, '_struct_conn')[0]
        self.assertEqual((link['ptnr1_label_asym_id'], link['ptnr2_label_asym_id']), ('Q', 'Q'))

    def test_rplchain(self):
        out = tool(cif_rplchain, WWPDB, ('B', 'Q'), scheme='label')
        self.assertEqual(column(out, 'label_asym_id')[15:17], ['Q', 'Q'])
        self.assertEqual(column(out, 'pdbx_label_asym_id', '_atom_site_anisotrop'),
                         ['A', 'A', 'Q'])
        self.assertEqual(rows(out, '_struct_conn')[0]['ptnr2_label_asym_id'], 'Q')
        self.assertEqual(tool(cif_rplchain, WWPDB, ('X', 'Y')), WWPDB)

    def test_b_and_occ(self):
        self.assertEqual(set(column(tool(cif_b, WWPDB, 5), 'B_iso_or_equiv')), {'5.00'})
        self.assertEqual(set(column(tool(cif_occ, WWPDB, 0.5), 'occupancy')), {'0.50'})

    def test_reres(self):
        out = tool(cif_reres, WWPDB, 1)
        numbers = column(out, 'auth_seq_id')
        self.assertEqual(numbers, ['1'] * 4 + ['2'] * 7 + ['3'] * 4 + ['4', '4', '5', '6'])
        self.assertEqual(column(out, 'pdbx_PDB_ins_code')[11:15], ['A'] * 4)  # kept
        nmr = column(tool(cif_reres, read('nmr_models.cif'), 7), 'auth_seq_id')
        self.assertEqual(nmr, ['7'] * 6)  # restarts in each model

    def test_reres_updates_links(self):
        text = WWPDB.replace('_atom_site_anisotrop.U[3][3] \n',
                             '_atom_site_anisotrop.U[3][3] \n_atom_site_anisotrop.pdbx_auth_seq_id \n'
                             '_atom_site_anisotrop.pdbx_auth_asym_id \n'
                             ).replace('0.3911 0.3602 0.4012 \n', '0.3911 0.3602 0.4012 10 A \n'
                             ).replace('0.3801 0.3699 0.3920 \n', '0.3801 0.3699 0.3920 10 A \n'
                             ).replace('0.4520 0.4401 0.4488 \n', '0.4520 0.4401 0.4488 201 A \n')
        out = tool(cif_reres, text, 100)
        self.assertEqual(column(out, 'pdbx_auth_seq_id', '_atom_site_anisotrop'),
                         ['100', '100', '103'])

    def test_shiftres(self):
        out = tool(cif_shiftres, WWPDB, 5)
        self.assertEqual(sorted(set(column(out, 'auth_seq_id')), key=int),
                         ['15', '16', '206', '306', '307'])
        out = tool(cif_shiftres, WWPDB, -1, scheme='label')
        self.assertEqual(column(out, 'label_seq_id'), column(WWPDB, 'label_seq_id')[:0] +
                         ['0'] * 4 + ['1'] * 7 + ['2'] * 4 + ['.'] * 4)
        link = rows(out, '_struct_conn')[0]
        self.assertEqual((link['ptnr1_label_seq_id'], link['ptnr2_label_seq_id']), ('1', '.'))

    def test_edits_change_only_their_column(self):
        out = tool(cif_b, WWPDB, 29.41)
        changed = [(a, b) for a, b in zip(WWPDB.splitlines(), out.splitlines()) if a != b]
        self.assertEqual(len(changed), 18)  # every atom but the one already at 29.41
        for old, new in changed:
            self.assertEqual(old.split()[:14], new.split()[:14])
            self.assertEqual(old.split()[15:], new.split()[15:])


class TestReports(unittest.TestCase):

    def test_wc(self):
        self.assertEqual(tool(cif_wc, WWPDB, ''), (
            'No. models:\t1\n'
            'No. chains:\t1\t(   1.0/model)\n'
            'No. residues:\t3\t(   3.0/model)\n'
            'No. atoms:\t15\t(  15.0/model)\n'
            'No. HETATM:\t4\n'
            'Multiple Occ.:\tTrue\n'
            'Res. Inserts:\tTrue\n'))
        self.assertEqual(tool(cif_wc, read('nmr_models.cif'), 'mc'), (
            'No. models:\t2\n\t->\t1,2\n'
            'No. chains:\t2\t(   1.0/model)\n\t->\tA,A\n'))
        self.assertEqual(tool(cif_wc, WWPDB, 'h'), 'No. HETATM:\t4\n\t->\tHOH,SO4\n')

    def test_tofasta(self):
        self.assertEqual(tool(cif_tofasta, WWPDB, False, False), '>1XYZ|A\nMSGXXX\n')
        self.assertEqual(tool(cif_tofasta, read('nucleic.cif'), True, False), '>3RNA|A\nGC\n')
        chains = tool(cif_tofasta, read('multichar_chains.cif'), False, False)
        self.assertEqual(chains, '>4BIG|A1,Axp,zz\nDKE\n')
        with self.assertRaises(_cif.CIFError):
            tool(cif_tofasta, WWPDB, False, True)

    def test_tofasta_seqres(self):
        text = WWPDB + ('loop_\n_pdbx_poly_seq_scheme.asym_id\n_pdbx_poly_seq_scheme.mon_id\n'
                        '_pdbx_poly_seq_scheme.pdb_strand_id\n'
                        + ''.join('A {} A\n'.format(n) for n in ['GLY', 'MET', 'SER', 'GLY', 'DA']))
        self.assertEqual(tool(cif_tofasta, text, True, True), '>1XYZ|A\nGMSGA\n')

    def test_tofasta_wraps(self):
        text = ('data_long\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.auth_asym_id\n'
                '_atom_site.auth_seq_id\n_atom_site.auth_comp_id\n')
        text += ''.join('ATOM {0} A {0} ALA\n'.format(i) for i in range(1, 131))
        self.assertEqual(tool(cif_tofasta, text, False, False),
                         '>long|A\n' + 'A' * 60 + '\n' + 'A' * 60 + '\n' + 'A' * 10 + '\n')


class SeekablePipe(io.TextIOWrapper):
    """A pipe that claims to be seekable, as pipes do on Windows."""

    def seekable(self):
        return True


def pipe_with(text):
    read_fd, write_fd = os.pipe()
    with os.fdopen(write_fd, 'w', newline='') as writer:  # small: fits the pipe buffer
        writer.write(text)
    return SeekablePipe(os.fdopen(read_fd, 'rb'), encoding='utf-8', newline='')


class TestSource(unittest.TestCase):
    """Tools that read their input twice must copy pipes, even 'seekable' ones."""

    def test_windows_like_pipe_is_copied(self):
        pipe = pipe_with(WWPDB)
        try:
            self.assertFalse(_atoms.rereadable(pipe))
            source = _atoms.Source(pipe)
            self.assertEqual(''.join(source.lines()), WWPDB)
            self.assertEqual(''.join(source.lines()), WWPDB)  # and again
        finally:
            pipe.close()

    def test_tool_reading_a_windows_like_pipe(self):
        pipe = pipe_with(WWPDB)
        try:
            out = io.StringIO()
            _cif.write(cif_selchain.run(pipe, ['C'], scheme='label'), out)
            self.assertEqual(ids(out.getvalue()), ['18', '19'])
        finally:
            pipe.close()

    def test_files_are_read_in_place(self):
        self.assertTrue(_atoms.rereadable(io.StringIO(WWPDB)))
        with _cif.open_input(os.path.join(DATA, 'wwpdb_like.cif')) as fh:
            self.assertTrue(_atoms.rereadable(fh))
        tmp = tempfile.mkdtemp()
        try:
            path = os.path.join(tmp, 'x.cif.gz')
            with gzip.open(path, 'wt', newline='') as fh:
                fh.write(WWPDB)
            with _cif.open_input(path) as fh:
                self.assertTrue(_atoms.rereadable(fh))
                source = _atoms.Source(fh)
                self.assertEqual(''.join(source.lines()), ''.join(source.lines()))
        finally:
            shutil.rmtree(tmp)


class FakeResponse(io.BytesIO):
    pass


class TestFetch(unittest.TestCase):

    def opener(self, responses):
        calls = []

        def opener(request, timeout):
            calls.append(request.full_url)
            data = responses[request.full_url]
            if isinstance(data, Exception):
                raise data
            return FakeResponse(data)
        return opener, calls

    def test_rcsb(self):
        url = 'https://files.rcsb.org/download/1xyz.cif.gz'
        opener, calls = self.opener({url: gzip.compress(WWPDB.encode())})
        self.assertEqual(''.join(cif_fetch.fetch('1XYZ', opener=opener)), WWPDB)
        self.assertEqual(calls, [url])

    def test_urls(self):
        self.assertEqual(cif_fetch.url_for('1ABC', True)[0],
                         'https://files.rcsb.org/download/1abc-assembly1.cif.gz')
        self.assertEqual(cif_fetch.url_for('pdb_00001abc')[0],
                         'https://files.rcsb.org/download/1abc.cif.gz')

    def test_alphafold(self):
        api = 'https://alphafold.ebi.ac.uk/api/prediction/P69905'
        model = 'https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_v6.cif'
        listing = b'[{"entryId": "AF-P69905-F1", "cifUrl": "' + model.encode() + b'"}]'
        opener, calls = self.opener({api: listing, model: WWPDB.encode()})
        self.assertEqual(''.join(cif_fetch.fetch('AF-P69905-F1', opener=opener)), WWPDB)
        self.assertEqual(calls, [api, model])

    def test_errors(self):
        url = 'https://files.rcsb.org/download/9zzz.cif.gz'
        opener, _ = self.opener({url: urllib.error.HTTPError(url, 404, 'Not Found', {}, None)})
        with self.assertRaises(cif_fetch.FetchError) as ctx:
            cif_fetch.fetch('9zzz', opener=opener)
        self.assertIn('404', str(ctx.exception))
        url = 'https://files.rcsb.org/download/1xyz.cif.gz'
        opener, _ = self.opener({url: gzip.compress(WWPDB.encode())[:200]})
        with self.assertRaises(_cif.CIFError):
            list(cif_fetch.fetch('1xyz', opener=opener))


def cli(name, *args, stdin=b''):
    return subprocess.run([sys.executable, '-m', 'ciftools.' + name] + list(args), input=stdin,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=ROOT)


class TestCommandLine(unittest.TestCase):
    PATH = os.path.join(DATA, 'wwpdb_like.cif')

    def test_help(self):
        proc = cli('cif_selchain', '--help')
        self.assertEqual(proc.returncode, 0)
        self.assertIn(b'Usage:', proc.stdout)

    def test_errors(self):
        path = self.PATH
        cases = [
            ('cif_selchain', [path], b'option is required'),
            ('cif_selchain', ['-A', 'missing.cif'], b'File not found'),
            ('cif_selchain', ['-A', '-B', path], b'Unexpected argument'),
            ('cif_delhetatm', ['-A', path], b'takes no options'),
            ('cif_delres', ['-5', path], b'Residue range must be'),
            ('cif_b', ['-x', path], b'invalid b-factor'),
            ('cif_wc', ['-x', path], b'not valid'),
            ('cif_selchain', ['-A', '--segid', path], b'Unexpected argument'),
        ]
        for name, args, message in cases:
            with self.subTest(tool=name, args=args):
                proc = cli(name, *args)
                self.assertEqual(proc.returncode, 1)
                self.assertIn(b'ERROR!!', proc.stderr)
                self.assertIn(message, proc.stderr)
                self.assertEqual(proc.stdout, b'')

    def test_malformed_input(self):
        proc = cli('cif_selchain', '-A',
                   stdin=b"data_x\nloop_\n_atom_site.id\n_atom_site.auth_asym_id\n1 'A\n")
        self.assertEqual(proc.returncode, 1)
        self.assertIn(b'line 5', proc.stderr)
        self.assertNotIn(b'Traceback', proc.stderr)

    def test_pipeline(self):
        with open(self.PATH, 'rb') as fh:
            data = fh.read()
        first = cli('cif_selchain', '--label', '-A', stdin=data)
        second = cli('cif_wc', '-a', stdin=first.stdout)
        self.assertEqual(second.returncode, 0)
        self.assertTrue(second.stdout.startswith(b'No. atoms:\t15\t'))

    def test_scheme_from_environment(self):
        with open(self.PATH, 'rb') as fh:
            data = fh.read()
        env = dict(os.environ, CIFTOOLS_SCHEME='label')
        proc = subprocess.run([sys.executable, '-m', 'ciftools.cif_selchain', '-C'], input=data,
                              stdout=subprocess.PIPE, cwd=ROOT, env=env)
        self.assertEqual(ids(proc.stdout.decode()), ['18', '19'])

    def test_gzip_file(self):
        tmp = tempfile.mkdtemp()
        try:
            path = os.path.join(tmp, 'x.cif.gz')
            with open(self.PATH, 'rb') as src, gzip.open(path, 'wb') as dst:
                dst.write(src.read())
            proc = cli('cif_selhetatm', path)
            self.assertEqual(proc.returncode, 0)
            self.assertEqual(ids(proc.stdout.decode()), ['16', '17', '18', '19'])
        finally:
            shutil.rmtree(tmp)

    def test_closed_pipe(self):
        """`cif_x big.cif | head -1` ends quietly."""
        tmp = tempfile.mkdtemp()
        try:
            path = os.path.join(tmp, 'big.cif')
            with open(path, 'w') as fh:
                fh.write('data_big\nloop_\n_atom_site.id\n_atom_site.B_iso_or_equiv\n')
                fh.writelines('{} 1.0\n'.format(i) for i in range(300000))
            proc = subprocess.Popen([sys.executable, '-m', 'ciftools.cif_b', '-2', path],
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=ROOT)
            proc.stdout.readline()
            proc.stdout.close()
            stderr = proc.stderr.read()
            proc.wait()
            proc.stderr.close()
            self.assertEqual(stderr, b'')
            self.assertEqual(proc.returncode, 0)
        finally:
            shutil.rmtree(tmp)


@unittest.skipUnless(gemmi, 'gemmi not installed (optional cross-check)')
class TestOutputIsValid(unittest.TestCase):
    """Every tool's output on every fixture is read the same way by gemmi."""

    CALLS = [
        (cif_selchain, (['A'],)), (cif_delchain, (['B'],)), (cif_selres, ({1, 11}, [])),
        (cif_delres, ((1, 10, 1),)), (cif_selresname, (['HOH'],)), (cif_delresname, (['HOH'],)),
        (cif_selatom, (['CA'],)), (cif_selelem, (['C'],)), (cif_selhetatm, ()),
        (cif_delhetatm, ()), (cif_selaltloc, (None,)), (cif_selmodel, ([1],)),
        (cif_chain, ('XY',)), (cif_rplchain, (('A', 'B'),)), (cif_b, (1.5,)), (cif_occ, (0.5,)),
        (cif_reres, (-3,)), (cif_shiftres, (1000,)),
    ]

    def test_fixtures(self):
        for name in sorted(os.listdir(DATA)):
            if not name.endswith('.cif') or name == 'empty_atom_site.cif':
                continue
            text = read(name)
            for module, args in self.CALLS:
                with self.subTest(fixture=name, tool=module.__name__):
                    try:
                        out = tool(module, text, *args)
                    except _cif.CIFError as err:  # a fixture without the needed column
                        self.assertIn('column', str(err))
                        continue
                    self.assertEqual(our_tables(out), gemmi_tables(out))


if __name__ == '__main__':
    unittest.main()
