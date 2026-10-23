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
Cross-check: each cif_X tool against the pdb_X tool it ports.

Each structure is converted to PDB format (row by row, so the atom order is
the same), both tools are run with the same command line, and the atoms
they output are compared. Skipped unless pdb-tools is installed.

The structures are the fixtures plus the real entries, when fetched;
entries that do not fit the PDB format (multi-character chains, more than
99,999 atoms) are skipped.
"""

import contextlib
import importlib
import io
import os
import shutil
import tempfile
import unittest

from support import DATA, _cif, read, read_path, real_entries  # noqa: E402

from ciftools import _atoms  # noqa: E402

try:
    import pdbtools  # noqa: F401
except ImportError:
    pdbtools = None

FIXTURES = ['wwpdb_like.cif', 'nmr_models.cif', 'nucleic.cif', 'column_order.cif', 'af3_like.cif']
NULLS = ('.', '?')

# Differences that follow from documented deviations (pdb-tools behaviour
# that cif-tools does not copy): {(tool, args, structure or None): reason}
KNOWN = {
    ('selaltloc', '', None): 'pdb_selaltloc groups residues by comparing numbers as text and '
                             'ignores insertion codes; cif_selaltloc chooses per residue',
    ('selaltloc', '-A', None): 'pdb_selaltloc -A drops atoms without altloc in residues that '
                               'have one; cif_selaltloc keeps them',
    ('tofasta', '', '1kx5'): 'pdb_tofasta writes DNA as X; cif_tofasta gives the bases',
    ('tofasta', '-multi', '1kx5'): 'pdb_tofasta writes DNA as X; cif_tofasta gives the bases',
}


def known(name, args, entry):
    args = ' '.join(args)
    return KNOWN.get((name, args, entry)) or KNOWN.get((name, args, None))


def to_pdb(text):
    """PDB-format text of a structure (auth identifiers), or None if it does not fit."""
    out, entry_id, models = [], None, []
    atoms = []
    for ev in _cif.parse(io.StringIO(text, newline=''), ('_atom_site', '_entry')):
        if isinstance(ev, _cif.Loop) and ev.category == '_atom_site':
            f = _atoms.Fields(ev, 'auth')
        elif isinstance(ev, _cif.Row) and ev.loop.category == '_entry':
            entry_id = ev.get('id')
        elif isinstance(ev, _cif.Row):
            atoms.append((f, ev.values))
    for f, v in atoms:
        model = v[f.model] if f.model is not None else '1'
        if model not in models:
            models.append(model)
    multi = len(models) > 1
    out.append('HEADER' + ' ' * 56 + (entry_id or 'XXXX')[:4] + '\n')
    current = None
    for f, v in atoms:
        model = v[f.model] if f.model is not None else '1'
        if multi and model != current:
            if current is not None:
                out.append('ENDMDL\n')
            out.append('MODEL     {:>4}\n'.format(model))
            current = model
        chain, name, resname = v[f.chain], v[f.atom], v[f.resname]
        element = v[f.element] if f.element is not None else ''
        serial = v[f.id] if f.id is not None else '0'
        resnum = _atoms.as_int(v[f.resnum])
        if (len(chain) != 1 or len(name) > 4 or len(resname) > 3 or len(element) > 2
                or resnum is None or not -999 <= resnum <= 9999 or int(serial) > 99999):
            return None
        altloc = v[f.altloc] if f.altloc is not None else '.'
        icode = v[f.icode] if f.icode is not None else '?'
        occ = v[f.occupancy] if f.occupancy is not None else '1.0'
        b = v[f.bfactor] if f.bfactor is not None else '0.0'
        x, y, z = (float(v[f.loop.index(c)]) for c in ('Cartn_x', 'Cartn_y', 'Cartn_z'))
        padded = ' ' + name.ljust(3) if len(name) < 4 and len(element) == 1 else name.ljust(4)
        out.append('{:<6}{:>5} {:4}{:1}{:>3} {:1}{:>4}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}'
                   '{:>6.2f}{:>6.2f}          {:>2}  \n'.format(
                       v[f.group], serial, padded, '' if altloc in NULLS else altloc, resname,
                       chain, resnum, '' if icode in NULLS else icode, x, y, z,
                       float(occ), float(b), element))
    if multi:
        out.append('ENDMDL\n')
    out.append('END\n')
    return ''.join(out)


def pdb_atoms(text):
    """Comparable atom tuples of PDB-format text."""
    atoms, model = [], '1'
    for line in text.splitlines():
        if line.startswith('MODEL'):
            model = str(int(line[10:14]))
        elif line.startswith(('ATOM', 'HETATM')):
            atoms.append((model, line[21], line[22:26].strip(), line[26].strip(),
                          line[17:20].strip(), line[12:16].strip(), line[16].strip(),
                          '{:.2f}'.format(float(line[54:60])), '{:.2f}'.format(float(line[60:66])),
                          line[30:38].strip(), line[38:46].strip(), line[46:54].strip()))
    return atoms


def cif_atoms(text):
    """The same tuples for mmCIF text (auth identifiers)."""
    atoms = []
    multi = None
    for ev in _cif.parse(io.StringIO(text, newline='')):
        if isinstance(ev, _cif.Loop):
            f = _atoms.Fields(ev, 'auth')
        elif isinstance(ev, _cif.Row):
            v = ev.values
            model = v[f.model] if f.model is not None else '1'
            icode = v[f.icode] if f.icode is not None else '?'
            altloc = v[f.altloc] if f.altloc is not None else '.'
            occ = v[f.occupancy] if f.occupancy is not None else '1.0'
            b = v[f.bfactor] if f.bfactor is not None else '0.0'
            coords = ['{:.3f}'.format(float(v[f.loop.index(c)]))
                      for c in ('Cartn_x', 'Cartn_y', 'Cartn_z')]
            atoms.append((model, v[f.chain], str(_atoms.as_int(v[f.resnum])),
                          '' if icode in NULLS else icode, v[f.resname], v[f.atom],
                          '' if altloc in NULLS else altloc, '{:.2f}'.format(float(occ)),
                          '{:.2f}'.format(float(b))) + tuple(coords))
            multi = multi or model != '1'
    return atoms


def run_pdb_tool(name, args, path):
    module = importlib.import_module('pdbtools.pdb_' + name)
    parsed = module.check_input(args + [path])
    parsed = parsed if isinstance(parsed, tuple) else (parsed,)
    fhandle, options = parsed[0], parsed[1:]
    out = io.StringIO()
    try:
        if name == 'wc':
            with contextlib.redirect_stdout(out):
                module.run(fhandle, *options)
        else:
            out.write(''.join(module.run(fhandle, *options)))
    except SystemExit:  # e.g. pdb_reres: 'Cannot set residue number above 9999'
        return None
    finally:
        if hasattr(fhandle, 'close'):  # pdb_selres hands back an itertools.chain
            fhandle.close()
    return out.getvalue()


def run_cif_tool(name, args, path):
    module = importlib.import_module('ciftools.cif_' + name)
    parsed = module.check_input(args + [path])
    fhandle, options = parsed[0], parsed[1:]
    try:
        result = module.run(fhandle, *options)
        if name in ('wc', 'tofasta'):
            return ''.join(result)
        out = io.StringIO()
        _cif.write(result, out)
        return out.getvalue()
    finally:
        fhandle.close()


def cases(text):
    """(tool, args) pairs to compare for a structure."""
    chains = list(dict.fromkeys(a[1] for a in cif_atoms(text)))
    return [
        ('selchain', ['-' + chains[0]]), ('delchain', ['-' + chains[0]]),
        ('selres', ['-1:20']), ('selres', ['-::3']), ('selres', ['-5,10:15']),
        ('delres', ['-5:30']), ('delres', ['-::4']),
        ('selresname', ['-ALA,GLY']), ('delresname', ['-HOH']),
        ('selatom', ['-CA,N']), ('selelem', ['-C,N']),
        ('selhetatm', []), ('delhetatm', []),
        ('selaltloc', []), ('selaltloc', ['-A']),
        ('selmodel', ['-1']), ('selmodel', ['-2']),
        ('chain', ['-Z']), ('rplchain', ['-{}:Z'.format(chains[0])]),
        ('b', ['-5.5']), ('occ', ['-0.5']),
        ('reres', ['-1']), ('reres', ['-100']), ('shiftres', ['-10']), ('shiftres', ['--5']),
        ('wc', []), ('wc', ['-mcroi']),
        ('tofasta', []), ('tofasta', ['-multi']),
    ]


def structures():
    found = [(name, read(name)) for name in FIXTURES]
    found += [(os.path.basename(p).split('.')[0], read_path(p)) for p in real_entries()]
    return found


def compare(name, args, pdb_path, cif_path):
    """None if both tools agree (or pdb-tools cannot run it), else a short
    description of the difference."""
    pdb_out = run_pdb_tool(name, args, pdb_path)
    if pdb_out is None:  # beyond a PDB-format limit; cif-tools has none
        return None
    cif_out = run_cif_tool(name, args, cif_path)
    if name == 'tofasta':  # the PDB header holds only 4 characters of the entry id
        pdb_out = '\n'.join(line.split('|', 1)[-1] for line in pdb_out.split('\n'))
        cif_out = '\n'.join(line.split('|', 1)[-1] for line in cif_out.split('\n'))
    if name in ('wc', 'tofasta'):
        if pdb_out != cif_out:
            return 'pdb: {!r}\n    cif: {!r}'.format(pdb_out[:300], cif_out[:300])
        return None
    expected, got = pdb_atoms(pdb_out), cif_atoms(cif_out)
    if expected == got:
        return None
    if sorted(expected) == sorted(got):
        return 'same atoms, different order'
    missing = [a for a in expected if a not in set(got)]
    extra = [a for a in got if a not in set(expected)]
    return '{} atoms from pdb-tools, {} from cif-tools; e.g. only pdb: {} only cif: {}'.format(
        len(expected), len(got), missing[:2], extra[:2])


@unittest.skipUnless(pdbtools, 'pdb-tools not installed')
class TestAgainstPdbTools(unittest.TestCase):

    def test_tools(self):
        tmp = tempfile.mkdtemp()
        report = []
        try:
            for entry, text in structures():
                pdb = to_pdb(text)
                if pdb is None:
                    continue
                pdb_path = os.path.join(tmp, entry + '.pdb')
                cif_path = os.path.join(tmp, entry + '.cif')
                with open(pdb_path, 'w') as fh:
                    fh.write(pdb)
                with open(cif_path, 'w', newline='') as fh:
                    fh.write(text)
                for name, args in cases(text):
                    with self.subTest(entry=entry, tool=name, args=args):
                        difference = compare(name, args, pdb_path, cif_path)
                        if difference and not known(name, args, entry):
                            report.append('{} {} {}: {}'.format(entry, name, ' '.join(args),
                                                                 difference))
                            self.fail(difference)
        finally:
            shutil.rmtree(tmp)
            if report and os.environ.get('CIFTOOLS_REPORT'):
                print('\n' + '\n'.join(report))


if __name__ == '__main__':
    unittest.main()
