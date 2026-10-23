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
Summarizes the contents of an mmCIF file, like the wc command in UNIX.

By default, this tool produces a general summary, but you can use several
options to produce focused but more detailed summaries:
    [m] - no. of models.
    [c] - no. of chains (plus per-model if multi-model file).
    [r] - no. of residues (plus per-model if multi-model file).
    [a] - no. of atoms (plus per-model if multi-model file).
    [h] - no. of HETATM (plus per-model if multi-model file).
    [o] - presence of disordered atoms (altloc).
    [i] - presence of insertion codes.

Counts follow pdb_wc: residues and atoms count ATOM records only, HETATM
records are counted separately, and the altloc and insertion-code flags
look at ATOM records. The lists printed with [a] and [h] give atom names
and HETATM residue names (pdb_wc prints padded atom names for both).

Usage:
    cif_wc [-<option>] <cif file>

Example:
    cif_wc 1abc.cif
    cif_wc -mc 1abc.cif

This program is part of the cif-tools suite, a port of pdb-tools to mmCIF.
"""

import sys

from ciftools import _atoms, _cif, _cli

OPTIONS = 'mcrahoi'
NULLS = ('.', '?')


def check_input(args):
    """Return (file handle, option letters, identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__)
    option = option or ''
    bad = sorted(set(option) - set(OPTIONS))
    if bad:
        _cli.fail(__doc__, 'The following options are not valid: {!r}'.format(''.join(bad)))
    return fhandle, option, scheme


def summarize(fhandle, scheme=None):
    """Collect the sets counted by cif_wc."""
    models, chains, residues, atoms, hetatm = set(), set(), set(), set(), set()
    has_altloc = has_icode = False
    for ev in _cif.parse(fhandle):
        if type(ev) is _cif.Row:
            values = ev.values
            model = values[f.model] if f.model is not None else None
            if model is not None:
                models.add(model)
            atom = identity(values)
            chains.add((model, atom[3]))
            if f.group is not None and values[f.group] == 'HETATM':
                hetatm.add((model,) + atom)
                continue
            residues.add((model, atom[2], atom[3], atom[4]))
            atoms.add((model,) + atom)
            if not has_altloc and atom[1] not in NULLS:
                has_altloc = True
            if not has_icode and atom[5] not in NULLS:
                has_icode = True
        elif type(ev) is _cif.Loop:
            f = _atoms.Fields(ev, scheme)
            f.need('chain', 'resnum', 'resname', 'atom')
            # (atom name, altloc, resname, chain, resnum, icode)
            identity = _atoms.getter(f.atom, f.altloc, f.resname, f.chain, f.resnum, f.icode)
    return {
        'models': models, 'chains': chains, 'residues': residues, 'atoms': atoms,
        'hetatm': hetatm, 'altloc': has_altloc, 'icode': has_icode,
    }


def run(fhandle, option='', scheme=None):
    """Return the summary lines."""
    s = summarize(fhandle, scheme)
    n_models = len(s['models']) or 1

    def per_model(n):
        return '{0:>6.1f}'.format(n / n_models)

    if not option:
        return [
            'No. models:\t{0}\n'.format(n_models),
            'No. chains:\t{0}\t({1}/model)\n'.format(len(s['chains']), per_model(len(s['chains']))),
            'No. residues:\t{0}\t({1}/model)\n'.format(
                len(s['residues']), per_model(len(s['residues']))),
            'No. atoms:\t{0}\t({1}/model)\n'.format(len(s['atoms']), per_model(len(s['atoms']))),
            'No. HETATM:\t{0}\n'.format(len(s['hetatm'])),
            'Multiple Occ.:\t{0}\n'.format(s['altloc']),
            'Res. Inserts:\t{0}\n'.format(s['icode']),
        ]

    lines = []
    if 'm' in option:
        lines.append('No. models:\t{0}\n'.format(n_models))
        if n_models > 1:
            ids = sorted(s['models'], key=lambda m: (_atoms.as_int(m) is None, _atoms.as_int(m), m))
            lines.append('\t->\t{0}\n'.format(','.join(ids)))
    if 'c' in option:
        lines.append('No. chains:\t{0}\t({1}/model)\n'.format(
            len(s['chains']), per_model(len(s['chains']))))
        lines.append('\t->\t{0}\n'.format(','.join(sorted(c[1] for c in s['chains']))))
    if 'r' in option:
        lines.append('No. residues:\t{0}\t({1}/model)\n'.format(
            len(s['residues']), per_model(len(s['residues']))))
        lines.append('\t->\t{0}\n'.format(','.join(sorted({r[1] for r in s['residues']}))))
    if 'a' in option:
        lines.append('No. atoms:\t{0}\t({1}/model)\n'.format(
            len(s['atoms']), per_model(len(s['atoms']))))
        lines.append('\t->\t{0}\n'.format(','.join(sorted({a[1] for a in s['atoms']}))))
    if 'h' in option:
        lines.append('No. HETATM:\t{0}\n'.format(len(s['hetatm'])))
        lines.append('\t->\t{0}\n'.format(','.join(sorted({h[3] for h in s['hetatm']}))))
    if 'o' in option:
        lines.append('Multiple Occ.:\t{0}\n'.format(s['altloc']))
    if 'i' in option:
        lines.append('Res. Inserts:\t{0}\n'.format(s['icode']))
    return lines


def main():
    fhandle, option, scheme = check_input(sys.argv[1:])
    _cli.output_text(_cli.lazy(run, fhandle, option, scheme), fhandle)


if __name__ == '__main__':
    main()
