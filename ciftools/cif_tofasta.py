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
Extracts the residue sequence of an mmCIF file in FASTA format.

Canonical amino acids and nucleotides (RNA and DNA) are represented by their
one-letter code while all others are represented by 'X'. As in pdb_tofasta,
the sequence is read from the atoms of the first model, ligands and waters
(HETATM) included, so pipe through cif_delhetatm to drop them.

The -multi option splits the different chains into different records in the
FASTA file. The -seqres option reads the full deposited sequence from
_pdbx_poly_seq_scheme instead, including residues without coordinates. The
two can be combined: -multi,seqres.

Usage:
    cif_tofasta [-multi] [-seqres] <cif file>

Example:
    cif_tofasta 1abc.cif
    cif_tofasta -multi 1abc.cif
    cif_tofasta -multi,seqres 1abc.cif

This program is part of the cif-tools suite, a port of pdb-tools to mmCIF.
"""

import sys

from ciftools import _atoms, _cif, _cli

THREE_TO_ONE = {
    # 20 canonical amino acids
    'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
    'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
    # RNA and DNA (pdb_tofasta maps DNA to X)
    'A': 'A', 'C': 'C', 'G': 'G', 'U': 'U', 'T': 'T',
    'DA': 'A', 'DC': 'C', 'DG': 'G', 'DT': 'T', 'DU': 'U',
}
WIDTH = 60
VALID = {'multi', 'seqres'}


def check_input(args):
    """Return (file handle, multi, seqres, identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__)
    flags = set(_atoms.parse_ids(option)) if option is not None else set()
    if option is not None and (not flags or flags - VALID):
        _cli.fail(__doc__, 'You provided an invalid option: {!r}'.format(option))
    return fhandle, 'multi' in flags, 'seqres' in flags, scheme


def _from_atoms(fhandle, scheme):
    """(structure id, [(chain, [residue names])]) from the first model's atoms."""
    entry_id, block = None, None
    chains, seen = [], set()
    first_model = None
    for ev in _cif.parse(fhandle, (_atoms.ATOM_SITE, '_entry')):
        if type(ev) is _cif.Row:
            if ev.loop.category.lower() == '_entry':
                entry_id = entry_id or ev.get('id')
                continue
            values = ev.values
            if f.model is not None:
                if first_model is None:
                    first_model = values[f.model]
                elif values[f.model] != first_model:
                    continue
            chain = values[f.chain]
            if not chains or chains[-1][0] != chain:
                chains.append((chain, []))
            uid = residue(values)
            if uid not in seen:
                seen.add(uid)
                chains[-1][1].append(values[f.resname])
        elif type(ev) is _cif.Loop and ev.category.lower() == _atoms.ATOM_SITE:
            f = _atoms.Fields(ev, scheme)
            f.need('chain', 'resnum', 'resname')
            residue = _atoms.getter(f.resname, f.chain, f.resnum, f.icode)
        elif type(ev) is _cif.DataBlock and block is None:
            block = ev.name
    return entry_id or block or 'PDB', chains


def _from_seqres(fhandle, scheme):
    """(structure id, [(chain, [residue names])]) from _pdbx_poly_seq_scheme."""
    entry_id, block = None, None
    chains, found = [], False
    column = 'pdb_strand_id' if _cif.get_scheme(scheme) == 'auth' else 'asym_id'
    for ev in _cif.parse(fhandle, ('_pdbx_poly_seq_scheme', '_entry')):
        if type(ev) is _cif.Row:
            if ev.loop.category.lower() == '_entry':
                entry_id = entry_id or ev.get('id')
                continue
            found = True
            chain = ev[column]
            if not chains or chains[-1][0] != chain:
                chains.append((chain, []))
            chains[-1][1].append(ev['mon_id'])
        elif type(ev) is _cif.DataBlock and block is None:
            block = ev.name
    if not found:
        raise _cif.CIFError('-seqres needs _pdbx_poly_seq_scheme, which this file does not have')
    return entry_id or block or 'PDB', chains


def _wrap(sequence):
    return [sequence[i:i + WIDTH] + '\n' for i in range(0, len(sequence), WIDTH)]


def run(fhandle, multi=False, seqres=False, scheme=None):
    """Return the FASTA lines."""
    reader = _from_seqres if seqres else _from_atoms
    struct_id, chains = reader(fhandle, scheme)
    if multi:
        lines = []
        for chain, names in chains:
            lines.append('>{}|{}\n'.format(struct_id, chain))
            lines.extend(_wrap(''.join(THREE_TO_ONE.get(n, 'X') for n in names)))
        return lines
    ids = sorted({chain for chain, _ in chains})
    label = ''.join(ids) if all(len(c) == 1 for c in ids) else ','.join(ids)
    sequence = ''.join(THREE_TO_ONE.get(n, 'X') for _, names in chains for n in names)
    return ['>{}|{}\n'.format(struct_id, label)] + _wrap(sequence)


def main():
    fhandle, multi, seqres, scheme = check_input(sys.argv[1:])
    _cli.output_text(_cli.lazy(run, fhandle, multi, seqres, scheme), fhandle)


if __name__ == '__main__':
    main()
