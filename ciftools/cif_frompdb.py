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
Converts a PDB file to mmCIF format.

Converts the coordinates (ATOM/HETATM records, models) into _atom_site,
plus _entry and, from CRYST1, _cell and _symmetry. Author identifiers
(auth_*) are those of the PDB file. Unlike pdb_tocif, which copies them into
the label_* columns, label identifiers are assigned as in PDB entries:

    label_asym_id    one per polymer chain, one per ligand, and one for the
                     waters of each chain (A, B, ..., Z, AA, BA, ...)
    label_entity_id  polymer chains with the same sequence share an entity;
                     each ligand type and water are entities of their own
    label_seq_id     1, 2, ... along each polymer chain; '.' elsewhere

A TER record ends a polymer: the residues of a chain up to a TER record
are one polymer, and the next TER record of the same chain ends another.
So a file whose chains are blank or repeated, such as MD output with one
TER record per molecule, gets one label_asym_id per molecule. In chains
without TER records, the polymer runs from the first to the last ATOM
record (so modified residues inside a chain belong to it).

Atoms with a blank chain get their label_asym_id as their author chain, or
the next free name when another chain of the file already has it (this
replaces pdb_chainbows | pdb_tocif).

For chains without missing residues this gives the label identifiers of the
PDB entry; but label_seq_id only counts the residues present in the file
(SEQRES is not used), so it differs from the entry after missing residues,
and PDB entries order ligand instances and group glycans (branched
entities) by their own conventions. Missing elements are assigned as by
cif_element. Atom ids run from 1 over the whole file; models keep their
MODEL numbers. ANISOU, CONECT and the other PDB records are not converted.

Usage:
    cif_frompdb <pdb file>

Example:
    cif_frompdb 1abc.pdb > 1abc.cif

This program is part of the cif-tools suite, a port of pdb-tools to mmCIF.
"""

import os
import sys

from ciftools import _atoms, _cif, _cli

WATERS = frozenset(('HOH', 'WAT', 'DOD', 'H2O', 'TIP', 'TIP3', 'SOL'))
COLUMNS = ['group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id', 'label_comp_id',
           'label_asym_id', 'label_entity_id', 'label_seq_id', 'pdbx_PDB_ins_code', 'Cartn_x',
           'Cartn_y', 'Cartn_z', 'occupancy', 'B_iso_or_equiv', 'pdbx_formal_charge',
           'auth_seq_id', 'auth_comp_id', 'auth_asym_id', 'auth_atom_id', 'pdbx_PDB_model_num']


def check_input(args):
    """Return (file handle, structure name)."""
    fhandle, _, _ = _cli.parse_args(args, __doc__, takes_option=False)
    path = next((a for a in args if not a.startswith('-')), None)
    name = os.path.splitext(os.path.basename(path))[0] if path else 'structure'
    return fhandle, name


def asym_name(n):
    """0 -> A, 25 -> Z, 26 -> AA, 27 -> BA, ... (as in PDB entries)."""
    letters = ''
    n += 1
    while n:
        n, rest = divmod(n - 1, 26)
        letters += chr(ord('A') + rest)
    return letters


def _charge(text):
    text = text.strip()
    if len(text) == 2 and text[0].isdigit() and text[1] in '+-':
        return ('-' if text[1] == '-' else '') + text[0]
    return '?'


def _atom(line, model):
    """The fields of an ATOM/HETATM record."""
    line = line.rstrip('\r\n').ljust(80)
    return {
        'group': line[0:6].strip(), 'name': line[12:16].strip(), 'altloc': line[16].strip(),
        'resname': line[17:20].strip(), 'chain': line[21].strip(), 'resnum': line[22:26].strip(),
        'icode': line[26].strip(), 'x': line[30:38].strip(), 'y': line[38:46].strip(),
        'z': line[46:54].strip(), 'occ': line[54:60].strip(), 'b': line[60:66].strip(),
        'element': line[76:78].strip(), 'charge': _charge(line[78:80]), 'model': model,
    }


def read_pdb(lines):
    """(atoms, id code or None, CRYST1 line or None, closed segments).

    Each atom gets a 'segment': the number of TER records of its chain seen
    before it in its model, so that a TER starts a new molecule even when
    chains are blank. Closed segments are the (chain, segment) pairs of the
    first model that end with a TER record.
    """
    atoms, idcode, cryst1 = [], None, None
    model = '1'
    ters, closed = {}, set()  # chain -> TER records so far in this model
    for line in lines:
        record = line[:6]
        if record.startswith('HEADER') and line[62:66].strip():
            idcode = line[62:66].strip()
        elif record.startswith('CRYST1'):
            cryst1 = line
        elif record.startswith('MODEL'):
            model = line[10:14].strip() or line[5:].strip() or '1'
            ters = {}
        elif record.startswith(('ATOM', 'HETATM')):
            atom = _atom(line, model)
            atom['segment'] = ters.get(atom['chain'], 0)
            atoms.append(atom)
        elif record.startswith('TER') and atoms and atoms[-1]['model'] == model:
            chain = atoms[-1]['chain']
            if model == atoms[0]['model']:
                closed.add((chain, ters.get(chain, 0)))
            ters[chain] = ters.get(chain, 0) + 1
    return atoms, idcode, cryst1, closed


def _residue(atom):
    """A residue position; residues modelled as two residue types (altlocs)
    share it."""
    return atom['chain'], atom['segment'], atom['resnum'], atom['icode']


def assign_labels(atoms, closed):
    """Add label_asym_id, label_entity_id and label_seq_id to the atoms, and
    give atoms of a blank chain an author chain (their label_asym_id, unless
    another chain of the file has that name)."""
    first_model = atoms[0]['model']
    position, names = {}, {}  # residue -> index in file order, first residue name
    first_atom, last_atom = {}, {}  # segment -> first / last index of an ATOM residue
    ter_chains = {chain for chain, _ in closed}
    for atom in atoms:
        if atom['model'] != first_model:
            continue
        key = _residue(atom)
        i = position.setdefault(key, len(position))
        names.setdefault(key, atom['resname'])
        if atom['group'] == 'ATOM':
            segment = key[:2]
            first_atom.setdefault(segment, i)
            last_atom[segment] = i

    def is_polymer(i, key):
        segment, resname = key[:2], names[key]
        if resname in WATERS or segment not in first_atom:
            return False
        if key[0] in ter_chains:
            return segment in closed  # a molecule ended by TER
        return first_atom[segment] <= i <= last_atom[segment]

    labels, asym_of, sequences, counters = {}, {}, {}, {}
    for key, i in position.items():
        if is_polymer(i, key):
            instance = ('polymer',) + key[:2]
        elif names[key] in WATERS:
            instance = ('water', key[0])
        else:
            instance = ('ligand', key)
        if instance not in asym_of:
            asym_of[instance] = asym_name(len(asym_of))
        if instance[0] == 'polymer':
            counters[instance] = counters.get(instance, 0) + 1
            sequences.setdefault(instance, []).append(names[key])
            labels[key] = (instance, str(counters[instance]))
        else:
            labels[key] = (instance, '.')

    # Entities: polymers by sequence, ligands by residue name, water
    entity_ids, entity_of = {}, {}
    for instance in asym_of:
        if instance[0] == 'polymer':
            ident = ('polymer', tuple(sequences[instance]))
        elif instance[0] == 'water':
            ident = ('water',)
        else:
            ident = ('ligand', names[instance[1]])
        entity_of[instance] = entity_ids.setdefault(ident, str(len(entity_ids) + 1))
    # Later models without the TER records of the first: match residues by
    # chain and number when that is unambiguous.
    by_number = {}
    for key, label in labels.items():
        number = (key[0],) + key[2:]
        by_number[number] = label if number not in by_number else None
    used = {atom['chain'] for atom in atoms if atom['chain']}
    blank_chains = {}  # label_asym_id -> author chain, for atoms of a blank chain
    for atom in atoms:
        key = _residue(atom)
        label = labels.get(key)
        if label is None and atom['model'] != first_model:
            label = by_number.get((key[0],) + key[2:])
        if label is None:  # a residue only in a later model
            atom['asym'], atom['entity'], atom['seq'] = '?', '?', '.'
        else:
            instance, seq = label
            atom['asym'], atom['entity'], atom['seq'] = asym_of[instance], entity_of[instance], seq
        if not atom['chain'] and atom['asym'] != '?':
            asym = atom['asym']
            if asym not in blank_chains:
                name, n = asym, 0
                while name in used:
                    name, n = asym_name(n), n + 1
                used.add(name)
                blank_chains[asym] = name
            atom['chain'] = blank_chains[asym]


def run(fhandle, name='structure'):
    """Yield the mmCIF text."""
    atoms, idcode, cryst1, closed = read_pdb(fhandle)
    if not atoms:
        raise _cif.CIFError('No ATOM/HETATM records to convert')
    assign_labels(atoms, closed)
    block = idcode or name
    yield '# Converted to mmCIF by cif-tools (cif_frompdb)\n'
    yield 'data_{}\n#\n'.format(block)
    yield '_entry.id {}\n#\n'.format(_cif.quote(block))
    if cryst1:
        cell = [cryst1[6:15], cryst1[15:24], cryst1[24:33], cryst1[33:40], cryst1[40:47], cryst1[47:54]]
        names = ['length_a', 'length_b', 'length_c', 'angle_alpha', 'angle_beta', 'angle_gamma']
        for item, value in zip(names, cell):
            yield '_cell.{:<12} {}\n'.format(item, value.strip() or '?')
        z = cryst1[66:70].strip()
        yield '_cell.Z_PDB        {}\n#\n'.format(z or '?')
        yield '_symmetry.space_group_name_H-M {}\n#\n'.format(_cif.quote(cryst1[55:66].strip() or '?'))
    loop = _cif.Loop.new('_atom_site', COLUMNS)
    rows = []
    for n, a in enumerate(atoms, 1):
        element = a['element'].upper() or _atoms.guess_element(
            a['name'], a['resname'], a['group'] == 'HETATM') or '?'
        rows.append(_cif.Row(loop, [
            a['group'] or 'ATOM', str(n), element, a['name'], a['altloc'] or '.', a['resname'],
            a['asym'], a['entity'], a['seq'], a['icode'] or '?', a['x'], a['y'], a['z'],
            a['occ'] or '?', a['b'] or '?', a['charge'], a['resnum'], a['resname'],
            a['chain'] or '?', a['name'], a['model']]))
    yield _cif._format_table(rows, loop)
    yield '#\n'


def main():
    fhandle, name = check_input(sys.argv[1:])
    _cli.output_text(_cli.lazy(run, fhandle, name), fhandle)


if __name__ == '__main__':
    main()
