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
Adds the core tables that an mmCIF file lacks, built from its atoms.

Files written by other programs often hold only _atom_site, sometimes
without label identifiers, while structure prediction programs and the
other mmcif-tools read the molecules from the core tables (mmcif_validate
warns about the ones missing). This tool adds, from the atoms, each of these
tables the file does not have; the tables it has are left as they are:

    _entity                  the molecules: polymer, non-polymer, water
    _entity_poly             the polymers: type, one-letter sequences
    _entity_poly_seq         their sequences
    _struct_asym             the chains (label_asym_id) of each molecule
    _pdbx_entity_nonpoly     the ligands and water
    _chem_comp               the residue types
    _pdbx_poly_seq_scheme    the residues of each chain along its sequence
    _pdbx_nonpoly_scheme     the chains and numbers of the ligands and waters

The sequences are those of the residues with atoms: residues missing from
the model (see mmcif_gap) cannot be known and are not added. Copies of a
molecule are one entity when their residue types agree at every position
they share, so chains modelled over different ranges share an entity whose
sequence has the residues of all of them. Residues are classified by name:
the standard amino acids and nucleotides, water (HOH, WAT, DOD, ...), and
the others by the polymer they sit in (a modified residue in a protein is
'L-peptide linking') or as ligands; glycans are ligands, one per residue.

When _atom_site has no label_entity_id, the label identifiers are assigned
first, as mmcif_frompdb does: one label_asym_id per polymer chain, per
ligand and for the waters of each chain, and label_seq_id along each
entity's sequence. A polymer runs from the first to the last ATOM record of
its author chain, with the HETATM residues bonded to its ends (a peptide or
phosphodiester bond: MSE, ACE, NH2, ...); in a chain without a name, a new
molecule starts where the author numbers start again. Copies are matched by
author number (they must share half of the shorter one), and the sequence
keeps the order of the chains, also where the author numbers are not in
order (an inserted domain). A file with label_entity_id keeps its label
identifiers. Elements are not assigned: see mmcif_element. Only the first
data block with atoms changes.

Usage:
    mmcif_rebuildmeta <cif file>

Example:
    mmcif_rebuildmeta model.cif > model_full.cif
    mmcif_keepcoord 1abc.cif | mmcif_rebuildmeta | mmcif_validate

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cif, _cli
from mmciftools.mmcif_delmeta import SCHEME_COLUMNS, _Present
from mmciftools.mmcif_frompdb import WATERS, asym_name

NULLS = ('.', '?')
AMINO_ACIDS = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H',
    'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V'}
DNA = {'DA': 'A', 'DC': 'C', 'DG': 'G', 'DT': 'T'}
RNA = {'A': 'A', 'C': 'C', 'G': 'G', 'U': 'U'}
# A polymer's type, and the type of the components it is made of (also its non-standard ones).
LINKING = {'polypeptide(L)': 'L-peptide linking', 'polydeoxyribonucleotide': 'DNA linking',
           'polyribonucleotide': 'RNA linking'}
HYBRID = 'polydeoxyribonucleotide/polyribonucleotide hybrid'
BEFORE_ATOMS = ('_entity', '_entity_poly', '_entity_poly_seq', '_struct_asym', '_pdbx_entity_nonpoly', '_chem_comp')
AFTER_ATOMS = ('_pdbx_poly_seq_scheme', '_pdbx_nonpoly_scheme')
CORE = BEFORE_ATOMS + AFTER_ATOMS
COLUMNS = {
    '_entity': ('id', 'type', 'pdbx_description', 'pdbx_number_of_molecules'),
    '_entity_poly': ('entity_id', 'type', 'nstd_linkage', 'nstd_monomer', 'pdbx_seq_one_letter_code',
                     'pdbx_seq_one_letter_code_can', 'pdbx_strand_id'),
    '_entity_poly_seq': ('entity_id', 'num', 'mon_id', 'hetero'),
    '_struct_asym': ('id', 'pdbx_blank_PDB_chainid_flag', 'pdbx_modified', 'entity_id', 'details'),
    '_pdbx_entity_nonpoly': ('entity_id', 'name', 'comp_id'),
    '_chem_comp': ('id', 'type', 'mon_nstd_flag', 'name', 'formula', 'formula_weight'),
    '_pdbx_poly_seq_scheme': SCHEME_COLUMNS,
    '_pdbx_nonpoly_scheme': ('asym_id', 'entity_id', 'mon_id', 'ndb_seq_num', 'pdb_seq_num', 'auth_seq_num',
                             'pdb_mon_id', 'auth_mon_id', 'pdb_strand_id', 'pdb_ins_code'),
}
LABELS = ('label_asym_id', 'label_entity_id', 'label_seq_id', 'label_comp_id', 'label_atom_id')
# The bonds that link residues in a polymer: (atom of a residue, atom of the next one, longest distance in A).
LINKS = (('C', 'N', 1.75), ("O3'", 'P', 1.8))
LINKING_ATOMS = frozenset(name for link in LINKS for name in link[:2])
# _atom_site columns read: the residue (author and label ids) and the record.
READ = ('auth_asym_id', 'auth_seq_id', 'pdbx_PDB_ins_code', 'auth_comp_id', 'label_comp_id', 'label_asym_id',
        'label_entity_id', 'label_seq_id', 'pdbx_PDB_model_num', 'group_PDB')


def check_input(args):
    """Return the file handle."""
    fhandle, _, _ = _cli.parse_args(args, __doc__, takes_option=False)
    return fhandle


def _value(value):
    return None if value is None or value in NULLS else value


class _Residues:
    """Reads the residues of _atom_site rows in file order: their author
    ids, their name, their label ids, and the molecule they belong to in a
    chain without a name (a new one where the author numbers start again)."""

    def __init__(self, loop):
        self.get = _atoms.getter(*(loop.index(c) if c in loop else None for c in READ))
        self.atom_name = loop.index('auth_atom_id') if 'auth_atom_id' in loop else (
            loop.index('label_atom_id') if 'label_atom_id' in loop else None)
        xyz = [loop.index(c) if c in loop else None for c in ('Cartn_x', 'Cartn_y', 'Cartn_z')]
        self.xyz = xyz if None not in xyz and self.atom_name is not None else None
        self.model = self.last = None

    def linking_atom(self, values):
        """(name, coordinates) of an atom that links residues in a polymer
        (N and C of a peptide bond, P and O3' of a phosphodiester), or None."""
        if self.xyz is None:
            return None
        name = values[self.atom_name]
        if name not in LINKING_ATOMS:
            return None
        try:
            return name, tuple(float(values[p]) for p in self.xyz)
        except ValueError:
            return None

    def residue(self, values):
        """(key, label ids, is an ATOM record, model) of a row; the key
        (author chain, molecule, number, insertion code, name, label chain,
        label number) finds the residue in every model."""
        (chain, number, icode, auth_comp, label_comp, label_asym, label_entity, label_seq, model,
         group) = (_value(v) for v in self.get(values))
        if model != self.model:
            self.model, self.numbers, self.molecules = model, {}, {}
        comp = label_comp or auth_comp or '?'
        if chain is None:
            chain = label_asym
        if chain is None:  # a chain without a name: a new molecule where the numbers start again
            n = _atoms.as_int(number)
            if n is not None:
                if n < self.numbers.get(chain, n):
                    self.molecules[chain] = self.molecules.get(chain, 0) + 1
                self.numbers[chain] = n
        key = (chain, self.molecules.get(chain, 0) if chain is None else 0, number, icode, comp, label_asym,
               label_seq)
        return key, (label_asym, label_entity, label_seq), group == 'ATOM', model


class _Block:
    """What the first data block with atoms holds: its residues and the core
    tables it has."""

    def __init__(self, index):
        self.index = index
        self.tables = set()   # core tables present
        self.kinds = {}       # _entity: id -> type
        self.struct_asym = None  # [(id, entity_id)]
        self.sequences = {}   # _entity_poly_seq: entity -> [(num, mon_id, hetero)]
        self.poly_types = {}  # _entity_poly: entity -> type
        self.residues = {}    # key -> [label ids, ATOM record, order, author ids] (first model)
        self.later = False    # residues in models after the first
        self.columns = None   # the _atom_site columns
        self.has_labels = False

    def read(self, row):
        category, get = row.loop.category.lower(), row.get
        if category == '_entity':
            self.kinds[get('id')] = (get('type') or '').lower()
        elif category == '_struct_asym':
            self.struct_asym.append((get('id'), get('entity_id')))
        elif category == '_entity_poly_seq':
            hetero = _value(get('hetero')) or 'n'
            self.sequences.setdefault(get('entity_id'), []).append((get('num'), get('mon_id'), hetero))
        elif category == '_entity_poly':
            self.poly_types[get('entity_id')] = get('type')


def _read(source):
    """The _Block of the first data block with atoms, or None."""
    block, residues, first_model = None, None, None
    for ev in _cif.parse(source.lines(), (_atoms.ATOM_SITE,) + CORE):
        etype = type(ev)
        if etype is _cif.DataBlock:
            if block is not None and block.columns is not None:
                break  # the first block with atoms is read
            block, first_model = _Block(0 if block is None else block.index + 1), None
            continue
        if etype is not _cif.Loop and etype is not _cif.Row:
            continue  # text between categories (a comment before data_ is not a block)
        if block is None:
            block = _Block(0)  # categories before any data_ line
        if etype is _cif.Row:
            if ev.loop.category.lower() != _atoms.ATOM_SITE:
                block.read(ev)
                continue
            key, labels, atom, model = residues.residue(ev.values)
            if first_model is None:
                first_model = model
            if model != first_model:
                block.later = True
                continue
            here = block.residues.get(key)
            if here is None:
                here = block.residues[key] = [labels, atom, len(block.residues), {}]
            elif atom:
                here[1] = True
            link = residues.linking_atom(ev.values)
            if link is not None:
                here[3].setdefault(*link)
        elif etype is _cif.Loop:
            category = ev.category.lower()
            if category != _atoms.ATOM_SITE:
                block.tables.add(category)
                if category == '_struct_asym':
                    block.struct_asym = []
                continue
            residues = _Residues(ev)
            if block.columns is None:
                block.columns = list(ev.columns)
                block.has_labels = 'label_entity_id' in ev
                if block.has_labels and not ('label_asym_id' in ev and 'label_seq_id' in ev):
                    raise _cif.CIFError('_atom_site has label_entity_id but no label_asym_id or label_seq_id: '
                                        'add them, or remove label_entity_id to have them all assigned')
    return block if block is not None and block.columns is not None else None


def _position(key, order):
    """Where a residue goes along a sequence: by author number and insertion
    code (by file order without a number)."""
    number = _atoms.as_int(key[2])
    return (0, number, key[3] or '') if number is not None else (1, order, key[3] or '')


def _merge(union, ranks, mine):
    """Add the residues of a copy to its entity's sequence: each one new to
    it goes right after the residue before it in the copy (so the sequence
    keeps the order of the chains, also where author numbers are not in
    order, as in fusion proteins)."""
    before, added = None, len(union)
    for p, names in mine.items():
        here = union.get(p)
        if here is None:
            added += 1
            here = union[p] = []
            ranks[p] = (ranks[before] if before is not None else (0,)) + (added,)
        here += [n for n in names if n not in here]
        before = p


def _linked(before, after):
    """Whether two residues (entries of _Block.residues) are bonded as in a
    polymer: the C of the first to the N of the second, or its O3' to the P."""
    for one, other, longest in LINKS:
        a, b = before[3].get(one), after[3].get(other)
        if a is not None and b is not None and sum((x - y) ** 2 for x, y in zip(a, b)) <= longest ** 2:
            return True
    return False


def _assign(block):
    """{residue key: (label_asym_id, label_entity_id, label_seq_id)} for a
    file without label ids: polymers, ligands and the waters of each chain
    get a label_asym_id each (in file order), polymers that are copies of
    one molecule an entity, and their residues their place along it."""
    residues = sorted(block.residues.items(), key=lambda item: item[1][2])
    first_atom, last_atom, members = {}, {}, {}  # molecule -> first / last order of an ATOM record; its residues
    for key, (_, atom, order, _) in residues:
        molecule = key[:2]
        members.setdefault(molecule, []).append((key, order))
        if atom and key[4] not in WATERS:
            first_atom.setdefault(molecule, order)
            last_atom[molecule] = order
    for molecule, first in first_atom.items():  # HETATM residues bonded to the ends (MSE, ACE, NH2, ...)
        chain = members[molecule]
        i = next(n for n, (_, order) in enumerate(chain) if order == first)
        while i > 0 and _linked(block.residues[chain[i - 1][0]], block.residues[chain[i][0]]):
            i -= 1
        first_atom[molecule] = chain[i][1]
        i = next(n for n, (_, order) in enumerate(chain) if order == last_atom[molecule])
        while i + 1 < len(chain) and _linked(block.residues[chain[i][0]], block.residues[chain[i + 1][0]]):
            i += 1
        last_atom[molecule] = chain[i][1]
    instance_of, maps = {}, {}  # key -> its instance; polymer instance -> {position: [names]}
    for key, (_, atom, order, _) in residues:
        molecule, name = key[:2], key[4]
        if name in WATERS:
            instance = ('water', key[0], key[1])
        elif molecule in first_atom and first_atom[molecule] <= order <= last_atom[molecule]:
            instance = ('polymer',) + molecule
            names = maps.setdefault(instance, {}).setdefault(_position(key, order), [])
            if name not in names:
                names.append(name)
        else:
            instance = ('ligand', key)
        instance_of[key] = instance
    asyms = {}
    for key, _ in residues:
        asyms.setdefault(instance_of[key], asym_name(len(asyms)))
    # entities: copies of a polymer whose names agree at every position they share (half of the
    # shorter one at least), each kind of ligand, and water
    polymers, entity_of, ids = [], {}, {}  # [(entity id, union of the copies' residues, their ranks)]
    for instance in asyms:
        if instance[0] == 'polymer':
            mine = maps[instance]
            for entity, union, ranks in polymers:
                shared = [p for p in mine if p in union]
                if shared and 2 * len(shared) >= min(len(mine), len(union)) and \
                        all(set(mine[p]) & set(union[p]) for p in shared):
                    _merge(union, ranks, mine)
                    entity_of[instance] = entity
                    break
            else:
                entity = ids.setdefault(('polymer', len(polymers)), str(len(ids) + 1))
                union = {p: list(names) for p, names in mine.items()}
                polymers.append((entity, union, {p: (i,) for i, p in enumerate(mine, 1)}))
                entity_of[instance] = entity
        elif instance[0] == 'water':
            entity_of[instance] = ids.setdefault(('water',), str(len(ids) + 1))
        else:
            entity_of[instance] = ids.setdefault(('ligand', instance[1][4]), str(len(ids) + 1))
    numbers = {}  # (entity, position) -> label_seq_id
    for entity, union, ranks in polymers:
        for n, p in enumerate(sorted(union, key=ranks.get), 1):
            numbers[(entity, p)] = str(n)
    labels = {}
    for key, (_, atom, order, _) in residues:
        instance = instance_of[key]
        entity = entity_of[instance]
        seq = numbers[(entity, _position(key, order))] if instance[0] == 'polymer' else '.'
        labels[key] = (asyms[instance], entity, seq)
    return labels


def _kind(names):
    """The type of a polymer made of these residue names."""
    amino = sum(n in AMINO_ACIDS for n in names)
    dna, rna = sum(n in DNA for n in names), sum(n in RNA for n in names)
    if amino and amino >= dna + rna:
        return 'polypeptide(L)'
    if dna and rna:
        return HYBRID
    if dna:
        return 'polydeoxyribonucleotide'
    if rna:
        return 'polyribonucleotide'
    return 'other'


def _one_letter(names, kind):
    """(pdbx_seq_one_letter_code, its canonical form), in lines of up to 80
    characters broken between residues."""
    codes, canonical = [], []
    for name in names:
        if name in AMINO_ACIDS and kind == 'polypeptide(L)':
            codes.append(AMINO_ACIDS[name])
        elif name in RNA and kind in ('polyribonucleotide', HYBRID):
            codes.append(RNA[name])
        else:
            codes.append('({})'.format(name))
        canonical.append(AMINO_ACIDS.get(name) if kind == 'polypeptide(L)' else DNA.get(name, RNA.get(name)))
    return _lines(codes), _lines([c or 'X' for c in canonical])


def _lines(codes):
    lines, line = [], ''
    for code in codes:
        if line and len(line) + len(code) > 80:
            lines.append(line)
            line = ''
        line += code
    return '\n'.join(lines + [line]) if line else '?'


def _tables(block, labels):
    """{category: text} of the core tables the block lacks."""
    residues = sorted(block.residues.items(), key=lambda item: item[1][2])
    chains, strands, sequences, ligands, waters = {}, {}, {}, [], {}  # in file order
    first_model = []
    for key, (ids, atom, order, _) in residues:
        asym, entity, seq = (_value(v) for v in (labels[key] if labels is not None else ids))
        if asym is None or entity is None:
            continue
        chain, number, icode, name = key[0], key[2], key[3], key[4]
        chains.setdefault(asym, (entity, chain))
        if chain is not None and chain not in strands.setdefault(entity, []):
            strands[entity].append(chain)
        if seq is not None:
            names = sequences.setdefault(entity, {}).setdefault(seq, [])
            if name not in names:
                names.append(name)
        else:
            ligands.append((asym, entity, name, number, chain, icode))
            if name in WATERS:
                waters[entity] = waters.get(entity, 0) + 1
        first_model.append((asym, name, seq, chain, name, number, icode, entity))

    kinds = dict(block.kinds)
    for asym, (entity, _) in chains.items():
        if entity not in kinds:
            if entity in sequences:
                kinds[entity] = 'polymer'
            elif all(lig[2] in WATERS for lig in ligands if lig[1] == entity):
                kinds[entity] = 'water'
            else:
                kinds[entity] = 'non-polymer'
    if block.sequences:  # the file's own sequences
        ordered = {e: [(num, mon, hetero) for num, mon, hetero in s] for e, s in block.sequences.items()}
    else:
        ordered = {}
        for entity, by_number in sequences.items():
            if kinds.get(entity) != 'polymer':
                continue
            rows = ordered.setdefault(entity, [])
            for num in sorted(by_number, key=lambda n: (_atoms.as_int(n) is None, _atoms.as_int(n) or 0, n)):
                names = by_number[num]
                for name in names:
                    rows.append((num, name, 'y' if len(names) > 1 else 'n'))
    poly_types = {}
    for entity, rows in ordered.items():
        poly_types[entity] = block.poly_types.get(entity) or _kind([mon for _, mon, _ in rows])
    entities = list(dict.fromkeys(entity for entity, _ in chains.values()))
    chains_of = {}
    for asym, (entity, _) in chains.items():
        chains_of.setdefault(entity, []).append(asym)

    out = {}
    if '_entity' not in block.tables:
        rows = []
        for entity in entities:
            count = waters.get(entity, 0) if kinds[entity] == 'water' else len(chains_of[entity])
            rows.append([entity, kinds[entity], 'water' if kinds[entity] == 'water' else '?', str(count)])
        out['_entity'] = rows
    polymers = [e for e in entities if e in ordered]
    if '_entity_poly' not in block.tables and polymers:
        rows = []
        for entity in polymers:
            kind = poly_types[entity]
            names = [mon for i, (num, mon, _) in enumerate(ordered[entity])
                     if i == 0 or ordered[entity][i - 1][0] != num]  # the first type at each position
            code, canonical = _one_letter(names, kind)
            standard = AMINO_ACIDS if kind == 'polypeptide(L)' else dict(DNA, **RNA)
            rows.append([entity, kind, 'no', 'yes' if any(n not in standard for n in names) else 'no', code,
                         canonical, ','.join(strands.get(entity, ())) or '?'])
        out['_entity_poly'] = rows
    if '_entity_poly_seq' not in block.tables and polymers:
        out['_entity_poly_seq'] = [[entity, num, mon, hetero] for entity in polymers
                                   for num, mon, hetero in ordered[entity]]
    if '_struct_asym' not in block.tables:
        out['_struct_asym'] = [[asym, 'N' if chain is not None else 'Y', 'N', entity, '?']
                               for asym, (entity, chain) in chains.items()]
    others = [e for e in entities if kinds[e] in ('non-polymer', 'water')]
    if '_pdbx_entity_nonpoly' not in block.tables and others:
        names = {}
        for asym, entity, name, number, chain, icode in ligands:
            names.setdefault(entity, name)
        out['_pdbx_entity_nonpoly'] = [[e, 'water' if kinds[e] == 'water' else '?', names[e]]
                                       for e in others if e in names]
    if '_chem_comp' not in block.tables:
        linking = {}
        for entity, rows in ordered.items():
            for _, mon, _ in rows:
                linking.setdefault(mon, LINKING.get(poly_types[entity], 'other'))
        components = {}
        for key, _ in residues:
            name = key[4]
            if name in components:
                continue
            if name in AMINO_ACIDS:
                kind = 'peptide linking' if name == 'GLY' else 'L-peptide linking'
            elif name in DNA:
                kind = 'DNA linking'
            elif name in RNA:
                kind = 'RNA linking'
            else:
                kind = linking.get(name, 'non-polymer')
            standard = name in AMINO_ACIDS or name in DNA or name in RNA
            components[name] = [name, kind, 'y' if standard else '.', 'WATER' if name in WATERS else '?', '?', '?']
        out['_chem_comp'] = [components[n] for n in sorted(components)]
    if '_pdbx_poly_seq_scheme' not in block.tables and ordered:
        info = _Present()
        info.first_model = first_model
        info.sequences = ordered
        info.chains = {asym: entity for asym, (entity, _) in chains.items()}
        info.asyms = set(info.chains)
        info.struct_asym = block.struct_asym
        rows = info.scheme()
        strand = SCHEME_COLUMNS.index('pdb_strand_id')
        for row in rows:
            if row[0] in chains and chains[row[0]][1] is None:
                row[strand] = '?'  # as the atoms: a chain without a name (not its label_asym_id)
        if rows:
            out['_pdbx_poly_seq_scheme'] = rows
    if '_pdbx_nonpoly_scheme' not in block.tables and ligands:
        rows, counters = [], {}
        for asym, entity, name, number, chain, icode in ligands:
            counters[asym] = counters.get(asym, 0) + 1
            rows.append([asym, entity, name, str(counters[asym]), number or '?', number or '?', name, name,
                         chain or '?', icode or '.'])
        out['_pdbx_nonpoly_scheme'] = rows
    return {category: _format(category, rows) for category, rows in out.items()}


def _format(category, rows):
    """The text of a new table: key-value pairs for one row, else a loop."""
    loop = _cif.Loop.new(category, COLUMNS[category])
    rows = [_cif.Row(loop, [v if v is not None else '?' for v in row]) for row in rows]
    if len(rows) == 1:
        return _cif._format_pairs(rows[0], loop) + '#\n'
    return _cif._format_table(rows, loop) + '#\n'


def run(fhandle, scheme=None):
    """Yield the events of the file, with the core tables it lacks added
    (and label ids assigned to its atoms when it has none)."""
    source = _atoms.as_source(fhandle)
    block = _read(source)
    if block is None:  # nothing to describe
        yield from _cif.parse(source.lines(), _atoms.ATOM_SITE)
        return
    labels = None if block.has_labels else _assign(block)
    tables = _tables(block, labels)
    blocks, loop, residues, done = -1, None, None, False
    for ev in _cif.parse(source.lines(), _atoms.ATOM_SITE):
        etype = type(ev)
        if etype is _cif.DataBlock:
            blocks += 1
        elif etype is _cif.Loop and not done:
            blocks = max(blocks, 0)
            if blocks == block.index:
                yield from (tables[c] for c in BEFORE_ATOMS if c in tables)
                if labels is not None:
                    added = [c for c in LABELS if c not in ev]
                    loop = ev.with_columns(ev.columns + added)
                    residues = _Residues(ev)
                    yield loop
                    continue
        elif etype is _cif.Row and loop is not None:
            key = residues.residue(ev.values)[0]
            asym, entity, seq = labels.get(key, ('?', '?', '.'))
            name = ev.values[residues.atom_name] if residues.atom_name is not None else '?'
            yield ev.moved(loop, {'label_asym_id': asym, 'label_entity_id': entity, 'label_seq_id': seq,
                                  'label_comp_id': key[4], 'label_atom_id': name})
            continue
        elif etype is _cif.LoopEnd and not done and blocks == block.index:
            done, loop = True, None
            yield ev
            yield from (tables[c] for c in AFTER_ATOMS if c in tables)
            continue
        yield ev


def main():
    fhandle = check_input(sys.argv[1:])
    _cli.output(_cli.lazy(run, fhandle), fhandle)


if __name__ == '__main__':
    main()
