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
Deletes the metadata that edits make stale, keeping the coordinates and the
metadata that stays true.

PDB entries carry dozens of categories that refer to residues by chain and
number: validation reports, secondary structure, sequence schemes, sites,
assemblies, ... The other tools do not update them, so after an edit they no
longer match the atoms. This tool keeps:

    _atom_site                 the atoms
    _atom_site_anisotrop       their anisotropic B-factors
    _struct_conn               the links (disulfides, metal sites, ...)
    _entry                     the entry id
    _cell, _symmetry           the unit cell and space group
    _exptl, _refine, _reflns,  the method and resolution
    _em_3d_reconstruction
    _pdbx_audit_revision_history  the release and revision dates
    _entity, _entity_poly,     the molecules, their sequences, and the
    _entity_poly_seq,          ligands
    _pdbx_entity_nonpoly
    _struct_asym               the chains (label_asym_id) of each molecule
    _pdbx_poly_seq_scheme      the residues of each chain along its sequence
    _pdbx_struct_mod_residue   the modified residues
    _chem_comp, _chem_comp_atom,  the chemical components, their atoms and
    _chem_comp_bond, _atom_type   bonds, and the elements

and deletes every other category. These are what structure prediction
programs (AlphaFold, ColabFold, OpenFold3, Boltz, Chai-1, Protenix, BioEmu)
and MD preparation programs (OpenMM, PDBFixer) read from a file.

The kept tables are brought in step with the atoms: rows of chains and
entities with no atoms left are removed, and so are the atom and bond lists
of components whose atoms were renamed. Modified residues are looked up in
the atoms by their label identifiers and get the current author numbers
(rows of residues no longer there are removed).
_entity.pdbx_number_of_molecules and _entity_poly.pdbx_strand_id are
recomputed, and _pdbx_poly_seq_scheme is rebuilt from the sequences and the
atoms (it is added when missing; residues without atoms get no author
number). So run it last in a pipeline, after the edits. Tables that cannot
be checked (no atoms, no label identifiers) are kept as they are.

With -all, only the atoms, anisotropic B-factors and links are kept
(mmcif_keepcoord also drops the anisotropic B-factors, as pdb_keepcoord
drops ANISOU records).

Usage:
    mmcif_delmeta [-all] <cif file>

Example:
    mmcif_delmeta 1abc.cif > 1abc_clean.cif
    mmcif_selchain -A 1abc.cif | mmcif_delentity -water | mmcif_delmeta > 1abc_A.cif
    mmcif_delmeta -all 1abc.cif > atoms_only.cif

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cif, _cli

COORDINATES = (_atoms.ATOM_SITE, _atoms.ANISOTROP, _atoms.STRUCT_CONN)
KEPT = ('_entry', '_cell', '_symmetry', '_exptl', '_refine', '_reflns', '_em_3d_reconstruction',
        '_pdbx_audit_revision_history', '_entity', '_entity_poly', '_entity_poly_seq', '_pdbx_entity_nonpoly',
        '_struct_asym', '_pdbx_struct_mod_residue', '_chem_comp', '_chem_comp_atom', '_chem_comp_bond',
        '_atom_type')
SCHEME = '_pdbx_poly_seq_scheme'  # rebuilt, not kept
MOD_RESIDUE = '_pdbx_struct_mod_residue'
COMPONENTS = ('_chem_comp_atom', '_chem_comp_bond')
BY_ENTITY = {'_entity': 'id', '_entity_poly': 'entity_id', '_entity_poly_seq': 'entity_id',
             '_pdbx_entity_nonpoly': 'entity_id'}
# The columns that identify a residue: (in _atom_site, in _pdbx_struct_mod_residue).
# Residues are tuples of these values, followed by the label_entity_id.
RESIDUE = (('label_asym_id', 'label_asym_id'), ('label_comp_id', 'label_comp_id'),
           ('label_seq_id', 'label_seq_id'), ('auth_asym_id', 'auth_asym_id'),
           ('auth_comp_id', 'auth_comp_id'), ('auth_seq_id', 'auth_seq_id'),
           ('pdbx_PDB_ins_code', 'PDB_ins_code'))
NAMES = ('label_comp_id', 'label_atom_id', 'auth_comp_id', 'auth_atom_id', 'type_symbol')
HYDROGENS = ('H', 'D')  # left out of the comparison with _chem_comp_atom (terminal and NMR hydrogens)
SCHEME_COLUMNS = ('asym_id', 'entity_id', 'seq_id', 'mon_id', 'ndb_seq_num', 'pdb_seq_num', 'auth_seq_num',
                  'pdb_mon_id', 'auth_mon_id', 'pdb_strand_id', 'pdb_ins_code', 'hetero')
NULLS = ('.', '?')


def check_input(args):
    """Return (file handle, delete all metadata, identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__)
    if option not in (None, 'all'):
        _cli.fail(__doc__, 'Invalid option: {!r} (the only option is -all)'.format('-' + option))
    return fhandle, option == 'all', scheme


def _plain(value):
    return '' if value is None or value in NULLS else value


class _Present:
    """What one data block holds: the chains, entities, residues and atom
    names of _atom_site, and the tables that are checked against them."""

    def __init__(self):
        self.asyms = None         # label_asym_id values (None: no such column)
        self.entities = None      # entity id -> set of label_asym_id
        self.molecules = {}       # entity id -> number of residues in the first model (for waters)
        self.strands = {}         # entity id -> author chains, in order
        self.chains = {}          # label_asym_id -> entity id, in order
        self.first_model = []     # residues of the first model, in order
        self.columns = None       # which RESIDUE columns _atom_site has
        self.names = set()        # NAMES values of the atoms
        self.struct_asym = None   # [(id, entity_id)]
        self.sequences = {}       # entity id -> [(num, mon_id, hetero)]
        self.component_atoms = {}  # comp_id -> atom names (_chem_comp_atom)
        self.bonded_atoms = {}    # comp_id -> atom names (_chem_comp_bond)
        self.had_scheme = False
        self.renamed = set()      # components whose atoms are not all in their definition
        self._matches = {}        # _pdbx_struct_mod_residue loop -> how to find its residues

    def read(self, row):
        """Note a row of a table other than _atom_site."""
        category = row.loop.category.lower()
        if category == '_struct_asym':
            self.struct_asym.append((row.get('id'), row.get('entity_id')))
        elif category == '_entity_poly_seq':
            hetero = row.get('hetero')
            self.sequences.setdefault(row.get('entity_id'), []).append(
                (row.get('num'), row.get('mon_id'), 'n' if hetero is None or hetero in NULLS else hetero))
        elif category == '_chem_comp_atom':
            self.component_atoms.setdefault(row.get('comp_id'), set()).add(row.get('atom_id'))
        elif category == '_chem_comp_bond':
            self.bonded_atoms.setdefault(row.get('comp_id'), set()).update(
                (row.get('atom_id_1'), row.get('atom_id_2')))

    def finish(self):
        for residue in set(self.first_model):
            self.molecules[residue[-1]] = self.molecules.get(residue[-1], 0) + 1
        used = {}
        for label_comp, label_atom, auth_comp, auth_atom, element in self.names:
            if element is not None and element.upper() in HYDROGENS:
                continue
            for comp, atom in ((label_comp, label_atom), (auth_comp, auth_atom)):
                if comp is not None and atom is not None:
                    used.setdefault(comp, set()).add(atom)
        for comp, atoms in used.items():
            defined = self.component_atoms.get(comp, self.bonded_atoms.get(comp))
            if defined is not None and not atoms <= defined:
                self.renamed.add(comp)

    def mod_residue(self, row):
        """Whether the residue of a _pdbx_struct_mod_residue row is in the
        atoms. It is found by its label identifiers when both tables have
        them (and its author identifiers are then updated from the atoms),
        else by the identifiers they share."""
        if self.columns is None:
            return True  # nothing to check against
        loop = row.loop
        if loop not in self._matches:
            by_label = None
            if self.columns[0] and self.columns[2] and 'label_asym_id' in loop and 'label_seq_id' in loop:
                by_label = {}
                for r in self.first_model:
                    by_label.setdefault((r[0], r[2]), {}).setdefault(r[1], r)
            pairs = [(loop.index(theirs), i) for i, (ours, theirs) in enumerate(RESIDUE)
                     if self.columns[i] and theirs in loop]
            residues = {tuple(_plain(r[i]) for _, i in pairs) for r in self.first_model}
            self._matches[loop] = (by_label, pairs, residues)
        by_label, pairs, residues = self._matches[loop]
        if by_label is None:
            return not pairs or tuple(_plain(row.values[p]) for p, _ in pairs) in residues
        here = by_label.get((row['label_asym_id'], row['label_seq_id']), {})
        comp = row.get('label_comp_id') if self.columns[1] else None
        residue = here.get(comp) if comp is not None and comp not in NULLS else next(iter(here.values()), None)
        if residue is None:
            return False
        for i in range(3, len(RESIDUE)):  # the author identifiers
            column = RESIDUE[i][1]
            if self.columns[i] and column in loop and _plain(row[column]) != _plain(residue[i]):
                row[column] = residue[i] if _plain(residue[i]) else '?'
        return True

    def can_rebuild(self):
        """Whether there are sequences, and atoms with label_asym_id and label_seq_id."""
        return bool(self.sequences) and self.columns is not None and self.columns[0] and self.columns[2]

    def scheme(self):
        """The rows of _pdbx_poly_seq_scheme, rebuilt from the sequences, the
        chains and the atoms of the first model."""
        first, others, strand = {}, {}, {}  # (chain, seq_id) -> its first residue, and the others
        for residue in self.first_model:
            asym, comp, seq, chain = residue[:4]
            here = first.setdefault((asym, seq), residue)
            if here[1] != comp:  # two residue types at one position
                others.setdefault((asym, seq, comp), residue)
            strand.setdefault(asym, chain)
        if self.struct_asym is not None:
            chains = [(asym, entity) for asym, entity in self.struct_asym if asym in self.asyms]
        else:
            chains = list(self.chains.items())
        rows = []
        for asym, entity in chains:
            chain = _plain(strand.get(asym)) or asym
            for num, mon, hetero in self.sequences.get(entity, ()):
                here = first.get((asym, num))  # gives its number to all residues at a position
                mine = here if here and here[1] == mon else others.get((asym, num, mon))
                # author identifiers, or the label ones where a file has none (as readers do)
                number = (_plain(here[5]) or _plain(here[2])) if here else ''
                icode = _plain(here[6]) if here else ''
                auth_number = (_plain(mine[5]) or _plain(mine[2])) if mine else ''
                name = (_plain(mine[4]) or mon) if mine else ''
                rows.append([asym, entity, num, mon, num, number or '.', auth_number or '?', name or '?',
                             name or '?', chain, icode or '.', hetero])
        return rows


def table(category, columns, rows):
    """Yield the lines of a new loop_ with aligned columns (a line at a time,
    for the large tables of large structures)."""
    for row in rows:
        row[:] = [_cif.quote(v) for v in row]
    widths = [max(len(row[i]) for row in rows) for i in range(len(columns))]
    yield 'loop_\n' + ''.join('{}.{}\n'.format(category, c) for c in columns)
    for row in rows:
        yield ' '.join(v.ljust(w) for v, w in zip(row, widths)).rstrip(' ') + '\n'


def present(source):
    """[_Present] for each data block."""
    blocks = []
    atoms = info = None
    for ev in _cif.parse(source.lines(), (_atoms.ATOM_SITE, '_struct_asym', '_entity_poly_seq',
                                          '_chem_comp_atom', '_chem_comp_bond', SCHEME)):
        etype = type(ev)
        if etype is _cif.Row:
            if ev.loop is not atoms:
                info.read(ev)
                continue
            v = ev.values
            names_add(name(v))
            key = residue(v)
            if key == last:
                continue  # the next atom of the same residue
            last = key
            if model_pos is not None:
                if first_model is None:
                    first_model = v[model_pos]
                in_first = v[model_pos] == first_model
            else:
                in_first = True
            asym = key[0]
            if asym_pos is not None:
                info.asyms.add(asym)
            if entity_pos is not None:
                entity = key[-1]
                info.entities.setdefault(entity, set()).add(asym)
                info.chains.setdefault(asym, entity)
                if auth_pos is not None:
                    strands = info.strands.setdefault(entity, [])
                    if v[auth_pos] not in strands:
                        strands.append(v[auth_pos])
            if in_first:
                info.first_model.append(tuple([shared(x, x) for x in key]))  # one copy of each value
        elif etype is _cif.DataBlock:
            blocks.append(_Present())
        elif etype is _cif.Loop:
            if not blocks:
                blocks.append(_Present())
            info = blocks[-1]
            category = ev.category.lower()
            if category == '_struct_asym':
                info.struct_asym = []
            elif category == SCHEME:
                info.had_scheme = True
            elif category == _atoms.ATOM_SITE:
                atoms = ev
                pos = {c: ev.index(c) if c in ev else None for c in (
                    'label_entity_id', 'auth_asym_id', 'pdbx_PDB_model_num')}
                entity_pos, auth_pos, model_pos = (pos['label_entity_id'], pos['auth_asym_id'],
                                                   pos['pdbx_PDB_model_num'])
                positions = [ev.index(c) if c in ev else None for c, _ in RESIDUE]
                asym_pos = positions[0]
                info.columns = [p is not None for p in positions]
                info.asyms = set() if asym_pos is not None else None
                info.entities = {} if entity_pos is not None else None
                residue = _atoms.getter(*positions, entity_pos)
                name = _atoms.getter(*(ev.index(c) if c in ev else None for c in NAMES))
                names_add, shared = info.names.add, {}.setdefault
                first_model = last = None
    for info in blocks:
        info.finish()
    return blocks


def _keep(ev, info):
    """Whether a row of a kept table still matches the atoms; recomputes the counts in _entity rows."""
    category = ev.loop.category.lower()
    if category == '_struct_asym':
        return info.asyms is None or ev.get('id') in info.asyms  # else a chain with no atoms left
    if category in COMPONENTS:
        return ev.get('comp_id') not in info.renamed
    if category == MOD_RESIDUE:
        return info.mod_residue(ev)
    if category in BY_ENTITY and info.entities is not None:
        entity = ev.get(BY_ENTITY[category])
        if entity not in info.entities:
            return False  # a molecule with no atoms left
        if category == '_entity' and 'pdbx_number_of_molecules' in ev:
            water = (ev.get('type') or '').lower() == 'water'
            count = str(info.molecules.get(entity, 0) if water else len(info.entities[entity]))
            if ev['pdbx_number_of_molecules'] not in NULLS and ev['pdbx_number_of_molecules'] != count:
                ev['pdbx_number_of_molecules'] = count
        if category == '_entity_poly' and 'pdbx_strand_id' in ev and info.strands.get(entity):
            strands = ','.join(info.strands[entity])
            if ev['pdbx_strand_id'] not in NULLS and ev['pdbx_strand_id'] != strands:
                ev['pdbx_strand_id'] = strands
    return True


def run(fhandle, everything=False, scheme=None):
    """Yield the events of the file without its metadata (all of it if
    `everything`, else all but KEPT, brought in step with the atoms)."""
    if everything:
        for ev in _cif.parse(fhandle, COORDINATES):
            if type(ev) is not _cif.Verbatim:  # every other category is Verbatim here
                yield ev
        return
    source = _atoms.as_source(fhandle)
    blocks = present(source)
    block, info, rebuild = -1, None, False
    for ev in _cif.parse(source.lines(), COORDINATES + KEPT + (SCHEME,)):
        etype = type(ev)
        if etype is _cif.Verbatim:
            continue
        if etype is _cif.DataBlock or (etype is _cif.Loop and block < 0):
            block += 1
            info = blocks[block] if block < len(blocks) else None
            rebuild = info is not None and info.can_rebuild()  # else the scheme is left as it is
        if etype is _cif.Row:
            category = ev.loop.category.lower()
            if (category == SCHEME and rebuild) or (info is not None and not _keep(ev, info)):
                continue
        elif etype is _cif.Loop and ev.category.lower() == SCHEME and rebuild:
            continue  # replaced at its LoopEnd
        elif etype is _cif.LoopEnd:
            category = ev.loop.category.lower()
            if category == SCHEME and rebuild:
                rows = info.scheme()
                if rows:
                    yield from table(SCHEME, SCHEME_COLUMNS, rows)
                if ev.text:
                    yield ev.text
                continue
            if category == '_entity_poly_seq' and rebuild and not info.had_scheme:
                rows = info.scheme()
                if rows:
                    yield ev
                    yield from table(SCHEME, SCHEME_COLUMNS, rows)
                    yield '#\n'
                    continue
        yield ev


def main():
    fhandle, everything, scheme = check_input(sys.argv[1:])
    _cli.output(_cli.lazy(run, fhandle, everything, scheme), fhandle)


if __name__ == '__main__':
    main()
