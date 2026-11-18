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
Validates the structure of an mmCIF file.

Checks the CIF syntax and the coordinates (_atom_site); it does not check
the file against the PDBx/mmCIF dictionary (item types, allowed values,
mandatory items). Errors:

    - CIF syntax errors (with their line number)
    - no _atom_site category, or required _atom_site columns missing (id,
      type_symbol, Cartn_x/y/z, and an atom name, residue name, chain and
      residue number in at least one identifier scheme)
    - coordinates, occupancies or B-factors that are not numbers
    - group_PDB other than ATOM / HETATM
    - residue numbers (auth_seq_id, label_seq_id) or model numbers that are
      not integers ('.' and '?' are allowed)
    - the same atom id twice in a data block
    - rows of _atom_site_anisotrop or _struct_conn that match no atom

Warnings: only one identifier scheme (auth_* or label_*) present,
occupancies outside 0-1, and metadata out of step with the atoms: rows of
the tables that name residues or atoms (_struct_conf,
_pdbx_poly_seq_scheme, _pdbx_validate_*, ...) whose author ids differ from
those of the atoms with the same label ids, or that name atoms or residues
the file does not have (tables of missing residues aside).

Warnings also for molecules the file does not describe, which structure
prediction programs read (files written by other programs often lack them):

    - a missing core table: _entity (the molecules), _struct_asym (the
      chains), _chem_comp (the residue types), and for polymers
      _entity_poly, _entity_poly_seq (the sequences) and
      _pdbx_poly_seq_scheme (the residues along them); or a missing
      _atom_site.label_entity_id column
    - an entity, a chain (label_asym_id), a polymer residue or a residue
      type of the atoms that these tables do not list

Prints one line per problem (the first 20 of each kind) and a summary.
The exit status is 1 if there are errors, 0 otherwise.

Usage:
    mmcif_validate <cif file>

Example:
    mmcif_validate 1abc.cif

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cif, _cli, _refs

SHOWN = 20  # messages printed per kind of problem
NULLS = ('.', '?')
NUMBERS = ('Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy', 'B_iso_or_equiv')
INTEGERS = ('auth_seq_id', 'label_seq_id', 'pdbx_PDB_model_num')
SCHEME_FIELDS = (('atom name', 'auth_atom_id', 'label_atom_id'),
                 ('residue name', 'auth_comp_id', 'label_comp_id'),
                 ('chain', 'auth_asym_id', 'label_asym_id'),
                 ('residue number', 'auth_seq_id', 'label_seq_id'))
# The tables that describe the molecules of the atoms, and what each one is missing tells.
CORE_TABLES = {
    '_entity': 'the molecules (entities) of the atoms are not described',
    '_struct_asym': 'the chains (label_asym_id) of the atoms are not described',
    '_chem_comp': 'the residue types are not described',
    '_entity_poly': 'the polymers are not described (type, one-letter sequence)',
    '_entity_poly_seq': 'the sequences of the polymers are not given',
    '_pdbx_poly_seq_scheme': 'the residues of the polymers are not placed along their sequences',
}
POLYMER_TABLES = ('_entity_poly', '_entity_poly_seq', '_pdbx_poly_seq_scheme')


class Molecules:
    """What the atoms of a data block name (entities, chains, polymer
    residues, residue types) and what its core tables describe."""

    def __init__(self):
        self.tables = set()  # the core tables present
        self.kinds = {}  # _entity: id -> type
        self.chains = {}  # _struct_asym: label_asym_id -> its entity ids (one, unless the id is repeated)
        self.sequence = set()  # _entity_poly_seq: (entity, num, mon_id)
        self.scheme = set()  # _pdbx_poly_seq_scheme: (asym_id, seq_id, mon_id)
        self.components = set()  # _chem_comp: id
        self.atoms = False
        self.no_entity_column = None  # line of an _atom_site loop without label_entity_id
        self.entities, self.atom_chains, self.residues, self.types = {}, {}, {}, {}  # -> first atom line
        self.positions = None
        self.last = None

    def table_row(self, category, row):
        get = row.get
        if category == '_entity':
            self.kinds[get('id')] = (get('type') or '').lower()
        elif category == '_struct_asym':
            self.chains.setdefault(get('id'), []).append(get('entity_id'))
        elif category == '_entity_poly_seq':
            self.sequence.add((get('entity_id'), get('num'), get('mon_id')))
        elif category == '_pdbx_poly_seq_scheme':
            self.scheme.add((get('asym_id'), get('seq_id'), get('mon_id')))
        elif category == '_chem_comp':
            self.components.add(get('id'))

    def atom_loop(self, loop):
        self.atoms = True
        self.last = None
        if 'label_entity_id' not in loop:
            self.no_entity_column = loop.lineno
        comp = 'label_comp_id' if 'label_comp_id' in loop else 'auth_comp_id'
        self.positions = tuple(loop.index(c) if c in loop else None
                               for c in ('label_entity_id', 'label_asym_id', 'label_seq_id', comp))

    def atom(self, values, line):
        key = tuple(None if p is None else values[p] for p in self.positions)
        if key == self.last:
            return
        self.last = key
        entity, chain, number, comp = (None if v in NULLS else v for v in key)
        if entity is not None:
            self.entities.setdefault(entity, line)
        if chain is not None:
            self.atom_chains.setdefault((chain, entity), line)
        if number is not None and chain is not None:
            self.residues.setdefault((entity, chain, number, comp), line)
        if comp is not None:
            self.types.setdefault(comp, line)

    def check(self, report):
        """Add the warnings of a data block with atoms."""
        if not self.atoms:
            return
        polymers = bool(self.residues) or 'polymer' in self.kinds.values()
        for table, why in CORE_TABLES.items():
            if table not in self.tables and (polymers or table not in POLYMER_TABLES):
                report.add('WARNING', 'core table', 'no {}: {}'.format(table, why))
        if self.no_entity_column is not None:
            report.add('WARNING', 'core table', 'line {}: _atom_site has no label_entity_id: the atoms are not '
                       'linked to their molecules'.format(self.no_entity_column))
        if '_entity' in self.tables:
            for entity, line in self.entities.items():
                if entity not in self.kinds:
                    report.add('WARNING', 'entity', 'line {}: entity {} of the atoms is not in _entity'.format(
                        line, entity))
        if '_struct_asym' in self.tables:
            for chain, entities in self.chains.items():
                if len(entities) > 1:
                    report.add('WARNING', 'chain', 'chain {} is listed {} times in _struct_asym (entities {}): a '
                               'chain (label_asym_id) names one molecule'.format(chain, len(entities),
                                                                                 ', '.join(dict.fromkeys(entities))))
            for (chain, entity), line in self.atom_chains.items():
                listed = self.chains.get(chain)
                if listed is None:
                    report.add('WARNING', 'chain', 'line {}: chain {} of the atoms is not in _struct_asym'.format(
                        line, chain))
                elif entity is not None and len(listed) == 1 and listed[0] not in NULLS and listed[0] != entity:
                    report.add('WARNING', 'chain', 'line {}: chain {} is entity {} in the atoms, {} in '
                               '_struct_asym'.format(line, chain, entity, listed[0]))
        for (entity, chain, number, comp), line in self.residues.items():
            if self.kinds.get(entity, 'polymer') != 'polymer':
                continue  # a branched molecule (glycan) or a ligand: no sequence
            if '_entity_poly_seq' in self.tables and entity is not None and \
                    (entity, number, comp) not in self.sequence:
                report.add('WARNING', 'sequence', 'line {}: residue {} {} of chain {} is not in the sequence of '
                           'entity {} (_entity_poly_seq)'.format(line, number, comp, chain, entity))
            if '_pdbx_poly_seq_scheme' in self.tables and (chain, number, comp) not in self.scheme:
                report.add('WARNING', 'scheme', 'line {}: residue {} {} of chain {} is not in '
                           '_pdbx_poly_seq_scheme'.format(line, number, comp, chain))
        if '_chem_comp' in self.tables:
            for comp, line in self.types.items():
                if comp not in self.components:
                    report.add('WARNING', 'component', 'line {}: residue type {} is not in _chem_comp'.format(
                        line, comp))


def check_input(args):
    """Return the file handle."""
    fhandle, _, _ = _cli.parse_args(args, __doc__, takes_option=False)
    return fhandle


class Report:
    """Problems found, grouped by kind."""

    def __init__(self):
        self.kinds = {}  # kind -> [count, level, [messages]]

    def add(self, level, kind, message):
        entry = self.kinds.setdefault(kind, [0, level, []])
        entry[0] += 1
        if len(entry[2]) < SHOWN:
            entry[2].append(message)

    def count(self, level):
        return sum(e[0] for e in self.kinds.values() if e[1] == level)

    def lines(self):
        out = []
        for kind, (count, level, messages) in self.kinds.items():
            for message in messages:
                out.append('{}: {}\n'.format(level, message))
            if count > len(messages):
                out.append('{}: ... and {} more ({})\n'.format(level, count - len(messages), kind))
        errors, warnings = self.count('ERROR'), self.count('WARNING')
        if errors:
            out.append('{} error(s), {} warning(s): the file is NOT valid.\n'.format(errors, warnings))
        else:
            out.append('No errors, {} warning(s): the file looks valid.\n'.format(warnings))
        return out


def _check_columns(loop, report):
    for column in ('id', 'type_symbol', 'Cartn_x', 'Cartn_y', 'Cartn_z'):
        if column not in loop:
            report.add('ERROR', 'missing column',
                       'line {}: _atom_site has no {} column'.format(loop.lineno, column))
    for name, auth, label in SCHEME_FIELDS:
        if auth not in loop and label not in loop:
            report.add('ERROR', 'missing column', 'line {}: _atom_site has no {} ({} or {})'.format(
                loop.lineno, name, auth, label))
        elif auth not in loop or label not in loop:
            report.add('WARNING', 'one scheme', 'line {}: _atom_site has only {} for the {}'.format(
                loop.lineno, auth if auth in loop else label, name))


def validate(fhandle):
    """Return a Report on the file."""
    report = Report()
    source = _atoms.as_source(fhandle)
    blocks, seen_atoms, molecules = 0, False, None
    try:
        for ev in _cif.parse(source.lines(), _cif.ALL):
            etype = type(ev)
            if etype is _cif.Row:
                if kind != _atoms.ATOM_SITE:
                    if kind in CORE_TABLES:
                        molecules.table_row(kind, ev)
                    continue
                values, line = ev.values, ev.lineno
                for column, pos in numbers:
                    value = values[pos]
                    if value not in NULLS:
                        try:
                            number = float(value)
                        except ValueError:
                            report.add('ERROR', 'not a number', 'line {}: {} is not a number: {!r}'.format(
                                line, column, value))
                            continue
                        if column == 'occupancy' and not 0.0 <= number <= 1.0:
                            report.add('WARNING', 'occupancy', 'line {}: occupancy outside 0-1: {}'.format(
                                line, value))
                for column, pos in integers:
                    value = values[pos]
                    if value not in NULLS and _atoms.as_int(value) is None:
                        report.add('ERROR', 'not an integer', 'line {}: {} is not an integer: {!r}'.format(
                            line, column, value))
                if group is not None and values[group] not in ('ATOM', 'HETATM'):
                    report.add('ERROR', 'group_PDB', 'line {}: group_PDB is not ATOM or HETATM: {!r}'.format(
                        line, values[group]))
                if id_pos is not None:
                    atom_id = values[id_pos]
                    if atom_id in ids:
                        report.add('ERROR', 'duplicate id', 'line {}: atom id {} used twice'.format(line, atom_id))
                    ids.add(atom_id)
                molecules.atom(values, line)
            elif etype is _cif.Loop:
                kind = ev.category.lower()
                if kind in CORE_TABLES:
                    molecules.tables.add(kind)
                if kind == _atoms.ATOM_SITE:
                    seen_atoms = True
                    molecules.atom_loop(ev)
                    _check_columns(ev, report)
                    numbers = [(c, ev.index(c)) for c in NUMBERS if c in ev]
                    integers = [(c, ev.index(c)) for c in INTEGERS if c in ev]
                    group = ev.index('group_PDB') if 'group_PDB' in ev else None
                    id_pos = ev.index('id') if 'id' in ev else None
            elif etype is _cif.LoopEnd:
                kind = None
            elif etype is _cif.DataBlock:
                if molecules is not None:
                    molecules.check(report)
                blocks += 1
                ids, kind, molecules = set(), None, Molecules()
    except _cif.CIFError as err:
        report.add('ERROR', 'syntax', 'syntax: {}'.format(err))
        return report
    if molecules is not None:
        molecules.check(report)
    if not blocks:
        report.add('ERROR', 'no data', 'no data_ block')
    elif not seen_atoms:
        report.add('ERROR', 'no atoms', 'no _atom_site category')
    else:
        for message in _atoms.reference_problems(source):
            report.add('ERROR', 'reference', message)
        for category, kind, message in _refs.problems(source):
            if category not in (_atoms.ANISOTROP, _atoms.STRUCT_CONN):  # errors above
                report.add('WARNING', 'metadata', message)
    return report


def run(fhandle, scheme=None):
    """Return (report lines, number of errors)."""
    report = validate(fhandle)
    return report.lines(), report.count('ERROR')


def main():
    fhandle = check_input(sys.argv[1:])
    lines, errors = run(fhandle)
    fhandle.close()
    sys.stdout.writelines(lines)
    sys.exit(1 if errors else 0)


if __name__ == '__main__':
    main()
