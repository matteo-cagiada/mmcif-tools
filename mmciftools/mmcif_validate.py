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
occupancies outside 0-1.

Prints one line per problem (the first 20 of each kind) and a summary.
The exit status is 1 if there are errors, 0 otherwise.

Usage:
    mmcif_validate <cif file>

Example:
    mmcif_validate 1abc.cif

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cif, _cli

SHOWN = 20  # messages printed per kind of problem
NULLS = ('.', '?')
NUMBERS = ('Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy', 'B_iso_or_equiv')
INTEGERS = ('auth_seq_id', 'label_seq_id', 'pdbx_PDB_model_num')
SCHEME_FIELDS = (('atom name', 'auth_atom_id', 'label_atom_id'),
                 ('residue name', 'auth_comp_id', 'label_comp_id'),
                 ('chain', 'auth_asym_id', 'label_asym_id'),
                 ('residue number', 'auth_seq_id', 'label_seq_id'))


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
    blocks, seen_atoms = 0, False
    try:
        for ev in _cif.parse(source.lines(), _cif.ALL):
            etype = type(ev)
            if etype is _cif.Row:
                if kind != _atoms.ATOM_SITE:
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
            elif etype is _cif.Loop:
                kind = ev.category.lower()
                if kind == _atoms.ATOM_SITE:
                    seen_atoms = True
                    _check_columns(ev, report)
                    numbers = [(c, ev.index(c)) for c in NUMBERS if c in ev]
                    integers = [(c, ev.index(c)) for c in INTEGERS if c in ev]
                    group = ev.index('group_PDB') if 'group_PDB' in ev else None
                    id_pos = ev.index('id') if 'id' in ev else None
            elif etype is _cif.LoopEnd:
                kind = None
            elif etype is _cif.DataBlock:
                blocks += 1
                ids, kind = set(), None
    except _cif.CIFError as err:
        report.add('ERROR', 'syntax', 'syntax: {}'.format(err))
        return report
    if not blocks:
        report.add('ERROR', 'no data', 'no data_ block')
    elif not seen_atoms:
        report.add('ERROR', 'no atoms', 'no _atom_site category')
    else:
        for message in _atoms.reference_problems(source):
            report.add('ERROR', 'reference', message)
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
