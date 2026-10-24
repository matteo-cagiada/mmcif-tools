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
Builds a biological assembly from the operations stored in an mmCIF file.

PDB entries describe their biological assemblies (_pdbx_struct_assembly) as
operations (_pdbx_struct_oper_list: a rotation matrix and a translation)
applied to label chains (_pdbx_struct_assembly_gen). This tool applies them
and writes the assembly's atoms in place of the asymmetric unit, offline and
for any file that carries these categories (mmcif_fetch -biounit downloads
the RCSB's prebuilt files instead). The option is the assembly id (default
1); mmcif_info lists the assemblies of a file.

Operation expressions follow the PDBx dictionary: 1, 1,2,5, (1-60), and
products such as (1-60)(61-88), whose operations are applied right to left.

Naming follows the RCSB's assembly files: the copy made by the identity
operation keeps its chain names, and every other copy gets the operation
added to its label and author chains (A-2; A-1x61 for a product). Atoms are
written chain by chain, each chain followed by its copies, in every model;
atom ids are renumbered from 1. Links (_struct_conn) are repeated in each
copy that holds both partners; links between copies or to symmetry mates
are not written. _atom_site_anisotrop is removed (its tensors would need
rotating). The other categories still describe the asymmetric unit.

Usage:
    mmcif_assembly [-<assembly id>] <cif file>

Example:
    mmcif_assembly 1abc.cif       # assembly 1
    mmcif_assembly -2 1abc.cif    # assembly 2

This program is part of the mmcif-tools suite.
"""

import itertools
import re
import sys

from mmciftools import _atoms, _cif, _cli

NULLS = ('.', '?')
IDENTITY = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)), (0.0, 0.0, 0.0)
_GROUP = re.compile(r'\(([^()]*)\)')


def check_input(args):
    """Return (file handle, assembly id)."""
    fhandle, option, _ = _cli.parse_args(args, __doc__)
    if option is not None and not option.strip():
        _cli.fail(__doc__, 'Invalid assembly id: {!r}'.format(option))
    return fhandle, option or '1'


def parse_expression(text):
    """'(1-3)(X0)' -> [['1', '2', '3'], ['X0']]; '1,2' -> [['1', '2']]."""
    text = ''.join(text.split())
    groups = _GROUP.findall(text) if text.startswith('(') else [text]
    if text.startswith('(') and ''.join('({})'.format(g) for g in groups) != text:
        raise _cif.CIFError('Invalid operation expression: {!r}'.format(text))
    out = []
    for group in groups:
        ids = []
        for part in group.split(','):
            first, dash, last = part.partition('-')
            if dash:
                try:
                    ids.extend(str(n) for n in range(int(first), int(last) + 1))
                except ValueError:
                    raise _cif.CIFError('Invalid operation range: {!r}'.format(part))
            elif part:
                ids.append(part)
        if not ids:
            raise _cif.CIFError('Invalid operation expression: {!r}'.format(text))
        out.append(ids)
    return out


def _operation(row):
    """(rotation, translation) of a _pdbx_struct_oper_list row."""
    def number(column, default):
        value = row.get(column)
        if value is None or value in NULLS:
            return default
        try:
            return float(value)
        except ValueError:
            raise _cif.CIFError('_pdbx_struct_oper_list {}: {} is not a number: {!r}'.format(
                row.get('id'), column, value))
    matrix = tuple(tuple(number('matrix[{}][{}]'.format(i, j), float(i == j)) for j in (1, 2, 3))
                   for i in (1, 2, 3))
    vector = tuple(number('vector[{}]'.format(i), 0.0) for i in (1, 2, 3))
    return matrix, vector


def _then(first, second):
    """The operation that applies `first`, then `second`."""
    (r1, t1), (r2, t2) = first, second
    matrix = tuple(tuple(sum(r2[i][k] * r1[k][j] for k in range(3)) for j in range(3))
                   for i in range(3))
    vector = tuple(sum(r2[i][k] * t1[k] for k in range(3)) + t2[i] for i in range(3))
    return matrix, vector


def _is_identity(op):
    return all(abs(a - b) < 1e-6 for row, ref in zip(op[0] + (op[1],), IDENTITY[0] + (IDENTITY[1],))
               for a, b in zip(row, ref))


def read_assembly(source, assembly):
    """[(label chains, [(operation ids, operation)])] for one assembly."""
    tables = {}
    for ev in _cif.parse(source.lines(), ('_pdbx_struct_assembly_gen', '_pdbx_struct_oper_list')):
        if type(ev) is _cif.Row:
            tables.setdefault(ev.loop.category.lower(), []).append(ev)
    gens = tables.get('_pdbx_struct_assembly_gen', [])
    if not gens:
        raise _cif.CIFError('No assemblies in the file (_pdbx_struct_assembly_gen)')
    known = list(dict.fromkeys(row.get('assembly_id') for row in gens))
    if assembly not in known:
        raise _cif.CIFError('No assembly {!r}; the file has {}'.format(assembly, ', '.join(known)))
    operations = {row.get('id'): _operation(row) for row in tables.get('_pdbx_struct_oper_list', [])}
    out = []
    for row in gens:
        if row.get('assembly_id') != assembly:
            continue
        combos = []
        for ids in itertools.product(*parse_expression(row.get('oper_expression') or '')):
            missing = [i for i in ids if i not in operations]
            if missing:
                raise _cif.CIFError('No operation {} in _pdbx_struct_oper_list'.format(missing[0]))
            op = IDENTITY
            for i in reversed(ids):  # right to left
                op = _then(op, operations[i])
            combos.append((ids, op))
        chains = [c for c in _atoms.parse_ids(row.get('asym_id_list') or '')]
        out.append((chains, combos))
    return out


def copies(generators):
    """[(label chains, operation ids, operation, {chain: name suffix})]; the
    first identity copy of a chain keeps its names (suffix '')."""
    kept, out = set(), []
    for chains, combos in generators:
        for ids, op in combos:
            suffixes = {}
            for chain in chains:
                if _is_identity(op) and chain not in kept:
                    kept.add(chain)
                    suffixes[chain] = ''
                else:
                    suffixes[chain] = '-' + 'x'.join(ids)
            out.append((chains, ids, op, suffixes))
    return out


def _coordinates(op, xyz):
    (r, t) = op
    return ['{:.3f}'.format(r[i][0] * xyz[0] + r[i][1] * xyz[1] + r[i][2] * xyz[2] + t[i])
            for i in range(3)]


def run(fhandle, assembly='1', scheme=None):
    """Yield the events of the file with the atoms of the assembly."""
    source = _atoms.as_source(fhandle)
    made = copies(read_assembly(source, assembly))
    kind, atoms = None, None
    categories = (_atoms.ATOM_SITE, _atoms.ANISOTROP, _atoms.STRUCT_CONN)
    for ev in _cif.parse(source.lines(), categories):
        etype = type(ev)
        if etype is _cif.Row:
            if kind == _atoms.ATOM_SITE:
                atoms.append(ev)
            elif kind == _atoms.STRUCT_CONN:
                for row in _link_copies(ev, made):
                    yield row
            continue  # anisotropic rows are dropped
        if etype is _cif.Loop:
            kind = ev.category.lower()
            if kind == _atoms.ATOM_SITE:
                atoms = []
        elif etype is _cif.LoopEnd and kind == _atoms.ATOM_SITE:
            for row in _atom_copies(ev.loop, atoms, made):
                yield row
            kind = None
        elif etype is _cif.LoopEnd:
            kind = None
        yield ev


def _atom_copies(loop, atoms, made):
    """The rows of the assembly, from the rows of the asymmetric unit."""
    for name in ('label_asym_id', 'Cartn_x', 'Cartn_y', 'Cartn_z'):
        if name not in loop:
            raise _cif.CIFError('_atom_site has no {} column (needed here)'.format(name), loop.lineno)
    label = loop.index('label_asym_id')
    xyz = [loop.index(c) for c in ('Cartn_x', 'Cartn_y', 'Cartn_z')]
    model = loop.index('pdbx_PDB_model_num') if 'pdbx_PDB_model_num' in loop else None
    auth_col = 'auth_asym_id' if 'auth_asym_id' in loop else None
    by_chain = {}  # (model, label chain) -> rows, in file order
    models = []
    for row in atoms:
        key = (row.values[model] if model is not None else None, row.values[label])
        if key[0] not in models:
            models.append(key[0])
        by_chain.setdefault(key, []).append(row)
    serial = 0
    for model_id in models:
        for chains, _, _, _ in _unique_generators(made):
            for chain in chains:
                rows = by_chain.get((model_id, chain), ())
                for _, ids, op, suffixes in (c for c in made if c[0] is chains):
                    suffix = suffixes[chain]
                    for row in rows:
                        try:
                            point = [float(row.values[p]) for p in xyz]
                        except ValueError:
                            raise _cif.CIFError('Coordinates are not numbers', row.lineno)
                        new = _cif.Row(loop, row.values, row.lead, row.text, row.lineno)
                        serial += 1
                        if 'id' in loop:
                            new['id'] = str(serial)
                        if suffix:
                            new['label_asym_id'] = chain + suffix
                            if auth_col is not None:
                                new[auth_col] = row.values[loop.index(auth_col)] + suffix
                        if not _is_identity(op):
                            for p, value in zip(xyz, _coordinates(op, point)):
                                new[loop.columns[p]] = value
                        yield new


def _unique_generators(made):
    """The distinct generator entries of `made`, in order."""
    seen = []
    for entry in made:
        if not any(entry[0] is s[0] for s in seen):
            seen.append(entry)
    return seen


def _link_copies(row, made):
    """The copies of a _struct_conn row (none if its partners are not in one
    copy, or if a partner is a symmetry mate)."""
    asyms = [row.get('ptnr{}_label_asym_id'.format(n)) for n in (1, 2)]
    symmetry = [row.get('ptnr{}_symmetry'.format(n)) for n in (1, 2)]
    if None in asyms or any(s not in (None, '1_555') + NULLS for s in symmetry):
        return
    for chains, _, _, suffixes in made:
        if asyms[0] not in chains or asyms[1] not in chains:
            continue
        new = _cif.Row(row.loop, row.values, row.lead, row.text, row.lineno)
        ends = [suffixes[a] for a in asyms]
        if 'id' in row.loop and (ends[0] or ends[1]):
            new['id'] = row['id'] + (ends[0] or ends[1])
        for n, suffix in zip((1, 2), ends):
            if suffix:
                for column in ('ptnr{}_label_asym_id'.format(n), 'ptnr{}_auth_asym_id'.format(n)):
                    if column in row.loop and row[column] not in NULLS:
                        new[column] = row[column] + suffix
        yield new


def main():
    fhandle, assembly = check_input(sys.argv[1:])
    _cli.output(_cli.lazy(run, fhandle, assembly), fhandle)


if __name__ == '__main__':
    main()
