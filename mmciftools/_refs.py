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
Where the metadata of an mmCIF file names chains, residues and atoms.

A table (category) names them through groups of columns: _struct_conf names
the first and last residue of a helix (beg_label_asym_id, beg_label_seq_id,
beg_auth_seq_id, ...), _pdbx_validate_close_contact two atoms by their author
ids, _ma_qa_metric_local one residue by its label ids. Each column of a group
has a role, named after the _atom_site column it matches:

    label_asym  label_seq  label_comp  label_atom    the label ids
    auth_asym   auth_seq   auth_comp   auth_atom     the author ids
    icode                                            pdbx_PDB_ins_code
    alt                                              label_alt_id

A group is compared with the atoms on the roles both have; '.' and '?' are
the same, alt ids are not compared, and an empty atom or chain name leaves
its role out (a group with no chain names nothing).
"""

import operator
import re

from mmciftools import _cif

ATOM_SITE = '_atom_site'
ATOM_COLUMNS = {
    'label_asym': 'label_asym_id', 'label_seq': 'label_seq_id', 'label_comp': 'label_comp_id',
    'label_atom': 'label_atom_id', 'auth_asym': 'auth_asym_id', 'auth_seq': 'auth_seq_id',
    'auth_comp': 'auth_comp_id', 'auth_atom': 'auth_atom_id', 'icode': 'pdbx_PDB_ins_code',
    'alt': 'label_alt_id',
}
ROLES = tuple(ATOM_COLUMNS)
RESIDUE_ROLES = ('label_asym', 'label_seq', 'label_comp', 'auth_asym', 'auth_seq', 'auth_comp', 'icode')
ATOM_ROLES = ('label_atom', 'auth_atom')
CHAIN_ROLES = ('label_asym', 'auth_asym')
LABEL_RESIDUE = ('label_asym', 'label_seq', 'label_comp')
AUTHOR_RESIDUE = ('auth_asym', 'auth_seq', 'icode')
# Roles whose null value is a value: no label number (a ligand or water), no insertion code.
NULL_IS_VALUE = ('label_seq', 'icode')
NULLS = ('.', '?')


def _norm(value):
    return '' if value in NULLS else value


class Group:
    """Columns of a table that together name a chain, a residue or an atom.

    `columns` maps a role to the column, or to several columns that hold the
    same value (the first one present is compared, all are updated). `part`
    is None, or 'beg' / 'end' for the ends of a range.
    """

    __slots__ = ('name', 'columns', 'part')

    def __init__(self, name, columns, part=None):
        self.name = name
        self.columns = {role: (cols,) if isinstance(cols, str) else tuple(cols) for role, cols in columns.items()}
        self.part = part


class Table:
    """The groups of a table. `observed`: whether what it names has atoms
    (False for the sequence schemes and the lists of missing residues).
    `lists` maps a column holding a comma-separated list of chains to the
    role of its chains."""

    __slots__ = ('category', 'groups', 'observed', 'lists')

    def __init__(self, category, groups, observed=True, lists=None):
        self.category = category
        self.groups = groups
        self.observed = observed
        self.lists = lists or {}


def _ids(name, prefix='', suffix='', icode=None, alt=None, label=True, auth=True, atom=True, part=None):
    """A group whose columns follow _atom_site: prefix + 'label_asym_id' + suffix, ..."""
    roles = (['label_asym', 'label_seq', 'label_comp'] if label else []) + \
        (['auth_asym', 'auth_seq', 'auth_comp'] if auth else [])
    if atom:
        roles += (['label_atom'] if label else []) + (['auth_atom'] if auth else [])
    columns = {role: prefix + ATOM_COLUMNS[role] + suffix for role in roles}
    if icode:
        columns['icode'] = icode
    if alt:
        columns['alt'] = alt
    return Group(name, columns, part)


def _numbered(n, **kw):
    """Group n of a table whose columns end in _1, _2, ... (_pdbx_validate_*)."""
    s = '_{}'.format(n)
    return _ids('atom {}'.format(n), '', s, icode='PDB_ins_code' + s, alt='label_alt_id' + s, **kw)


def _tables():
    t = []

    def add(category, *groups, **kw):
        t.append(Table(category, list(groups), **kw))

    # Atoms and links (also kept in step by _atoms.filter_atoms and _atoms.rename)
    add('_atom_site_anisotrop', _ids('atom', 'pdbx_', icode='pdbx_PDB_ins_code', alt='pdbx_label_alt_id'))
    add('_struct_conn', *[_ids('partner {}'.format(n), 'ptnr{}_'.format(n), icode='pdbx_ptnr{}_PDB_ins_code'.format(n),
                               alt='pdbx_ptnr{}_label_alt_id'.format(n)) for n in (1, 2)],
        _ids('partner 3', 'pdbx_ptnr3_', icode='pdbx_ptnr3_PDB_ins_code', alt='pdbx_ptnr3_label_alt_id', auth=False))
    add('_pdbx_struct_conn_angle', *[_ids('partner {}'.format(n), 'ptnr{}_'.format(n),
                                          icode='ptnr{}_PDB_ins_code'.format(n), alt='ptnr{}_label_alt_id'.format(n))
                                     for n in (1, 2, 3)])
    # Secondary structure, cis peptides, sites
    for category in ('_struct_conf', '_struct_sheet_range'):
        add(category, *[_ids(part, part + '_', icode='pdbx_{}_PDB_ins_code'.format(part), atom=False, part=part)
                        for part in ('beg', 'end')])
    add('_pdbx_struct_sheet_hbond', *[_ids('range {}'.format(n), 'range_{}_'.format(n),
                                           icode='range_{}_PDB_ins_code'.format(n)) for n in (1, 2)])
    add('_struct_mon_prot_cis', _ids('residue 1', icode='pdbx_PDB_ins_code', alt='label_alt_id', atom=False),
        _ids('residue 2', 'pdbx_', '_2', icode='pdbx_PDB_ins_code_2', atom=False))
    add('_struct_site', Group('residue', {'auth_asym': 'pdbx_auth_asym_id', 'auth_seq': 'pdbx_auth_seq_id',
                                          'auth_comp': 'pdbx_auth_comp_id', 'icode': 'pdbx_auth_ins_code'}))
    add('_struct_site_gen', _ids('residue', icode='pdbx_auth_ins_code', alt='label_alt_id'))
    add('_pdbx_refine_tls_group', *[_ids(part, part + '_', atom=False, part=part) for part in ('beg', 'end')])
    # Modified residues, features, coordination, symmetry
    add('_pdbx_struct_mod_residue', _ids('residue', icode='PDB_ins_code', atom=False))
    add('_pdbx_modification_feature', _ids('residue', icode='PDB_ins_code', alt='label_alt_id', atom=False),
        _ids('modified residue', 'modified_residue_', icode='modified_residue_PDB_ins_code',
             alt='modified_residue_label_alt_id', atom=False))
    for category in ('_pdbx_nonpoly_atom_coordination', '_pdbx_nonpoly_atom_coordination_sphere',
                     '_pdbx_nonpoly_atom_coordination_sphere_order'):
        add(category, _ids('atom', icode='PDB_ins_code', alt='label_alt_id'))
    add('_pdbx_struct_special_symmetry', _ids('residue', icode='PDB_ins_code', atom=False))
    add('_pdbx_entity_instance_feature', Group('residue', {
        'label_asym': 'asym_id', 'label_seq': 'seq_num', 'label_comp': 'comp_id',
        'auth_asym': 'auth_asym_id', 'auth_seq': 'auth_seq_num', 'auth_comp': 'auth_comp_id'}))
    # Nucleic-acid base pairs
    add('_ndb_struct_na_base_pair', *[_ids(b, b + '_', icode=b + '_PDB_ins_code', atom=False) for b in 'ij'])
    add('_ndb_struct_na_base_pair_step', *[_ids('{} {}'.format(b, n), b + '_', '_{}'.format(n),
                                                icode='{}_PDB_ins_code_{}'.format(b, n), atom=False)
                                           for b in 'ij' for n in (1, 2)])
    # Validation reports: author ids only
    for category in ('_pdbx_validate_close_contact', '_pdbx_validate_symm_contact', '_pdbx_validate_rmsd_bond',
                     '_pdbx_validate_polymer_linkage'):
        add(category, _numbered(1, label=False), _numbered(2, label=False))
    add('_pdbx_validate_rmsd_angle', *[_numbered(n, label=False) for n in (1, 2, 3)])
    add('_pdbx_validate_peptide_omega', _numbered(1, label=False, atom=False), _numbered(2, label=False, atom=False))
    add('_pdbx_validate_chiral', _ids('atom', icode='PDB_ins_code', alt='label_alt_id', label=False))
    for category in ('_pdbx_validate_torsion', '_pdbx_validate_planes', '_pdbx_validate_main_chain_plane'):
        add(category, _ids('residue', icode='PDB_ins_code', alt='label_alt_id', label=False, atom=False))
    # Model quality (ModelCIF: AlphaFold DB, ...)
    add('_ma_qa_metric_local', _ids('residue', auth=False, atom=False))
    # Sequences: what they name may have no atoms
    add('_pdbx_poly_seq_scheme', Group('residue', {
        'label_asym': 'asym_id', 'label_seq': 'seq_id', 'label_comp': 'mon_id', 'auth_asym': 'pdb_strand_id',
        'auth_seq': ('pdb_seq_num', 'auth_seq_num'), 'auth_comp': ('pdb_mon_id', 'auth_mon_id'),
        'icode': 'pdb_ins_code'}), observed=False)
    add('_pdbx_nonpoly_scheme', Group('residue', {
        'label_asym': 'asym_id', 'label_comp': 'mon_id', 'auth_asym': 'pdb_strand_id',
        'auth_seq': ('pdb_seq_num', 'auth_seq_num'), 'auth_comp': ('pdb_mon_id', 'auth_mon_id'),
        'icode': 'pdb_ins_code'}), observed=False)
    add('_pdbx_branch_scheme', Group('residue', {
        'label_asym': 'asym_id', 'label_comp': 'mon_id', 'auth_asym': ('pdb_asym_id', 'auth_asym_id'),
        'auth_seq': ('pdb_seq_num', 'auth_seq_num'), 'auth_comp': ('pdb_mon_id', 'auth_mon_id')}), observed=False)
    add('_pdbx_unobs_or_zero_occ_residues', _ids('residue', icode='PDB_ins_code', atom=False), observed=False)
    add('_pdbx_unobs_or_zero_occ_atoms', _ids('atom', icode='PDB_ins_code', alt='label_alt_id'), observed=False)
    add('_struct_ref_seq', *[Group(part, {
        'auth_asym': 'pdbx_strand_id', 'label_seq': 'seq_align_' + part,
        'auth_seq': 'pdbx_auth_seq_align_' + part, 'icode': 'pdbx_seq_align_{}_ins_code'.format(part)}, part)
        for part in ('beg', 'end')], observed=False)
    add('_struct_ref_seq_dif', Group('residue', {
        'auth_asym': 'pdbx_pdb_strand_id', 'label_seq': 'seq_num', 'label_comp': 'mon_id',
        'auth_seq': 'pdbx_auth_seq_num', 'icode': 'pdbx_pdb_ins_code'}), observed=False)
    # Chains
    add('_struct_asym', Group('chain', {'label_asym': 'id'}), observed=False)
    add('_pdbx_molecule', Group('chain', {'label_asym': 'asym_id'}))
    add('_ma_target_entity_instance', Group('chain', {'label_asym': 'asym_id'}))
    add('_pdbx_struct_assembly_gen', lists={'asym_id_list': 'label_asym'})
    add('_entity_poly', lists={'pdbx_strand_id': 'auth_asym'})
    return {table.category: table for table in t}


TABLES = _tables()
CATEGORIES = tuple(TABLES)


class Bound:
    """A Group bound to the columns of one loop: the roles it has there."""

    __slots__ = ('group', 'roles', 'positions', 'updates')

    def __init__(self, group, loop):
        self.group = group
        self.roles, self.positions, self.updates = [], [], {}
        for role, columns in group.columns.items():
            present = [loop.index(c) for c in columns if c in loop]
            if present:
                self.updates[role] = present
                if role != 'alt':
                    self.roles.append(role)
                    self.positions.append(present[0])

    def named(self, values):
        """{role: value} of what the row names here (empty roles left out)."""
        out = {}
        for role, pos in zip(self.roles, self.positions):
            value = _norm(values[pos])
            if value or role in NULL_IS_VALUE:
                out[role] = value
        if not any(role in out for role in CHAIN_ROLES):
            return {}
        return out


def bind(loop):
    """[Bound] of a loop of a listed table (empty for other tables)."""
    table = TABLES.get(loop.category.lower())
    if table is None:
        return []
    bound = [Bound(g, loop) for g in table.groups]
    return [b for b in bound if b.roles]


def atom_roles(loop):
    """{role: position} of the _atom_site columns a loop has."""
    return {role: loop.index(col) for role, col in ATOM_COLUMNS.items() if col in loop and role != 'alt'}


def _key(named, roles):
    return tuple(named[r] for r in roles)


def _level(roles):
    return 'atom' if any(r in ATOM_ROLES for r in roles) else 'residue' if len(set(roles) - set(CHAIN_ROLES)) else 'chain'


_SHAPES = {}


def _shape(named, have):
    """The roles of `named` that the atoms have, in a fixed order (None: nothing to compare)."""
    cached = (tuple(named), tuple(have))
    if cached not in _SHAPES:
        roles = tuple(r for r in ROLES if r in named and r in have)
        _SHAPES[cached] = roles if any(r in roles for r in CHAIN_ROLES) else None
    return _SHAPES[cached]


def _describe(named):
    return ' '.join(named[r] or '.' for r in ROLES if r in named)


def look(source):
    """What the listed tables name that the atoms have: {shape: {key}} (and
    the atom roles, and the author ids found by label ids)."""
    have = None
    for ev in _cif.parse(source.lines(), ATOM_SITE):
        if type(ev) is _cif.Loop:
            have = atom_roles(ev)
            break
    if have is None:
        return None

    # What the tables name, by shape (the roles compared).
    wanted, by_label = {}, {}
    for named, table, bound, row in _references(source):
        shape = _shape(named, have)
        if shape is None:
            continue
        wanted.setdefault(shape, set()).add(_key(named, shape))
        label = tuple(r for r in shape if r in LABEL_RESIDUE)
        if 'label_asym' in label and len(label) > 1 and any(r in AUTHOR_RESIDUE for r in shape):
            by_label.setdefault(label, set()).add(_key(named, label))

    found = {shape: set() for shape in wanted}
    authors = {label: {} for label in by_label}  # label shape -> label key -> author ids of the atoms
    atom_shapes = [s for s in wanted if _level(s) == 'atom']
    other_shapes = [s for s in wanted if _level(s) != 'atom']
    author_roles = [r for r in AUTHOR_RESIDUE if r in have]
    getters = {}
    for ev in _cif.parse(source.lines(), ATOM_SITE):
        if type(ev) is _cif.Row:
            values = ev.values
            residue = residue_of(values)
            if residue != last:
                last = residue
                for shape in other_shapes:
                    key = tuple(_norm(v) for v in getters[shape](values))
                    if key in wanted[shape]:
                        found[shape].add(key)
                for label in by_label:
                    key = tuple(_norm(v) for v in getters[label](values))
                    if key in by_label[label]:
                        authors[label].setdefault(key, set()).add(tuple(_norm(v) for v in author_of(values)))
            for shape in atom_shapes:
                key = tuple(_norm(v) for v in getters[shape](values))
                if key in wanted[shape]:
                    found[shape].add(key)
        elif type(ev) is _cif.Loop:
            positions = atom_roles(ev)
            if set(positions) != set(have):
                break  # another data block
            for shape in list(wanted) + list(by_label):
                getters[shape] = _getter([positions[r] for r in shape])
            residue_of = _getter([positions[r] for r in RESIDUE_ROLES if r in positions])
            author_of = _getter([positions[r] for r in author_roles])
            last = None
        elif type(ev) is _cif.LoopEnd and ev.loop.category.lower() == ATOM_SITE:
            break
    return have, found, authors, author_roles


def problems(source, before=None):
    """Rows of the listed tables that are out of step with the atoms.

    Returns [(category, kind, message)], kind being 'missing' (it names a
    residue or atom with no atoms, in a table whose residues have atoms) or
    'author ids' (found by its label ids, its author ids differ from the
    atoms'). With `before`, the look() of the file the source was made from,
    also 'removed': it names something that had atoms there and has none
    here. Only the first data block with atoms is looked at.
    """
    seen = look(source)
    if seen is None:
        return []
    have, found, authors, author_roles = seen
    out = []
    for named, table, bound, row in _references(source):
        shape = _shape(named, have)
        if shape is None:
            continue
        where = 'line {}: {} {}'.format(row.lineno, table.category, bound.group.name)
        key = _key(named, shape)
        label = tuple(r for r in shape if r in LABEL_RESIDUE)
        if label in authors:
            atoms = authors[label].get(_key(named, label))
            if atoms:
                mine = tuple(named.get(r) for r in author_roles)
                if not any(all(m is None or m == a for m, a in zip(mine, theirs)) for theirs in atoms):
                    theirs = sorted(atoms)[0]
                    out.append((table.category, 'author ids', '{}: residue {} (label ids) has author ids {}, its atoms {}'.format(
                        where, ' '.join(named[r] or '.' for r in label),
                        ' '.join(m for m in mine if m),
                        ' '.join(a for m, a in zip(mine, theirs) if m is not None and a))))
                    continue
        if key in found[shape]:
            continue
        if before is not None and key in before[1].get(shape, ()):
            out.append((table.category, 'removed', '{}: {} {} has no atoms left'.format(
                where, _level(shape), _describe(named))))
        elif table.observed:
            out.append((table.category, 'missing', '{}: {} {} has no atoms'.format(
                where, _level(shape), _describe(named))))
    return out


def _references(source):
    """(named, table, bound group, row) for each group of each row of the listed tables."""
    for ev in _cif.parse(source.lines(), CATEGORIES):
        if type(ev) is _cif.Row:
            for b in bound:
                named = b.named(ev.values)
                if named:
                    yield named, table, b, ev
        elif type(ev) is _cif.Loop:
            table = TABLES[ev.category.lower()]
            bound = bind(ev)


def _getter(positions):
    if len(positions) == 1:
        p = positions[0]
        return lambda values: (values[p],)
    if not positions:
        return lambda values: ()
    return operator.itemgetter(*positions)


# ---------------------------------------------------------------------------
# Renames
# ---------------------------------------------------------------------------

def scheme_columns(loop, scheme):
    """Where a loop of a listed table names chains and residues in one scheme
    (auth or label), for renames: ([(chain positions, number positions,
    icode position, name position)], [positions of lists of chains])."""
    table = TABLES.get(loop.category.lower())
    if table is None:
        return [], []
    roles = {'auth': ('auth_asym', 'auth_seq', 'auth_comp'), 'label': ('label_asym', 'label_seq', 'label_comp')}[scheme]
    out = []
    for b in bind(loop):
        chains, numbers = b.updates.get(roles[0], []), b.updates.get(roles[1], [])
        if not chains and not numbers:
            continue
        icode = b.updates['icode'][0] if scheme == 'auth' and 'icode' in b.updates else None
        other = 'label_comp' if scheme == 'auth' else 'auth_comp'  # the same name, if the scheme has none
        name = b.updates.get(roles[2], b.updates.get(other, [None]))[0]
        out.append((chains, numbers, icode, name))
    lists = [loop.index(column) for column, role in table.lists.items() if role == roles[0] and column in loop]
    return out, lists


def rename_list(value, chain):
    """A comma-separated list of chains with chain(old) -> new (None: keep) applied."""
    if value in NULLS:
        return value
    names = value.split(',')
    new = [chain(n.strip()) for n in names]
    return ','.join(n if m is None else m for n, m in zip(names, new))


# ---------------------------------------------------------------------------
# Removals
# ---------------------------------------------------------------------------

# Tables of a molecule (entity): their rows go with the molecule.
ENTITY_TABLES = {
    '_entity': 'id', '_entity_poly': 'entity_id', '_entity_poly_seq': 'entity_id',
    '_pdbx_entity_nonpoly': 'entity_id', '_entity_src_gen': 'entity_id', '_entity_src_nat': 'entity_id',
    '_pdbx_entity_src_syn': 'entity_id', '_entity_name_com': 'entity_id', '_entity_name_sys': 'entity_id',
    '_entity_keywords': 'entity_id', '_pdbx_entity_branch': 'entity_id', '_pdbx_entity_branch_list': 'entity_id',
    '_pdbx_entity_branch_link': 'entity_id', '_pdbx_entity_branch_descriptor': 'entity_id',
    '_struct_ref': 'entity_id', '_ma_target_entity': 'entity_id',
}
# Rows that go with a removed row of another table: (table, key columns, child table, its key columns).
CHILDREN = [
    ('_struct_sheet_range', ('sheet_id', 'id'), '_struct_sheet_order', ('sheet_id', 'range_id_1')),
    ('_struct_sheet_range', ('sheet_id', 'id'), '_struct_sheet_order', ('sheet_id', 'range_id_2')),
    ('_struct_sheet_range', ('sheet_id', 'id'), '_pdbx_struct_sheet_hbond', ('sheet_id', 'range_id_1')),
    ('_struct_sheet_range', ('sheet_id', 'id'), '_pdbx_struct_sheet_hbond', ('sheet_id', 'range_id_2')),
    ('_struct_site', ('id',), '_struct_site_gen', ('site_id',)),
    ('_struct_site', ('id',), '_struct_site_keywords', ('site_id',)),
    ('_struct_sheet', ('id',), '_struct_sheet_order', ('sheet_id',)),
    ('_struct_sheet', ('id',), '_pdbx_struct_sheet_hbond', ('sheet_id',)),
    ('_pdbx_struct_assembly', ('id',), '_pdbx_struct_assembly_prop', ('biol_id',)),
    ('_pdbx_struct_assembly', ('id',), '_pdbx_struct_assembly_auth_evidence', ('assembly_id',)),
]
# Rows that go when every row of another table that belongs to them went: (table, key columns, other table, its
# key columns, column counting them or None).
EMPTIED = [
    ('_struct_sheet', ('id',), '_struct_sheet_range', ('sheet_id',), 'number_strands'),
    ('_pdbx_struct_assembly', ('id',), '_pdbx_struct_assembly_gen', ('assembly_id',), None),
]
ONE_LETTER = ('pdbx_seq_one_letter_code', 'pdbx_seq_one_letter_code_can')
RANGE_SHIFT = {'_struct_ref_seq': ('db_align_beg', 'db_align_end')}  # numbers that move with the ends


def _int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


class Removal:
    """What a tool removes, found by a first look at the file, and how the
    rows of the other tables follow.

    `make_owner(loop)` is called once per _atom_site loop and returns a
    function that takes a row's values and gives the atom's owner: True or
    False for a tool that keeps some atoms, the chain for a tool that splits
    the file. It is called on every row in file order. view(owner) then
    tells what goes when only that owner's atoms are kept: a residue goes
    when none of its atoms is kept, a chain or a molecule when none of its
    residues is. Only what had atoms in the input counts: an unobserved
    residue does not go.
    """

    def __init__(self, source, make_owner):
        self.source = source
        self.have = None
        for ev in _cif.parse(source.lines(), ATOM_SITE):
            if type(ev) is _cif.Loop:
                self.have = atom_roles(ev)
                break
        wanted, ranged = {}, set()
        if self.have is not None:
            for named, table, bound, row in _references(source):
                shape = _shape(named, self.have)
                if shape is not None:
                    wanted.setdefault(shape, set()).add(_key(named, shape))
                    if bound.group.part:
                        ranged.add(_chain_of(named))
        self._look(make_owner, wanted, ranged)

    def _look(self, make_owner, wanted, ranged):
        """The first look at the atoms: the owners of what the tables name."""
        sets = {}

        def add(owners, owner):  # one frozenset per combination of owners
            new = owners | {owner} if owners else frozenset((owner,))
            return sets.setdefault(new, new)

        owned = {shape: {} for shape in wanted}  # shape -> {key: owners}
        atom_shapes = [s for s in wanted if _level(s) == 'atom']
        other_shapes = [s for s in wanted if _level(s) != 'atom']
        # an atom is looked up only in a residue that a table names atoms of
        residue_part = {s: tuple(r for r in s if r not in ATOM_ROLES) for s in atom_shapes}
        named_residues = {s: {tuple(v for v, r in zip(k, s) if r not in ATOM_ROLES) for k in wanted[s]}
                          for s in atom_shapes}
        chains = {role: {} for role in CHAIN_ROLES}  # role -> {name: owners}
        entities, positions = {}, {}  # entity, (entity, label number) -> owners
        molecules = {}  # entity -> {owner: label chains}, with None for all owners
        waters = {}  # entity -> {owner: residues of the first model}, with None for all
        residues = {chain: {} for chain in ranged}  # (role, chain) -> {label number: [(owner, ids)]}
        loop = None
        for ev in _cif.parse(self.source.lines(), ATOM_SITE):
            if type(ev) is not _cif.Row:
                if type(ev) is _cif.Loop:
                    owner_of = make_owner(ev)
                continue
            values = ev.values
            owner = owner_of(values)
            if ev.loop is not loop:
                loop = ev.loop
                where = atom_roles(loop)
                getters = {shape: _getter([where[r] for r in shape]) if all(r in where for r in shape)
                           else (lambda values: None) for shape in list(wanted) + list(residue_part.values())}
                residue_of = _getter([where[r] for r in RESIDUE_ROLES if r in where])
                chain_pos = [(r, where[r]) for r in CHAIN_ROLES if r in where]
                id_pos = [(r, where[r]) for r in ('label_seq', 'label_comp', 'auth_seq', 'auth_comp', 'icode')
                          if r in where]
                entity_pos = loop.index('label_entity_id') if 'label_entity_id' in loop else None
                seq_pos, label_pos = where.get('label_seq'), where.get('label_asym')
                model_pos = loop.index('pdbx_PDB_model_num') if 'pdbx_PDB_model_num' in loop else None
                first_model = values[model_pos] if model_pos is not None else None
                last = None
            residue = residue_of(values)
            if residue != last:
                last, marked = residue, set()
                here = []  # (owners by key, key) of the tables' names for this residue
                for shape in other_shapes:
                    key = _normed(getters[shape](values))
                    if key in wanted[shape]:
                        here.append((owned[shape], key))
                active = [s for s in atom_shapes if _normed(getters[residue_part[s]](values)) in named_residues[s]]
                names = [(chains[role], _norm(values[pos])) for role, pos in chain_pos]
                entity = values[entity_pos] if entity_pos is not None else None
                number = _int(values[seq_pos]) if seq_pos is not None else None
                first = model_pos is None or values[model_pos] == first_model
                trims = [residues[(role, name)] for role, name in
                         ((r, _norm(values[p])) for r, p in chain_pos) if (role, name) in residues]
                if entity is not None and first:
                    count = waters.setdefault(entity, {})
                    count[None] = count.get(None, 0) + 1
                if entity is not None and label_pos is not None:
                    molecules.setdefault(entity, {}).setdefault(None, set()).add(values[label_pos])
            if owner not in marked:  # the first atom of this residue with this owner
                marked.add(owner)
                for by_key, key in here:
                    by_key[key] = add(by_key.get(key), owner)
                for by_name, name in names:
                    by_name[name] = add(by_name.get(name), owner)
                if entity is not None:
                    entities[entity] = add(entities.get(entity), owner)
                    if number is not None:
                        positions[(entity, number)] = add(positions.get((entity, number)), owner)
                    if label_pos is not None:
                        molecules[entity].setdefault(owner, set()).add(values[label_pos])
                    if first:
                        waters[entity][owner] = waters[entity].get(owner, 0) + 1
                if number is not None and trims:
                    ids = {role: values[pos] for role, pos in id_pos}
                    for trim in trims:
                        trim.setdefault(number, []).append((owner, ids))
            for shape in active:
                key = _normed(getters[shape](values))
                if key in wanted[shape]:
                    by_key = owned[shape]
                    by_key[key] = add(by_key.get(key), owner)
        self.owned, self.chains, self.entities, self.positions = owned, chains, entities, positions
        self.molecules, self.waters, self.residues = molecules, waters, residues
        self.sets = list(sets)

    def view(self, owner):
        """A View: what goes when only `owner`'s atoms are kept."""
        return self.views([owner])[0]

    def views(self, owners):
        """Views for several owners, with one more reading of the file."""
        views = [View(self, owner) for owner in owners]
        _parents(self, [v for v in views if v.active])
        return views

    def candidates(self, ev):
        """The owners whose view may keep a row (None: any)."""
        category = ev.loop.category.lower()
        out = None
        if category in ENTITY_TABLES:
            entity = ev.get(ENTITY_TABLES[category])
            owners = self.entities.get(entity)
            if category == '_entity_poly_seq':
                owners = self.positions.get((entity, _int(ev.get('num'))), owners)
            out = owners
        table = TABLES.get(category)
        if table is not None:
            for b in bind(ev.loop):
                role, name = _chain_of(b.named(ev.values))
                owners = self.chains[role].get(name) if role else None
                if owners is not None:
                    out = owners if out is None else out & owners
        return out


class View:
    """What goes when only one owner's atoms are kept (Removal.view), and
    how a row follows (row())."""

    def __init__(self, removal, owner):
        self.removal, self.owner, self.have = removal, owner, removal.have
        self.active = any(owner not in owners for owners in removal.sets)
        self.gone, self.left, self.sequences, self._loops, self._residues = {}, {}, {}, {}, {}
        self.categories = FOLLOWING

    def _gone(self, owners):
        return owners is not None and self.owner not in owners

    def row(self, ev):
        """Whether to keep a row of one of `categories` (edited in place when kept)."""
        return self._fate(ev)

    def _fate(self, ev, decide_only=False):
        removal = self.removal
        category = ev.loop.category.lower()
        if category in ENTITY_TABLES:
            entity = ev.get(ENTITY_TABLES[category])
            if self._gone(removal.entities.get(entity)):
                return False
            if category == '_entity_poly_seq' and self._gone(removal.positions.get((entity, _int(ev.get('num'))))):
                return False
            if not decide_only:
                if category == '_entity':
                    self._count(ev, entity)
                elif category == '_entity_poly' and entity in self.sequences:
                    self._one_letter(ev, entity)
        if not decide_only:
            for parent, cols, child, child_cols in CHILDREN:
                if child == category and tuple(ev.get(c) for c in child_cols) in self.gone.get((parent, cols), ()):
                    return False
            for table, cols, other, other_cols, count in EMPTIED:
                if table == category:
                    left = self.left.get((other, other_cols), {}).get(tuple(ev.get(c) for c in cols))
                    if left == 0:
                        return False
                    if left is not None and count and count in ev and _int(ev[count]) not in (None, left):
                        ev[count] = str(left)
        table = TABLES.get(category)
        if table is None:
            return True
        bound = self._loops.get(ev.loop)
        if bound is None:
            bound = self._loops[ev.loop] = (bind(ev.loop), [(ev.loop.index(c), role) for c, role in table.lists.items()
                                                            if c in ev.loop])
        groups, lists = bound
        values = ev.values
        ends = {}
        for b in groups:
            named = b.named(values)
            if not named:
                continue
            role, name = _chain_of(named)
            if self._gone(removal.chains[role].get(name)):
                return False
            shape = _shape(named, self.have)
            if shape is not None and self._gone(removal.owned.get(shape, {}).get(_key(named, shape))):
                if not b.group.part:
                    return False
                ends[b.group.part] = (b, named, True)
            elif b.group.part:
                ends[b.group.part] = (b, named, False)
        if any(removed for _, _, removed in ends.values()) and not self._trim(ev, ends, decide_only):
            return False
        for pos, role in lists:
            value = values[pos]
            if value in NULLS:
                continue
            names = [n for n in value.split(',') if not self._gone(removal.chains[role].get(n.strip()))]
            if not names:
                return False
            if not decide_only and len(names) != len(value.split(',')):
                ev[ev.loop.columns[pos]] = ','.join(names)
        return True

    def _kept_residues(self, chain):
        """[(label number, ids)] of the kept residues of a chain named by a range."""
        if chain not in self._residues:
            by_number = self.removal.residues.get(chain, {})
            self._residues[chain] = [(n, next(ids for o, ids in by_number[n] if o == self.owner))
                                     for n in sorted(by_number) if any(o == self.owner for o, _ in by_number[n])]
        return self._residues[chain]

    def _trim(self, ev, ends, decide_only):
        """Move the removed end(s) of a range to the nearest residue left in it;
        False if there is none."""
        if 'beg' not in ends or 'end' not in ends:
            return False
        (b_bound, b_named, b_removed), (e_bound, e_named, e_removed) = ends['beg'], ends['end']
        beg, end = _int(b_named.get('label_seq')), _int(e_named.get('label_seq'))
        if beg is None or end is None:
            return False
        inside = [(n, ids) for n, ids in self._kept_residues(_chain_of(b_named)) if beg <= n <= end]
        if not inside:
            return False
        new = {}
        if b_removed:
            new['beg'] = inside[0]
        if e_removed:
            new['end'] = inside[-1]
        if (new['beg'][0] if b_removed else beg) > (new['end'][0] if e_removed else end):
            return False
        if decide_only:
            return True
        columns = ev.loop.columns
        shifts = RANGE_SHIFT.get(ev.loop.category.lower(), ())
        for part, (number, ids) in new.items():
            bound = b_bound if part == 'beg' else e_bound
            for role, positions in bound.updates.items():
                if role in ids:
                    for pos in positions:
                        if ev.values[pos] not in NULLS or role == 'icode':
                            ev[columns[pos]] = ids[role]
            delta = number - (beg if part == 'beg' else end)
            if shifts and delta:
                column = shifts[0 if part == 'beg' else 1]
                old = _int(ev.get(column))
                if old is not None:
                    ev[column] = str(old + delta)
        return True

    def _count(self, ev, entity):
        """_entity.pdbx_number_of_molecules of a molecule that lost copies."""
        if 'pdbx_number_of_molecules' not in ev or ev['pdbx_number_of_molecules'] in NULLS:
            return
        if (ev.get('type') or '').lower() == 'water':
            by_owner = self.removal.waters.get(entity, {})
            count, total = by_owner.get(self.owner, 0), by_owner.get(None, 0)
        else:
            by_owner = self.removal.molecules.get(entity, {})
            count, total = len(by_owner.get(self.owner, ())), len(by_owner.get(None, ()))
        if count != total and ev['pdbx_number_of_molecules'] != str(count):
            ev['pdbx_number_of_molecules'] = str(count)

    def _one_letter(self, ev, entity):
        """The one-letter sequences without the removed positions."""
        numbers = self.sequences[entity]
        for column in ONE_LETTER:
            value = ev.get(column)
            if value is None or value in NULLS:
                continue
            letters = re.findall(r'\([^)]*\)|[^\s()]', value)
            if len(letters) != len(numbers):
                continue  # not one letter per position: left as it is
            positions = self.removal.positions
            kept = ''.join(l for l, n in zip(letters, numbers) if not self._gone(positions.get((entity, n))))
            ev[column] = '\n'.join(kept[i:i + 80] for i in range(0, len(kept), 80)) or '?'


def _parents(removal, views):
    """For each view, before writing: which rows of the parent tables of
    CHILDREN go, and how many rows of the tables counted by EMPTIED are left
    (they may come before the rows that decide it); and the positions of the
    molecules that lose some (for the one-letter codes). One reading of the
    file for all the views."""
    if not views:
        return
    parents = {(p, cols) for p, cols, _, _ in CHILDREN}
    counted = {(c, cols) for _, _, c, cols, _ in EMPTIED}
    tables = {t for t, _ in parents | counted}
    entities = {e for e, n in removal.positions}
    if any(view._gone(owners) for view in views for owners in set(removal.positions.values())):
        tables.add('_entity_poly_seq')
    numbers = {}  # entity -> its positions, in order
    for ev in _cif.parse(removal.source.lines(), tuple(tables)):
        if type(ev) is not _cif.Row:
            continue
        category = ev.loop.category.lower()
        if category == '_entity_poly_seq':
            entity, number = ev.get('entity_id'), _int(ev.get('num'))
            if entity in entities and number is not None:
                here = numbers.setdefault(entity, [])
                if not here or here[-1] != number:  # one per position (two residue types at one place)
                    here.append(number)
            continue
        candidates = removal.candidates(ev)
        for view in views:
            fate = (candidates is None or view.owner in candidates) and view._fate(ev, decide_only=True)
            for table, cols in parents:
                if table == category and not fate:
                    view.gone.setdefault((table, cols), set()).add(tuple(ev.get(c) for c in cols))
            for table, cols in counted:
                if table == category:
                    left = view.left.setdefault((table, cols), {})
                    key = tuple(ev.get(c) for c in cols)
                    left[key] = left.get(key, 0) + (1 if fate else 0)
    for view in views:
        for table, cols, other, other_cols, _ in EMPTIED:  # emptied rows take their children with them
            for key, left in view.left.get((other, other_cols), {}).items():
                if not left:
                    view.gone.setdefault((table, cols), set()).add(key)
        for entity, here in numbers.items():
            if any(view._gone(removal.positions.get((entity, n))) for n in here):
                view.sequences[entity] = here


FOLLOWING = tuple(dict.fromkeys(
    list(TABLES) + list(ENTITY_TABLES) + [t for parent, _, child, _ in CHILDREN for t in (parent, child)] +
    [t for parent, _, child, _, _ in EMPTIED for t in (parent, child)]))


def _chain_of(named):
    """(role, name) of the chain a group names: its label chain, else its author chain."""
    for role in CHAIN_ROLES:
        if role in named:
            return role, named[role]
    return None, None


def _normed(values):
    return None if values is None else tuple('' if v in NULLS else v for v in values)


# ---------------------------------------------------------------------------
# Names
# ---------------------------------------------------------------------------

# Residue names outside the listed tables that follow a renamed residue: the sequences of the molecules.
SEQUENCE_NAMES = {'_entity_poly_seq': ['mon_id'], '_pdbx_entity_nonpoly': ['comp_id'],
                  '_pdbx_entity_branch_list': ['comp_id']}


def residue_name_columns():
    """{table: [columns]} of the residue names of the listed tables and of the sequences."""
    out = {ATOM_SITE: ['auth_comp_id', 'label_comp_id']}
    for category, table in TABLES.items():
        columns = [c for g in table.groups for role in ('label_comp', 'auth_comp') for c in g.columns.get(role, ())]
        if columns:
            out[category] = list(dict.fromkeys(columns))
    out.update(SEQUENCE_NAMES)
    return out


class AtomNames:
    """New names of renamed atoms, for the rows of the listed tables that
    name them. `renamed` is [(atom values, new name)] of an _atom_site loop."""

    def __init__(self, loop, renamed):
        self.where = atom_roles(loop)
        self.renamed = renamed
        self._by_shape = {}

    def follow(self, ev):
        """Rename the atoms a row of a listed table names (in place)."""
        values = ev.values  # as read: the groups are looked up before any change
        for b in bind(ev.loop):
            if not any(role in b.roles for role in ATOM_ROLES):
                continue
            named = b.named(values)
            shape = _shape(named, self.where)
            if shape is None or not any(role in shape for role in ATOM_ROLES):
                continue
            if shape not in self._by_shape:
                get = _getter([self.where[r] for r in shape])
                self._by_shape[shape] = {_normed(get(v)): new for v, new in self.renamed}
            new = self._by_shape[shape].get(_key(named, shape))
            if new is None:
                continue
            for role in ATOM_ROLES:
                for pos in b.updates.get(role, ()):
                    if ev.values[pos] not in NULLS and ev.values[pos] != new:
                        ev[ev.loop.columns[pos]] = new
