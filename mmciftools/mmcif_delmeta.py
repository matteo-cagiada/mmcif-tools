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
    _entity, _entity_poly,     the molecules and their sequences
    _entity_poly_seq
    _struct_asym               the chains (label_asym_id) of each molecule
    _chem_comp, _atom_type     the chemical components and elements

and deletes every other category. These are what AlphaFold and ColabFold
read from template files, for example. The kept tables are brought in step
with the atoms: rows of chains (_struct_asym) and entities that have no
atoms left are removed, and _entity.pdbx_number_of_molecules and
_entity_poly.pdbx_strand_id are recomputed. So run it last in a pipeline,
after the edits.

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
        '_pdbx_audit_revision_history', '_entity', '_entity_poly', '_entity_poly_seq', '_struct_asym',
        '_chem_comp', '_atom_type')
BY_ENTITY = {'_entity': 'id', '_entity_poly': 'entity_id', '_entity_poly_seq': 'entity_id'}
NULLS = ('.', '?')


def check_input(args):
    """Return (file handle, delete all metadata, identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__)
    if option not in (None, 'all'):
        _cli.fail(__doc__, 'Invalid option: {!r} (the only option is -all)'.format('-' + option))
    return fhandle, option == 'all', scheme


class _Present:
    """The chains and entities that have atoms in one data block."""

    def __init__(self):
        self.asyms = None     # label_asym_id values (None: no such column)
        self.entities = None  # entity id -> set of label_asym_id
        self.residues = {}    # entity id -> residues of the first model (for waters)
        self.strands = {}     # entity id -> author chains, in order


def present(source):
    """[_Present] for each data block, from _atom_site."""
    blocks = []
    for ev in _cif.parse(source.lines(), _atoms.ATOM_SITE):
        etype = type(ev)
        if etype is _cif.Row:
            v = ev.values
            asym = v[asym_pos] if asym_pos is not None else None
            if asym_pos is not None:
                info.asyms.add(asym)
            if entity_pos is not None:
                entity = v[entity_pos]
                info.entities.setdefault(entity, set()).add(asym)
                if model_pos is None or v[model_pos] == first_model[0]:
                    info.residues.setdefault(entity, set()).add(residue(v))
                if auth_pos is not None:
                    strands = info.strands.setdefault(entity, [])
                    if v[auth_pos] not in strands:
                        strands.append(v[auth_pos])
            if model_pos is not None and first_model[0] is None:
                first_model[0] = v[model_pos]
        elif etype is _cif.DataBlock:
            blocks.append(_Present())
        elif etype is _cif.Loop:
            if not blocks:
                blocks.append(_Present())
            info = blocks[-1]
            pos = {c: ev.index(c) if c in ev else None for c in (
                'label_asym_id', 'label_entity_id', 'auth_asym_id', 'pdbx_PDB_model_num')}
            asym_pos, entity_pos = pos['label_asym_id'], pos['label_entity_id']
            auth_pos, model_pos = pos['auth_asym_id'], pos['pdbx_PDB_model_num']
            info.asyms = set() if asym_pos is not None else None
            info.entities = {} if entity_pos is not None else None
            fields = _atoms.Fields(ev, 'auth')
            residue = _atoms.getter(asym_pos, fields.chain, fields.resnum, fields.icode)
            first_model = [None]
    return blocks


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
    block = -1
    for ev in _cif.parse(source.lines(), COORDINATES + KEPT):
        etype = type(ev)
        if etype is _cif.Verbatim:
            continue
        if etype is _cif.DataBlock:
            block += 1
        elif etype is _cif.Row:
            info = blocks[block] if 0 <= block < len(blocks) else None
            category = ev.loop.category.lower()
            if info is not None and category == '_struct_asym' and info.asyms is not None:
                if ev.get('id') not in info.asyms:
                    continue  # a chain with no atoms left
            elif info is not None and category in BY_ENTITY and info.entities is not None:
                entity = ev.get(BY_ENTITY[category])
                if entity not in info.entities:
                    continue  # a molecule with no atoms left
                if category == '_entity' and 'pdbx_number_of_molecules' in ev:
                    water = (ev.get('type') or '').lower() == 'water'
                    count = str(len(info.residues.get(entity, ())) if water else len(info.entities[entity]))
                    if ev['pdbx_number_of_molecules'] not in NULLS and ev['pdbx_number_of_molecules'] != count:
                        ev['pdbx_number_of_molecules'] = count
                if category == '_entity_poly' and 'pdbx_strand_id' in ev and info.strands.get(entity):
                    strands = ','.join(info.strands[entity])
                    if ev['pdbx_strand_id'] not in NULLS and ev['pdbx_strand_id'] != strands:
                        ev['pdbx_strand_id'] = strands
        yield ev


def main():
    fhandle, everything, scheme = check_input(sys.argv[1:])
    _cli.output(_cli.lazy(run, fhandle, everything, scheme), fhandle)


if __name__ == '__main__':
    main()
