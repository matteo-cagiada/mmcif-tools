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
Selects atoms by entity: by entity number (label_entity_id) or by type.

In mmCIF every molecule belongs to an entity (_entity), whose type is one
of:

    polymer       protein and nucleic-acid chains
    non-polymer   ligands and ions
    branched      oligosaccharides (glycans)
    macrolide     macrolides
    water         water

Unlike mmcif_selhetatm, this keeps modified residues with their chain: they
are HETATM records inside a polymer (MSE in proteins, PSU and 2MG in tRNA,
...). Numbers and types can be mixed (-polymer,5); types are
case-insensitive. Files without an _entity category (some non-PDB files) are
classified from their atoms: atoms with a label_seq_id are polymer, HOH,
WAT and DOD are water, and the rest non-polymer. Rows of
_atom_site_anisotrop and _struct_conn of removed atoms are removed too.

Usage:
    mmcif_selentity -<entity id or type>[,...] <cif file>

Example:
    mmcif_selentity -polymer 1abc.cif           # chains only, with modified residues
    mmcif_selentity -polymer,branched 1abc.cif  # chains and glycans
    mmcif_selentity -1,3 1abc.cif               # entities 1 and 3

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cif, _cli

TYPES = ('polymer', 'non-polymer', 'branched', 'macrolide', 'water')
WATERS = frozenset(('HOH', 'WAT', 'DOD', 'H2O'))


def parse_entities(option, doc):
    """'-polymer,3' -> ({'3'}, {'polymer'}) (fails with the usage when empty)."""
    ids, types = set(), set()
    for item in _atoms.parse_ids(option):
        if item.lower() in TYPES:
            types.add(item.lower())
        else:
            ids.add(item)
    if not ids and not types:
        _cli.fail(doc, 'You must provide at least one entity number or type')
    return ids, types


def check_input(args):
    """Return (file handle, entity ids, entity types)."""
    fhandle, option, _ = _cli.parse_args(args, __doc__, needs_option=True)
    ids, types = parse_entities(option, __doc__)
    return fhandle, ids, types


def entity_types(source):
    """{entity id: type} from the _entity category."""
    types = {}
    for ev in _cif.parse(source.lines(), '_entity'):
        if type(ev) is _cif.Row and 'id' in ev and 'type' in ev:
            types[ev['id']] = ev['type'].lower()
    return types


def matcher(source, ids, types):
    """make_keep-style function: Fields -> function(values) -> True if the
    atom belongs to one of the entities."""
    known = entity_types(source) if types else {}

    def make_match(fields):
        loop = fields.loop
        entity = loop.index('label_entity_id') if 'label_entity_id' in loop else None
        if ids and entity is None:
            raise _cif.CIFError('_atom_site has no label_entity_id column (needed to select '
                                'entities by number)', loop.lineno)
        seq = loop.index('label_seq_id') if 'label_seq_id' in loop else None
        resname = fields.resname

        def kind(values):
            found = known.get(values[entity]) if entity is not None else None
            if found:
                return found
            if seq is not None and _atoms.as_int(values[seq]) is not None:
                return 'polymer'
            if resname is not None and values[resname].upper() in WATERS:
                return 'water'
            return 'non-polymer'

        def match(values):
            if entity is not None and values[entity] in ids:
                return True
            return bool(types) and kind(values) in types
        return match
    return make_match


def run(fhandle, ids, types, scheme=None):
    """Yield the events of the file with only the atoms of the entities."""
    source = _atoms.as_source(fhandle)
    return _atoms.filter_atoms(source, matcher(source, ids, types), scheme)


def main():
    fhandle, ids, types = check_input(sys.argv[1:])
    _cli.output(_cli.lazy(run, fhandle, ids, types), fhandle)


if __name__ == '__main__':
    main()
