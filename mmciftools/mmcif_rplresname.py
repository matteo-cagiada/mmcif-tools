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
Replaces one residue name by another in an mmCIF file.

Residue names are the same chemical identity in both identifier schemes,
and programs differ in which one they read, so both auth_comp_id and
label_comp_id are changed (whatever --auth/--label says), together with the
residue names of the other tables that name residues (_struct_conn,
_pdbx_poly_seq_scheme, _struct_conf, ...) and of the sequences
(_entity_poly_seq, _pdbx_entity_nonpoly), and of the definitions of the
residues (_chem_comp, _chem_comp_atom, ...): the atoms keep their
definition under the new name. If the new name is already defined, the
old name's definition is removed instead (the new name keeps its own).
Names are matched exactly (case-sensitive) and may be longer than three
characters.

Usage:
    mmcif_rplresname -<from>:<to> <cif file>

Example:
    mmcif_rplresname -HIP:HIS 1abc.cif   # changes all HIP residues to HIS

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cif, _cli, _refs

COLUMNS = _refs.residue_name_columns()
# The definitions of the residues, by the column naming the residue.
DEFINITIONS = {'_chem_comp': 'id', '_chem_comp_atom': 'comp_id', '_chem_comp_bond': 'comp_id',
               '_pdbx_chem_comp_descriptor': 'comp_id', '_pdbx_chem_comp_identifier': 'comp_id',
               '_pdbx_chem_comp_synonyms': 'comp_id', '_pdbx_chem_comp_feature': 'comp_id',
               '_pdbx_chem_comp_audit': 'comp_id', '_pdbx_chem_comp_pcm': 'comp_id'}


def check_input(args):
    """Return (file handle, (from, to), identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__, needs_option=True)
    if option.count(':') != 1:
        _cli.fail(__doc__, 'Invalid option value: {!r}'.format(option))
    old, new = option.split(':')
    if not old or not new or len(old.split()) != 1 or len(new.split()) != 1 \
            or old != old.strip() or new != new.strip():
        _cli.fail(__doc__, 'Residue names must be given as -<from>:<to>: {!r}'.format(option))
    return fhandle, (old, new), scheme


def _defined(source, name):
    """True if _chem_comp defines a residue `name`."""
    for ev in _cif.parse(source.lines(), ('_chem_comp',)):
        if type(ev) is _cif.Row and ev.get('id') == name:
            return True
    return False


def run(fhandle, names, scheme=None):
    """Yield the events of the file with residue `from` renamed to `to`."""
    old, new = names
    source = _atoms.as_source(fhandle)
    drop = old != new and _defined(source, new)  # else the definition of `old` is renamed
    columns = dict(COLUMNS, **{c: [] if drop else [k] for c, k in DEFINITIONS.items()})
    for ev in _atoms.replace_values(source.lines(), columns, {old: new}):
        if drop and type(ev) is _cif.Row:
            key = DEFINITIONS.get(ev.loop.category.lower())
            if key is not None and ev.get(key) == old:
                continue  # the new name keeps its own definition
        yield ev


def main():
    fhandle, names, scheme = check_input(sys.argv[1:])
    _cli.output(run(fhandle, names, scheme), fhandle)


if __name__ == '__main__':
    main()
