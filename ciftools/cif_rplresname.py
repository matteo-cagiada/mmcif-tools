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
Replaces one residue name by another in an mmCIF file.

Residue names are the same chemical identity in both identifier schemes,
and programs differ in which one they read, so both auth_comp_id and
label_comp_id are changed (whatever --auth/--label says), together with the
residue names in _atom_site_anisotrop and _struct_conn. Names are matched
exactly (case-sensitive) and may be longer than three characters. Other
categories that list residue names (_entity_poly_seq, _pdbx_poly_seq_scheme,
_chem_comp, ...) are not changed.

Usage:
    cif_rplresname -<from>:<to> <cif file>

Example:
    cif_rplresname -HIP:HIS 1abc.cif   # changes all HIP residues to HIS

This program is part of the cif-tools suite, a port of pdb-tools to mmCIF.
"""

import sys

from ciftools import _atoms, _cli

COLUMNS = {
    _atoms.ATOM_SITE: ['auth_comp_id', 'label_comp_id'],
    _atoms.ANISOTROP: ['pdbx_auth_comp_id', 'pdbx_label_comp_id'],
    _atoms.STRUCT_CONN: ['ptnr1_auth_comp_id', 'ptnr1_label_comp_id',
                         'ptnr2_auth_comp_id', 'ptnr2_label_comp_id'],
}


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


def run(fhandle, names, scheme=None):
    """Yield the events of the file with residue `from` renamed to `to`."""
    old, new = names
    return _atoms.replace_values(fhandle, COLUMNS, {old: new})


def main():
    fhandle, names, scheme = check_input(sys.argv[1:])
    _cli.output(run(fhandle, names, scheme), fhandle)


if __name__ == '__main__':
    main()
