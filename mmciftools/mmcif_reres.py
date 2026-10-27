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
Renumbers the residues of an mmCIF file starting from a given number
(default 1).

As in pdb_reres, numbering continues across chains and restarts in each
model; a new residue starts whenever the residue name, chain, number or
insertion code changes from one atom to the next (so 52 and 52A become two
numbers), and insertion codes are kept. There is no upper limit.

Changes auth_seq_id by default, in _atom_site and in the matching columns of
_atom_site_anisotrop and _struct_conn. With --label it changes label_seq_id
of polymer residues instead; that number links atoms to the sequence in
_entity_poly_seq, so a warning is printed.

Usage:
    mmcif_reres -<number> [--label] <cif file>

Example:
    mmcif_reres -10 1abc.cif   # renumbers from 10
    mmcif_reres --1 1abc.cif   # renumbers from -1

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cif, _cli

NULLS = ('.', '?')
LABEL_WARNING = ('WARNING: renumbering label_seq_id breaks the link between atoms and the '
                 'sequence in _entity_poly_seq.\n')


def check_input(args):
    """Return (file handle, first number, identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__)
    try:
        start = 1 if option is None else int(option)
    except ValueError:
        _cli.fail(__doc__, 'You provided an invalid residue number: {!r}'.format(option))
    return fhandle, start, scheme


def _numberer(fields, start):
    """A function giving each atom row its new residue number (None: keep)."""
    fields.need('resnum')
    model, resnum = fields.model, fields.resnum
    uid_of = _atoms.getter(fields.resname, fields.chain, fields.resnum, fields.icode)
    label = fields.scheme == 'label'
    prev_model, prev, count = object(), object(), start - 1

    def number(values):
        nonlocal prev_model, prev, count
        if model is not None and values[model] != prev_model:
            prev_model, prev, count = values[model], object(), start - 1
        if label and values[resnum] in NULLS:
            return None
        uid = uid_of(values)
        if uid != prev:
            prev = uid
            count += 1
        return str(count)
    return number


def _reference_key(chain, resnum, icode, resname=None):
    return chain, resnum, None if icode in NULLS else icode, resname


def run(fhandle, start, scheme=None):
    """Yield the events of the file with residues renumbered from start."""
    source = _atoms.as_source(fhandle)
    scheme = _cif.get_scheme(scheme)

    # First look: new numbers of the residues of the first model, for the
    # rows of _struct_conn (before _atom_site) and _atom_site_anisotrop.
    new_numbers = {}
    first_model = None
    for ev in _cif.parse(source.lines(), _atoms.ATOM_SITE):
        if type(ev) is _cif.Row:
            values = ev.values
            if fields.model is not None:
                if first_model is None:
                    first_model = values[fields.model]
                elif values[fields.model] != first_model:
                    break
            num = number(values)
            if num is not None:
                chain, resnum, icode, resname = where(values)
                # Keyed with and without the residue name: a residue modelled
                # as two residue types gets two numbers (as in pdb_reres).
                new_numbers[_reference_key(chain, resnum, icode, resname)] = num
                new_numbers[_reference_key(chain, resnum, icode)] = num
        elif type(ev) is _cif.Loop:
            fields = _atoms.Fields(ev, scheme)
            number = _numberer(fields, start)
            where = _atoms.getter(fields.chain, fields.resnum, fields.icode, fields.resname)

    def ref_resnum(chain, resnum, icode, resname):
        num = new_numbers.get(_reference_key(chain, resnum, icode, resname))
        return num if num is not None else new_numbers.get(_reference_key(chain, resnum, icode))

    # Anisotropic rows follow their atom: models may hold different residues,
    # and then the numbers of model 1 do not fit the others.
    has_aniso = any(type(ev) is _cif.Loop for ev in _cif.parse(source.lines(), _atoms.ANISOTROP))
    return _atoms.rename(source.lines(), scheme,
                         atom_resnum=lambda fields: _numberer(fields, start),
                         ref_resnum=ref_resnum, aniso_by_id=has_aniso)


def main():
    fhandle, start, scheme = check_input(sys.argv[1:])
    if scheme == 'label':
        sys.stderr.write(LABEL_WARNING)
    _cli.output(run(fhandle, start, scheme), fhandle)


if __name__ == '__main__':
    main()
