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
Shifts the residue numbers of an mmCIF file by a constant value.

Changes auth_seq_id by default, in _atom_site and in the matching columns of
_atom_site_anisotrop and _struct_conn; insertion codes are kept. With
--label it shifts label_seq_id of polymer residues instead; that number
links atoms to the sequence in _entity_poly_seq, so a warning is printed.
There is no upper or lower limit.

Usage:
    cif_shiftres -<number> [--label] <cif file>

Example:
    cif_shiftres -10 1abc.cif   # adds 10 to the original numbering
    cif_shiftres --5 1abc.cif   # subtracts 5 from the original numbering

This program is part of the cif-tools suite, a port of pdb-tools to mmCIF.
"""

import sys

from ciftools import _atoms, _cli

LABEL_WARNING = ('WARNING: shifting label_seq_id breaks the link between atoms and the '
                 'sequence in _entity_poly_seq.\n')


def check_input(args):
    """Return (file handle, shift, identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__)
    try:
        shift = 0 if option is None else int(option)
    except ValueError:
        _cli.fail(__doc__, 'You provided an invalid number: {!r}'.format(option))
    return fhandle, shift, scheme


def run(fhandle, shift, scheme=None):
    """Yield the events of the file with residue numbers shifted by shift."""
    def shifted(resnum):
        num = _atoms.as_int(resnum)
        return None if num is None else str(num + shift)

    def atom_resnum(fields):
        resnum = fields.resnum
        return lambda values: shifted(values[resnum])

    return _atoms.rename(fhandle, scheme, atom_resnum=atom_resnum,
                         ref_resnum=lambda chain, resnum, icode, resname: shifted(resnum))


def main():
    fhandle, shift, scheme = check_input(sys.argv[1:])
    if scheme == 'label':
        sys.stderr.write(LABEL_WARNING)
    _cli.output(run(fhandle, shift, scheme), fhandle)


if __name__ == '__main__':
    main()
