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
Deletes a range of residues from an mmCIF file.

The range option has three components: start, end, and step. Start and end
are optional and if omitted the range will start at the first residue or
end at the last, respectively. As in pdb_delres, the start and end values
are residue numbers (matched in every chain and model), while the step
counts residue occurrences in file order, regardless of their number.
Unlike pdb_delres, a range may start and end at the same residue (-5:5).

Residue numbers are auth_seq_id, or label_seq_id with --label. Rows of
_atom_site_anisotrop and _struct_conn that refer to removed atoms are
removed too.

Usage:
    mmcif_delres -[resid]:[resid]:[step] [--label] <cif file>

Example:
    mmcif_delres -1:10 1abc.cif     # deletes residues 1 to 10
    mmcif_delres -1: 1abc.cif       # deletes residues 1 to END
    mmcif_delres -:5 1abc.cif       # deletes residues from START to 5
    mmcif_delres -::5 1abc.cif      # deletes every 5th residue
    mmcif_delres -1:10:5 1abc.cif   # deletes every 5th residue from 1 to 10

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cli


def check_input(args):
    """Return (file handle, (start, end, step), identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__, needs_option=True)
    if ':' not in option:
        _cli.fail(__doc__, _atoms.RANGE_HELP)
    try:
        resrange = _atoms.parse_range(option)
    except ValueError as err:
        _cli.fail(__doc__, err)
    return fhandle, resrange, scheme


def run(fhandle, resrange, scheme=None):
    """Yield the events of the file without the residues in the range."""
    start, end, step = resrange

    def make_keep(fields):
        fields.need('chain', 'resnum')
        chain, resnum = fields.chain, fields.resnum
        prev, count = None, -1

        def keep(values):
            nonlocal prev, count
            key = (values[chain], values[resnum])
            if key != prev:
                prev = key
                count += 1
            num = _atoms.as_int(values[resnum])
            return num is None or not (start <= num <= end and count % step == 0)
        return keep

    return _atoms.filter_atoms(fhandle, make_keep, scheme)


def main():
    fhandle, resrange, scheme = check_input(sys.argv[1:])
    _cli.output(run(fhandle, resrange, scheme), fhandle)


if __name__ == '__main__':
    main()
