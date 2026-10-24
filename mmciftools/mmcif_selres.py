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
Selects residues by their number, piecewise or in a range.

The range option has three components: start, end, and step. Start and end
are optional and if omitted the range will start at the first residue or
end at the last, respectively. As in pdb_selres, residues are matched by
number only (in every chain and model, with or without insertion code), and
the step counts residue occurrences in file order.

Residue numbers are auth_seq_id, or label_seq_id with --label. Rows of
_atom_site_anisotrop and _struct_conn that refer to removed atoms are
removed too.

Usage:
    mmcif_selres -[resid]:[resid]:[step] [--label] <cif file>

Example:
    mmcif_selres -1,2,4,6 1abc.cif     # extracts residues 1, 2, 4 and 6
    mmcif_selres -1:10 1abc.cif        # extracts residues 1 to 10
    mmcif_selres -1:10,20:30 1abc.cif  # extracts residues 1 to 10 and 20 to 30
    mmcif_selres -1: 1abc.cif          # extracts residues 1 to END
    mmcif_selres -:5 1abc.cif          # extracts residues from START to 5
    mmcif_selres -::5 1abc.cif         # extracts every 5th residue
    mmcif_selres -1:10:5 1abc.cif      # extracts every 5th residue from 1 to 10
    mmcif_selres --5:5 1abc.cif        # extracts residues -5 to 5

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cli


def parse_selection(option):
    """'1,4:6,10::2' -> (set of single numbers, list of (start, end, step))."""
    singles, ranges = set(), []
    for entry in option.split(','):
        if ':' in entry:
            ranges.append(_atoms.parse_range(entry))
        else:
            singles.add(_atoms.parse_number(entry))
    return singles, ranges


def check_input(args):
    """Return (file handle, single numbers, ranges, identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__)
    try:
        singles, ranges = parse_selection('::' if option is None else option)
    except ValueError as err:
        _cli.fail(__doc__, err)
    return fhandle, singles, ranges, scheme


def run(fhandle, singles, ranges, scheme=None):
    """Yield the events of the file, keeping only the selected residues."""
    source = _atoms.as_source(fhandle)
    selected = set(singles)
    plain = [(start, end) for start, end, step in ranges if step == 1]
    stepped = [r for r in ranges if r[2] != 1]
    if stepped:
        numbers = _atoms.residue_numbers(source, scheme)
        for start, end, step in stepped:
            selected.update([n for n in numbers if start <= n <= end][::step])

    def make_keep(fields):
        fields.need('resnum')
        resnum = fields.resnum

        def keep(values):
            num = _atoms.as_int(values[resnum])
            if num is None:
                return False
            if num in selected:
                return True
            for start, end in plain:
                if start <= num <= end:
                    return True
            return False
        return keep

    return _atoms.filter_atoms(source, make_keep, scheme)


def main():
    fhandle, singles, ranges, scheme = check_input(sys.argv[1:])
    _cli.output(run(fhandle, singles, ranges, scheme), fhandle)


if __name__ == '__main__':
    main()
