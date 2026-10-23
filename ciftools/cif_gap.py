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
Finds gaps between consecutive protein residues in an mmCIF file.

Detects gaps both by a distance criterion (consecutive alpha carbons more
than 4.0 A apart) and by discontinuous residue numbering, within each chain
and model. Only applies to protein residues (ATOM records), as pdb_gap. The
output lines are those of pdb_gap; unlike pdb_gap, 52 followed by 52A is
not a gap, and the second alpha carbon of a residue with alternate
locations is not compared with the first.

Chains and residue numbers are auth identifiers (label ones with --label).

Usage:
    cif_gap <cif file>

Example:
    cif_gap 1abc.cif

This program is part of the cif-tools suite, a port of pdb-tools to mmCIF.
"""

import sys

from ciftools import _atoms, _cif, _cli

NULLS = ('.', '?')
MAX_CA_DISTANCE = 4.0
FMT_DISTANCE = '{0}:{1}{2} < {3:7.2f}A > {4}:{5}{6}\n'
FMT_SEQUENCE = '{0}:{1}{2} < Seq. Gap > {3}:{4}{5}\n'


def check_input(args):
    """Return (file handle, identifier scheme)."""
    fhandle, _, scheme = _cli.parse_args(args, __doc__, takes_option=False)
    return fhandle, scheme


def run(fhandle, scheme=None):
    """Yield the report lines."""
    gaps = 0
    prev = None  # (model, chain, resname, number, icode, x, y, z)
    for ev in _cif.parse(fhandle):
        if type(ev) is _cif.Row:
            values = ev.values
            if (values[f.atom] != 'CA' or (f.group is not None and values[f.group] != 'ATOM')
                    or (f.element is not None and values[f.element].upper() != 'C')):
                continue
            number = _atoms.as_int(values[f.resnum])
            if number is None:
                continue
            model = values[f.model] if f.model is not None else None
            icode = values[f.icode] if f.icode is not None else '?'
            icode = '' if icode in NULLS else icode
            here = (model, values[f.chain], values[f.resname], number, icode,
                    float(values[x]), float(values[y]), float(values[z]))
            if prev is not None and prev[:2] == here[:2]:
                if prev[3:5] == here[3:5]:
                    continue  # another alpha carbon of the same residue (altloc)
                distance = sum((a - b) ** 2 for a, b in zip(prev[5:], here[5:])) ** 0.5
                left = (prev[1], prev[2], '{}{}'.format(prev[3], prev[4]))
                right = (here[1], here[2], '{}{}'.format(here[3], here[4]))
                if distance > MAX_CA_DISTANCE:
                    gaps += 1
                    yield FMT_DISTANCE.format(*(left + (distance,) + right))
                elif here[3] not in (prev[3], prev[3] + 1):
                    gaps += 1
                    yield FMT_SEQUENCE.format(*(left + right))
            prev = here
        elif type(ev) is _cif.Loop:
            f = _atoms.Fields(ev, scheme)
            f.need('atom', 'chain', 'resnum', 'resname')
            x, y, z = (ev.index(c) for c in ('Cartn_x', 'Cartn_y', 'Cartn_z'))
            prev = None
    yield 'Found {} gap(s) in the structure\n'.format(gaps)


def main():
    fhandle, scheme = check_input(sys.argv[1:])
    _cli.output_text(_cli.lazy(run, fhandle, scheme), fhandle)


if __name__ == '__main__':
    main()
