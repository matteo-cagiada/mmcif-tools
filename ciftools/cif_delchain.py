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
Deletes all atoms of one or more chains from an mmCIF file.

Chains are matched on auth_asym_id (the chain shown by viewers), or on
label_asym_id with --label. Chain identifiers are case-sensitive (unlike
pdb_delchain, which upper-cases them) and may be longer than one character.
Rows of _atom_site_anisotrop and _struct_conn that refer to removed atoms
are removed too.

Usage:
    cif_delchain -<chain id>[,<chain id>...] [--label] <cif file>

Example:
    cif_delchain -A 1abc.cif     # removes chain A
    cif_delchain -A,B 1abc.cif   # removes chains A and B

This program is part of the cif-tools suite, a port of pdb-tools to mmCIF.
"""

import sys

from ciftools import _atoms, _cli


def check_input(args):
    """Return (file handle, chain ids, identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__, needs_option=True)
    chains = _atoms.parse_ids(option)
    if not chains:
        _cli.fail(__doc__, 'You must provide at least one chain identifier')
    return fhandle, chains, scheme


def run(fhandle, chains, scheme=None):
    """Yield the events of the file without the atoms of the given chains."""
    chains = set(chains)

    def make_keep(fields):
        fields.need('chain')
        chain = fields.chain
        return lambda values: values[chain] not in chains

    return _atoms.filter_atoms(fhandle, make_keep, scheme)


def main():
    fhandle, chains, scheme = check_input(sys.argv[1:])
    _cli.output(run(fhandle, chains, scheme), fhandle)


if __name__ == '__main__':
    main()
