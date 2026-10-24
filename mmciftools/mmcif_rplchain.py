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
Replaces one chain identifier by another in an mmCIF file.

Changes auth_asym_id (label_asym_id with --label), in _atom_site and in the
matching columns of _atom_site_anisotrop and _struct_conn. Identifiers are
case-sensitive and may be longer than one character; mmCIF has no blank
chain, so both must be given.

Usage:
    mmcif_rplchain -<from>:<to> [--label] <cif file>

Example:
    mmcif_rplchain -A:B 1abc.cif   # replaces chain A by chain B

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cli


def check_input(args):
    """Return (file handle, (from, to), identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__, needs_option=True)
    if option.count(':') != 1:
        _cli.fail(__doc__, 'Invalid option value: {!r}'.format(option))
    old, new = (part.strip() for part in option.split(':'))
    if not old or not new:
        _cli.fail(__doc__, 'You must provide both chain identifiers, e.g. -A:B')
    if len(old.split()) != 1 or len(new.split()) != 1:
        _cli.fail(__doc__, 'Invalid chain identifier in {!r}'.format(option))
    return fhandle, (old, new), scheme


def run(fhandle, chains, scheme=None):
    """Yield the events of the file with chain `from` renamed to `to`."""
    old, new = chains
    return _atoms.rename(fhandle, scheme, chain=lambda chain: new if chain == old else None)


def main():
    fhandle, chains, scheme = check_input(sys.argv[1:])
    _cli.output(run(fhandle, chains, scheme), fhandle)


if __name__ == '__main__':
    main()
