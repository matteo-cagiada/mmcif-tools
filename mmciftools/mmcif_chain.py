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
Sets the chain identifier of all atoms in an mmCIF file.

Changes auth_asym_id (label_asym_id with --label), in _atom_site and in the
matching columns of _atom_site_anisotrop and _struct_conn. mmCIF has no
blank chain, so unlike pdb_chain an identifier must be given; it may be
longer than one character.

Usage:
    mmcif_chain -<chain id> [--label] <cif file>

Example:
    mmcif_chain -C 1abc.cif

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cli


def check_input(args):
    """Return (file handle, chain id, identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__, needs_option=True)
    if not option or len(option.split()) != 1 or option != option.strip():
        _cli.fail(__doc__, 'Invalid chain identifier: {!r}'.format(option))
    return fhandle, option, scheme


def run(fhandle, chain_id, scheme=None):
    """Yield the events of the file with every chain renamed to chain_id."""
    return _atoms.rename(fhandle, scheme, chain=lambda old: chain_id)


def main():
    fhandle, chain_id, scheme = check_input(sys.argv[1:])
    _cli.output(run(fhandle, chain_id, scheme), fhandle)


if __name__ == '__main__':
    main()
