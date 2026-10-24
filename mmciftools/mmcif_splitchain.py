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
Splits an mmCIF file into several, each containing one chain.

As pdb_splitchain, writes <name>_<chain>.cif files in the current directory
(<name> is the input file name without its extension, or 'splitchains' when
reading stdin), overwriting existing ones, and prints nothing. Unlike
pdb_splitchain, each file is a complete mmCIF file: the metadata
categories are kept, the data block is renamed <block>_<chain>, and rows of
_atom_site_anisotrop and _struct_conn of other chains are removed.

Chains are auth_asym_id (label_asym_id with --label).

Usage:
    mmcif_splitchain <cif file>

Example:
    mmcif_splitchain 1abc.cif   # writes 1abc_A.cif, 1abc_B.cif, ...

This program is part of the mmcif-tools suite.
"""

import os
import sys

from mmciftools import _atoms, _cif, _cli


def check_input(args):
    """Return (file handle, output name base, identifier scheme)."""
    fhandle, _, scheme = _cli.parse_args(args, __doc__, takes_option=False)
    path = next((a for a in args if not a.startswith('-')), None)
    return fhandle, _atoms.base_name(path, 'splitchains'), scheme


def chains_of(source, scheme=None):
    """Chain ids in order of first appearance."""
    chains = {}
    for ev in _cif.parse(source.lines()):
        if type(ev) is _cif.Row:
            chains.setdefault(ev.values[chain], None)
        elif type(ev) is _cif.Loop:
            fields = _atoms.Fields(ev, scheme)
            fields.need('chain')
            chain = fields.chain
    return list(chains)


def run(fhandle, base, scheme=None, outdir='.'):
    """Write one file per chain; return their paths."""
    source = _atoms.as_source(fhandle)
    paths = []
    for chain_id in chains_of(source, scheme):
        def make_keep(fields, chain_id=chain_id):
            chain = fields.chain
            return lambda values: values[chain] == chain_id
        path = os.path.join(outdir, '{}_{}.cif'.format(base, chain_id))
        events = _atoms.filter_atoms(source, make_keep, scheme)
        _cli.write_file(path, _atoms.rename_block(events, '_' + chain_id))
        paths.append(path)
    return paths


def main():
    fhandle, base, scheme = check_input(sys.argv[1:])
    try:
        run(fhandle, base, scheme)
    except _cif.CIFError as err:
        sys.stderr.write('ERROR!! {}\n'.format(err))
        sys.exit(1)
    finally:
        fhandle.close()
    sys.exit(0)


if __name__ == '__main__':
    main()
