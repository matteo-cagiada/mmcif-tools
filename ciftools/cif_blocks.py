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
Lists or extracts the data blocks of a CIF file.

Without an option, prints the name of each data block (data_<name>), one
per line. With block names, comma-separated, writes those blocks (in file
order) and nothing else. Names are case-insensitive, as in CIF; reading
stops once every requested block has been written, so extracting one
component from the whole Chemical Component Dictionary is quick.

Usage:
    cif_blocks [-<name>[,<name>...]] <cif file>

Example:
    cif_blocks components.cif              # lists the blocks
    cif_blocks -ATP,GTP components.cif     # writes the ATP and GTP blocks

This program is part of the cif-tools suite, a port of pdb-tools to mmCIF.
"""

import sys

from ciftools import _atoms, _cif, _cli


def check_input(args):
    """Return (file handle, block names or None)."""
    fhandle, option, _ = _cli.parse_args(args, __doc__)
    names = None
    if option is not None:
        names = _atoms.parse_ids(option)
        if not names:
            _cli.fail(__doc__, 'You must provide at least one block name')
    return fhandle, names


def run(fhandle, names=None, scheme=None):
    """Yield the block names (as lines) or, given names, the events of those
    blocks."""
    events = _cif.parse(fhandle, ())
    if names is None:
        for ev in events:
            if type(ev) is _cif.DataBlock:
                yield ev.name + '\n'
        return
    missing = {n.lower(): n for n in names}
    inside = False
    for ev in events:
        if type(ev) is _cif.DataBlock:
            if not missing:
                return  # every block written
            inside = missing.pop(ev.name.lower(), None) is not None
        if inside:
            yield ev
    if missing:
        raise _cif.CIFError('No data block named {}'.format(', '.join(missing.values())))


def main():
    fhandle, names = check_input(sys.argv[1:])
    _cli.output(run(fhandle, names), fhandle)


if __name__ == '__main__':
    main()
