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
Deletes whole categories from an mmCIF file.

Categories are given by name, comma-separated; a name ending in '*' removes
every category that starts with it (quote it, so that the shell does not
expand it). Names are case-insensitive, and the leading '_' is optional.
Everything else is copied unchanged.

Nothing is checked or updated in the categories that remain: e.g. removing
_atom_site leaves _atom_site_anisotrop without its atoms.

Usage:
    mmcif_delcat -<category>[,<category>...] <cif file>

Example:
    mmcif_delcat -_atom_site_anisotrop 1abc.cif
    mmcif_delcat '-_pdbx_validate*' 1abc.cif      # PDB validation reports
    mmcif_delcat -_struct_conf,_struct_sheet_range 1abc.cif

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cif, _cli


def parse_categories(option, doc):
    """'_a,b*' -> ['_a', '_b*'] (fails with the usage on bad names)."""
    categories = []
    for name in _atoms.parse_ids(option):
        name = name if name.startswith('_') else '_' + name
        if name.rstrip('*') in ('', '_') or '*' in name.rstrip('*') or '.' in name:
            _cli.fail(doc, 'Invalid category name: {!r}'.format(name))
        categories.append(name)
    if not categories:
        _cli.fail(doc, 'You must provide at least one category')
    return categories


def check_input(args):
    """Return (file handle, categories)."""
    fhandle, option, _ = _cli.parse_args(args, __doc__, needs_option=True)
    return fhandle, parse_categories(option, __doc__)


def run(fhandle, categories, scheme=None):
    """Yield the events of the file without the given categories."""
    for ev in _cif.parse(fhandle, categories):
        if type(ev) in (_cif.Verbatim, _cif.DataBlock):
            yield ev


def main():
    fhandle, categories = check_input(sys.argv[1:])
    _cli.output(run(fhandle, categories), fhandle)


if __name__ == '__main__':
    main()
