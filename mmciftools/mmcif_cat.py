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
Prints one category of an mmCIF file as a table (tab-separated values).

The first line holds the column names (without the category prefix), then
each row of the category follows on a line of its own, so the output can go
to awk, sort, a spreadsheet or pandas. Columns can be chosen and ordered
after a dot (-_atom_site.id,Cartn_x). A category written as key-value pairs
(_cell.length_a 50.0, ...) is printed as one row.

Values are printed without their CIF quotes; nulls stay '.' and '?'. Tabs,
newlines and backslashes inside a value are written as \\t, \\n and \\\\.
In a file with several data blocks, the rows of every block are printed,
and the header is repeated only when the columns change. Category and
column names are case-insensitive, as in CIF.

Usage:
    mmcif_cat -<category>[.<column>,<column>...] <cif file>

Example:
    mmcif_cat -_atom_site 1abc.cif
    mmcif_cat -_atom_site.auth_asym_id,auth_seq_id,Cartn_x 1abc.cif
    mmcif_cat -_chem_comp components.cif

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _cif, _cli

_ESCAPES = str.maketrans({'\\': '\\\\', '\t': '\\t', '\n': '\\n', '\r': '\\r'})


def check_input(args):
    """Return (file handle, category, columns or None)."""
    fhandle, option, _ = _cli.parse_args(args, __doc__, needs_option=True)
    category, _, columns = option.partition('.')
    category = category if category.startswith('_') else '_' + category
    if category == '_' or '*' in category:
        _cli.fail(__doc__, 'Invalid category name: {!r}'.format(option))
    columns = [c.strip() for c in columns.split(',') if c.strip()] if columns else None
    return fhandle, category, columns


def _tsv(values):
    return '\t'.join(v.translate(_ESCAPES) for v in values) + '\n'


def run(fhandle, category, columns=None, scheme=None):
    """Yield the lines of the table."""
    header, found = None, False
    for ev in _cif.parse(fhandle, category):
        if type(ev) is _cif.Row:
            yield _tsv(pick(ev.values))
        elif type(ev) is _cif.Loop:
            found = True
            if columns is None:
                names = ev.columns
                pick = tuple
            else:
                missing = [c for c in columns if c not in ev]
                if missing:
                    raise _cif.CIFError('{} has no column {}'.format(
                        ev.category, ', '.join(missing)), ev.lineno)
                names = columns
                positions = [ev.index(c) for c in columns]
                pick = lambda values: [values[p] for p in positions]  # noqa: E731
            if header is None or [n.lower() for n in names] != header:
                header = [n.lower() for n in names]
                yield _tsv(names)
    if not found:
        raise _cif.CIFError('No {} category in the input'.format(category))


def main():
    fhandle, category, columns = check_input(sys.argv[1:])
    _cli.output_text(_cli.lazy(run, fhandle, category, columns), fhandle)


if __name__ == '__main__':
    main()
