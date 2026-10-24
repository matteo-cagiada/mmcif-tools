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
Keeps only some categories of an mmCIF file.

Categories are given by name, comma-separated; a name ending in '*' keeps
every category that starts with it (quote it, so that the shell does not
expand it). Names are case-insensitive, and the leading '_' is optional.
The data_ line of each block is kept; everything else (other categories,
comments) is removed. mmcif_keepcoord is the same as
mmcif_keepcat -_atom_site,_struct_conn.

Nothing is checked or added: e.g. keeping _atom_site_anisotrop without
_atom_site leaves it without its atoms.

Usage:
    mmcif_keepcat -<category>[,<category>...] <cif file>

Example:
    mmcif_keepcat -_atom_site 1abc.cif
    mmcif_keepcat '-_entity*,_struct_asym' 1abc.cif

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _cif, _cli
from mmciftools.mmcif_delcat import parse_categories


def check_input(args):
    """Return (file handle, categories)."""
    fhandle, option, _ = _cli.parse_args(args, __doc__, needs_option=True)
    return fhandle, parse_categories(option, __doc__)


def run(fhandle, categories, scheme=None):
    """Yield the data_ lines and the given categories of the file."""
    for ev in _cif.parse(fhandle, categories):
        if type(ev) is not _cif.Verbatim:  # every other category is Verbatim here
            yield ev


def main():
    fhandle, categories = check_input(sys.argv[1:])
    _cli.output(run(fhandle, categories), fhandle)


if __name__ == '__main__':
    main()
