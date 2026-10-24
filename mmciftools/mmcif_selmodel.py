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
Extracts one or more models from an mmCIF file.

Models are read from _atom_site.pdbx_PDB_model_num. As in pdb_selmodel, a
file with a single model (or without model numbers) is returned whole,
whatever model is asked for. Rows of _atom_site_anisotrop and _struct_conn
that refer to removed atoms are removed too.

Usage:
    mmcif_selmodel -<model id>[,<model id>...] <cif file>

Example:
    mmcif_selmodel -1 1abc.cif     # selects model 1
    mmcif_selmodel -1,3 1abc.cif   # selects models 1 and 3

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cif, _cli


def check_input(args):
    """Return (file handle, model numbers, identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__, needs_option=True)
    models = []
    for model in _atoms.parse_ids(option):
        try:
            models.append(int(model))
        except ValueError:
            _cli.fail(__doc__, 'Model identifier is invalid: {!r}'.format(model))
    if not models:
        _cli.fail(__doc__, 'You must provide at least one model identifier')
    return fhandle, models, scheme


def _has_several_models(source):
    """True if the atoms belong to more than one model."""
    seen = set()
    for ev in _cif.parse(source.lines()):
        if type(ev) is _cif.Row:
            if position is None:
                return False
            seen.add(ev.values[position])
            if len(seen) > 1:
                return True
        elif type(ev) is _cif.Loop:
            position = _atoms.Fields(ev).model
    return False


def run(fhandle, models, scheme=None):
    """Yield the events of the file, keeping only atoms of the given models."""
    models = set(models)
    source = _atoms.as_source(fhandle)
    if not _has_several_models(source):
        return _cif.parse(source.lines())  # nothing to select from, as in pdb_selmodel

    def make_keep(fields):
        model = fields.model
        return lambda values: _atoms.as_int(values[model]) in models

    return _atoms.filter_atoms(source, make_keep, scheme)


def main():
    fhandle, models, scheme = check_input(sys.argv[1:])
    _cli.output(run(fhandle, models, scheme), fhandle)


if __name__ == '__main__':
    main()
