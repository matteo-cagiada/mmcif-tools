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
Splits an mmCIF file into several, each containing one model.

As pdb_splitmodel, writes <name>_<model>.cif files in the current directory
(<name> is the input file name without its extension, or 'splitmodels'
when reading stdin), named by the model numbers as written, and prints
nothing. Unlike pdb_splitmodel, which writes nothing for a file without
MODEL records, a single-model file gives one file; and each file is a
complete mmCIF file: the metadata categories are kept, the data block is
renamed <block>_<model>, and rows of _atom_site_anisotrop and _struct_conn
of removed atoms are removed.

Usage:
    cif_splitmodel <cif file>

Example:
    cif_splitmodel 1abc.cif   # writes 1abc_1.cif, 1abc_2.cif, ...

This program is part of the cif-tools suite, a port of pdb-tools to mmCIF.
"""

import os
import sys

from ciftools import _atoms, _cif, _cli


def check_input(args):
    """Return (file handle, output name base, identifier scheme)."""
    fhandle, _, scheme = _cli.parse_args(args, __doc__, takes_option=False)
    path = next((a for a in args if not a.startswith('-')), None)
    return fhandle, _atoms.base_name(path, 'splitmodels'), scheme


def models_of(source):
    """Model numbers in order of first appearance ('1' if there are none)."""
    models = {}
    for ev in _cif.parse(source.lines()):
        if type(ev) is _cif.Row:
            models.setdefault(ev.values[model] if model is not None else '1', None)
        elif type(ev) is _cif.Loop:
            model = _atoms.Fields(ev).model
    return list(models)


def run(fhandle, base, scheme=None, outdir='.'):
    """Write one file per model; return their paths."""
    source = _atoms.as_source(fhandle)
    paths = []
    for model_id in models_of(source):
        def make_keep(fields, model_id=model_id):
            model = fields.model
            if model is None:
                return lambda values: True
            return lambda values: values[model] == model_id
        path = os.path.join(outdir, '{}_{}.cif'.format(base, model_id))
        events = _atoms.filter_atoms(source, make_keep, scheme)
        _cli.write_file(path, _atoms.rename_block(events, '_' + model_id))
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
