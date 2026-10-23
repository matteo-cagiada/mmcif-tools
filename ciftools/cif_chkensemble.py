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
Checks that all models of a multi-model mmCIF file have the same
composition.

Composition is defined as same atoms/residues/chains: an atom is its chain,
residue name, number, insertion code, atom name and altloc (auth
identifiers, or label ones with --label). Unlike pdb_chkensemble, the atom
id is not part of it (mmCIF atom ids are unique across models, so every
ensemble would differ). Differences are written to stderr as by
pdb_chkensemble, one atom per line; the exit status is 1 if models differ.

Usage:
    cif_chkensemble <cif file>

Example:
    cif_chkensemble 1abc.cif

This program is part of the cif-tools suite, a port of pdb-tools to mmCIF.
"""

import itertools
import sys

from ciftools import _atoms, _cif, _cli

NULLS = ('.', '?')


def check_input(args):
    """Return (file handle, identifier scheme)."""
    fhandle, _, scheme = _cli.parse_args(args, __doc__, takes_option=False)
    return fhandle, scheme


def _label(atom):
    chain, resname, resnum, icode, name, altloc = atom
    text = '{}:{}{}{}:{}'.format(chain, resname, resnum, '' if icode in NULLS else icode, name)
    return text if altloc in NULLS else text + ':' + altloc


def compositions(fhandle, scheme=None):
    """{model: set of atoms}, models in order of appearance."""
    models = {}
    for ev in _cif.parse(fhandle):
        if type(ev) is _cif.Row:
            values = ev.values
            model = values[f.model] if f.model is not None else '1'
            models.setdefault(model, set()).add(identity(values))
        elif type(ev) is _cif.Loop:
            f = _atoms.Fields(ev, scheme)
            f.need('chain', 'resnum', 'resname', 'atom')
            identity = _atoms.getter(f.chain, f.resname, f.resnum, f.icode, f.atom, f.altloc)
    return models


def run(fhandle, scheme=None):
    """Return (stdout lines, stderr lines, exit status)."""
    models = compositions(fhandle, scheme)
    errors = []
    for (i, atoms_i), (j, atoms_j) in itertools.combinations(models.items(), 2):
        if atoms_i == atoms_j:
            continue
        errors.append('Models {} and {} differ:\n'.format(i, j))
        for model, only in ((i, atoms_i - atoms_j), (j, atoms_j - atoms_i)):
            if only:
                errors.append('Atoms in model {} only:\n'.format(model))
                errors.extend(_label(atom) + '\n' for atom in sorted(only, key=_label))
    if errors:
        return [], errors, 1
    return ['Ensemble of {} models *seems* OK\n'.format(len(models))], [], 0


def main():
    fhandle, scheme = check_input(sys.argv[1:])
    try:
        out, err, status = run(fhandle, scheme)
    except _cif.CIFError as error:
        sys.stderr.write('ERROR!! {}\n'.format(error))
        sys.exit(1)
    finally:
        fhandle.close()
    sys.stdout.writelines(out)
    sys.stderr.writelines(err)
    sys.exit(status)


if __name__ == '__main__':
    main()
