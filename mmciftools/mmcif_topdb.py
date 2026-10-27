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
Converts the coordinates of an mmCIF file to the PDB format.

As pdb_fromcif, writes the ATOM/HETATM records of _atom_site (with MODEL /
ENDMDL records when there are several models) and END; other categories,
and data blocks after the first, are not converted. Chains, residue numbers and names are auth identifiers,
falling back to label ones when a file lacks them (with --label, label ones
are used). Atom serials restart at 1 in each model. Unlike pdb_fromcif, a
TER record follows the last residue of each polymer (atoms with a
label_seq_id, or ATOM records in files without one) and takes a serial
number, as the PDB format specifies.

Will not convert a structure that does not fit the PDB format: chains of
more than one character, more than 99,999 atoms in a model, residue numbers
outside -999..9999, residue names longer than 3 characters or atom names
longer than 4; nothing is written then. Unlike pdb_fromcif, residue names
and elements are right-justified and MODEL numbers placed as the PDB
format specifies, formal charges are written as 1-, 2+, ..., and missing
occupancies / B-factors are written as 1.00 / 0.00.

Usage:
    mmcif_topdb <cif file>

Example:
    mmcif_topdb 1abc.cif > 1abc.pdb

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cif, _cli

NULLS = ('.', '?')
ATOM_LINE = ('{:<6}{:>5} {:4}{:1}{:>3} {:1}{:>4}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}'
             '{:>6.2f}{:>6.2f}          {:>2}{:2}\n')
TER_LINE = 'TER   {:>5}      {:>3} {:1}{:>4}{:1}'


def check_input(args):
    """Return (file handle, identifier scheme)."""
    fhandle, _, scheme = _cli.parse_args(args, __doc__, takes_option=False)
    return fhandle, scheme


def _charge(value):
    num = _atoms.as_int(value)
    if not num:
        return ''
    return '{}{}'.format(abs(num), '-' if num < 0 else '+')


def _number(value, default):
    return float(value) if value not in NULLS else default


def run(fhandle, scheme=None):
    """Return the PDB lines (raises CIFError if the structure does not fit)."""
    models, current, serial = [], None, 0
    polymer = None  # (label_asym_id, TER fields) of the polymer being written

    def end_polymer():
        nonlocal serial, polymer
        serial += 1
        models[-1][1].append('{:<80}\n'.format(TER_LINE.format(serial, *polymer[1])))
        polymer = None

    for ev in _cif.parse(fhandle):
        if type(ev) is _cif.Row:
            v = ev.values
            model = v[f.model] if f.model is not None else '1'
            if model != current:
                if polymer is not None:
                    end_polymer()
                models.append((model, []))
                current, serial = model, 0
            if seq_pos is not None:
                key = (v[asym_pos] if asym_pos is not None else v[f.chain]
                       ) if _atoms.as_int(v[seq_pos]) is not None else None
            else:
                key = v[f.chain] if f.group is not None and v[f.group] == 'ATOM' else None
            if polymer is not None and key != polymer[0]:
                end_polymer()
            serial += 1
            chain, name, resname = v[f.chain], v[f.atom], v[f.resname]
            resnum = _atoms.as_int(v[f.resnum])
            element = v[f.element] if f.element is not None else ''
            element = '' if element in NULLS else element.upper()
            where = 'line {}'.format(ev.lineno)
            if serial > 99999:  # TER records count too
                raise _cif.CIFError('Number of atoms exceeds PDB format limit in model {}'.format(model))
            if len(chain) > 1:
                raise _cif.CIFError('Chain IDs is too large: {!r} ({})'.format(chain, where))
            if resnum is None or not -999 <= resnum <= 9999:
                raise _cif.CIFError('Residue number does not fit the PDB format: {!r} ({})'.format(
                    v[f.resnum], where))
            if len(resname) > 3 or len(name) > 4:
                raise _cif.CIFError('Name does not fit the PDB format: {} {} ({})'.format(
                    resname, name, where))
            altloc = v[f.altloc] if f.altloc is not None else '.'
            icode = v[f.icode] if f.icode is not None else '?'
            padded = ' ' + name if len(name) < 4 and name[:1].isalpha() and len(element) < 2 else name
            try:
                coords = [float(v[p]) for p in xyz]
                occupancy = _number(v[f.occupancy], 1.0) if f.occupancy is not None else 1.0
                bfactor = _number(v[f.bfactor], 0.0) if f.bfactor is not None else 0.0
            except ValueError:
                raise _cif.CIFError('Coordinates, occupancy or B-factor are not numbers ({})'.format(where))
            group = v[f.group] if f.group is not None else 'ATOM'
            charge = _charge(v[charge_pos]) if charge_pos is not None else ''
            icode = '' if icode in NULLS else icode[:1]
            models[-1][1].append(ATOM_LINE.format(
                group, serial, padded, '' if altloc in NULLS else altloc[:1], resname, chain,
                resnum, icode, coords[0], coords[1], coords[2], occupancy, bfactor, element,
                charge))
            if key is not None:
                polymer = (key, (resname, chain, resnum, icode))
        elif type(ev) is _cif.Loop:
            if models:
                break  # only the first _atom_site of the file
            f = _atoms.Fields(ev, scheme)
            f.need('chain', 'resnum', 'resname', 'atom')
            xyz = [ev.index(c) for c in ('Cartn_x', 'Cartn_y', 'Cartn_z')]
            charge_pos = ev.index('pdbx_formal_charge') if 'pdbx_formal_charge' in ev else None
            seq_pos = ev.index('label_seq_id') if 'label_seq_id' in ev else None
            asym_pos = ev.index('label_asym_id') if 'label_asym_id' in ev else None
    if polymer is not None:
        end_polymer()
    if not models:
        raise _cif.CIFError('No atoms (_atom_site) to convert')
    lines = []
    for model, atoms in models:
        if len(models) > 1:
            lines.append('MODEL     {:>4}\n'.format(model))
        lines.extend(atoms)
        if len(models) > 1:
            lines.append('ENDMDL\n')
    lines.append('{:<80}\n'.format('END'))
    return lines


def main():
    fhandle, scheme = check_input(sys.argv[1:])
    _cli.output_text(_cli.lazy(run, fhandle, scheme), fhandle)


if __name__ == '__main__':
    main()
