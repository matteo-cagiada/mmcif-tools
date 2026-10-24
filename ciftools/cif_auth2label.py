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
Copies author chains onto the label chains of polymers (auth_asym_id ->
label_asym_id) in an mmCIF file.

For programs that read label_asym_id, e.g. after renaming chains with
cif_rplchain. Only chains are copied, and only for polymers (atoms with a
label_seq_id):

  - Ligands and waters keep their label chains. In mmCIF each of them is a
    molecule of its own, while in PDB entries they share the author chain
    of their protein, so copying would merge them.
  - Residue numbers are never copied: label_seq_id counts 1, 2, ... along
    the deposited sequence, which author numbers (with their offsets, gaps
    and insertion codes) do not. cif_label2auth -res goes the other way.

Nothing is written if an author chain holds more than one polymer, or if a
new label chain is already the label chain of another molecule. Rows of
_atom_site_anisotrop and _struct_conn follow; other categories that mention
label chains (_struct_asym, _pdbx_poly_seq_scheme, ...) are not updated, as
with cif_rplchain --label.

Usage:
    cif_auth2label <cif file>

Example:
    cif_rplchain -A:H 1abc.cif | cif_auth2label

This program is part of the cif-tools suite, a port of pdb-tools to mmCIF.
"""

import sys

from ciftools import _atoms, _cif, _cli


def check_input(args):
    """Return the file handle."""
    fhandle, _, _ = _cli.parse_args(args, __doc__, takes_option=False)
    return fhandle


def chain_map(source):
    """{label chain of a polymer: its author chain}, for the chains that
    change (CIFError when that would merge molecules)."""
    polymers, others = {}, set()  # label chain -> author chains; non-polymer label chains
    for ev in _cif.parse(source.lines(), _atoms.ATOM_SITE):
        if type(ev) is _cif.Row:
            values = ev.values
            if is_polymer(values):
                polymers.setdefault(values[label], set()).add(values[auth])
            else:
                others.add(values[label])
        elif type(ev) is _cif.Loop:
            for name in ('label_asym_id', 'auth_asym_id'):
                if name not in ev:
                    raise _cif.CIFError('_atom_site has no {} column (needed here)'.format(name),
                                        ev.lineno)
            label, auth = ev.index('label_asym_id'), ev.index('auth_asym_id')
            if 'label_seq_id' in ev:
                seq = ev.index('label_seq_id')
                is_polymer = lambda values: _atoms.as_int(values[seq]) is not None  # noqa: E731
            elif 'group_PDB' in ev:
                group = ev.index('group_PDB')
                is_polymer = lambda values: values[group] == 'ATOM'  # noqa: E731
            else:
                raise _cif.CIFError('_atom_site has no label_seq_id column (needed here)',
                                    ev.lineno)

    owners = {}  # author chain -> polymer label chain
    for chain, authors in sorted(polymers.items()):
        if len(authors) > 1:
            raise _cif.CIFError('Polymer {} has atoms in author chains {}'.format(
                chain, ', '.join(sorted(authors))))
        author = authors.pop()
        if author in owners:
            raise _cif.CIFError('Author chain {} holds two polymers ({} and {})'.format(
                author, owners[author], chain))
        owners[author] = chain
        if author in others - set(polymers):
            raise _cif.CIFError('Label chain {} is already used by a ligand or water'.format(
                author))
    return {chain: author for author, chain in owners.items() if author != chain}


def run(fhandle, scheme=None):
    """Yield the events of the file with author chains copied onto the label
    chains of polymers."""
    source = _atoms.as_source(fhandle)
    mapping = chain_map(source)
    return _atoms.rename(source.lines(), 'label', chain=mapping.get)


def main():
    fhandle = check_input(sys.argv[1:])
    _cli.output(_cli.lazy(run, fhandle), fhandle)


if __name__ == '__main__':
    main()
