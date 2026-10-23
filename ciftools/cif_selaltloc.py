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
Selects alternate locations (_atom_site.label_alt_id) for the entire file.

By default, for each residue the alternate location with the highest
occupancy is kept and the others are removed. With an option (e.g. -A),
location A is kept wherever a residue has it; atoms that have other
locations but no A (e.g. only B and C) are kept as they are. Kept atoms
lose their altloc label ('.'); occupancies are not changed.

The choice is made per residue, so a residue modelled as two different
residue types (e.g. SER/THR) is never mixed. Atoms keep their order, and
rows of _atom_site_anisotrop of removed atoms are removed too.

Usage:
    cif_selaltloc [-<altloc>] <cif file>

Example:
    cif_selaltloc 1abc.cif      # picks locations with highest occupancy
    cif_selaltloc -A 1abc.cif   # picks alternate locations labelled 'A'

This program is part of the cif-tools suite, a port of pdb-tools to mmCIF.
"""

import sys

from ciftools import _atoms, _cif, _cli

NULLS = ('.', '?')


def check_input(args):
    """Return (file handle, altloc or None, identifier scheme)."""
    fhandle, option, scheme = _cli.parse_args(args, __doc__)
    if option is not None and (not option or len(option.split()) != 1):
        _cli.fail(__doc__, 'Invalid alternate location identifier: {!r}'.format(option))
    return fhandle, option, scheme


def _occupancy(value):
    try:
        return float(value)
    except ValueError:
        return 0.0


def run(fhandle, altloc=None, scheme=None):
    """Yield the events of the file with one alternate location kept."""
    removed, blanked = set(), set()

    def flush(group):
        locs = []
        for row in group:
            loc = row.values[alt]
            if loc not in NULLS and loc not in locs:
                locs.append(loc)
        if not locs:
            return group
        if altloc is None:
            best = {}
            for row in group:
                loc = row.values[alt]
                if loc not in NULLS:
                    occ = _occupancy(row.values[occ_pos]) if occ_pos is not None else 0.0
                    best[loc] = max(best.get(loc, occ), occ)
            chosen = max(locs, key=lambda loc: best[loc])  # first one on ties
        elif altloc in locs:
            chosen = altloc
        else:
            return group
        # Atoms (by name) that have a conformer at the chosen location
        has_chosen = {row.values[name] for row in group if row.values[alt] == chosen}
        out = []
        for row in group:
            loc = row.values[alt]
            if loc in NULLS:
                out.append(row)
            elif loc == chosen:
                row[alt_col] = '.'
                if id_pos is not None:
                    blanked.add(row.values[id_pos])
                out.append(row)
            elif altloc is not None and row.values[name] not in has_chosen:
                out.append(row)  # e.g. -A and an atom with only B and C
            elif id_pos is not None:
                removed.add(row.values[id_pos])
        return out

    kind = None
    group, group_key = [], None
    for ev in _cif.parse(fhandle, (_atoms.ATOM_SITE, _atoms.ANISOTROP)):
        etype = type(ev)
        if etype is _cif.Row and kind == _atoms.ATOM_SITE:
            key = residue(ev.values)
            if key != group_key:
                for row in flush(group):
                    yield row
                group, group_key = [], key
            group.append(ev)
            continue
        if group:
            for row in flush(group):
                yield row
            group, group_key = [], None
        if etype is _cif.Row and kind == _atoms.ANISOTROP and aniso_id is not None:
            atom_id = ev.values[aniso_id]
            if atom_id in removed:
                continue
            if aniso_alt is not None and atom_id in blanked:
                ev[ev.loop.columns[aniso_alt]] = '.'
        elif etype is _cif.Loop:
            kind = ev.category.lower()
            if kind == _atoms.ATOM_SITE:
                fields = _atoms.Fields(ev, scheme)
                fields.need('altloc', 'atom')
                alt, name, occ_pos, id_pos = fields.altloc, fields.atom, fields.occupancy, fields.id
                alt_col = ev.columns[alt]
                residue = _atoms.getter(fields.model, fields.chain, fields.resnum, fields.icode)
            elif kind == _atoms.ANISOTROP:
                aniso_id = ev.index('id') if 'id' in ev else None
                aniso_alt = ev.index('pdbx_label_alt_id') if 'pdbx_label_alt_id' in ev else None
        elif etype is _cif.LoopEnd:
            kind = None
        yield ev


def main():
    fhandle, altloc, scheme = check_input(sys.argv[1:])
    _cli.output(run(fhandle, altloc, scheme), fhandle)


if __name__ == '__main__':
    main()
