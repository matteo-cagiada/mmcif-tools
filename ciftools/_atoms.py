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
Helpers for tools that select or edit atoms (rows of _atom_site).

Removing or renaming atoms also affects the rows of other categories that
point at them. The helpers here keep two of those in step:

    _atom_site_anisotrop   rows of removed atoms are removed; renamed chains
                           and residues are renamed
    _struct_conn           links to a removed atom are removed; renamed
                           chains and residues are renamed

Other categories that mention chains or residues (_struct_conf,
_pdbx_poly_seq_scheme, ...) are left as they are.

_struct_conn normally comes before _atom_site, so deciding which links to
drop needs a first look at the atoms. Tools therefore read their input
through a Source, which can be read more than once (stdin is copied to a
temporary file, kept in memory while small).
"""

import io
import operator
import os
import stat
import tempfile

from ciftools import _cif

ATOM_SITE = '_atom_site'
ANISOTROP = '_atom_site_anisotrop'
STRUCT_CONN = '_struct_conn'


class Fields:
    """Positions of the atom identifiers in an _atom_site loop.

    chain, resnum, resname and atom follow the identifier scheme (auth_* or
    label_*, falling back to the other one when a file lacks a column);
    icode is only used with the auth scheme. Absent columns are None.
    """

    __slots__ = ('loop', 'scheme', 'id', 'group', 'model', 'chain', 'resnum', 'icode',
                 'resname', 'atom', 'altloc', 'element', 'occupancy', 'bfactor')

    def __init__(self, loop, scheme=None):
        self.loop = loop
        self.scheme = scheme = _cif.get_scheme(scheme)

        def pos(name):
            return loop.index(name) if name in loop else None

        def scheme_pos(field):
            try:
                return loop.index(_cif.column(loop, field, scheme))
            except _cif.CIFError:
                return None

        self.id = pos('id')
        self.group = pos('group_PDB')
        self.model = pos('pdbx_PDB_model_num')
        self.chain = scheme_pos('chain')
        self.resnum = scheme_pos('resnum')
        self.icode = pos('pdbx_PDB_ins_code') if scheme == 'auth' else None
        self.resname = scheme_pos('resname')
        self.atom = scheme_pos('atom')
        self.altloc = pos('label_alt_id')
        self.element = pos('type_symbol')
        self.occupancy = pos('occupancy')
        self.bfactor = pos('B_iso_or_equiv')

    def need(self, *names):
        """Raise CIFError unless the loop has the given fields."""
        for name in names:
            if getattr(self, name) is None:
                raise _cif.CIFError('{} has no {} column (needed here)'.format(
                    self.loop.category, name), self.loop.lineno)


def getter(*positions):
    """Fast tuple-of-values getter; None positions give None."""
    present = [p for p in positions if p is not None]
    if len(present) == len(positions):
        if len(positions) == 1:
            p = positions[0]
            return lambda values: (values[p],)
        return operator.itemgetter(*positions)
    return lambda values: tuple(None if p is None else values[p] for p in positions)


def parse_ids(option):
    """'A,B, C' -> ['A', 'B', 'C'] (order kept, blanks dropped)."""
    return [o.strip() for o in option.split(',') if o.strip()]


def as_int(value):
    """int(value), or None for nulls and other non-numbers."""
    try:
        return int(value)
    except ValueError:
        return None


def parse_number(text):
    """A residue number given on the command line (ValueError if invalid)."""
    try:
        return int(text)
    except ValueError:
        raise ValueError('Not a valid number: {!r}'.format(text))


RANGE_HELP = ("Residue range must be in 'a:z:s' where a and z are optional (default to "
              "first residue and last respectively), and s is an optional step value "
              "(to return every s-th residue).")


def parse_range(text):
    """'a:z:s' -> (start, end, step); a missing start/end is -inf/+inf, step 1."""
    parts = text.split(':')
    if len(parts) not in (2, 3):
        raise ValueError(RANGE_HELP)
    start, end = parts[0].strip(), parts[1].strip()
    step = parts[2].strip() if len(parts) == 3 else ''
    start = parse_number(start) if start else float('-inf')
    end = parse_number(end) if end else float('inf')
    step = parse_number(step) if step else 1
    if step < 1:
        raise ValueError('Step value must be a positive number: {!r}'.format(step))
    if start > end:
        raise ValueError('Start ({}) cannot be larger than end ({})'.format(start, end))
    return start, end, step


def residue_numbers(source, scheme=None):
    """Residue numbers in file order, one per residue occurrence (a new
    occurrence starts whenever the chain or residue number changes)."""
    numbers = []
    prev = None
    for ev in _cif.parse(source.lines(), ATOM_SITE):
        if type(ev) is _cif.Row:
            values = ev.values
            key = (values[chain], values[resnum])
            if key != prev:
                prev = key
                num = as_int(values[resnum])
                if num is not None:
                    numbers.append(num)
        elif type(ev) is _cif.Loop:
            fields = Fields(ev, scheme)
            fields.need('chain', 'resnum')
            chain, resnum = fields.chain, fields.resnum
    return numbers


def rereadable(fhandle):
    """True if fhandle can be rewound and read again.

    That is a regular file (plain or gzip) or an in-memory stream. seekable()
    alone is not enough: on Windows a pipe reports itself as seekable, but
    seeking back does not bring its data back.
    """
    try:
        if not fhandle.seekable():
            return False
    except (AttributeError, ValueError):
        return False
    try:
        fd = fhandle.fileno()
    except (AttributeError, OSError, io.UnsupportedOperation):
        return True  # e.g. io.StringIO
    try:
        return stat.S_ISREG(os.fstat(fd).st_mode)
    except OSError:
        return False


class Source:
    """The input of a tool, readable more than once.

    A regular file is re-read in place; anything else (stdin, a pipe) is
    first copied to a temporary file.
    """

    def __init__(self, fhandle):
        if rereadable(fhandle):
            self.fh = fhandle
            self.start = fhandle.tell()
        else:
            spool = tempfile.SpooledTemporaryFile(
                max_size=32 * 1024 * 1024, mode='w+', encoding=_cif.ENCODING,
                errors=_cif.ERRORS, newline='')
            for line in fhandle:
                spool.write(line)
            self.fh = spool
            self.start = 0

    def lines(self):
        """An iterator over the lines, from the beginning."""
        self.fh.seek(self.start)
        return iter(self.fh.readline, '')


def as_source(fhandle):
    return fhandle if isinstance(fhandle, Source) else Source(fhandle)


# ---------------------------------------------------------------------------
# Selecting atoms
# ---------------------------------------------------------------------------

# (_atom_site column, _struct_conn column template) pairs identifying a link
# partner; the columns present in both categories are used.
_LINK_KEYS = [
    ('label_asym_id', 'ptnr{}_label_asym_id'),
    ('label_comp_id', 'ptnr{}_label_comp_id'),
    ('label_seq_id', 'ptnr{}_label_seq_id'),
    ('label_atom_id', 'ptnr{}_label_atom_id'),
    ('auth_asym_id', 'ptnr{}_auth_asym_id'),
    ('auth_seq_id', 'ptnr{}_auth_seq_id'),
    ('pdbx_PDB_ins_code', 'pdbx_ptnr{}_PDB_ins_code'),
]


def _link_columns(atom_loop, conn_loop):
    return [(a, p) for a, p in _LINK_KEYS
            if a in atom_loop and p.format(1) in conn_loop and p.format(2) in conn_loop]


def _partner_getters(conn_loop, columns):
    first = getter(*[conn_loop.index(p.format(1)) for _, p in columns])
    second = getter(*[conn_loop.index(p.format(2)) for _, p in columns])
    return first, second


def _dangling_links(source, make_keep, scheme):
    """Identities of link partners that exist in the input but not in the output.

    Returns (dangling identities, columns used), or None when there is
    nothing to check.
    """
    conn_loop, partners = None, []
    for ev in _cif.parse(source.lines(), STRUCT_CONN):  # _atom_site is not tokenized here
        if type(ev) is _cif.Row:
            partners.append(ev)
        elif type(ev) is _cif.Loop:
            conn_loop = ev
    if not partners:
        return None

    seen, kept = set(), set()
    columns = None
    for ev in _cif.parse(source.lines(), ATOM_SITE):
        if type(ev) is _cif.Loop:
            columns = _link_columns(ev, conn_loop)
            if not columns:
                return None
            atom_key = getter(*[ev.index(a) for a, _ in columns])
            first, second = _partner_getters(conn_loop, columns)
            wanted = {first(r.values) for r in partners} | {second(r.values) for r in partners}
            keep = make_keep(Fields(ev, scheme))
        elif type(ev) is _cif.Row:
            keeping = keep(ev.values)  # on every row: keep() may count residues
            key = atom_key(ev.values)
            if key in wanted:
                seen.add(key)
                if keeping:
                    kept.add(key)
    if columns is None:
        return None
    return seen - kept, columns


def filter_atoms(fhandle, make_keep, scheme=None):
    """Yield the events of the input, keeping only the atoms that pass.

    `make_keep(fields)` is called once per _atom_site loop (and once more
    for the first look at the file) with its Fields, and returns a function
    that takes a row's values and says whether to keep the atom; it is
    called on every row in file order, so it may keep state. Rows of
    _atom_site_anisotrop and _struct_conn that refer to removed atoms are
    removed as well.
    """
    source = as_source(fhandle)
    links = _dangling_links(source, make_keep, scheme)
    dangling = links[0] if links else None

    removed = set()  # ids of removed atoms, for _atom_site_anisotrop
    kind = None
    for ev in _cif.parse(source.lines(), (ATOM_SITE, ANISOTROP, STRUCT_CONN)):
        etype = type(ev)
        if etype is _cif.Row:
            if kind == ATOM_SITE:
                if keep(ev.values):
                    yield ev
                elif id_pos is not None:
                    removed.add(ev.values[id_pos])
                continue
            if kind == ANISOTROP:
                if aniso_id is not None and ev.values[aniso_id] in removed:
                    continue
            elif kind == STRUCT_CONN and dangling:
                if first(ev.values) in dangling or second(ev.values) in dangling:
                    continue
            yield ev
            continue
        if etype is _cif.Loop:
            kind = ev.category.lower()
            if kind == ATOM_SITE:
                fields = Fields(ev, scheme)
                keep = make_keep(fields)
                id_pos = fields.id
            elif kind == ANISOTROP:
                aniso_id = ev.index('id') if 'id' in ev else None
            elif kind == STRUCT_CONN and dangling:
                first, second = _partner_getters(ev, links[1])
        elif etype is _cif.LoopEnd:
            kind = None
        yield ev


# ---------------------------------------------------------------------------
# Renaming chains and residues
# ---------------------------------------------------------------------------

def _reference_columns(loop, scheme):
    """(chain, resnum, icode, resname) column positions in
    _atom_site_anisotrop or _struct_conn rows, as a list (one entry per
    partner)."""
    def pos(name):
        return loop.index(name) if name in loop else None

    kind = loop.category.lower()
    auth = scheme == 'auth'
    if kind == ANISOTROP:
        if auth:
            return [(pos('pdbx_auth_asym_id'), pos('pdbx_auth_seq_id'), pos('pdbx_PDB_ins_code'),
                     pos('pdbx_auth_comp_id'))]
        return [(pos('pdbx_label_asym_id'), pos('pdbx_label_seq_id'), None,
                 pos('pdbx_label_comp_id'))]
    out = []
    for n in (1, 2):
        if auth:
            out.append((pos('ptnr{}_auth_asym_id'.format(n)), pos('ptnr{}_auth_seq_id'.format(n)),
                        pos('pdbx_ptnr{}_PDB_ins_code'.format(n)),
                        pos('ptnr{}_auth_comp_id'.format(n))))
        else:
            out.append((pos('ptnr{}_label_asym_id'.format(n)),
                        pos('ptnr{}_label_seq_id'.format(n)), None,
                        pos('ptnr{}_label_comp_id'.format(n))))
    return out


def rename(fhandle, scheme=None, chain=None, atom_resnum=None, ref_resnum=None, start=None):
    """Yield the events of the input with chains and/or residues renamed.

    chain(old_chain) -> new chain (or None to keep it)
    atom_resnum(fields) -> function(values) -> new residue number (or None);
        called once per _atom_site loop, may keep state between rows
    ref_resnum(chain, resnum, icode, resname) -> new residue number (or
        None), used for _atom_site_anisotrop and _struct_conn rows (with the
        old names; None for a column the category lacks)
    start(fields), if given, is called at each _atom_site loop (e.g. to
        check that the needed columns exist).
    """
    scheme = _cif.get_scheme(scheme)
    kind = None
    for ev in _cif.parse(fhandle, (ATOM_SITE, ANISOTROP, STRUCT_CONN)):
        etype = type(ev)
        if etype is _cif.Row:
            if kind == ATOM_SITE:
                values = ev.values
                if new_resnum is not None:
                    num = new_resnum(values)
                    if num is not None and num != values[f_resnum]:
                        ev[resnum_col] = num
                if chain is not None:
                    new = chain(values[f_chain])
                    if new is not None and new != values[f_chain]:
                        ev[chain_col] = new
            elif kind is not None:
                for c_pos, r_pos, i_pos, n_pos in refs:
                    values = ev.values
                    old_chain = values[c_pos] if c_pos is not None else None
                    if ref_resnum is not None and r_pos is not None:
                        icode = values[i_pos] if i_pos is not None else None
                        resname = values[n_pos] if n_pos is not None else None
                        num = ref_resnum(old_chain, values[r_pos], icode, resname)
                        if num is not None and num != values[r_pos]:
                            ev[ev.loop.columns[r_pos]] = num
                    if chain is not None and c_pos is not None:
                        new = chain(old_chain)
                        if new is not None and new != old_chain:
                            ev[ev.loop.columns[c_pos]] = new
            yield ev
            continue
        if etype is _cif.Loop:
            kind = ev.category.lower()
            if kind == ATOM_SITE:
                fields = Fields(ev, scheme)
                if start is not None:
                    start(fields)
                if chain is not None:
                    fields.need('chain')
                    f_chain = fields.chain
                    chain_col = ev.columns[f_chain]
                new_resnum = None
                if atom_resnum is not None:
                    fields.need('resnum')
                    f_resnum = fields.resnum
                    resnum_col = ev.columns[f_resnum]
                    new_resnum = atom_resnum(fields)
            else:
                refs = _reference_columns(ev, scheme)
        elif etype is _cif.LoopEnd:
            kind = None
        yield ev
