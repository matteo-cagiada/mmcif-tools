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
Re-formats an mmCIF file in a uniform layout.

Every category is rewritten: a '#' line between categories, key-value
categories with their values aligned, loops with each column padded to its
widest value, values quoted only when needed (text fields on lines of
their own) and '\\n' line endings. Comments between categories are kept
(the '#' separators are rewritten); comments inside a category are
dropped. The data itself is not changed. Use it after a pipeline of other
tools, which keep the original layout and only rewrite what they change.

Unlike pdb_tidy, which fixes PDB-format records (TER, END, line padding),
there is nothing to repair in mmCIF, so this tool only changes the layout.

Usage:
    mmcif_tidy <cif file>

Example:
    mmcif_tidy 1abc.cif

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _atoms, _cif, _cli


def check_input(args):
    """Return (file handle, identifier scheme)."""
    fhandle, _, scheme = _cli.parse_args(args, __doc__, takes_option=False)
    return fhandle, scheme


def _widths(source):
    """Column widths of every loop, in order of appearance (None for
    key-value categories)."""
    widths, current = [], None
    for ev in _cif.parse(source.lines(), _cif.ALL):
        etype = type(ev)
        if etype is _cif.Row:
            if current is not None:
                for i, value in enumerate(ev.values):
                    tok = _cif.quote(value)
                    if tok[0] != ';' and len(tok) > current[i]:
                        current[i] = len(tok)
        elif etype is _cif.Loop:
            current = None if ev.pairs else [0] * len(ev.columns)
            widths.append(current)
    return widths


def _comments(text):
    """The comment lines of a Verbatim chunk worth keeping."""
    out = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped == '#':
            continue
        out.append(line.rstrip() + '\n')
    return out


def run(fhandle, scheme=None):
    """Yield the lines of the re-formatted file."""
    source = _atoms.as_source(fhandle)
    widths = _widths(source)
    n = -1
    in_block = False
    for ev in _cif.parse(source.lines(), _cif.ALL):
        etype = type(ev)
        if etype is _cif.Row:
            if loop.pairs:
                width = max(len(t) for t in loop.tags) + 1
                for tag, value in zip(loop.tags, ev.values):
                    tok = _cif.quote(value)
                    if tok[0] == ';':
                        yield tag + '\n' + tok + '\n'
                    else:
                        yield tag.ljust(width) + tok + '\n'
                continue
            line = ''
            for i, value in enumerate(ev.values):
                tok = _cif.quote(value)
                if tok[0] == ';':
                    if line:
                        yield line.rstrip() + '\n'
                        line = ''
                    yield tok + '\n'
                else:
                    line += tok.ljust(width_of[i]) + ' '
            if line:
                yield line.rstrip() + '\n'
        elif etype is _cif.Loop:
            n += 1
            loop, width_of = ev, widths[n]
            yield '#\n'
            if not loop.pairs:
                yield 'loop_\n'
                for tag in loop.tags:
                    yield tag + '\n'
        elif etype is _cif.DataBlock:
            if in_block:
                yield '#\n'
            yield 'data_{}\n'.format(ev.name)
            in_block = True
        elif etype is _cif.Verbatim:
            for line in _comments(ev.text):
                yield line
    if in_block:
        yield '#\n'


def main():
    fhandle, scheme = check_input(sys.argv[1:])
    _cli.output_text(_cli.lazy(run, fhandle, scheme), fhandle)


if __name__ == '__main__':
    main()
