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
Shared mmCIF streaming layer for cif-tools (standard library only).

A CIF file is read as a flat stream of events:

    DataBlock  a `data_` line
    Verbatim   any text the tool did not ask to parse (copied as-is)
    Loop       start of a requested category (columns + raw header text)
    Row        one record of that category
    LoopEnd    end of the category (trailing comments / separators)

Only the categories a tool asks for (by default `_atom_site`) are tokenized;
everything else is passed through untouched. Each Row keeps its raw text, so
rows that a tool does not modify are written back byte-for-byte, and in
modified rows only the changed values are rewritten (the spaces after a
value absorb a change of length, so columns stay aligned).

A typical tool:

    from ciftools import _cif

    def run(fhandle, chains):
        for ev in _cif.parse(fhandle):
            if isinstance(ev, _cif.Row) and ev['auth_asym_id'] not in chains:
                continue
            yield ev

    _cif.write(run(fhandle, {'A'}), _cif.stdout())

Run as a script, this module round-trips a file through the parser and
writer, which must reproduce it exactly:

    python -m ciftools._cif [--all] [file.cif]
"""

import gzip
import io
import os
import re
import sys

__all__ = [
    'ALL', 'CIFError', 'DataBlock', 'Loop', 'LoopEnd', 'Row', 'Verbatim',
    'Writer', 'column', 'get_scheme', 'open_input', 'parse', 'quote',
    'stdin', 'stdout', 'write',
]

# Pass as `categories` to parse every category in the file.
ALL = '*'

_NEW = object()  # Row._orig of rows created by a tool

DEFAULT_CATEGORIES = ('_atom_site',)

ENCODING = 'utf-8'
ERRORS = 'surrogateescape'  # lets arbitrary bytes survive a round-trip

_VERBATIM_CHUNK = 1000  # lines buffered before a Verbatim event is emitted

# Token kinds
VALUE, TAG, KEYWORD = 'value', 'tag', 'keyword'

# One token per match: comment, single-quoted, double-quoted or bare.
# A quoted string ends at the first matching quote followed by whitespace
# (CIF 1.1), so embedded quotes such as 'N,N'-dimethyl' are allowed.
_TOKEN_RE = re.compile(r"""
    \s*(?:
        (?P<comment>\#.*)
      | '(?P<sq>.*?)'(?=\s|$)
      | "(?P<dq>.*?)"(?=\s|$)
      | (?P<bare>\S+)
    )""", re.VERBOSE)

_WHITESPACE = ' \t\r\n\f\v'
_SPECIAL_START = frozenset('_#$\'"[];')
_HAS_WHITESPACE = re.compile(r'\s').search
_BARE_RUN = re.compile(r'\S+').finditer
# Values that never need quotes: no whitespace, quotes or special first
# character, and no '_' (so they cannot be a keyword such as data_x).
_PLAIN = re.compile(r'[A-Za-z0-9.+\-][A-Za-z0-9.+\-*:,()/]*\Z').match
_DQUOTE_BREAK = re.compile(r'"\s').search


class CIFError(Exception):
    """Raised on malformed or unsupported CIF input."""

    def __init__(self, message, lineno=None):
        self.message = message
        self.lineno = lineno
        if lineno is not None:
            message = 'line {}: {}'.format(lineno, message)
        super().__init__(message)


# --------------------------------------------------------------------------
# Events
# --------------------------------------------------------------------------
class Verbatim:
    """Raw text that is copied to the output unchanged."""

    __slots__ = ('text',)

    def __init__(self, text):
        self.text = text

    def __repr__(self):
        return 'Verbatim({!r})'.format(self.text[:40])


class DataBlock:
    """A `data_<name>` line."""

    __slots__ = ('name', 'text')

    def __init__(self, name, text):
        self.name = name
        self.text = text

    def __repr__(self):
        return 'DataBlock({!r})'.format(self.name)


class Loop:
    """Start of a parsed category.

    `columns` are item names without the category prefix (`'auth_asym_id'`
    for `_atom_site.auth_asym_id`). Lookups are case-insensitive, as in CIF.
    `pairs` is True when the category was written as key-value pairs
    instead of a `loop_`.
    """

    __slots__ = ('category', 'columns', 'text', 'pairs', 'lineno', '_index')

    def __init__(self, category, columns, text='', pairs=False, lineno=None):
        self.category = category
        self.columns = list(columns)
        self.text = text
        self.pairs = pairs
        self.lineno = lineno
        self._index = {}
        for i, name in reversed(list(enumerate(self.columns))):
            self._index[name] = i
            self._index[name.lower()] = i

    def index(self, column):
        """Position of `column` in each row (KeyError if absent)."""
        try:
            return self._index[column]
        except KeyError:
            return self._index[column.lower()]

    def __contains__(self, column):
        return column in self._index or column.lower() in self._index

    @property
    def tags(self):
        return ['{}.{}'.format(self.category, c) for c in self.columns]

    def __repr__(self):
        return 'Loop({!r}, {} columns)'.format(self.category, len(self.columns))


class Row:
    """One record of a parsed category.

    Read and modify values by column name: `row['B_iso_or_equiv'] = '10.00'`.
    Values are plain strings; CIF nulls are the strings '.' and '?'.

    The values of a row read from a file are a tuple, so they cannot be
    changed behind the writer's back (`row.values[3] = 'x'` raises
    TypeError); the first item assignment turns them into a list.
    """

    # _orig: None while untouched, the original values tuple once modified,
    # or _NEW for a row that was never read from a file.
    __slots__ = ('loop', 'values', 'lead', 'text', 'lineno', '_orig')

    def __init__(self, loop, values, lead='', text=None, lineno=None):
        self.loop = loop
        self.lead = lead  # whitespace/comments preceding the first value
        self.lineno = lineno
        if text is None:
            self.values = list(values)
            self.text = ''
            self._orig = _NEW
        else:
            self.values = values if type(values) is tuple else tuple(values)
            self.text = text
            self._orig = None

    def __getitem__(self, column):
        try:
            return self.values[self.loop._index[column]]
        except KeyError:
            return self.values[self.loop.index(column)]

    def __setitem__(self, column, value):
        if not isinstance(value, str):
            value = str(value)
        i = self.loop.index(column)
        if self._orig is None:
            self._orig = self.values
            self.values = list(self.values)
        self.values[i] = value

    def __contains__(self, column):
        return column in self.loop

    def get(self, column, default=None):
        try:
            return self[column]
        except KeyError:
            return default

    def as_dict(self):
        return dict(zip(self.loop.columns, self.values))

    @property
    def modified(self):
        return bool(self._changed())

    def _changed(self):
        """Sorted indices of the values that differ from the file (all if new)."""
        orig = self._orig
        if orig is None:
            return []
        values = self.values
        if orig is _NEW or len(values) != len(orig):
            return list(range(len(values)))
        return [i for i, (new, old) in enumerate(zip(values, orig)) if new is not old and new != old]

    def __repr__(self):
        return 'Row({!r})'.format(self.values)


class LoopEnd:
    """End of a parsed category; `text` holds trailing comments/blank lines.

    `nrows` is the number of rows the reader found (None for events a tool
    created); the writer uses it to keep the header of a loop that was
    already empty in the input.
    """

    __slots__ = ('loop', 'text', 'nrows')

    def __init__(self, loop, text='', nrows=None):
        self.loop = loop
        self.text = text
        self.nrows = nrows

    def __repr__(self):
        return 'LoopEnd({!r})'.format(self.loop.category)


# --------------------------------------------------------------------------
# Tokenizer
# --------------------------------------------------------------------------
def _is_keyword(token):
    low = token.lower()
    return low.startswith(('data_', 'save_')) or low in ('loop_', 'global_', 'stop_')


def _is_blank(text):
    """True if `text` holds only whitespace and/or a comment."""
    text = text.lstrip()
    return not text or text[0] == '#'


def _category(tag):
    """Category of a tag (or of a stripped line starting with one)."""
    tag = tag.split(None, 1)[0]
    dot = tag.find('.')
    return tag[:dot] if dot > 0 else tag


def _tokens(line, lineno=None):
    """Yield (start, end, value, kind) for each token of a single line.

    `line` must not open a text field (a `;` in column 0); those are handled
    by the callers. Comments end the line. `lineno` is used in errors.
    """
    for m in _TOKEN_RE.finditer(line):
        group = m.lastgroup
        if group == 'comment':
            return
        if group == 'bare':
            tok = m.group('bare')
            first = tok[0]
            if first == '_':
                kind = TAG
            elif first == "'" or first == '"':
                raise CIFError('unterminated quoted string {!r}'.format(tok), lineno)
            elif _is_keyword(tok):
                kind = KEYWORD
            else:
                kind = VALUE
            yield m.start('bare'), m.end('bare'), tok, kind
        else:
            yield m.start(group) - 1, m.end(group) + 1, m.group(group), VALUE


def _unquote_tokens(parts):
    """Unquote whitespace-free tokens such as 'x' or "O5'".

    Returns None if a token starts with a quote but does not end with the
    same quote (a quoted value containing whitespace, or an error); the
    caller then uses the full tokenizer.
    """
    out = []
    for tok in parts:
        first = tok[0]
        if first == "'" or first == '"':
            if len(tok) < 2 or tok[-1] != first:
                return None
            out.append(tok[1:-1])
        else:
            out.append(tok)
    return out


def _strip_eol(text):
    if text.endswith('\r\n'):
        return text[:-2]
    if text.endswith('\n') or text.endswith('\r'):
        return text[:-1]
    return text


def _split_lines(text):
    """Split on '\\n' only, keeping line endings (unlike str.splitlines)."""
    parts = text.split('\n')
    lines = [p + '\n' for p in parts[:-1]]
    if parts[-1]:
        lines.append(parts[-1])
    return lines


def _token_spans(text):
    """(start, end) offsets of every token in the raw text of a row."""
    if ("'" not in text and '"' not in text and '#' not in text
            and text[:1] != ';' and '\n;' not in text):
        return [m.span() for m in _BARE_RUN(text)]  # nothing but bare values
    spans = []
    lines = _split_lines(text)
    pos = 0  # offset of lines[i] in text
    i = 0
    while i < len(lines):
        line = lines[i]
        base = pos
        if line[:1] == ';':  # text field: spans up to the closing ';'
            j = i + 1
            end = pos + len(line)
            while j < len(lines) and lines[j][:1] != ';':
                end += len(lines[j])
                j += 1
            if j == len(lines):
                raise CIFError('unterminated text field')
            spans.append((pos, end + 1))
            line = lines[j][1:]
            base = end + 1
            pos = end
            i = j
        for start, end, _, _ in _tokens(line):
            spans.append((base + start, base + end))
        pos += len(lines[i])
        i += 1
    return spans


# --------------------------------------------------------------------------
# Reader
# --------------------------------------------------------------------------
class _Lines:
    """Line iterator with push-back and line counting."""

    __slots__ = ('_it', '_back', 'lineno')

    def __init__(self, fhandle):
        self._it = iter(fhandle)
        self._back = []
        self.lineno = 0

    def next(self):
        if self._back:
            self.lineno += 1
            return self._back.pop()
        line = next(self._it, None)
        if line is not None:
            self.lineno += 1
        return line

    def push(self, line):
        self._back.append(line)
        self.lineno -= 1


def parse(fhandle, categories=DEFAULT_CATEGORIES):
    """Stream the events of a CIF file.

    `fhandle` is any iterable of lines (keep line endings, e.g. open the file
    with `newline=''` or use `open_input`). `categories` lists the categories
    to parse into Loop/Row/LoopEnd events, e.g. ('_atom_site',
    '_atom_site_anisotrop'); pass ALL to parse every category.
    """
    if isinstance(categories, str):
        categories = (categories,)
    wanted = set()
    for cat in categories:
        cat = cat.lower()
        wanted.add(cat if cat == ALL or cat.startswith('_') else '_' + cat)
    return _Reader(fhandle, wanted).events()


class _Reader:

    def __init__(self, fhandle, wanted):
        self.lines = _Lines(fhandle)
        self.everything = ALL in wanted
        self.wanted = wanted

    def _wanted(self, category):
        return self.everything or category.lower() in self.wanted

    def _text_field(self, first):
        """Read a text field opened by `first`; return (lines, closing_line, value)."""
        start = self.lines.lineno
        parts = [first]
        while True:
            line = self.lines.next()
            if line is None:
                raise CIFError('unterminated text field', start)
            if line[:1] == ';':
                return parts, line, _strip_eol(''.join(parts)[1:])
            parts.append(line)

    def events(self):
        lines = self.lines
        buf = []
        while True:
            line = lines.next()
            if line is None:
                break

            if line[:1] == ';':  # text field of an unparsed category
                parts, closing, _ = self._text_field(line)
                buf.extend(parts)
                buf.append(closing)
                continue

            s = line.lstrip()
            first = s[:1]
            if first == '_':
                if self._wanted(_category(s)):
                    if buf:
                        yield Verbatim(''.join(buf))
                        buf = []
                    for ev in self._pairs(line):
                        yield ev
                    continue

            elif first and first in 'lLdD':
                keyword = s.split(None, 1)[0].lower()
                if keyword == 'loop_':
                    hold, tag_line = self._loop_lookahead(line, s)
                    if tag_line is None:
                        buf.extend(hold)
                        continue
                    if buf:
                        yield Verbatim(''.join(buf))
                        buf = []
                    for ev in self._loop(hold, tag_line):
                        yield ev
                    continue
                if keyword.startswith('data_'):
                    if buf:
                        yield Verbatim(''.join(buf))
                        buf = []
                    yield DataBlock(s.split(None, 1)[0][5:], line)
                    continue

            buf.append(line)
            if len(buf) >= _VERBATIM_CHUNK:
                yield Verbatim(''.join(buf))
                buf = []

        if buf:
            yield Verbatim(''.join(buf))

    def _loop_lookahead(self, line, s):
        """After `loop_`, find the first tag; return (held_lines, tag_line).

        tag_line is None when the loop is not one we parse; the held lines
        are then copied verbatim.
        """
        lines = self.lines
        hold = [line]
        rest = s[5:].strip()
        if rest and rest[0] != '#':
            if rest[0] == '_' and self._wanted(_category(rest)):
                raise CIFError('tags on the same line as loop_ are not supported',
                               lines.lineno)
            return hold, None
        while True:
            nxt = lines.next()
            if nxt is None:
                return hold, None
            ns = nxt.lstrip()
            if not ns or ns[0] == '#':
                hold.append(nxt)
                continue
            if ns[0] == '_' and self._wanted(_category(ns)):
                return hold, nxt
            lines.push(nxt)
            return hold, None

    def _loop(self, hold, line):
        lines = self.lines
        lineno = lines.lineno
        category = _category(line.lstrip())
        cat_low = category.lower()
        prefix_len = len(category) + 1

        # Header: tag lines (plus blank/comment lines between them)
        header = list(hold)
        columns = []
        while line is not None:
            s = line.lstrip()
            if s[:1] == '_':
                if _category(s).lower() != cat_low:
                    break  # a loop without rows, followed by another category
                for _, _, tok, kind in _tokens(line, lines.lineno):
                    if kind != TAG:
                        raise CIFError('values on the same line as loop tags '
                                       'are not supported', lines.lineno)
                    if _category(tok).lower() != cat_low:
                        raise CIFError('mixed categories in one loop', lines.lineno)
                    columns.append(tok[prefix_len:])
                header.append(line)
            elif not s or s[0] == '#':
                header.append(line)
            else:
                break
            line = lines.next()

        loop = Loop(category, columns, ''.join(header), lineno=lineno)
        yield loop

        ncol = len(columns)
        nrows = 0
        vals = []
        lead = []   # text before the first value of the current row
        body = []   # raw text of the current row
        row_line = None

        def feed(seg):
            """Tokenize part of a line; return the rows it completes."""
            nonlocal vals, lead, body, row_line
            done = []
            pos = 0
            for start, end, tok, kind in _tokens(seg, lines.lineno):
                if kind != VALUE:
                    raise CIFError('unexpected {!r} inside loop {}'.format(tok, category),
                                   lines.lineno)
                if vals:
                    body.append(seg[pos:end])
                else:
                    lead.append(seg[pos:start])
                    body.append(seg[start:end])
                    row_line = lines.lineno
                vals.append(tok)
                pos = end
                if len(vals) == ncol:
                    tail = seg[end:]
                    if _is_blank(tail):
                        body.append(tail)
                        pos = len(seg)
                    done.append(Row(loop, vals, ''.join(lead), ''.join(body), row_line))
                    vals, lead, body = [], [], []
            tail = seg[pos:]
            if tail:
                (body if vals else lead).append(tail)
            return done

        while line is not None:
            if line[:1] == ';':
                if not vals:
                    row_line = lines.lineno
                parts, closing, value = self._text_field(line)
                body.extend(parts)
                body.append(';')
                vals.append(value)
                rest = closing[1:]
                if len(vals) == ncol:
                    if _is_blank(rest):
                        body.append(rest)
                        rest = ''
                    nrows += 1
                    yield Row(loop, vals, ''.join(lead), ''.join(body), row_line)
                    vals, lead, body = [], [], []
                for row in feed(rest):
                    nrows += 1
                    yield row
                line = lines.next()
                continue

            s = line.lstrip()
            first = s[:1]
            if first == '_' or (first and first in 'dDlLsSgG'
                                and _is_keyword(s.split(None, 1)[0])):
                if vals:
                    raise CIFError('incomplete row in {}: {} of {} values'.format(
                        category, len(vals), ncol), row_line)
                lines.push(line)
                break

            # Fast path: one complete row on one line, without comments. If
            # it splits into exactly ncol tokens, no value contains
            # whitespace, so a quoted value is a single token such as "O5'".
            if not vals and '#' not in line and ' _' not in line and '\t_' not in line:
                parts = line.split()
                if len(parts) == ncol and ("'" in line or '"' in line):
                    parts = _unquote_tokens(parts)
                if parts is not None and len(parts) == ncol:
                    nrows += 1
                    indent = len(line) - len(s)
                    if indent or lead:
                        lead.append(line[:indent])
                        yield Row(loop, parts, ''.join(lead), line[indent:], lines.lineno)
                        lead = []
                    else:
                        yield Row(loop, parts, '', line, lines.lineno)
                    line = lines.next()
                    continue
                if parts == []:
                    lead.append(line)
                    line = lines.next()
                    continue

            for row in feed(line):
                nrows += 1
                yield row
            line = lines.next()

        if vals:
            raise CIFError('incomplete row in {}: {} of {} values'.format(
                category, len(vals), ncol), row_line)
        yield LoopEnd(loop, ''.join(lead), nrows)

    def _pairs(self, line):
        """Parse a category written as `_cat.item value` pairs into one Row."""
        lines = self.lines
        lineno = lines.lineno
        category = _category(line.lstrip())
        cat_low = category.lower()
        prefix_len = len(category) + 1
        columns, vals = [], []
        raw, pending = [], []  # pending: blank/comment lines not yet attributed
        expect = False  # a tag is waiting for its value

        while line is not None:
            if line[:1] == ';':
                if not expect:
                    raise CIFError('unexpected text field in {}'.format(category), lines.lineno)
                parts, closing, value = self._text_field(line)
                if not _is_blank(closing[1:]):
                    raise CIFError('content after a closing ";" is not supported here',
                                   lines.lineno)
                raw.extend(pending)
                raw.extend(parts)
                raw.append(closing)
                pending = []
                vals.append(value)
                expect = False
            else:
                s = line.lstrip()
                if not s or s[0] == '#':
                    (raw if expect else pending).append(line)
                    line = lines.next()
                    continue
                if not expect and (s[0] != '_' or _category(s).lower() != cat_low):
                    break
                for _, _, tok, kind in _tokens(line, lines.lineno):
                    if kind == TAG:
                        if expect:
                            raise CIFError('no value for {}.{}'.format(category, columns[-1]),
                                           lines.lineno)
                        if _category(tok).lower() != cat_low:
                            raise CIFError('mixed categories on one line', lines.lineno)
                        columns.append(tok[prefix_len:])
                        expect = True
                    elif kind == VALUE and expect:
                        vals.append(tok)
                        expect = False
                    else:
                        raise CIFError('unexpected {!r} in {}'.format(tok, category),
                                       lines.lineno)
                raw.extend(pending)
                raw.append(line)
                pending = []
            line = lines.next()

        if expect:
            raise CIFError('no value for {}.{}'.format(category, columns[-1]), lines.lineno)
        if line is not None:
            lines.push(line)
        loop = Loop(category, columns, '', pairs=True, lineno=lineno)
        yield loop
        yield Row(loop, vals, '', ''.join(raw), lineno)
        yield LoopEnd(loop, ''.join(pending), 1)


# --------------------------------------------------------------------------
# Writer
# --------------------------------------------------------------------------
def quote(value, eol='\n'):
    """Format a value as a CIF token, quoting only when needed.

    '.' and '?' are written bare (CIF null values).
    """
    if not isinstance(value, str):
        value = str(value)
    if _PLAIN(value):  # numbers, names, '.', ...: the common case
        return value
    if value == '?':
        return value
    if not value:
        return "''"
    if '\n' in value or '\r' in value:
        return _text_field(value, eol)
    if not (value[0] in _SPECIAL_START or "'" in value or '"' in value
            or _HAS_WHITESPACE(value) or _is_keyword(value)):
        return value
    if "'" not in value:
        return "'" + value + "'"
    if not _DQUOTE_BREAK(value):
        return '"' + value + '"'
    return _text_field(value, eol)


def _text_field(value, eol):
    if '\n;' in value or '\r;' in value:
        raise CIFError('cannot write a value with a line starting with ";"')
    return ';' + value + eol + ';'


def _eol_of(text):
    """Line ending used by `text` ('' if it does not end with one)."""
    if text.endswith('\r\n'):
        return '\r\n'
    if text.endswith('\n'):
        return '\n'
    return ''


def _values_for(row, loop):
    """Row values in the column order of `loop`."""
    if row.loop is loop:
        return row.values
    return [row.get(c, '?') for c in loop.columns]


def _absorb(gap, delta):
    """Shrink or grow a run of spaces by `delta`, keeping at least one."""
    if not gap or gap.strip(' '):  # tabs or a line break: leave it alone
        return gap
    return gap[:max(1, len(gap) - delta)] if delta > 0 else gap + ' ' * -delta


def _splice_bare(text, values, changed):
    """Fast _splice for rows whose tokens contain no whitespace; None if not
    applicable.

    Then token i is found with str.split(None, i). Changes are applied from
    the last one backwards, which leaves the earlier tokens where they were.
    """
    if '#' in text or text[0] in _WHITESPACE or text[0] == ';' or '\n;' in text:
        return None
    if ("'" in text or '"' in text) and len(text.split()) != len(values):
        return None  # a quoted value contains whitespace
    for i in reversed(changed):
        tok = quote(values[i])
        if tok[0] == ';':
            return None
        rest = text.split(None, i)[i] if i else text  # text from token i on
        start = len(text) - len(rest)
        old_len = len(rest.split(None, 1)[0])
        end = start + old_len
        after = text[end:]
        delta = len(tok) - old_len
        if delta:  # let the spaces after the token absorb the length change
            stripped = after.lstrip(' ')
            gap = len(after) - len(stripped)
            if gap and stripped and stripped[0] not in _WHITESPACE:
                after = ' ' * (max(1, gap - delta) if delta > 0 else gap - delta) + stripped
        text = text[:start] + tok + after
    return text


def _splice(text, spans, values, changed):
    """Rewrite `text`, replacing only the tokens listed in `changed`.

    Everything else (separators, tabs, comments, line breaks) is kept. When a
    token changes length, the run of spaces after it absorbs the difference
    so that the following columns stay aligned.
    """
    eol = '\r\n' if '\r\n' in text else '\n'
    out = []
    prev = 0  # text[prev:] has not been copied yet
    for i in changed:
        start, end = spans[i]
        gap = text[prev:start]
        tok = quote(values[i], eol)
        if tok[0] == ';':  # a text field must start a line
            if gap:
                if not gap.endswith('\n'):
                    gap = gap.rstrip(' \t') + eol
            elif out and not out[-1].endswith('\n'):
                out[-1] = out[-1].rstrip(' \t') + eol
        out.append(gap)
        out.append(tok)
        prev = end
        if tok[0] == ';' or '\n' in text[start:end] or i + 1 == len(spans):
            continue
        delta = len(tok) - (end - start)
        if delta:  # let the spaces after the token absorb the length change
            following = spans[i + 1][0]
            out.append(_absorb(text[end:following], delta))
            prev = following
    out.append(text[prev:])
    return ''.join(out)


def _fresh_row(values):
    """Format the values of a row that has no raw text to start from."""
    out = []
    for value in values:
        tok = quote(value)
        if tok[0] == ';':
            if out and out[-1][-1:] != '\n':
                out.append('\n')
            out.append(tok + '\n')
        else:
            if out and out[-1][-1:] != '\n':
                out.append(' ')
            out.append(tok)
    if out and out[-1][-1:] != '\n':
        out.append('\n')
    return ''.join(out)


def _format_row(row, loop):
    """Re-format a row; only the values that changed are rewritten."""
    values = _values_for(row, loop)
    text = row.text
    if row.loop is loop and text:  # a row read from the file
        orig = row._orig
        if orig is None:
            return text
        if len(orig) == len(values):
            changed = row._changed()
            if not changed:
                return text
            out = _splice_bare(text, values, changed)
            if out is not None:
                return out
            spans = _token_spans(text)
            if len(spans) == len(values):
                return _splice(text, spans, values, changed)
    return _fresh_row(values)


def _format_pairs(row, loop):
    """Re-format a key-value category; only the changed values are rewritten."""
    values = _values_for(row, loop)
    if row.loop is loop and row.text:  # read from the file
        changed = row._changed()
        if not changed:
            return row.text
        spans = _token_spans(row.text)  # tag, value, tag, value, ...
        if len(spans) == 2 * len(values):
            return _splice(row.text, spans[1::2], values, changed)
    width = max(len(t) for t in loop.tags)
    out = []
    for tag, value in zip(loop.tags, values):
        tok = quote(value)
        if tok[0] == ';':
            out.append('{}\n{}\n'.format(tag, tok))
        else:
            out.append('{} {}\n'.format(tag.ljust(width), tok))
    return ''.join(out)


def _format_table(rows, loop):
    """Write rows as a new `loop_` with aligned columns."""
    table = [[quote(v) for v in _values_for(r, loop)] for r in rows]
    widths = [0] * len(loop.columns)
    for toks in table:
        for i, tok in enumerate(toks):
            if tok[0] != ';':
                widths[i] = max(widths[i], len(tok))
    out = ['loop_\n']
    out.extend(t + '\n' for t in loop.tags)
    for toks in table:
        line = ''
        for i, tok in enumerate(toks):
            if tok[0] == ';':  # text fields go on lines of their own
                if line:
                    out.append(line.rstrip(' ') + '\n')
                    line = ''
                out.append(tok + '\n')
            else:
                line += tok.ljust(widths[i]) + ' '
        if line:
            out.append(line.rstrip(' ') + '\n')
    return ''.join(out)


class Writer:
    """Serialize events back to CIF text.

    Unmodified rows are written from their raw text; in modified rows only
    the changed values are rewritten. A loop whose rows were all removed is
    dropped entirely (header included).
    """

    def __init__(self, out):
        self.out = out
        self._last = '\n'
        self._loop = None
        self._nrows = 0
        self._pairs = None

    def _w(self, text):
        if text:
            self.out.write(text)
            self._last = text[-1]

    def _newline(self):
        if self._last != '\n':
            self._w('\n')

    def write(self, ev):
        if isinstance(ev, Row):
            self._row(ev)
        elif isinstance(ev, (Verbatim, DataBlock)):
            self._newline()
            self._w(ev.text)
        elif isinstance(ev, Loop):
            self._loop = ev
            self._nrows = 0
            self._pairs = [] if ev.pairs else None
        elif isinstance(ev, LoopEnd):
            self._end(ev)
        elif isinstance(ev, str):
            self._newline()
            self._w(ev)
        else:
            raise TypeError('cannot write {!r}'.format(ev))

    def _row(self, row):
        loop = self._loop
        if loop is None:
            raise CIFError('row written outside of a loop')
        if self._pairs is not None:
            self._pairs.append(row)
            return
        if not self._nrows:
            self._newline()
            self._w(loop.text)
        self._nrows += 1

        if row._orig is None and row.loop is loop:
            text = row.text  # untouched
        else:
            text = _format_row(row, loop)
        if row.lead:
            if self._last not in _WHITESPACE and row.lead[0] not in _WHITESPACE:
                self._w('\n')
            self._w(row.lead)
        if text[:1] == ';':
            self._newline()
        elif self._last not in _WHITESPACE and text[:1] not in _WHITESPACE:
            self._w('\n')
        self._w(text)

    def _end(self, ev):
        rows = self._pairs
        if (rows is None and not self._nrows and ev.nrows == 0
                and self._loop is not None and ev.loop is self._loop):
            # Empty in the input already (not emptied by a tool): keep as-is.
            self._newline()
            self._w(self._loop.text)
        if rows:
            loop = self._loop
            self._newline()
            if len(rows) > 1:
                self._w(_format_table(rows, loop))
            else:
                self._w(_format_pairs(rows[0], loop))
        self._loop = None
        self._pairs = None
        if ev.text:
            self._newline()
            self._w(ev.text)

    def close(self):
        """Flush a category left open (an event stream missing its LoopEnd)."""
        if self._loop is not None:
            self._end(LoopEnd(self._loop))


def write(events, out):
    """Write an event stream to the text stream `out`."""
    writer = Writer(out)
    for ev in events:
        writer.write(ev)
    writer.close()


# --------------------------------------------------------------------------
# Identifier schemes (auth_* vs label_*)
# --------------------------------------------------------------------------
SCHEMES = ('auth', 'label')
DEFAULT_SCHEME = 'auth'

_FIELDS = {
    'chain': ('auth_asym_id', 'label_asym_id'),
    'resnum': ('auth_seq_id', 'label_seq_id'),
    'resname': ('auth_comp_id', 'label_comp_id'),
    'atom': ('auth_atom_id', 'label_atom_id'),
}


def get_scheme(name=None):
    """Resolve the identifier scheme: explicit name, $CIFTOOLS_SCHEME, or 'auth'."""
    scheme = (name or os.environ.get('CIFTOOLS_SCHEME') or DEFAULT_SCHEME).lower()
    if scheme not in SCHEMES:
        raise ValueError('unknown identifier scheme {!r} (use auth or label)'.format(scheme))
    return scheme


def column(loop, field, scheme=None):
    """Name of the `_atom_site` column holding `field` ('chain', 'resnum',
    'resname' or 'atom') in the given scheme.

    Falls back to the other scheme when the file lacks the preferred column
    (e.g. a minimal CIF with label_* only).
    """
    auth, label = _FIELDS[field]
    first, second = (auth, label) if get_scheme(scheme) == 'auth' else (label, auth)
    if first in loop:
        return first
    if second in loop:
        return second
    raise CIFError('{} has neither {} nor {}'.format(loop.category, first, second))


# --------------------------------------------------------------------------
# I/O helpers
# --------------------------------------------------------------------------
def open_input(path):
    """Open a CIF file (plain or gzip) for streaming, keeping line endings."""
    if path == '-':
        return stdin()
    with open(path, 'rb') as fh:
        magic = fh.read(2)
    if magic == b'\x1f\x8b':
        return gzip.open(path, 'rt', encoding=ENCODING, errors=ERRORS, newline='')
    return open(path, 'r', encoding=ENCODING, errors=ERRORS, newline='')


def _reconfigure(stream):
    try:
        stream.reconfigure(encoding=ENCODING, errors=ERRORS, newline='')
    except (AttributeError, ValueError, io.UnsupportedOperation):
        pass
    return stream


def stdin():
    """sys.stdin set up to keep line endings (byte-exact round-trips)."""
    return _reconfigure(sys.stdin)


def stdout():
    """sys.stdout set up not to translate line endings."""
    return _reconfigure(sys.stdout)


def _main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(
        prog='python -m ciftools._cif',
        description='Round-trip a CIF file through the cif-tools parser and writer.')
    ap.add_argument('file', nargs='?', default='-', help='CIF file (default: stdin)')
    ap.add_argument('-c', '--category', action='append',
                    help='category to parse (repeatable; default: _atom_site)')
    ap.add_argument('--all', action='store_true', help='parse every category')
    args = ap.parse_args(argv)

    categories = ALL if args.all else (args.category or DEFAULT_CATEGORIES)
    fhandle = open_input(args.file)
    out = stdout()
    try:
        write(parse(fhandle, categories), out)
        out.flush()
    except CIFError as err:
        sys.stderr.write('ERROR: {}\n'.format(err))
        sys.exit(1)
    except BrokenPipeError:
        pass
    finally:
        fhandle.close()


if __name__ == '__main__':
    _main()
