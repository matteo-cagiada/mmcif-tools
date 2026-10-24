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
Command-line plumbing shared by the mmcif_* tools.

It follows pdb-tools' conventions: a tool takes at most one option (which
starts with '-') and at most one file; without a file it reads stdin, and
without either it prints its usage. Errors go to stderr as `ERROR!! ...`
followed by the usage, with exit code 1.

On top of that every tool accepts `--auth` / `--label` to choose which
identifiers (chain, residue number, ...) it uses (the default is auth, or
$MMCIFTOOLS_SCHEME), and `--help`. Input files may be gzip-compressed.
"""

import os
import sys

from mmciftools import _cif


class UsageError(Exception):
    """A command-line mistake; reported with the tool's usage text."""


def fail(doc, message):
    """Print an error and the usage text, then exit with status 1."""
    sys.stderr.write('ERROR!! {}\n'.format(message))
    sys.stderr.write(doc)
    sys.exit(1)


def parse_args(args, doc, takes_option=True, needs_option=False):
    """Parse `[-option] [file]` plus `--auth`/`--label`.

    Returns (fhandle, option, scheme), where `option` is the text after the
    leading '-' (or None) and `scheme` is 'auth' or 'label'.
    """
    scheme = None
    rest = []
    for arg in args:
        if arg == '--help':  # not -h: that is an option (mmcif_wc -h, chain h, ...)
            sys.stdout.write(doc)
            sys.exit(0)
        if arg in ('--auth', '--label'):
            scheme = arg[2:]
        else:
            rest.append(arg)
    try:
        scheme = _cif.get_scheme(scheme)
    except ValueError as err:
        fail(doc, err)

    option, path = None, None
    for arg in rest:
        if arg.startswith('-') and len(arg) > 1 and option is None and path is None:
            option = arg[1:]
        elif path is None and not (arg.startswith('-') and len(arg) > 1):
            path = arg
        else:
            fail(doc, 'Unexpected argument: {!r}'.format(arg))

    if option is not None and not takes_option:
        fail(doc, 'This tool takes no options: {!r}'.format('-' + option))
    if option is None and needs_option:
        fail(doc, 'An option is required.')

    if path is None:
        if sys.stdin.isatty():
            if option is None and not needs_option and not rest:
                sys.stderr.write(doc)
                sys.exit(1)
            fail(doc, 'No data to process!')
        return _cif.stdin(), option, scheme
    if not os.path.isfile(path):
        fail(doc, 'File not found or not readable: {!r}'.format(path))
    return _cif.open_input(path), option, scheme


def parse_files(args, doc, minimum=1, takes_option=False):
    """Parse `[-option] file [file ...]` plus `--auth`/`--label`, for tools
    that read several files (they do not read stdin).

    Returns (paths, option, scheme).
    """
    scheme, option, paths = None, None, []
    for arg in args:
        if arg == '--help':
            sys.stdout.write(doc)
            sys.exit(0)
        if arg in ('--auth', '--label'):
            scheme = arg[2:]
        elif arg.startswith('-') and len(arg) > 1 and option is None and not paths:
            option = arg[1:]
        else:
            paths.append(arg)
    try:
        scheme = _cif.get_scheme(scheme)
    except ValueError as err:
        fail(doc, err)
    if option is not None and not takes_option:
        fail(doc, 'This tool takes no options: {!r}'.format('-' + option))
    if not paths and not args:
        sys.stderr.write(doc)
        sys.exit(1)
    if len(paths) < minimum:
        fail(doc, 'This tool needs at least {} file{}.'.format(minimum, 's' if minimum > 1 else ''))
    for path in paths:
        if not os.path.isfile(path):
            fail(doc, 'File not found or not readable: {!r}'.format(path))
    return paths, option, scheme


def write_file(path, events):
    """Write an event stream to a new file (used by the split tools)."""
    with open(path, 'w', encoding=_cif.ENCODING, errors=_cif.ERRORS, newline='') as out:
        _cif.write(events, out)


def output(events, fhandle=None):
    """Write an event stream to stdout and exit.

    CIFErrors in the input are reported (exit 1); a closed pipe, e.g. from
    `| head`, ends the tool quietly.
    """
    out = _cif.stdout()
    try:
        _cif.write(events, out)
        out.flush()
    except _cif.CIFError as err:
        sys.stderr.write('ERROR!! {}\n'.format(err))
        sys.exit(1)
    except OSError as err:
        if not _cif.pipe_closed(err):
            raise
        _cif.silence_stdout()
    finally:
        if fhandle is not None:
            fhandle.close()
    sys.exit(0)


def lazy(func, *args):
    """Iterate over func(*args), calling it only when iteration starts (so
    that output_text reports the errors it raises)."""
    for item in func(*args):
        yield item


def output_text(lines, fhandle=None):
    """Like output(), for tools that print plain text (mmcif_wc, mmcif_tofasta)."""
    out = _cif.stdout()
    try:
        for line in lines:
            out.write(line)
        out.flush()
    except _cif.CIFError as err:
        sys.stderr.write('ERROR!! {}\n'.format(err))
        sys.exit(1)
    except OSError as err:
        if not _cif.pipe_closed(err):
            raise
        _cif.silence_stdout()
    finally:
        if fhandle is not None:
            fhandle.close()
    sys.exit(0)
