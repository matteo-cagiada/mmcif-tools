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
Checks CIF files with gemmi, an independent CIF reader.

This script is not part of mmcif-tools and is not installed with it. It only
runs the gemmi program, which you install yourself, e.g. with
`pip install gemmi-program` or from https://github.com/project-gemmi/gemmi.

For each file it runs `gemmi validate`, which checks the CIF syntax and, with
--dict, also checks the file against a DDL2 dictionary such as the PDBx/mmCIF
dictionary (mmcif_pdbx_v50.dic, from https://mmcif.wwpdb.org): item types,
allowed values, mandatory items and unique keys. --parents adds gemmi's check
that the parent of each child item is present (e.g. that every label_asym_id
of _atom_site is in _struct_asym); files from the PDB often fail it already,
so compare with the file you started from. gemmi's
messages are printed as they are, followed by one OK / FAILED line per file.
The exit status is 0 if every file passes, 1 if one fails, and 2 if gemmi
cannot be run.

Usage:
    python tools/check_with_gemmi.py [--gemmi PATH] [--dict FILE [--parents]] FILE [FILE ...]

Example:
    python tools/check_with_gemmi.py out.cif
    python tools/check_with_gemmi.py --gemmi ~/bin/gemmi --dict mmcif_pdbx_v50.dic out.cif
    mmcif_selchain -A 1abc.cif | python tools/check_with_gemmi.py -

The gemmi program is found with --gemmi, else $GEMMI, else `gemmi` on the
PATH. Use '-' to check standard input (e.g. the end of a pipeline).
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile


def find_gemmi(path=None):
    """Path of the gemmi program, or None."""
    path = path or os.environ.get('GEMMI') or 'gemmi'
    if os.path.dirname(path):
        path = os.path.expanduser(path)
        return path if os.path.isfile(path) and os.access(path, os.X_OK) else None
    return shutil.which(path)


def check(gemmi, path, dictionary=None, parents=False):
    """Run `gemmi validate` on one file; return (passed, output)."""
    command = [gemmi, 'validate']
    if dictionary:
        command += ['-d', dictionary] + (['-p'] if parents else [])
    proc = subprocess.run(command + [path], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          universal_newlines=True)
    return proc.returncode == 0, proc.stdout


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Check CIF files with the gemmi program (installed separately).',
        epilog='See the top of this script for details.')
    parser.add_argument('files', nargs='+', metavar='FILE', help="CIF file(s); '-' for stdin")
    parser.add_argument('--gemmi', metavar='PATH',
                        help='the gemmi program (default: $GEMMI, else gemmi on the PATH)')
    parser.add_argument('--dict', metavar='FILE', dest='dictionary',
                        help='also check against this DDL2 dictionary (e.g. mmcif_pdbx_v50.dic)')
    parser.add_argument('--parents', action='store_true',
                        help='with --dict, also check that parent items are present')
    args = parser.parse_args(argv)
    if args.parents and not args.dictionary:
        parser.error('--parents needs --dict')

    gemmi = find_gemmi(args.gemmi)
    if gemmi is None:
        wanted = args.gemmi or os.environ.get('GEMMI') or 'gemmi'
        sys.stderr.write('ERROR!! gemmi program not found ({}). Install it (e.g. pip install '
                         'gemmi-program) and/or pass its path with --gemmi.\n'.format(wanted))
        return 2
    if args.dictionary and not os.path.isfile(args.dictionary):
        sys.stderr.write('ERROR!! Dictionary not found: {}\n'.format(args.dictionary))
        return 2

    tmpdir = None
    failed = 0
    try:
        for path in args.files:
            name = path
            if path == '-':
                tmpdir = tmpdir or tempfile.mkdtemp()
                path = os.path.join(tmpdir, 'stdin.cif')
                with open(path, 'wb') as out:
                    shutil.copyfileobj(sys.stdin.buffer, out)
                name = '<stdin>'
            elif not os.path.isfile(path):
                sys.stdout.write('FAILED  {}: file not found\n'.format(name))
                failed += 1
                continue
            try:
                passed, output = check(gemmi, path, args.dictionary, args.parents)
            except OSError as err:
                sys.stderr.write('ERROR!! Cannot run {}: {}\n'.format(gemmi, err))
                return 2
            if output:
                sys.stdout.write(output if output.endswith('\n') else output + '\n')
            sys.stdout.write('{:<8}{}\n'.format('OK' if passed else 'FAILED', name))
            failed += not passed
    finally:
        if tmpdir:
            shutil.rmtree(tmpdir)
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
