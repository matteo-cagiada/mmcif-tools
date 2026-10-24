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
Downloads a structure in mmCIF format and writes it to standard output.

PDB entries (4-character codes, or extended ids such as pdb_00001abc) come
from the RCSB; -biounit downloads the first biological assembly instead of
the deposited coordinates. AlphaFold DB models (AF-<UniProt accession>, e.g.
AF-P69905-F1) are looked up through the AlphaFold DB API, so the latest
model version is used.

Usage:
    mmcif_fetch [-biounit] <pdb code | AlphaFold DB id>

Example:
    mmcif_fetch 1brs > 1brs.cif            # deposited coordinates, all 6 chains
    mmcif_fetch -biounit 1brs > 1brs.cif   # first biological assembly, 2 chains
    mmcif_fetch AF-P69905-F1 > hba.cif     # AlphaFold DB model

This program is part of the mmcif-tools suite.
"""

import gzip
import io
import json
import re
import sys
import urllib.error
import urllib.request

from mmciftools import __version__, _cif, _cli

RCSB = 'https://files.rcsb.org/download/{}.cif.gz'
RCSB_ASSEMBLY = 'https://files.rcsb.org/download/{}-assembly1.cif.gz'
AFDB_API = 'https://alphafold.ebi.ac.uk/api/prediction/{}'
TIMEOUT = 60

PDB_CODE = re.compile(r'[0-9][0-9A-Za-z]{3}\Z')
PDB_EXTENDED = re.compile(r'pdb_[0-9A-Za-z]{8}\Z', re.IGNORECASE)
AFDB_ID = re.compile(r'AF-([0-9A-Z]{6,10})(-F[0-9]+)?\Z', re.IGNORECASE)


class FetchError(Exception):
    pass


def check_input(args):
    """Return (identifier, biounit)."""
    if '--help' in args:
        sys.stdout.write(__doc__)
        sys.exit(0)
    if len(args) == 1:
        biounit, code = False, args[0]
    elif len(args) == 2:
        if args[0] != '-biounit':
            _cli.fail(__doc__, 'Invalid option: {!r}'.format(args[0]))
        biounit, code = True, args[1]
    else:
        sys.stderr.write(__doc__)
        sys.exit(1)
    if not (PDB_CODE.match(code) or PDB_EXTENDED.match(code) or AFDB_ID.match(code)):
        _cli.fail(__doc__, 'Invalid PDB code or AlphaFold DB id: {!r}'.format(code))
    if biounit and AFDB_ID.match(code):
        _cli.fail(__doc__, '-biounit only applies to PDB entries')
    return code, biounit


def _open(url, opener):
    request = urllib.request.Request(url, headers={'User-Agent': 'mmcif-tools/' + __version__})
    try:
        return opener(request, timeout=TIMEOUT)
    except urllib.error.HTTPError as err:
        raise FetchError('Could not fetch {}: HTTP {} {}'.format(url, err.code, err.reason))
    except urllib.error.URLError as err:
        raise FetchError('Could not fetch {}: {}'.format(url, err.reason))


def url_for(code, biounit=False, opener=urllib.request.urlopen):
    """Return (url, gzipped) of the mmCIF file for an identifier."""
    match = AFDB_ID.match(code)
    if match:
        accession = match.group(1).upper()
        with _open(AFDB_API.format(accession), opener) as response:
            entries = json.loads(response.read().decode('utf-8'))
        wanted = code.upper() if match.group(2) else None
        for entry in entries:
            if wanted is None or entry.get('entryId', '').upper() == wanted:
                return entry['cifUrl'], False
        raise FetchError('No AlphaFold DB model {} for {}'.format(code, accession))
    if PDB_EXTENDED.match(code) and code[4:8] == '0000':
        code = code[8:]  # pdb_00001abc -> 1abc
    return (RCSB_ASSEMBLY if biounit else RCSB).format(code.lower()), True


def _checked(stream, url):
    """Yield the lines of stream, reporting a broken download as a CIFError."""
    try:
        for line in stream:
            yield line
    except (OSError, EOFError) as err:  # e.g. a truncated gzip stream
        raise _cif.CIFError('Download of {} failed: {}'.format(url, err))
    finally:
        stream.close()


def fetch(code, biounit=False, opener=urllib.request.urlopen):
    """Return an iterator over the lines of the mmCIF file."""
    url, gzipped = url_for(code, biounit, opener)
    response = _open(url, opener)
    raw = gzip.GzipFile(fileobj=response) if gzipped else response
    text = io.TextIOWrapper(raw, encoding=_cif.ENCODING, errors=_cif.ERRORS, newline='')
    return _checked(text, url)


def main():
    code, biounit = check_input(sys.argv[1:])
    try:
        lines = fetch(code, biounit)
    except (FetchError, ValueError) as err:
        sys.stderr.write('ERROR!! {}\n'.format(err))
        sys.exit(1)
    _cli.output_text(lines)


if __name__ == '__main__':
    main()
