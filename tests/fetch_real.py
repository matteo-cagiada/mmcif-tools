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
Download the real entries listed in tests/real_set.txt into tests/real/.

PDB ids come from the RCSB; AlphaFold DB ids (AF-...) from the AlphaFold
DB, compressed here. Files already present are skipped. Standard library
only.

Usage:
    python tests/fetch_real.py
"""

import gzip
import os
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from mmciftools.mmcif_fetch import url_for  # noqa: E402

MANIFEST = os.path.join(HERE, 'real_set.txt')
TARGET = os.path.join(HERE, 'real')


def entries():
    with open(MANIFEST) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith('#'):
                yield line.split(None, 1)[0].upper()


def main():
    os.makedirs(TARGET, exist_ok=True)
    for pdb_id in entries():
        path = os.path.join(TARGET, pdb_id.lower() + '.cif.gz')
        if os.path.exists(path):
            print('have  {}'.format(pdb_id))
            continue
        tmp = path + '.part'
        try:
            url, gzipped = url_for(pdb_id)
            with urllib.request.urlopen(url, timeout=60) as resp, open(tmp, 'wb') as out:
                data = resp.read()
                out.write(data if gzipped else gzip.compress(data))
        except Exception as err:  # network errors, FetchError
            sys.stderr.write('ERROR fetching {}: {}\n'.format(pdb_id, err))
            if os.path.exists(tmp):
                os.remove(tmp)
            continue
        os.replace(tmp, path)
        print('got   {} ({:.2f} MB)'.format(pdb_id, os.path.getsize(path) / 1e6))


if __name__ == '__main__':
    main()
