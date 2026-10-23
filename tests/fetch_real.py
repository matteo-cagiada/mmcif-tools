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
Download the real PDB entries listed in tests/real_set.txt into tests/real/.

Files already present are skipped. Standard library only.

Usage:
    python tests/fetch_real.py
"""

import os
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
MANIFEST = os.path.join(HERE, 'real_set.txt')
TARGET = os.path.join(HERE, 'real')
URL = 'https://files.rcsb.org/download/{}.cif.gz'


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
            with urllib.request.urlopen(URL.format(pdb_id), timeout=60) as resp, \
                    open(tmp, 'wb') as out:
                out.write(resp.read())
        except OSError as err:
            sys.stderr.write('ERROR fetching {}: {}\n'.format(pdb_id, err))
            if os.path.exists(tmp):
                os.remove(tmp)
            continue
        os.replace(tmp, path)
        print('got   {} ({:.2f} MB)'.format(pdb_id, os.path.getsize(path) / 1e6))


if __name__ == '__main__':
    main()
