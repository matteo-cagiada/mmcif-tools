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
Deletes insertion codes and shifts the residue numbering to make room.

The same tool as cif_fixinsert, under the name of pdb_delinsertion (which
pdb-tools keeps as a deprecated alias of pdb_fixinsert). See cif_fixinsert
for the details.

Usage:
    cif_delinsertion [-<option>] <cif file>

Example:
    cif_delinsertion 1abc.cif            # deletes ALL insertion codes
    cif_delinsertion -A9,B12 1abc.cif    # deletes insertion codes of residue 9
                                         # of chain A and 12 of chain B

This program is part of the cif-tools suite, a port of pdb-tools to mmCIF.
"""

from ciftools.cif_fixinsert import check_input, run  # noqa: F401
from ciftools import cif_fixinsert


def main():
    cif_fixinsert.main(__doc__)


if __name__ == '__main__':
    main()
