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
Deletes insertion codes and shifts the residue numbering to make room.

The same tool as mmcif_fixinsert, under the name of pdb_delinsertion (which
pdb-tools keeps as a deprecated alias of pdb_fixinsert). See mmcif_fixinsert
for the details.

Usage:
    mmcif_delinsertion [-<option>] <cif file>

Example:
    mmcif_delinsertion 1abc.cif            # deletes ALL insertion codes
    mmcif_delinsertion -A9,B12 1abc.cif    # deletes insertion codes of residue 9
                                           # of chain A and 12 of chain B

This program is part of the mmcif-tools suite.
"""

from mmciftools.mmcif_fixinsert import check_input, run  # noqa: F401
from mmciftools import mmcif_fixinsert


def main():
    mmcif_fixinsert.main(__doc__)


if __name__ == '__main__':
    main()
