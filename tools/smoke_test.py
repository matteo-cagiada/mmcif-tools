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
Checks an installation of mmcif-tools: every command answers --help, and a few
pipelines give the expected results on a small structure written here.

Run it after `python -m pip install .`; it uses the installed commands (no
network, no other files). The exit status is 1 if a check fails.

Usage:
    python tools/smoke_test.py
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

STRUCTURE = """data_TEST
_entry.id TEST
loop_
_entity.id
_entity.type
1 polymer
2 non-polymer
3 water
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   1 N N   . MET A 1 1 ? 0.000 0.000 0.000 1.00 10.00 10 MET A N   1
ATOM   2 C CA  . MET A 1 1 ? 1.458 0.000 0.000 1.00 10.00 10 MET A CA  1
ATOM   3 N N   . GLY A 1 2 ? 2.000 1.200 0.000 1.00 10.00 11 GLY A N   1
ATOM   4 C CA  . GLY A 1 2 ? 3.400 1.400 0.000 1.00 10.00 11 GLY A CA  1
ATOM   5 N N   . SER A 1 3 ? 4.000 2.600 0.000 1.00 10.00 11 SER A N   1
ATOM   6 O OG  . SER A 1 3 ? 5.000 3.000 0.000 1.00 10.00 11 SER A OG  1
HETATM 7 S S   . SO4 B 2 . ? 8.000 0.000 0.000 1.00 20.00 201 SO4 A S   1
HETATM 8 O O   . HOH C 3 . ? 9.000 9.000 9.000 1.00 30.00 301 HOH A O   1
"""
# residue 11 appears twice (GLY, SER): mmcif_fixinsert below renumbers SER

CHECKS = [
    # (commands, expected text in the output)
    ([['mmcif_wc']], 'No. atoms:\t6'),
    ([['mmcif_delentity', '-non-polymer,water'], ['mmcif_wc', '-a']], 'No. atoms:\t6'),
    ([['mmcif_selchain', '-A'], ['mmcif_selatom', '-CA'], ['mmcif_tofasta']], '>TEST|A\nMG\n'),
    ([['mmcif_reres', '-1'], ['mmcif_cat', '-_atom_site.auth_seq_id']], 'auth_seq_id\n1\n1\n2\n2\n3\n3\n'),
    ([['mmcif_fixinsert'], ['mmcif_cat', '-_atom_site.auth_seq_id']], '11\n12\n12\n201'),
    ([['mmcif_info']], 'Entity:\t2\tnon-polymer'),
    ([['mmcif_delmeta'], ['mmcif_info']], 'Entity:\t1\tpolymer'),
    ([['mmcif_topdb']], 'TER       7      SER A  11'),
    ([['mmcif_tidy'], ['mmcif_validate']], 'the file looks valid'),
]


def commands():
    """The mmcif_* commands listed in pyproject.toml (or installed next to python)."""
    path = os.path.join(ROOT, 'pyproject.toml')
    if os.path.exists(path):
        with open(path, encoding='utf-8') as fh:
            table = fh.read().split('[project.scripts]', 1)[1].split('\n[', 1)[0]
        return re.findall(r'^(mmcif_\w+)\s*=', table, flags=re.M)
    return sorted({n.split('.')[0] for n in os.listdir(os.path.dirname(sys.executable)) if n.startswith('mmcif_')})


def run(argv, text):
    exe = shutil.which(argv[0])
    if exe is None:
        return 127, '', '{} is not installed (run: python -m pip install .)'.format(argv[0])
    proc = subprocess.run([exe] + argv[1:], input=text.encode(), stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    return proc.returncode, proc.stdout.decode(), proc.stderr.decode()


def main():
    failed = 0
    names = commands()
    for name in names:
        code, out, err = run([name, '--help'], '')
        if code != 0 or 'Usage:' not in out:
            print('FAILED  {} --help {}'.format(name, err.strip()))
            failed += 1
    print('{} commands answer --help'.format(len(names) - failed))
    tmp = tempfile.mkdtemp()
    try:
        for steps, expected in CHECKS:
            text, label = STRUCTURE, ' | '.join(' '.join(s) for s in steps)
            for step in steps:
                code, text, err = run(step, text)
                if code != 0 and step[0] != 'mmcif_validate':
                    break
            if expected not in text:
                print('FAILED  {}: expected {!r}, got {!r} {}'.format(label, expected, text[-200:], err[-200:]))
                failed += 1
            else:
                print('ok      {}'.format(label))
    finally:
        shutil.rmtree(tmp)
    print('{} failure(s)'.format(failed))
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
