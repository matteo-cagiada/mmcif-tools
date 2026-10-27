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
Writes docs/tools.md, the tool reference, from the tools' --help texts.

Every console script of pyproject.toml must be in one of the GROUPS below.
CI checks that docs/tools.md is up to date (--check), so run this after
changing a help text or adding a tool.

Usage:
    python tools/make_docs.py           # writes docs/tools.md
    python tools/make_docs.py --check   # exit 1 if docs/tools.md is out of date
"""

import importlib
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT = os.path.join(ROOT, 'docs', 'tools.md')
sys.path.insert(0, ROOT)

GROUPS = [
    ('Fetch / convert', ['fetch', 'topdb', 'frompdb']),
    ('Select / delete', ['selchain', 'delchain', 'selres', 'delres', 'selresname', 'delresname',
                         'selatom', 'selelem', 'delelem', 'selhetatm', 'delhetatm', 'selentity',
                         'delentity', 'selaltloc', 'selmodel', 'head', 'delmeta', 'keepcoord', 'intersect']),
    ('Edit', ['chain', 'rplchain', 'reres', 'shiftres', 'fixinsert', 'delinsertion', 'reatom',
              'rplresname', 'uniqname', 'element', 'b', 'occ']),
    ('Identifiers', ['label2auth', 'auth2label']),
    ('Models, files and assemblies', ['splitchain', 'splitmodel', 'mkensemble', 'chkensemble',
                                      'merge', 'assembly']),
    ('Categories and data blocks', ['cat', 'delcat', 'keepcat', 'blocks']),
    ('Inspect / format', ['info', 'wc', 'tofasta', 'gap', 'sort', 'tidy', 'validate']),
]


def console_scripts():
    """Tool names of the [project.scripts] table of pyproject.toml."""
    with open(os.path.join(ROOT, 'pyproject.toml'), encoding='utf-8') as fh:
        text = fh.read()
    table = text.split('[project.scripts]', 1)[1].split('\n[', 1)[0]
    return re.findall(r'^(mmcif_\w+)\s*=', table, flags=re.M)


def render():
    """The text of docs/tools.md."""
    from mmciftools import __version__
    scripts = console_scripts()
    grouped = ['mmcif_' + t for _, tools in GROUPS for t in tools]
    missing = sorted(set(scripts) - set(grouped))
    unknown = sorted(set(grouped) - set(scripts))
    if missing or unknown:
        raise SystemExit('make_docs.py: GROUPS is out of step with pyproject.toml '
                         '(missing: {}; unknown: {})'.format(missing, unknown))

    out = ['# mmcif-tools: tool reference\n\n',
           'Version {}, {} tools. Generated from the tools\' `--help` texts by '
           '`python tools/make_docs.py`; do not edit by hand.\n\n'.format(__version__, len(scripts)),
           '| Group | Tools |\n|---|---|\n']
    for group, tools in GROUPS:
        links = ', '.join('[`mmcif_{0}`](#mmcif_{0})'.format(t) for t in tools)
        out.append('| {} | {} |\n'.format(group, links))
    for group, tools in GROUPS:
        out.append('\n## {}\n'.format(group))
        for tool in tools:
            name = 'mmcif_' + tool
            doc = importlib.import_module('mmciftools.' + name).__doc__.strip('\n')
            summary = doc.strip().split('\n\n', 1)[0].replace('\n', ' ')
            out.append('\n### {}\n\n{}\n\n```text\n{}\n```\n'.format(name, summary, doc.rstrip()))
    return ''.join(out)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    text = render()
    if argv == ['--check']:
        try:
            with open(OUTPUT, encoding='utf-8') as fh:
                current = fh.read()
        except OSError:
            current = None
        if current != text:
            sys.stderr.write('docs/tools.md is out of date: run python tools/make_docs.py\n')
            return 1
        return 0
    if argv:
        sys.stderr.write(__doc__)
        return 1
    os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
    with open(OUTPUT, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write(text)
    print('wrote {}'.format(os.path.relpath(OUTPUT, ROOT)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
