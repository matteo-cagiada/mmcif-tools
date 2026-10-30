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
Writes the pages of the documentation website (docs/, built by zensical.toml):

- docs/tools.md, the tool reference, from the tools' --help texts. Every
  console script of pyproject.toml must be in one of the GROUPS below.
- the other pages from the sections of README.md (see PAGES) and from
  CHANGELOG.md, so that the README stays the one place to edit that text.

CI checks that the pages are up to date (--check), so run this after
changing a help text, the README or the changelog, or adding a tool.

Usage:
    python tools/make_docs.py           # writes the pages in docs/
    python tools/make_docs.py --check   # exit 1 if a page is out of date
"""

import importlib
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS = os.path.join(ROOT, 'docs')
OUTPUT = os.path.join(DOCS, 'tools.md')
REPOSITORY = 'https://github.com/matteo-cagiada/mmcif-tools'
SITE = 'https://matteo-cagiada.github.io/mmcif-tools/'  # the website: the main page and app/ (web/), and docs/
DOCUMENTATION = SITE + 'docs/'  # the pages written here (zensical.toml)
sys.path.insert(0, ROOT)

# The pages made from README.md: (file, title, the README sections it holds). A section is a "## " heading of
# the README; '' is the text before the first one. Every section must be on a page. A section named like its
# page goes under the page title, its headings one level up; the others keep their heading.
PAGES = [
    ('index.md', None, ['', 'Tools']),
    ('installation.md', 'Installation', ['Installation', 'Extra tools']),
    ('usage.md', 'Usage examples and tips', ['Usage examples and tips']),
    ('mmcif.md', 'Things to know about mmCIF', ['Things to know about mmCIF']),
    ('pdb-tools.md', 'Coming from pdb-tools', ['Coming from pdb-tools']),
    ('about.md', 'Citation and license', ['Citation and changes', 'Acknowledgements', 'License']),
]
GENERATED = ('<!-- Generated from {} by tools/make_docs.py: edit {} and run python tools/make_docs.py. -->\n\n')

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


def read(name):
    with open(os.path.join(ROOT, name), encoding='utf-8') as fh:
        return fh.read()


def code_lines(text):
    """(line, inside a fenced code block) for each line of text; the fences count as code."""
    fenced = False
    for line in text.splitlines(True):
        if line.lstrip().startswith('```'):
            fenced = not fenced
            yield line, True
        else:
            yield line, fenced


def slug(heading):
    """The anchor of a heading, as Zensical (and GitHub, for the README's headings) makes it."""
    text = re.sub(r'[^\w\s-]', '', heading).strip().lower()
    return re.sub(r'[-\s]+', '-', text)


def readme_sections(text):
    """{heading: text} of the "## " sections of the README; '' is the text before the first one."""
    sections, current = {'': []}, ''
    for line, code in code_lines(text):
        match = None if code else re.match(r'## (.+)\n', line)
        if match:
            current = match.group(1).strip()
            sections[current] = []
        else:
            sections[current].append(line)
    return {name: ''.join(lines).strip('\n') + '\n' for name, lines in sections.items()}


def promote(text):
    """text with its headings one level up ("### " becomes "## ")."""
    return ''.join(line[1:] if not code and re.match(r'#{3,} ', line) else line for line, code in code_lines(text))


def relink(text, anchors, page):
    """text with the README's links made to work on the site: other pages, or files of the repository."""
    def link(match):
        url = match.group(1)
        if url.startswith(DOCUMENTATION):  # a page of the documentation (so that the build checks it): its file
            path, sep, anchor = url[len(DOCUMENTATION):].partition('#')
            return '](' + (path.strip('/') or 'index') + '.md' + sep + anchor + ')'
        if re.match(r'[a-z]+:', url):
            return match.group(0)
        if url.startswith('#'):
            if url[1:] not in anchors:
                raise SystemExit('make_docs.py: no heading for the link {} of the README'.format(url))
            target = anchors[url[1:]]
            return '](' + (target[len(page):] if target.startswith(page + '#') else target) + ')'
        path, sep, anchor = url.partition('#')
        if path.startswith('docs/'):
            path = path[len('docs/'):]
        elif path in ('README.md', 'CHANGELOG.md'):
            path = {'README.md': 'index.md', 'CHANGELOG.md': 'changelog.md'}[path]
        else:
            return '](' + REPOSITORY + '/blob/main/' + url + ')'
        return '](' + path + sep + anchor + ')'
    return ''.join(line if code else re.sub(r'\]\(([^)\s]+)\)', link, line) for line, code in code_lines(text))


def site_pages():
    """{file name: text} of the pages made from README.md and CHANGELOG.md."""
    sections = readme_sections(read('README.md'))
    placed = [s for _, _, names in PAGES for s in names]
    missing = [s for s in sections if s not in placed]
    unknown = [s for s in placed if s not in sections]
    if missing or unknown:
        raise SystemExit('make_docs.py: PAGES is out of step with the sections of README.md '
                         '(missing: {}; unknown: {})'.format(missing, unknown))
    anchors = {}  # the README's anchors, as links on the site
    for page, title, names in PAGES:
        for section in filter(None, names):
            anchors[slug(section)] = page if section == title else '{}#{}'.format(page, slug(section))
            for line, code in code_lines(sections[section]):
                match = None if code else re.match(r'#{3,} (.+)\n', line)
                if match:
                    anchors[slug(match.group(1))] = '{}#{}'.format(page, slug(match.group(1)))
    pages = {}
    for page, title, names in PAGES:
        parts = [GENERATED.format('README.md', 'README.md')] + (['# {}\n\n'.format(title)] if title else [])
        for section in names:
            text = sections[section]
            if section == title:
                text = promote(text)
            elif section:
                text = '## {}\n\n{}'.format(section, text)
            parts.append(relink(text, anchors, page) + '\n')
        pages[page] = ''.join(parts).rstrip('\n') + '\n'
    pages['changelog.md'] = (GENERATED.format('CHANGELOG.md', 'CHANGELOG.md') +
                             relink(read('CHANGELOG.md'), anchors, 'changelog.md'))
    return pages


def outputs():
    """{path: text} of every page this script writes."""
    pages = {OUTPUT: render()}
    pages.update((os.path.join(DOCS, name), text) for name, text in site_pages().items())
    return pages


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    pages = outputs()
    if argv == ['--check']:
        stale = []
        for path, text in pages.items():
            try:
                with open(path, encoding='utf-8') as fh:
                    current = fh.read()
            except OSError:
                current = None
            if current != text:
                stale.append(os.path.relpath(path, ROOT).replace(os.sep, '/'))
        if stale:
            sys.stderr.write('{} out of date: run python tools/make_docs.py\n'.format(', '.join(stale)))
            return 1
        return 0
    if argv:
        sys.stderr.write(__doc__)
        return 1
    os.makedirs(DOCS, exist_ok=True)
    for path, text in pages.items():
        with open(path, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write(text)
        print('wrote {}'.format(os.path.relpath(path, ROOT)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
