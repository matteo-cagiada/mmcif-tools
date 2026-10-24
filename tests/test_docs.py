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
"""The documentation is in step with the tools."""

import importlib.util
import os
import subprocess
import sys
import unittest

from support import ROOT  # noqa: E402

from mmciftools import __version__  # noqa: E402


def make_docs():
    path = os.path.join(ROOT, 'tools', 'make_docs.py')
    spec = importlib.util.spec_from_file_location('make_docs', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@unittest.skipUnless(os.path.isfile(os.path.join(ROOT, 'tools', 'make_docs.py')),
                     'tools/ is not part of the source distribution')
class TestDocs(unittest.TestCase):

    def test_tool_reference_is_up_to_date(self):
        docs = make_docs()
        with open(docs.OUTPUT, encoding='utf-8') as fh:
            self.assertEqual(fh.read(), docs.render(),
                             'docs/tools.md is out of date: run python tools/make_docs.py')

    def test_readme_lists_every_tool(self):
        with open(os.path.join(ROOT, 'README.md'), encoding='utf-8') as fh:
            readme = fh.read()
        missing = [t for t in make_docs().console_scripts() if '`{}`'.format(t) not in readme]
        self.assertEqual(missing, [])

    def test_every_tool_answers_help(self):
        for tool in make_docs().console_scripts():
            with self.subTest(tool=tool):
                proc = subprocess.run([sys.executable, '-m', 'mmciftools.' + tool, '--help'],
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=ROOT)
                self.assertEqual(proc.returncode, 0)
                self.assertIn(b'Usage:', proc.stdout)
                self.assertIn(b'This program is part of the mmcif-tools suite.', proc.stdout)

    def test_versions_agree(self):
        with open(os.path.join(ROOT, 'CITATION.cff'), encoding='utf-8') as fh:
            self.assertIn('version: {}\n'.format(__version__), fh.read())
        with open(os.path.join(ROOT, 'CHANGELOG.md'), encoding='utf-8') as fh:
            self.assertIn('## {} '.format(__version__), fh.read())


if __name__ == '__main__':
    unittest.main()
