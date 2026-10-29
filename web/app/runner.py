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
Runs mmcif-tools for the web page (index.html), inside Pyodide.

Each tool runs through its main(), as on the command line: sys.argv holds
its arguments, stdout goes to a file and stderr is captured. The steps of a
pipeline pass files instead of pipes: each step reads the output file of
the one before.

The same code runs on CPython with mmcif-tools installed, which is how
README.md compares the two and how CI checks it (the best of N runs is
reported; the exit status is 1 if a step failed):

    python runner.py [--repeat N] '<tool>[:<option>] | <tool> ...' <file> [<file> ...]

Example:
    python runner.py 'selchain:A | delhetatm | tidy' 1abc.cif
"""

import importlib
import io
import json
import os
import pkgutil
import re
import shutil
import sys
import time
import traceback
import zipfile

# The groups of the tool menu: those of the tool reference (GROUPS in
# tools/make_docs.py; a test checks that they agree). A tool missing here is
# listed under 'Other', and one the installed release does not have yet is left out.
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

# Tools that are not "one mmCIF file in, one mmCIF file out".
UNAVAILABLE = {
    'fetch': 'downloads files, which Python cannot do inside the browser: get the file from '
             'the RCSB or the AlphaFold DB and open it here',
}
MULTI_INPUT = {'intersect', 'merge', 'mkensemble'}  # read every opened file: first step only
PDB_INPUT = {'frompdb'}                             # reads a PDB file: first step only
SPLIT = {'splitchain', 'splitmodel'}                # write several files: last step, zipped
NOT_CIF = {'cat': '.tsv', 'chkensemble': '.txt', 'gap': '.txt', 'info': '.txt',
           'tofasta': '.fasta', 'topdb': '.pdb', 'validate': '.txt', 'wc': '.txt'}  # last step only

MAX_STDERR = 100000  # characters of stderr kept per step
INPUT_EXTENSIONS = ('.cif', '.mmcif', '.pdb', '.ent')


def _module(tool):
    return importlib.import_module('mmciftools.mmcif_' + tool)


def _usage(doc):
    """The line after 'Usage:' in a help text."""
    lines = doc.splitlines()
    for i, line in enumerate(lines[:-1]):
        if line.strip() == 'Usage:':
            return lines[i + 1].strip()
    return ''


def _option(usage):
    """(kind, placeholder, label) of the option in a usage line: kind is
    'none', 'optional' or 'required', and label is True if the usage shows
    --label."""
    words, label = [], False
    for word in usage.split()[1:]:
        if word in ('[--label]', '--label'):
            label = True
        elif word.startswith('<'):
            break
        else:
            words.append(word)
    text = ' '.join(words)
    if not text:
        return 'none', '', label
    if not text.startswith('['):
        return 'required', text[1:], label
    # [-<option>], [-<name>[,<name>...]], [-multi] [-seqres], [-chain | -res]
    choices, depth = [], 0
    for char in text:
        if char == '[':
            depth += 1
            if depth == 1:
                choices.append('')
                continue
        elif char == ']':
            depth -= 1
            if depth == 0:
                continue
        if depth > 0:
            choices[-1] += char
    return 'optional', re.sub(r'(^|\| )-', r'\1', ' | '.join(choices)), label


def _describe(tool, group):
    doc = _module(tool).__doc__ or ''
    paragraphs = doc.strip().split('\n\n')
    usage = _usage(doc)
    option, placeholder, label = _option(usage)
    return {
        'name': tool,
        'group': group,
        'summary': ' '.join(paragraphs[0].split()) if paragraphs else '',
        'usage': usage,
        'help': doc.strip('\n'),
        'option': option,
        'placeholder': placeholder,
        'label': label,
        'unavailable': UNAVAILABLE.get(tool),
        'multi_input': tool in MULTI_INPUT,
        'first_only': tool in MULTI_INPUT or tool in PDB_INPUT,
        'last_only': tool in SPLIT or tool in NOT_CIF,
        'split': tool in SPLIT,
        'ext': '.zip' if tool in SPLIT else NOT_CIF.get(tool, '.cif'),
    }


def tools():
    """The tools, in the order of the tool menu, as a JSON list."""
    import mmciftools
    installed = sorted(m.name[6:] for m in pkgutil.iter_modules(mmciftools.__path__)
                       if m.name.startswith('mmcif_'))
    out = [_describe(tool, group) for group, names in GROUPS for tool in names if tool in installed]
    known = {tool for _, names in GROUPS for tool in names}
    out.extend(_describe(tool, 'Other') for tool in installed if tool not in known)
    return json.dumps(out)


def versions():
    """Versions of mmcif-tools and Python, as JSON."""
    import mmciftools
    return json.dumps({'mmciftools': mmciftools.__version__, 'python': sys.version.split()[0]})


def check(steps, n_inputs):
    """Why a pipeline cannot run (None if it can)."""
    if not steps:
        return 'The pipeline is empty.'
    if not n_inputs:
        return 'Open a file first.'
    for i, step in enumerate(steps):
        tool = step['tool']
        name = 'mmcif_' + tool
        if tool in UNAVAILABLE:
            return '{} is not available here: it {}.'.format(name, UNAVAILABLE[tool])
        if i > 0 and (tool in MULTI_INPUT or tool in PDB_INPUT):
            return '{} can only be the first step.'.format(name)
        if i < len(steps) - 1 and (tool in SPLIT or tool in NOT_CIF):
            return '{} does not write mmCIF, so it can only be the last step.'.format(name)
    return None


def stem(path):
    """An input file name without its extensions (1abc.cif.gz -> 1abc)."""
    name = os.path.basename(path)
    if name.lower().endswith('.gz'):
        name = name[:-3]
    for ext in INPUT_EXTENSIONS:
        if name.lower().endswith(ext):
            return name[:-len(ext)]
    return name


def run_step(tool, option, label, inputs, out_path, cwd):
    """Run mmcif_<tool> on the input files, as on the command line, in the
    directory `cwd` with stdout written to `out_path`.

    Returns a dict with the arguments, the exit status, the run time in
    seconds and stderr.
    """
    main = _module(tool).main
    args = ['mmcif_' + tool]
    if option:
        args.append('-' + option)
    if label:
        args.append('--label')
    args.extend(inputs)

    saved = sys.argv, sys.stdout, sys.stderr, os.getcwd()
    err = io.StringIO()
    status = 0
    out = open(out_path, 'w', encoding='utf-8', errors='surrogateescape', newline='')
    try:
        sys.argv, sys.stdout, sys.stderr = list(args), out, err
        os.chdir(cwd)  # the split tools write their files there
        start = time.perf_counter()
        try:
            main()
        except SystemExit as exc:
            if isinstance(exc.code, str):
                err.write(exc.code + '\n')
                status = 1
            else:
                status = exc.code or 0
        except MemoryError:
            err.write('ERROR!! Out of memory: this file is too large for the browser '
                      '(install mmcif-tools and run it on your computer instead).\n')
            status = 1
        except Exception:  # a bug: show the traceback, as Python would
            err.write(traceback.format_exc())
            status = 1
        seconds = time.perf_counter() - start
    finally:
        sys.argv, sys.stdout, sys.stderr = saved[:3]
        os.chdir(saved[3])
        out.close()

    text = err.getvalue()
    for path in inputs:  # messages name the files as the user knows them
        text = text.replace(path, os.path.basename(path))
    if len(text) > MAX_STDERR:
        text = text[:MAX_STDERR] + '\n[... {:,} more characters]\n'.format(len(text) - MAX_STDERR)
    return {
        'tool': tool,
        'args': [os.path.basename(a) if a in inputs else a for a in args],
        'status': status,
        'seconds': seconds,
        'stderr': text,
    }


def _zip(directory, path):
    """Zip the files of a directory; returns [(name, size)]."""
    names = sorted(os.listdir(directory))
    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as zf:
        for name in names:
            zf.write(os.path.join(directory, name), name)
    return [(name, os.path.getsize(os.path.join(directory, name))) for name in names]


def run_pipeline(spec, progress=None):
    """Run a pipeline. `spec` is a JSON string:

        {"inputs": [paths], "steps": [{"tool": ..., "option": ..., "label": ...}],
         "workdir": directory for the outputs (emptied first)}

    A first step that reads several files gets every input, others only the
    first one. `progress(i)` is called before step i. Returns a JSON string
    with the result of each step, and the output file when the last step
    ran (also when it failed: mmcif_validate exits with 1 on errors).
    """
    spec = json.loads(spec)
    steps, inputs, workdir = spec['steps'], spec['inputs'], spec['workdir']
    problem = check(steps, len(inputs))
    if problem:
        return json.dumps({'error': problem})

    shutil.rmtree(workdir, ignore_errors=True)
    base = stem(inputs[0])
    current = inputs if steps[0]['tool'] in MULTI_INPUT else inputs[:1]
    results, output, files, previous = [], None, None, None
    for i, step in enumerate(steps):
        if progress is not None:
            progress(i)
        tool = step['tool']
        stepdir = os.path.join(workdir, str(i))
        cwd = os.path.join(stepdir, 'files')
        os.makedirs(cwd)
        out_path = os.path.join(stepdir, base + NOT_CIF.get(tool, '.cif'))
        result = run_step(tool, step.get('option') or '', step.get('label'), current, out_path, cwd)
        result['size'] = os.path.getsize(out_path)
        if previous is not None:
            shutil.rmtree(previous)  # frees the memory of the intermediate file
        previous = stepdir
        current = [out_path]
        if tool in SPLIT and result['status'] == 0:
            out_path = os.path.join(stepdir, '{}_{}.zip'.format(base, tool))
            files = _zip(cwd, out_path)
            result['size'] = os.path.getsize(out_path)
        results.append(result)
        if i == len(steps) - 1:
            output = out_path
        elif result['status'] != 0:
            break
    name = '{}_{}{}'.format(base, steps[-1]['tool'], os.path.splitext(output)[1]) if output else None
    return json.dumps({'steps': results, 'output': output, 'name': name, 'files': files,
                       'seconds': sum(r['seconds'] for r in results)})


def _main(argv):
    repeat = 1
    if len(argv) > 2 and argv[0] == '--repeat':
        repeat, argv = int(argv[1]), argv[2:]
    if len(argv) < 2:
        sys.stderr.write(__doc__)
        sys.exit(1)
    steps = []
    for part in argv[0].split('|'):
        tool, _, option = part.strip().partition(':')
        steps.append({'tool': tool[6:] if tool.startswith('mmcif_') else tool, 'option': option})
    workdir = os.path.join(os.environ.get('TMPDIR', '/tmp'), 'mmcif-web-runner')
    spec = json.dumps({'inputs': [os.path.abspath(p) for p in argv[1:]], 'steps': steps,
                       'workdir': workdir})
    best = None
    for _ in range(repeat):
        result = json.loads(run_pipeline(spec))
        if 'error' in result:
            sys.exit('ERROR!! ' + result['error'])
        if best is None or result['seconds'] < best['seconds']:
            best = result
    for step in best['steps']:
        sys.stdout.write('{:<40} exit {}  {:8.3f} s  {:>12,} bytes\n'.format(
            ' '.join(step['args']), step['status'], step['seconds'], step['size']))
        sys.stdout.write(step['stderr'])
    sys.stdout.write('{:<40}         {:8.3f} s  -> {}\n'.format('total', best['seconds'], best['output']))
    sys.exit(1 if any(step['status'] for step in best['steps']) else 0)


if __name__ == '__main__':
    _main(sys.argv[1:])
