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
Summarizes the metadata of an mmCIF file: title, method, resolution, ...

Prints one line per item, as 'Item:<tab>value', for the items the file
has (each data block in turn):

    Entry           _entry.id (else the data block name)
    Title           _struct.title
    Keywords        _struct_keywords.pdbx_keywords
    Method          _exptl.method (several for hybrid methods)
    Resolution      _refine.ls_d_res_high, else _em_3d_reconstruction.
                    resolution, else _reflns.d_resolution_high (in A)
    R-work, R-free  _refine.ls_R_factor_R_work / ls_R_factor_R_free
    Deposited       _pdbx_database_status.recvd_initial_deposition_date
    Released        first _pdbx_audit_revision_history.revision_date
    Revised         last revision date, when there are several
    Organisms       scientific names in _entity_src_gen, _entity_src_nat
                    and _pdbx_entity_src_syn
    Model quality   global scores of models (_ma_qa_metric_global, e.g.
                    the pLDDT of AlphaFold models)
    Entity          one line per entity: id, type, copies, chains,
                    residues and name
    Assembly        one line per biological assembly (see mmcif_assembly):
                    id, number of chains and description

Usage:
    mmcif_info <cif file>

Example:
    mmcif_info 1abc.cif
    mmcif_fetch 4hhb | mmcif_info

This program is part of the mmcif-tools suite.
"""

import sys

from mmciftools import _cif, _cli

NULLS = ('.', '?')
CATEGORIES = (
    '_entry', '_struct', '_struct_keywords', '_exptl', '_refine', '_em_3d_reconstruction',
    '_reflns', '_pdbx_database_status', '_pdbx_audit_revision_history', '_entity_src_gen',
    '_entity_src_nat', '_pdbx_entity_src_syn', '_ma_qa_metric', '_ma_qa_metric_global', '_entity',
    '_entity_poly', '_entity_poly_seq', '_struct_asym', '_pdbx_struct_assembly')


def check_input(args):
    """Return the file handle."""
    fhandle, _, _ = _cli.parse_args(args, __doc__, takes_option=False)
    return fhandle


def _clean(value):
    """A value on one line (None for nulls)."""
    if value is None or value in NULLS:
        return None
    return ' '.join(value.split())


def _first(tables, category, *columns):
    """The first non-null value of the columns, in that order."""
    for column in columns:
        for row in tables.get(category, ()):
            value = _clean(row.get(column))
            if value is not None:
                return value
    return None


def _values(tables, category, column):
    """The distinct non-null values of a column, in file order."""
    out = []
    for row in tables.get(category, ()):
        value = _clean(row.get(column))
        if value is not None and value not in out:
            out.append(value)
    return out


def _plural(word, names):
    return '{}{} {}'.format(word, '' if len(names) == 1 else 's', ', '.join(names))


def summary(tables, block):
    """The lines for one data block."""
    lines = []

    def add(item, value):
        if value:
            lines.append('{}:\t{}\n'.format(item, value))

    add('Entry', _first(tables, '_entry', 'id') or block)
    add('Title', _first(tables, '_struct', 'title'))
    add('Keywords', _first(tables, '_struct_keywords', 'pdbx_keywords'))
    add('Method', ', '.join(_values(tables, '_exptl', 'method')))
    resolution = (_first(tables, '_refine', 'ls_d_res_high')
                  or _first(tables, '_em_3d_reconstruction', 'resolution')
                  or _first(tables, '_reflns', 'd_resolution_high'))
    add('Resolution', resolution and resolution + ' A')
    add('R-work', _first(tables, '_refine', 'ls_R_factor_R_work'))
    add('R-free', _first(tables, '_refine', 'ls_R_factor_R_free'))
    add('Deposited', _first(tables, '_pdbx_database_status', 'recvd_initial_deposition_date'))
    dates = _values(tables, '_pdbx_audit_revision_history', 'revision_date')
    add('Released', dates[0] if dates else None)
    add('Revised', dates[-1] if len(dates) > 1 else None)
    organisms = []
    for category, column in (('_entity_src_gen', 'pdbx_gene_src_scientific_name'),
                             ('_entity_src_nat', 'pdbx_organism_scientific'),
                             ('_pdbx_entity_src_syn', 'organism_scientific')):
        organisms += [o for o in _values(tables, category, column) if o not in organisms]
    add('Organisms', ', '.join(organisms))

    metrics = {row.get('id'): row for row in tables.get('_ma_qa_metric', ())}
    for row in tables.get('_ma_qa_metric_global', ()):
        metric = metrics.get(row.get('metric_id'))
        name = _clean(metric.get('name')) if metric is not None else None
        add('Model quality', ' '.join(v for v in (name or 'score', _clean(row.get('metric_value')))
                                      if v))

    chains, residues = {}, {}
    for row in tables.get('_entity_poly', ()):
        strands = _clean(row.get('pdbx_strand_id'))
        if strands:
            chains[row.get('entity_id')] = strands.replace(',', ', ')
    for row in tables.get('_entity_poly_seq', ()):
        residues[row.get('entity_id')] = residues.get(row.get('entity_id'), 0) + 1
    asyms = {}
    for row in tables.get('_struct_asym', ()):
        asyms.setdefault(row.get('entity_id'), []).append(row.get('id'))
    for row in tables.get('_entity', ()):
        entity = row.get('id')
        copies = _clean(row.get('pdbx_number_of_molecules'))
        parts = [entity, _clean(row.get('type')) or '?',
                 '{} cop{}'.format(copies, 'y' if copies == '1' else 'ies') if copies else None,
                 _plural('chain', chains[entity].split(', ')) if entity in chains
                 else _plural('label chain', asyms[entity]) if entity in asyms else None,
                 '{} residues'.format(residues[entity]) if entity in residues else None,
                 _clean(row.get('pdbx_description'))]
        add('Entity', '\t'.join(p for p in parts if p))

    for row in tables.get('_pdbx_struct_assembly', ()):
        count = _clean(row.get('oligomeric_count'))
        parts = [row.get('id'), '{} chain{}'.format(count, '' if count == '1' else 's') if count else None,
                 _clean(row.get('oligomeric_details')), _clean(row.get('details'))]
        add('Assembly', '\t'.join(p for p in parts if p))
    return lines


def run(fhandle, scheme=None):
    """Yield the summary lines."""
    tables, block = None, None
    for ev in _cif.parse(fhandle, CATEGORIES):
        if type(ev) is _cif.Row:
            tables.setdefault(ev.loop.category.lower(), []).append(ev)
        elif type(ev) is _cif.DataBlock:
            if tables is not None:
                for line in summary(tables, block):
                    yield line
                yield '\n'
            tables, block = {}, ev.name
    if tables is None:
        raise _cif.CIFError('No data block in the input')
    for line in summary(tables, block):
        yield line


def main():
    fhandle = check_input(sys.argv[1:])
    _cli.output_text(_cli.lazy(run, fhandle), fhandle)


if __name__ == '__main__':
    main()
