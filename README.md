# cif-tools

A swiss army knife for mmCIF files: small, single-purpose command-line tools that read and write mmCIF and chain
together with Unix pipes. It is a port of [pdb-tools](https://www.bonvinlab.org/pdb-tools/) to the format the PDB
has used since 2014.

```bash
cif_selchain -A,B 1abc.cif | cif_delhetatm | cif_reres -1 | cif_tidy > clean.cif   # the goal (phase 1+)
```

**Status: phase 0 (foundations).** The shared mmCIF layer and its tests are done. The `cif_*` tools themselves
come next. The full roadmap, the pdb-tools → cif-tools mapping, and the design decisions are in
[PLAN.md](PLAN.md).

## Design in one paragraph

[`ciftools/_cif.py`](ciftools/_cif.py) is a streaming, **standard-library-only** mmCIF reader and writer.
- It parses only the categories a tool asks for (by default `_atom_site`) and passes everything else through
  untouched.
- Rows a tool doesn't modify are written back **byte-for-byte**.
- In modified rows, only the changed values are rewritten. The spaces after a value absorb any change in length,
  so columns stay aligned and diffs stay minimal.
- Tools select residues and chains by **author identifiers (`auth_*`)** by default. You can switch to `label_*`
  per call or with `CIFTOOLS_SCHEME=label`.

## Install (development)

```bash
python -m pip install -e .
```

There are no runtime dependencies. `pip install -e ".[test]"` also installs gemmi, which is used *only* as an
optional, independent cross-check in the tests.

## Tests

```bash
python -m unittest discover -s tests -v
```

The suite passes without gemmi; the gemmi cross-checks are then skipped. The phase-0 criterion is that every
fixture round-trips exactly:

```bash
for f in tests/data/*.cif; do python -m ciftools._cif < "$f" | cmp - "$f" && echo "ok $f"; done
```

Benchmark on a synthetic 1M-atom structure, generated in memory:

```bash
python tests/bench_cif.py
```

## License

Apache-2.0, as pdb-tools. See [LICENSE](LICENSE) and [NOTICE](NOTICE) for credits to the pdb-tools authors.
