# cif-tools

A swiss army knife for mmCIF files: small, single-purpose command-line tools that read and write mmCIF and chain
together with Unix pipes. It is a port of [pdb-tools](https://www.bonvinlab.org/pdb-tools/) to the format the PDB
has used since 2014.

```bash
cif_fetch 1brs | cif_selchain -A,D | cif_delhetatm | cif_reres -1 > clean.cif
```

**Status: phase 1 (the 21 core tools) done.** The roadmap, the pdb-tools → cif-tools mapping and the design
decisions are in [PLAN.md](PLAN.md).

## Tools

| Tool | What it does |
|---|---|
| `cif_fetch` | Download an entry from the RCSB (or `-biounit`), or an AlphaFold DB model (`AF-P69905-F1`) |
| `cif_selchain` / `cif_delchain` | Keep / remove chains (`-A,B`) |
| `cif_selres` / `cif_delres` | Keep / remove residues by number (`-1:10,20:30`, `-::5`) |
| `cif_selresname` / `cif_delresname` | Keep / remove residues by name (`-HOH`) |
| `cif_selatom` | Keep atoms by name (`-CA,C,N,O`) |
| `cif_selelem` | Keep atoms by element (`-C,N`) |
| `cif_selhetatm` / `cif_delhetatm` | Keep / remove HETATM records (ligands, ions, waters) |
| `cif_selaltloc` | Keep one alternate location (highest occupancy, or `-A`) |
| `cif_selmodel` | Keep models (`-1,3`) |
| `cif_chain` / `cif_rplchain` | Set all chains (`-C`) / rename one chain (`-A:B`) |
| `cif_reres` / `cif_shiftres` | Renumber residues from a number / shift them by a number |
| `cif_b` / `cif_occ` | Set B-factors / occupancies |
| `cif_wc` | Summary: models, chains, residues, atoms, ... |
| `cif_tofasta` | Sequence in FASTA format (`-multi`, `-seqres`) |

Each tool prints its usage with `--help`. The options follow pdb-tools: one option (starting with `-`) and one
file. When no file is given, the tool reads stdin. Negative numbers need a double dash (`cif_shiftres --5`). Input
files may be gzip-compressed.

## How it differs from pdb-tools

- **Identifiers.** mmCIF has two sets of names: author names (`auth_*`, what papers and viewers show) and the
  PDB's own names (`label_*`). Tools use author names by default. To use `label_*` instead, pass `--label`, or set
  `CIFTOOLS_SCHEME=label` for every tool.
- **No fixed columns.** Chain identifiers can be longer than one character, and there are no limits on atom or
  residue numbers. Rows a tool does not change are written back byte-for-byte. In changed rows, only the edited
  values are rewritten, and the other columns keep their alignment.
- **Related tables stay in step.** Removing atoms also removes their rows in `_atom_site_anisotrop` and
  `_struct_conn`. Renaming chains or renumbering residues updates the same two tables. Other categories that
  mention chains or residues are left untouched: `_struct_conf`, `_pdbx_poly_seq_scheme`, ...
- **No segment IDs.** mmCIF has no segment identifier, so the segment tools (`pdb_seg`, `pdb_selseg`, ...) have no
  equivalent.

Where pdb-tools has a bug, the tool says so in its `--help` and does not copy it. The deviations are listed in
PLAN.md §9.

## Install (development)

```bash
python -m pip install -e .
```

There are no runtime dependencies. `pip install -e ".[test]"` adds gemmi and pdb-tools, which are used only in the
tests, as independent cross-checks.

## Tests

```bash
python -m unittest discover -s tests -v
```

- **Default run:** the hand-written fixtures in `tests/data/`. It passes with the standard library alone; the
  gemmi and pdb-tools cross-checks are skipped when those packages are missing.
- **Real PDB entries:** 10 entries (about 10 MB, listed in `tests/real_set.txt`) are included in the run once they
  have been fetched:

  ```bash
  python tests/fetch_real.py
  ```

- **Benchmark:** a synthetic 1M-atom structure, generated in memory:

  ```bash
  python tests/bench_cif.py
  ```

## License

Apache-2.0, as pdb-tools. See [LICENSE](LICENSE) and [NOTICE](NOTICE) for credits to the pdb-tools authors.
