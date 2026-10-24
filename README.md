# cif-tools

A swiss army knife for mmCIF files: small, single-purpose command-line tools that read and write mmCIF and chain
together with Unix pipes. It is a port of [pdb-tools](https://www.bonvinlab.org/pdb-tools/) to the format the PDB
has used since 2014.

```bash
cif_fetch 1brs | cif_selchain -A,D | cif_delhetatm | cif_reres -1 > clean.cif
```

**Status: phases 1–3 done: 48 tools.** They cover every pdb-tools tool that has an mmCIF meaning, plus tools that
only make sense in mmCIF: categories, data blocks, and copying one identifier scheme onto the other. The
roadmap, the pdb-tools → cif-tools mapping and the design decisions are in [PLAN.md](PLAN.md).

## Tools

| Group | Tool | What it does |
|---|---|---|
| Fetch / convert | `cif_fetch` | Download an entry from the RCSB (or `-biounit`), or an AlphaFold DB model (`AF-P69905-F1`) |
| | `cif_topdb` / `cif_frompdb` | Convert to / from the PDB format (`cif_frompdb` assigns label identifiers) |
| Select / delete | `cif_selchain` / `cif_delchain` | Keep / remove chains (`-A,B`) |
| | `cif_selres` / `cif_delres` | Keep / remove residues by number (`-1:10,20:30`, `-::5`) |
| | `cif_selresname` / `cif_delresname` | Keep / remove residues by name (`-HOH`) |
| | `cif_selatom` | Keep atoms by name (`-CA,C,N,O`) |
| | `cif_selelem` / `cif_delelem` | Keep / remove atoms by element (`-C,N`, `-H`) |
| | `cif_selhetatm` / `cif_delhetatm` | Keep / remove HETATM records (ligands, ions, waters) |
| | `cif_selaltloc` | Keep one alternate location (highest occupancy, or `-A`) |
| | `cif_selmodel` | Keep models (`-1,3`) |
| | `cif_head` | Keep the first N atoms (`-100`) |
| | `cif_keepcoord` | Keep only the coordinates (`_atom_site`, `_struct_conn`) |
| | `cif_intersect` | Keep the atoms of the first file found in all the others |
| Edit | `cif_chain` / `cif_rplchain` | Set all chains (`-C`) / rename one chain (`-A:B`) |
| | `cif_reres` / `cif_shiftres` | Renumber residues from a number / shift them by a number |
| | `cif_fixinsert` (`cif_delinsertion`) | Remove insertion codes, shifting the numbering |
| | `cif_reatom` | Renumber atom ids |
| | `cif_rplresname` | Rename a residue type (`-HIP:HIS`) |
| | `cif_uniqname` | Rename ligand atoms C1, C2, O1, … |
| | `cif_element` | Assign elements from atom names (`-all` to redo all) |
| | `cif_b` / `cif_occ` | Set B-factors / occupancies |
| Models / files | `cif_splitchain` / `cif_splitmodel` | Write one file per chain / model |
| | `cif_mkensemble` / `cif_chkensemble` | Build a multi-model file / check its models match |
| | `cif_merge` | Merge the atoms of several files into one |
| Identifiers | `cif_label2auth` | Copy label chains and residue numbers onto the author ones (`-chain`, `-res`) |
| | `cif_auth2label` | Copy author chains onto the label chains of polymers |
| Categories / blocks | `cif_cat` | Print one category as a table (TSV): `-_atom_site.id,Cartn_x` |
| | `cif_delcat` / `cif_keepcat` | Remove / keep whole categories (`-_struct_conf`, `'-_pdbx_validate*'`) |
| | `cif_blocks` | List the data blocks, or extract some (`-ATP`) |
| Inspect / format | `cif_wc` | Summary: models, chains, residues, atoms, … |
| | `cif_tofasta` | Sequence in FASTA format (`-multi`, `-seqres`) |
| | `cif_gap` | Report chain breaks |
| | `cif_sort` | Sort atoms by chain and residue |
| | `cif_tidy` | Rewrite in a uniform, aligned layout |
| | `cif_validate` | Check syntax, coordinates and consistency |

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
  mention chains or residues are left untouched: `_struct_conf`, `_pdbx_poly_seq_scheme`, ... (a tool to bring
  them back in step, `cif_prune`, is planned after the first release).
- **No segment IDs, no TER records.** mmCIF has no segment identifier, so the segment tools (`pdb_seg`,
  `pdb_selseg`, `pdb_splitseg`) have no equivalent: `label_asym_id` names each molecule (`cif_selchain --label`,
  `cif_splitchain --label`). Chains are explicit, so there is no `pdb_chainbows`: `cif_frompdb` uses the TER records
  of a PDB file to separate molecules and names blank chains, and `cif_label2auth -chain` gives each molecule its
  own author chain (also the job of `pdb_chainxseg` / `pdb_segxchain`).
- **Names in both schemes.** Residue and atom names are the same chemical identity in both schemes, so
  `cif_rplresname` and `cif_uniqname` change both name columns.

Where pdb-tools has a bug, the tool says so in its `--help` and does not copy it. The deviations are listed in
PLAN.md §9 and §10.

## Install (development)

```bash
python -m pip install -e .
```

There are no runtime dependencies. During development, `pip install -e ".[test]"` adds gemmi and pdb-tools, which
are used only in the tests, as independent cross-checks. The gemmi checks will be removed from the tests before
the first release (PLAN.md §7, decision 12).

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
