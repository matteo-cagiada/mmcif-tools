# cif-tools: porting pdb-tools to mmCIF

Draft plan, 2026-10-23.

## 1. Goal

Rebuild [pdb-tools](https://www.bonvinlab.org/pdb-tools/) (v2.6.1, 48 scripts, ~9k lines) as **cif-tools**: the
same small, single-purpose, pipeable command-line tools, but reading and writing mmCIF instead of PDB.

```
cif_selchain -A,B 1abc.cif | cif_delhetatm | cif_reres -1 | cif_tidy > clean.cif
```

What we keep from pdb-tools:

- **One script = one job**, named `cif_<verb>`, flags in the same style (`-A`, `-10`, `-1:50`).
- **Unix streaming**: read from a file or stdin, write to stdout, so the tools chain with pipes.
- **No dependencies** in the core: standard-library Python only, so it installs anywhere.
- **The same file layout**: `check_input(args)` → `run(fhandle, ...)` generator → `main()`. Contributors who know
  pdb-tools can read our code straight away.

## 2. The core difference: CIF isn't a fixed-column format

In pdb-tools every script works on fixed columns (`line[21]` is the chain, and so on). mmCIF has no fixed columns,
so that model doesn't carry over:

| Topic | PDB (pdb-tools) | mmCIF (cif-tools) |
|---|---|---|
| Atom record layout | Fixed columns | `loop_` of `_atom_site.*`; column order comes from the header, whitespace-separated tokens, with quoting (`'...'`, `"..."`, `;...;`) |
| Identifiers | One chain/resnum | Two sets: `auth_*` (author, what PDB shows) and `label_*` (canonical); `label_seq_id` is `.` for HETATM/waters |
| Chain id | 1 char | Any length (`label_asym_id`, `auth_asym_id`) |
| Models | `MODEL`/`ENDMDL` | `pdbx_PDB_model_num` column |
| Segment id | cols 73–76 | No equivalent (closest is `label_asym_id`) |
| TER | Record | Doesn't exist (chains are explicit) |
| Anisotropic B | `ANISOU` lines | `_atom_site_anisotrop` loop, keyed by atom id |
| Connectivity | `CONECT` | `_struct_conn` (and `_chem_comp_bond`) |
| Other metadata | Mostly ignorable | Cross-referenced categories (`_entity`, `_struct_asym`, `_entity_poly_seq`, `_pdbx_poly_seq_scheme`, ...) that can become inconsistent when we edit atoms |
| Size limits | 99,999 atoms, 9,999 residues | None |

So the main engineering work is **one small, shared mmCIF streaming layer**. After that, most tools are short.

## 3. Architecture

### 3.1 Shared library `ciftools/_cif.py` (stdlib only)

pdb-tools copies its helpers into every script. We share one module instead, because the CIF tokenizer is too
complex to duplicate 48 times.

1. **Tokenizer**: splits CIF lines into tokens and handles bare values, `'quoted'`/`"quoted"` values (including
   embedded quotes, as CIF 1.1 allows), `;` multi-line text fields, `#` comments, and `.`/`?` null values.
2. **Streaming block reader**: yields a sequence of *chunks*:
   - `Verbatim(lines)`: anything we don't touch, passed through **byte-for-byte**.
   - `Loop(category, columns, rows)`, only for the categories a tool asks for (by default `_atom_site`).
     Rows are yielded lazily, so memory stays flat even on large assemblies.
   - Also handles multiple `data_` blocks and single-row (key–value, non-loop) categories.
3. **Row object**: a light dict-like view (`row['auth_asym_id']`) over the token list, with typed helpers
   (`row.resnum`, `row.chain(scheme)`).
4. **Writer**: re-emits a loop with aligned columns and correct quoting (adds quotes only when a value needs them,
   e.g. atom names like `"O5'"`). Unmodified rows are written back as they were, so diffs stay small.
5. **Identifier scheme switch** (`--label` / `--auth`): one global convention for which identifiers a tool uses
   for selection and editing. **Default: `auth_*`**, because that matches what users see in PyMOL/ChimeraX and in
   papers, and it matches pdb-tools behaviour. Tools that *edit* an identifier edit only the chosen scheme and
   leave the other untouched; there is deliberately no `--both` (see §7, decision 3).
6. **Consistency helpers**: after deleting atoms, drop orphaned rows in dependent loops (`_atom_site_anisotrop`
   by atom id, `_struct_conn` by partner atoms). Optionally, prune `_struct_asym`/`_entity`/`_pdbx_poly_seq_scheme`.
   This is opt-in (`--prune`) in v1, so that by default we never silently rewrite metadata.

Optional speed path: if `gemmi` is installed, a few heavy tools (`cif_tidy`, `cif_validate`) may use it. The
stdlib path must always work.

### 3.2 Packaging

- `pyproject.toml`, package `ciftools`, console entry points `cif_*` (pdb-tools uses the same mechanism).
- Python ≥ 3.9, no runtime deps. Optional extra: `cif-tools[gemmi]`.
- Licence: pdb-tools is Apache-2.0, so we can derive from it if we keep the licence and NOTICE/attribution to the
  Bonvin lab authors. (Check if we want to coordinate with them upstream, see §7.)

### 3.3 Testing

- Fixtures: a small set of real CIFs covering the edge cases: multi-model NMR (e.g. 2K39-like), insertion codes,
  altlocs, HETATM + waters, multi-character chain ids (large assembly), nucleic acids with `'` in atom names,
  `_atom_site_anisotrop`, `_struct_conn`, an AlphaFold model, and a file whose `_atom_site` columns are in an
  unusual order.
- Per tool: unit tests mirroring pdb-tools' `tests/` layout (input fixture → expected output).
- **Cross-check harness**: for each ported tool, run `cif_X` on the CIF and `pdb_X` on the PDB version of the same
  entry, then compare atom sets/coordinates with gemmi (as a test-only dependency). This automatically catches
  semantic drift from the original tools.
- Round-trip test: a tool that does nothing (`cif_tidy --noop`) must reproduce the input byte-for-byte.
- Validity: every output must parse with gemmi, and optionally with the PDBx dictionary check.

## 4. Tool mapping (all 48 pdb-tools)

**Legend:** ✅ direct port (same semantics, different I/O) · 🔁 port with changed semantics · ❌ drop / replace ·
🆕 new CIF-specific tool.

### Editing columns
| pdb-tools | cif-tools | Notes |
|---|---|---|
| pdb_b | ✅ cif_b | `B_iso_or_equiv` |
| pdb_occ | ✅ cif_occ | `occupancy` |
| pdb_chain | ✅ cif_chain | Sets the chain for all atoms. Chain ids are no longer limited to 1 char |
| pdb_rplchain | ✅ cif_rplchain | `-A:B` syntax; multi-char ids allowed |
| pdb_rplresname | ✅ cif_rplresname | `auth_comp_id` + `label_comp_id` (comp id is the same in both) |
| pdb_element | ✅ cif_element | `type_symbol`; reuse pdb-tools' atom-name → element heuristics |
| pdb_reatom | ✅ cif_reatom | `id`; must also renumber references in `_atom_site_anisotrop` / `_struct_conn` |
| pdb_reres | 🔁 cif_reres | Renumbers `auth_seq_id` by default. `label_seq_id` is canonical, so renumbering it needs `--label` and a warning |
| pdb_shiftres | 🔁 cif_shiftres | Same caveat as reres |
| pdb_uniqname | ✅ cif_uniqname | HETATM atom names |
| pdb_delinsertion | ✅ cif_delinsertion | `pdbx_PDB_ins_code` → `?`, then shift numbering |
| pdb_fixinsert | ✅ cif_fixinsert | Same |
| pdb_seg | ❌ | No segment id in mmCIF. Suggest using `cif_chain --label` |
| pdb_chainxseg / pdb_segxchain | 🔁 cif_swapids | Replaces both: swap/copy `auth_asym_id` ↔ `label_asym_id` (the closest CIF analogue) |
| pdb_chainbows | 🔁 cif_chainbows | No TER records. Instead, split on chain breaks (C–N distance), the way `cif_gap` does |

### Selecting / deleting
| pdb-tools | cif-tools | Notes |
|---|---|---|
| pdb_selchain / pdb_delchain | ✅ cif_selchain / cif_delchain | Comma-separated, multi-char ids |
| pdb_selres / pdb_delres | ✅ cif_selres / cif_delres | Range syntax `1:10,20:`; uses the chosen scheme |
| pdb_selresname / pdb_delresname | ✅ | `comp_id` |
| pdb_selatom | ✅ cif_selatom | `auth_atom_id` / `label_atom_id` |
| pdb_selelem / pdb_delelem | ✅ | `type_symbol` |
| pdb_selhetatm / pdb_delhetatm | ✅ | `group_PDB` column |
| pdb_selaltloc | ✅ cif_selaltloc | `label_alt_id`; same "highest occupancy" default |
| pdb_selmodel | ✅ cif_selmodel | `pdbx_PDB_model_num` |
| pdb_selseg | ❌ | No segment id (see `cif_selchain --label`) |
| pdb_keepcoord | 🔁 cif_keepcoord | Keep only `_atom_site` (+ optional `_atom_site_anisotrop`, `_cell`, `_symmetry`), drop all other categories |
| pdb_head | ✅ cif_head | First N atom rows, with metadata kept |

### Splitting / merging
| pdb-tools | cif-tools | Notes |
|---|---|---|
| pdb_splitchain / pdb_splitmodel | ✅ | Each output gets its own `data_` block name (`<id>_A`) |
| pdb_splitseg | ❌ | |
| pdb_merge | 🔁 cif_merge | Harder than in PDB: must reconcile duplicate atom `id`s, chain clashes, and `_entity`/`_struct_asym` rows. v1 merges `_atom_site` plus minimal metadata, renumbers atoms, and warns on chain clashes |
| pdb_mkensemble | ✅ cif_mkensemble | Set `pdbx_PDB_model_num` per input |
| pdb_chkensemble | ✅ cif_chkensemble | |
| pdb_intersect | ✅ cif_intersect | Atom identity key = (model, chain, resnum, icode, resname, atom, altloc) |

### Inspecting / validating / formatting
| pdb-tools | cif-tools | Notes |
|---|---|---|
| pdb_wc | ✅ cif_wc | Plus a count of entities / categories |
| pdb_gap | ✅ cif_gap | Can also use `label_seq_id` jumps, which is more reliable than distances |
| pdb_tofasta | 🔁 cif_tofasta | Option to read the true sequence from `_entity_poly(_seq)` (includes unmodelled residues) instead of atoms. This is a real improvement over PDB |
| pdb_sort | ✅ cif_sort | Sort `_atom_site` rows by model/chain/resnum; re-sort the anisotrop loop to match |
| pdb_tidy | 🔁 cif_tidy | Normalise: column alignment, quoting, `#` category separators, `?`/`.` usage, sequential ids |
| pdb_validate | 🔁 cif_validate | Syntax check (tokenizer), required `_atom_site` columns, value types, and optionally the PDBx/mmCIF dictionary via gemmi |

### Fetching / conversion
| pdb-tools | cif-tools | Notes |
|---|---|---|
| pdb_fetch | ✅ cif_fetch | RCSB / PDBe `.cif(.gz)`, plus `-biounit` (assembly CIFs), plus AlphaFold DB (`AF-<uniprot>`) |
| pdb_tocif | ✅ cif_frompdb | Rewrite pdb-tools' "rudimentary" converter to output both auth and label ids |
| pdb_fromcif | ✅ cif_topdb | Fails with a clear message when the structure doesn't fit PDB (>99,999 atoms, multi-char chains) unless `--force` |

### New CIF-specific tools (🆕, phase 3, optional)
- `cif_cat <category>`: print one category as TSV (the CIF version of `grep ATOM`).
- `cif_delcat` / `cif_keepcat`: drop/keep whole categories.
- `cif_blocks`: list or select `data_` blocks.
- `cif_label2auth` / `cif_auth2label`: copy one identifier scheme onto the other.
- `cif_prune`: re-synchronise `_entity`/`_struct_asym`/`_pdbx_poly_seq_scheme` after heavy editing.

**Totals:** 48 originals → ~40 ported (incl. 4 merged/renamed), 4 dropped (segment-based), plus ~6 new.

## 5. Phases

| Phase | Scope | Done when |
|---|---|---|
| **0. Foundations** | Repo, `pyproject`, CI, fixture set, `_cif.py` tokenizer/reader/writer, round-trip test | `cat x.cif \| python -m ciftools._cif` reproduces input byte-for-byte on all fixtures |
| **1. Core 20** | The tools we actually use daily: selchain, delchain, selres, delres, selresname, delresname, selatom, selelem, selhetatm, delhetatm, selaltloc, selmodel, chain, rplchain, b, occ, reres, shiftres, wc, tofasta, fetch | Each passes unit tests + cross-check vs pdb-tools |
| **2. The rest** | Remaining ports: split/merge/ensemble, sort, tidy, validate, gap, element, reatom, uniqname, insertion tools, intersect, keepcoord, head, conversions | All ported tools in §4 implemented |
| **3. CIF extras** | 🆕 tools, `--prune`, optional gemmi path | |
| **4. Release** | Docs site (same style as the pdb-tools page: one entry per tool with `-h` text), PyPI, conda-forge, migration table pdb_X → cif_X | `pip install cif-tools` works |

Before phase 1, one useful step: **grep our own lab scripts/pipelines for `pdb_` calls** to rank the tools by how
much we use them, and adjust the "Core 20" list.

## 6. Risks

- **Quoting edge cases** in the writer (atom names with `'`, `"`, spaces; values starting with `_`, `#`, `$`,
  `[`, `;`, or reserved words like `data_`, `loop_`). Mitigation: a fuzz test that writes random values and
  re-parses them.
- **Silent metadata inconsistency** after edits (e.g. deleting chain B but `_struct_asym` still lists it).
  Mitigation: document it, `--prune`, `cif_validate` warns.
- **auth vs label confusion** for users. Mitigation: a global default, a `-h` note in every tool, and an
  `CIFTOOLS_SCHEME` env var.
- **Performance** of pure-Python tokenizing on huge assemblies (ribosomes, viral capsids). Mitigation: only fully
  tokenize the categories a tool needs, and benchmark early on a ~1M-atom file.
- **Non-wwPDB CIFs** (output of AlphaFold, Rosetta, ChimeraX, PyMOL) with missing or extra columns. Mitigation:
  have fixtures from each source; tools degrade gracefully when an optional column is absent.

## 7. Decisions

Decided (2026-10-23):

1. **Default identifier scheme: `auth_*`.** Switch per call with `--label`, or globally with `CIFTOOLS_SCHEME=label`.
2. **Core is strictly stdlib.** gemmi is used only by the tests (cross-checks) and, later, as an optional extra.
3. **Edits change one identifier scheme, never both.** One author chain usually spans several label chains
   (protein, ligands, waters), and `label_seq_id` is tied to the deposited sequence, so editing both at once would
   merge separate instances or break the sequence link. Copying one scheme onto the other is left to the planned
   `cif_label2auth` / `cif_auth2label` tools.
4. **Exception for names (phase 2):** `cif_rplresname` and `cif_uniqname` change *both* name columns (residue
   names: `auth_comp_id` and `label_comp_id`; atom names: `auth_atom_id` and `label_atom_id`). Names are the same
   chemical identity in both schemes, and programs differ in which one they read, so changing only one would be
   invisible in some of them.
5. **`cif_chainbows` and `cif_swapids` move to phase 3** (TER-based chain splitting has no mmCIF meaning, and
   swapping ids overlaps with `cif_label2auth` / `cif_auth2label`).

Still open:

3. **Naming**: `cif_*` (proposed) or `mmcif_*`?
4. **Upstream**: build an independent package, or propose it to the Bonvin lab as a sister project / `pdb-tools`
   v3 subpackage?
5. **Which of our existing pipelines** need to migrate first (to define the phase-1 list)?

## 8. Phase 0 status (2026-10-23): done

**What exists:**
- `ciftools/_cif.py`: the shared streaming layer, stdlib only.
  - **Events:** `DataBlock`, `Verbatim`, `Loop`, `Row`, `LoopEnd`.
  - **Functions:** `parse()`, `write()`/`Writer`, `quote()`, `get_scheme()`/`column()` (auth by default), and
    `open_input()` (handles gzip and CRLF).
  - **CLI:** `python -m ciftools._cif` round-trips a file.
- `tests/test_cif.py`: 47 tests, all passing; the 2 gemmi tests are skipped when gemmi is missing.
- `tests/data/`: 11 hand-written fixtures, one edge case each.
- `tests/bench_cif.py`: the benchmark.

**How the writer behaves.**
- Unmodified rows are written from their raw text.
- In modified rows only the changed values are spliced in. The spaces after a value absorb any change in length,
  and tabs, comments and line breaks are kept.
- Editing one item of a key–value category (`_cell.length_a 53.6`) changes only that line.
- A loop that a tool empties is dropped. A loop that was already empty in the input is kept.

**Tested:**
- Every fixture is byte-exact in both modes (only `_atom_site` parsed, or every category parsed), with CRLF line
  endings, without a final newline, gzipped, and through the CLI's stdin/stdout.
- Tokenizer edge cases.
- A 5k-string fuzz of `quote()`.
- Edits: alignment, quoting, new rows, rows mapped from another loop, deleting rows that share a line, text fields.
- A stress test that edits *every* value of every row of every fixture, then parses the result back. When gemmi is
  installed, it also checks that gemmi reads the same values as we do from every fixture and from the edited
  output.

**One-off check on real files** (local, not part of the suite; 1,827 files including two 450 MB CCD dumps):
- All files round-trip byte-exact except one with an empty `_atom_site` loop.
- That case is now fixed and covered by `empty_atom_site.cif`.

**Benchmark** (synthetic 1M atoms, 100 MB, MacBook):

| Operation | Time |
|---|---|
| Copy round-trip | 1.4 s |
| Edit 1 column in every row | 5.6 s |
| Edit 3 columns in every row | 7.1 s |
| Peak memory (100k or 1M atoms) | ~0.01 MB, so streaming works |

- For comparison, pdb-tools' `pdb_b` takes 0.36 s on a 1M-atom PDB, because it only slices fixed columns.
- Editing is about 4 µs per row above copying, mostly the `re.split` in `_splice_bare`. If phase-1 tools feel slow
  on large assemblies, that is the first thing to optimise.

**Known unsupported layouts.** These raise a `CIFError` with a line number; they are legal CIF, but none was seen
in the real-file check:
- tags on the same line as `loop_` in a parsed category;
- a tag after values on the same line;
- more content after a closing `;` inside a key–value category.

**Not done:** the CI workflow (`.github/workflows/tests.yml`) has not run yet, because there is no GitHub remote.

## 9. Phase 1 status (2026-10-23): done

**The 21 core tools:**

| Group | Tools |
|---|---|
| Select / delete | `cif_selchain`, `cif_delchain`, `cif_selres`, `cif_delres`, `cif_selresname`, `cif_delresname`, `cif_selatom`, `cif_selelem`, `cif_selhetatm`, `cif_delhetatm`, `cif_selaltloc`, `cif_selmodel` |
| Edit | `cif_chain`, `cif_rplchain`, `cif_b`, `cif_occ`, `cif_reres`, `cif_shiftres` |
| Inspect | `cif_wc`, `cif_tofasta` |
| Fetch | `cif_fetch` |

All 21 are registered as console scripts in `pyproject.toml`. Each tool has pdb-tools' layout: the docstring is
its `--help` text, then `check_input`, `run` and `main`.

**Shared modules:**
- `ciftools/_cli.py`: pdb-tools-style argument handling. It reads `[-option] [file]` or stdin, takes
  `--auth`/`--label` and `--help`, writes `ERROR!! ...` plus the usage with exit code 1, and exits quietly on a
  closed pipe.
- `ciftools/_atoms.py`:
  - `Fields`, the identifier columns for the chosen scheme;
  - `filter_atoms()`, for the select/delete tools;
  - `rename()`, for chain and residue edits;
  - `Source`, which lets a tool read its input twice. stdin is copied to a temporary file, kept in memory up to
    32 MB.

**Phase-0 problems addressed:**
1. **Speed.**
   - Edits splice a value in using `str.split(None, i)`, and `quote()` has a shortcut for plain values. On the
     synthetic 1M-atom benchmark, a one-column edit went from 5.6 to 4.1 s, and three columns from 7.1 to 6.1 s.
     A copy stays at about 1.5 s, and peak memory stays at about 0.01 MB.
   - Rows holding quoted values without spaces (e.g. `"O5'"` in every nucleic acid) now take the fast path. That
     made `cif_delhetatm` on the 80S ribosome (4V6X, 237k atoms) go from 2.5 to 1.3 s.
2. **Related tables.** Select/delete tools remove `_atom_site_anisotrop` rows (matched by atom id) and
   `_struct_conn` links whose partner atom was removed. A link counts as dangling only if its partner existed in
   the input. Chain and residue edits update the same two tables.
3. **The API trap.** A row read from a file holds its values as a tuple, so `row.values[i] = x` raises an error
   instead of being silently lost.
4. **Real files in the tests.** See below.

Still open:
- the CI workflow, which is waiting for a push;
- adding or removing columns (phase 2);
- the three unsupported layouts;
- a quoted `'.'` being read the same as a null.

**Real-entry tests** (`tests/test_real.py`; the set is listed in `tests/real_set.txt` and fetched with
`tests/fetch_real.py`):
- **The set:** 1EJG, 4HHB, 1D3Z, 1PPB, 1EHZ, 1KX5, 2RH1, 6VXX, 1JJ2 and 4V6X. Together they cover a tiny
  high-resolution structure with alternate locations and anisotropic B, a tetramer with hemes, an NMR ensemble,
  insertion codes, tRNA, a protein–DNA complex, a membrane protein, a cryo-EM glycoprotein, and two ribosomes. One
  ribosome is mmCIF-only.
- **The layer:** a byte-exact round trip; gemmi reading identical values in every category; editing every value
  of every row.
- **Every tool on every entry:**
  - the output re-reads identically, and gemmi reads it the same way (entries under 30k atoms);
  - select X + delete X give back every atom;
  - edit tools keep the atom count;
  - no tool creates an anisotropic row or link that does not match an atom.
- **Bugs this found and fixed:**
  - In the first pass, the keep-test was only run on link partners, which broke `cif_delres -::N`, a tool that
    counts residues (dangling links on 4V6X and 6VXX).
  - `cif_reres` mismatched the anisotropic rows of residues modelled as two residue types (1EJG PRO/SER and
    LEU/ILE).
  - Quoted atom names were slow.

**Cross-check against pdb-tools** (`tests/test_vs_pdbtools.py`): each structure is converted to PDB format, the
same command line goes to `pdb_X` and `cif_X`, and the atoms (or the text, for wc/tofasta) are compared.

Over 14 structures (5 fixtures + 9 real entries; 4V6X does not fit the PDB format), 29 command lines each:

| Outcome | Count | Cases |
|---|---|---|
| Identical | 395 of 406 | Atoms, order, occupancies, B and coordinates; `wc` and `tofasta` text |
| Documented deviations | 8 | `selaltloc` on the fixture, 1EJG and 1PPB; `tofasta` on 1KX5 (DNA) |
| Not comparable | 3 | `reres` / `shiftres` on 1JJ2, where pdb-tools stops at residue 9999 |

The cross-check also showed that `cif_selmodel` should return a single-model file whole, as `pdb_selmodel`
does. That is fixed.

**Deliberate deviations from pdb-tools** (each is noted in the tool's `--help`):
- **Chain identifiers.** They are case-sensitive and may be longer than one character (`pdb_delchain` upper-cases
  them). mmCIF has no blank chain, so `cif_chain` needs an identifier.
- **`cif_selaltloc`.** It chooses one location per residue, keeps atoms without an altloc, and keeps the atom
  order.
  - `pdb_selaltloc` compares residue numbers as text and ignores insertion codes. On 1PPB it loses
    insertion-code residues; on 1EJG it drops backbone atoms with `-A`.
- **`cif_selhetatm`.** It keeps the metadata categories; `pdb_selhetatm` drops every non-HETATM record.
- **`cif_selelem`.** It matches elements ignoring case.
- **`cif_delres`.** It accepts a one-residue range (`-5:5`).
- **`cif_tofasta`.** It maps DNA (DA, DC, DG, DT) to bases, where `pdb_tofasta` gives X. It reads the first model
  only, never writes empty records, and adds `-seqres`.
- **`cif_wc`.** The `-a` / `-h` lists show atom names and HETATM residue names. (`pdb_wc -h` lists the first 3
  characters of atom names.)
- **`cif_fetch`.**
  - It exits 1 on a failed download, where `pdb_fetch` exits 0.
  - It accepts extended PDB ids and AlphaFold DB ids, and `-biounit` uses `<id>-assembly1.cif.gz`.
- **Residue numbers.** There are no limits (pdb-tools stops at 9999). `--label` renumbering warns that it breaks
  the link to `_entity_poly_seq`.
- **Help.** It is `--help` only, because `-h` is a real option (`cif_wc -h`, chain `h`).

**Known limitations:**
- These categories are not updated by any tool: `_struct_conf`, `_struct_sheet_range`, `_pdbx_poly_seq_scheme`,
  `_pdbx_struct_assembly_gen`, `_struct_site_gen`, `_pdbx_struct_conn_angle`, ...
- `cif_selaltloc` leaves `_struct_conn.pdbx_ptnr*_label_alt_id` as it was.
- Anisotropic rows placed *before* `_atom_site` are not pruned. No real file has been seen with that order.
- `cif_fetch` was checked live on one RCSB entry (byte-identical to the copy fetched earlier). The AlphaFold DB
  path is tested with a fake network only.
- With the real entries fetched and gemmi + pdb-tools installed, the full suite (93 tests) takes about 6.5 minutes.
  Most of that is the edit-every-value and gemmi checks on the two ribosomes and the pdb-tools cross-check.
  Without the real entries, it takes about 1 second.

## 10. Phase 2 status (2026-10-23): done

**The 21 remaining ports** (42 tools in total; `cif_chainbows` and `cif_swapids` moved to phase 3, decision 5):

| Group | Tools |
|---|---|
| Single file | `cif_delelem`, `cif_keepcoord`, `cif_head`, `cif_element`, `cif_rplresname`, `cif_fixinsert` (+ alias `cif_delinsertion`), `cif_reatom`, `cif_uniqname`, `cif_gap`, `cif_sort` |
| Models / several files | `cif_splitchain`, `cif_splitmodel`, `cif_mkensemble`, `cif_chkensemble`, `cif_merge`, `cif_intersect` |
| Format | `cif_tidy`, `cif_validate`, `cif_topdb`, `cif_frompdb` |

**Library additions:**
- **Adding columns.** `Loop.new()` / `Loop.with_columns()` and `Row.moved()`: rows can be put into a loop with
  more columns. When the columns are only appended, each row keeps its original text and the new values are added
  at the end of the line.
- **Several inputs and file output.** `_cli.parse_files()` (several input files) and `_cli.write_file()` (the
  split tools).
- **In `_atoms`:**
  - `guess_element()`: element rules for mmCIF names, which have no PDB column alignment;
  - `replace_values()`;
  - `reference_problems()`: anisotropic rows and links that match no atom;
  - `base_name()`, `atom_columns()`, `rename_block()`.

**Tested:**
- `tests/test_tools_phase2.py`: 26 unit and command-line tests.
- **Real entries** (`tests/test_real.py`):
  - every new single-file tool runs on every entry, with the phase-1 consistency checks;
  - `cif_element -all` gives back the deposited `type_symbol` for **all 403,053 atoms** of the 10 entries (hemes,
    ions, glycans, lipids, nucleic acids);
  - `cif_topdb` → `cif_frompdb` gives back every author identifier on the 9 entries that fit PDB. It also gives
    back every label identifier on the 5 entries without missing residues;
  - 1D3Z split into its 10 models and rebuilt with `cif_mkensemble` passes `cif_chkensemble`;
  - `cif_validate` finds no errors in any entry;
  - `cif_tidy` keeps every value and is stable (running it twice changes nothing).
- **Cross-check against pdb-tools** (40 command lines × 14 structures = 630 comparisons): 606 identical, 15
  documented deviations, 9 not comparable. The 9 are `pdb_sort` refusing multi-model files, and the 9999-residue
  limit on 1JJ2.

**Deliberate deviations** (each in the tool's `--help`):
- **`cif_element`**:
  - It fills only missing elements by default (`-all` re-assigns), because deposited elements are authoritative.
  - mmCIF names have no column alignment, so it uses its own rules (ion = residue name; SE in MSE; two-letter
    elements in ligands; else the first letter).
- **`cif_fixinsert`.** The shift restarts at each chain instance (`label_asym_id`, the TER equivalent) and in each
  model. It accepts multi-character chains and negative numbers.
- **`cif_reatom`.** Ids continue across models (mmCIF ids are unique); ANISOU ids follow.
- **`cif_uniqname`.** It also updates the atom names in `_atom_site_anisotrop` and `_struct_conn`.
- **`cif_gap`.** `52 → 52A` is not a gap, and alternate-location Cα atoms are not compared.
- **`cif_sort`.** Multi-model files are sorted model by model, nothing is removed, and anisotropic rows follow
  their atoms.
- **Split tools.** Each output is a complete mmCIF file (metadata kept, block renamed). `cif_splitmodel` writes a
  single-model file too.
- **`cif_mkensemble` / `cif_merge`.**
  - The `_atom_site` loop has the union of the columns, and atom ids are renumbered to stay unique.
  - `cif_merge` merges into one data block and warns on chain clashes; `pdb_merge` concatenates.
- **`cif_chkensemble`.** The atom id is not part of an atom's identity, otherwise every ensemble would differ.
- **`cif_intersect`.** It keeps every model, and the first file's metadata.
- **`cif_head` / `cif_keepcoord`.** They keep a complete mmCIF file.
- **`cif_topdb`.** It writes the PDB format as specified: right-justified names and elements, MODEL columns
  11–14, charges as `1-`. It refuses structures that do not fit PDB, with no partial output.
- **`cif_frompdb`.** It assigns real label identifiers, plus `_entry`, `_cell` and `_symmetry` from CRYST1.
- **`cif_tidy` / `cif_validate`.** Both are mmCIF-specific: `pdb_tidy` and `pdb_validate` fix and check PDB
  columns, which mmCIF does not have.

**Known limitations:**
- **`cif_frompdb`:**
  - `label_seq_id` counts only the residues present in the file, because SEQRES is not used. It differs from the
    PDB entry after missing residues (1JJ2, 1KX5, 6VXX).
  - PDB entries order ligand instances and group glycans (branched entities) by their own conventions (2RH1,
    6VXX).
- **`cif_merge`.** It does not reconcile `_entity` / `_struct_asym`, and it drops the anisotropic and link rows of
  the files after the first.
- **`cif_topdb`.** It writes no TER, CONECT or header records.
- **`cif_validate`.** It does not check the file against the PDBx/mmCIF dictionary.
- **`cif_sort`.** It holds the atom rows in memory, as raw text.
- **Run time.** With the real entries fetched and gemmi + pdb-tools installed, the full suite (130 tests) takes
  about 7.5 minutes. Without the real entries, it takes a few seconds.
