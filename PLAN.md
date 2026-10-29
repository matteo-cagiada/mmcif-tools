# mmcif-tools: plan and design record

Started 2026-10-23 as "cif-tools", a port of pdb-tools; renamed mmcif-tools (commands `mmcif_*`, package
`mmciftools`) and made an independent project in phase 4 (§7, decisions 13-14). Sections 1-11 are kept as
written at the time, with the names updated.

## 1. Goal

Rebuild [pdb-tools](https://www.bonvinlab.org/pdb-tools/) (v2.6.1, 48 scripts, ~9k lines) as **mmcif-tools**: the
same small, single-purpose, pipeable command-line tools, but reading and writing mmCIF instead of PDB.

```
mmcif_selchain -A,B 1abc.cif | mmcif_delhetatm | mmcif_reres -1 | mmcif_tidy > clean.cif
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

| Topic | PDB (pdb-tools) | mmCIF (mmcif-tools) |
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

### 3.1 Shared library `mmciftools/_cif.py` (stdlib only)

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

~~Optional speed path: if `gemmi` is installed, a few heavy tools (`mmcif_tidy`, `mmcif_validate`) may use it.~~
Dropped in phase 3 (§7, decision 10): every tool is stdlib only.

### 3.2 Packaging

- `pyproject.toml`, package `mmciftools`, console entry points `mmcif_*` (pdb-tools uses the same mechanism).
- Python ≥ 3.9, no runtime deps and no gemmi extra: gemmi is not part of the release (§7, decision 12).
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
- Round-trip test: a tool that does nothing (`mmcif_tidy --noop`) must reproduce the input byte-for-byte.
- Validity: every output must parse with gemmi, and optionally with the PDBx dictionary check.

## 4. Tool mapping (all 48 pdb-tools)

**Legend:** ✅ direct port (same semantics, different I/O) · 🔁 port with changed semantics · ❌ drop / replace ·
🆕 new CIF-specific tool.

### Editing columns
| pdb-tools | mmcif-tools | Notes |
|---|---|---|
| pdb_b | ✅ mmcif_b | `B_iso_or_equiv` |
| pdb_occ | ✅ mmcif_occ | `occupancy` |
| pdb_chain | ✅ mmcif_chain | Sets the chain for all atoms. Chain ids are no longer limited to 1 char |
| pdb_rplchain | ✅ mmcif_rplchain | `-A:B` syntax; multi-char ids allowed |
| pdb_rplresname | ✅ mmcif_rplresname | `auth_comp_id` + `label_comp_id` (comp id is the same in both) |
| pdb_element | ✅ mmcif_element | `type_symbol`; reuse pdb-tools' atom-name → element heuristics |
| pdb_reatom | ✅ mmcif_reatom | `id`; must also renumber references in `_atom_site_anisotrop` / `_struct_conn` |
| pdb_reres | 🔁 mmcif_reres | Renumbers `auth_seq_id` by default. `label_seq_id` is canonical, so renumbering it needs `--label` and a warning |
| pdb_shiftres | 🔁 mmcif_shiftres | Same caveat as reres |
| pdb_uniqname | ✅ mmcif_uniqname | HETATM atom names |
| pdb_delinsertion | ✅ mmcif_delinsertion | `pdbx_PDB_ins_code` → `?`, then shift numbering |
| pdb_fixinsert | ✅ mmcif_fixinsert | Same |
| pdb_seg | ❌ | No segment id in mmCIF; `label_asym_id` names each molecule |
| pdb_chainxseg / pdb_segxchain | ❌ (was `mmcif_swapids`) | Replaced by 🆕 `mmcif_label2auth` / `mmcif_auth2label`, which copy one scheme onto the other (§7, decision 7) |
| pdb_chainbows | ❌ | No TER records in mmCIF. `mmcif_frompdb` splits molecules at TER and names blank chains; `mmcif_label2auth -chain` gives each molecule its own chain (§7, decision 6) |

### Selecting / deleting
| pdb-tools | mmcif-tools | Notes |
|---|---|---|
| pdb_selchain / pdb_delchain | ✅ mmcif_selchain / mmcif_delchain | Comma-separated, multi-char ids |
| pdb_selres / pdb_delres | ✅ mmcif_selres / mmcif_delres | Range syntax `1:10,20:`; uses the chosen scheme |
| pdb_selresname / pdb_delresname | ✅ | `comp_id` |
| pdb_selatom | ✅ mmcif_selatom | `auth_atom_id` / `label_atom_id` |
| pdb_selelem / pdb_delelem | ✅ | `type_symbol` |
| pdb_selhetatm / pdb_delhetatm | ✅ | `group_PDB` column |
| pdb_selaltloc | ✅ mmcif_selaltloc | `label_alt_id`; same "highest occupancy" default |
| pdb_selmodel | ✅ mmcif_selmodel | `pdbx_PDB_model_num` |
| pdb_selseg | ❌ | No segment id (see `mmcif_selchain --label`) |
| pdb_keepcoord | 🔁 mmcif_keepcoord | Keep only `_atom_site` (+ optional `_atom_site_anisotrop`, `_cell`, `_symmetry`), drop all other categories |
| pdb_head | ✅ mmcif_head | First N atom rows, with metadata kept |

### Splitting / merging
| pdb-tools | mmcif-tools | Notes |
|---|---|---|
| pdb_splitchain / pdb_splitmodel | ✅ | Each output gets its own `data_` block name (`<id>_A`) |
| pdb_splitseg | ❌ | |
| pdb_merge | 🔁 mmcif_merge | Harder than in PDB: must reconcile duplicate atom `id`s, chain clashes, and `_entity`/`_struct_asym` rows. v1 merges `_atom_site` plus minimal metadata, renumbers atoms, and warns on chain clashes |
| pdb_mkensemble | ✅ mmcif_mkensemble | Set `pdbx_PDB_model_num` per input |
| pdb_chkensemble | ✅ mmcif_chkensemble | |
| pdb_intersect | ✅ mmcif_intersect | Atom identity key = (model, chain, resnum, icode, resname, atom, altloc) |

### Inspecting / validating / formatting
| pdb-tools | mmcif-tools | Notes |
|---|---|---|
| pdb_wc | ✅ mmcif_wc | Plus a count of entities / categories |
| pdb_gap | ✅ mmcif_gap | Can also use `label_seq_id` jumps, which is more reliable than distances |
| pdb_tofasta | 🔁 mmcif_tofasta | Option to read the true sequence from `_entity_poly(_seq)` (includes unmodelled residues) instead of atoms. This is a real improvement over PDB |
| pdb_sort | ✅ mmcif_sort | Sort `_atom_site` rows by model/chain/resnum; re-sort the anisotrop loop to match |
| pdb_tidy | 🔁 mmcif_tidy | Normalise: column alignment, quoting, `#` category separators, `?`/`.` usage, sequential ids |
| pdb_validate | 🔁 mmcif_validate | Syntax check (tokenizer), required `_atom_site` columns, value types, and optionally the PDBx/mmCIF dictionary via gemmi |

### Fetching / conversion
| pdb-tools | mmcif-tools | Notes |
|---|---|---|
| pdb_fetch | ✅ mmcif_fetch | RCSB / PDBe `.cif(.gz)`, plus `-biounit` (assembly CIFs), plus AlphaFold DB (`AF-<uniprot>`) |
| pdb_tocif | ✅ mmcif_frompdb | Rewrite pdb-tools' "rudimentary" converter to output both auth and label ids |
| pdb_fromcif | ✅ mmcif_topdb | Fails with a clear message when the structure doesn't fit PDB (>99,999 atoms, multi-char chains) unless `--force` |

### New CIF-specific tools (🆕, phase 3)
- `mmcif_cat -<category>`: print one category as TSV (the CIF version of `grep ATOM`).
- `mmcif_delcat` / `mmcif_keepcat`: drop/keep whole categories.
- `mmcif_blocks`: list or select `data_` blocks.
- `mmcif_label2auth` / `mmcif_auth2label`: copy one identifier scheme onto the other.
- `mmcif_prune`: re-synchronise `_entity`/`_struct_asym`/`_pdbx_poly_seq_scheme` after heavy editing. **Postponed
  until after the first release** (§7, decision 9).

**Totals:** 48 originals → 42 ported (41 tools plus the `mmcif_delinsertion` alias), 6 dropped (the three segment
tools, chainbows and the two chain/segment swaps), plus 6 new: 48 tools.

## 5. Phases

| Phase | Scope | Done when |
|---|---|---|
| **0. Foundations** | Repo, `pyproject`, CI, fixture set, `_cif.py` tokenizer/reader/writer, round-trip test | `cat x.cif \| python -m mmciftools._cif` reproduces input byte-for-byte on all fixtures |
| **1. Core 20** | The tools we actually use daily: selchain, delchain, selres, delres, selresname, delresname, selatom, selelem, selhetatm, delhetatm, selaltloc, selmodel, chain, rplchain, b, occ, reres, shiftres, wc, tofasta, fetch | Each passes unit tests + cross-check vs pdb-tools |
| **2. The rest** | Remaining ports: split/merge/ensemble, sort, tidy, validate, gap, element, reatom, uniqname, insertion tools, intersect, keepcoord, head, conversions | All ported tools in §4 implemented |
| **3. CIF extras** | 🆕 tools (without `mmcif_prune`), the `mmcif_frompdb` TER fix, leftovers from phases 1–2 | Done (§11) |
| **4. Test release 0.1.0** | Rename to mmcif-tools; gemmi out of the package and tests (separate script); tool reference, CHANGELOG, CITATION; packaging metadata; CI up to Python 3.15; entity, assembly and info tools | Done (§12) |
| **5. Public release** | Feedback from testers; user guide (GitHub Pages); PyPI (TestPyPI first, trusted publishing on a tag), then conda-forge | `pip install mmcif-tools` works |
| **6. After the first release** | `mmcif_prune` (§7, decision 9) | |

Before phase 1, one useful step: **grep our own lab scripts/pipelines for `pdb_` calls** to rank the tools by how
much we use them, and adjust the "Core 20" list.

## 6. Risks

- **Quoting edge cases** in the writer (atom names with `'`, `"`, spaces; values starting with `_`, `#`, `$`,
  `[`, `;`, or reserved words like `data_`, `loop_`). Mitigation: a fuzz test that writes random values and
  re-parses them.
- **Silent metadata inconsistency** after edits (e.g. deleting chain B but `_struct_asym` still lists it).
  Mitigation: document it, `--prune`, `mmcif_validate` warns.
- **auth vs label confusion** for users. Mitigation: a global default, a `-h` note in every tool, and an
  `MMCIFTOOLS_SCHEME` env var.
- **Performance** of pure-Python tokenizing on huge assemblies (ribosomes, viral capsids). Mitigation: only fully
  tokenize the categories a tool needs, and benchmark early on a ~1M-atom file.
- **Non-wwPDB CIFs** (output of AlphaFold, Rosetta, ChimeraX, PyMOL) with missing or extra columns. Mitigation:
  have fixtures from each source; tools degrade gracefully when an optional column is absent.

## 7. Decisions

Decided (2026-10-23):

1. **Default identifier scheme: `auth_*`.** Switch per call with `--label`, or globally with `MMCIFTOOLS_SCHEME=label`.
2. **Core is strictly stdlib.** gemmi is used only by the tests (cross-checks); revised by decision 12: not as
   an optional extra, and not in the released test suite.
3. **Edits change one identifier scheme, never both.** One author chain usually spans several label chains
   (protein, ligands, waters), and `label_seq_id` is tied to the deposited sequence, so editing both at once would
   merge separate instances or break the sequence link. Copying one scheme onto the other is left to the planned
   `mmcif_label2auth` / `mmcif_auth2label` tools.
4. **Exception for names (phase 2):** `mmcif_rplresname` and `mmcif_uniqname` change *both* name columns (residue
   names: `auth_comp_id` and `label_comp_id`; atom names: `auth_atom_id` and `label_atom_id`). Names are the same
   chemical identity in both schemes, and programs differ in which one they read, so changing only one would be
   invisible in some of them.
5. **`mmcif_chainbows` and `mmcif_swapids` move to phase 3** (TER-based chain splitting has no mmCIF meaning, and
   swapping ids overlaps with `mmcif_label2auth` / `mmcif_auth2label`).

Decided for phase 3 (2026-10-24):

6. **No `mmcif_chainbows`.** TER records only separate molecules in PDB files, and mmCIF names every molecule. So
   `mmcif_frompdb` now starts a new molecule at each TER record (also in blank or repeated chains) and gives atoms
   of a blank chain their label chain as author chain. For mmCIF files, `mmcif_label2auth -chain` gives each
   molecule its own author chain.
7. **No `mmcif_swapids`.** A swap would put author values into the label columns. `mmcif_auth2label` is kept, for
   chains only: it copies the author chain onto the label chain of each *polymer*, leaves ligands and waters on
   their own label chains, and refuses to merge molecules. Its `--help` says why residue numbers are not copied.
8. **`mmcif_label2auth` copies everything by default** (chains and residue numbers); `-chain` or `-res` copies only
   one. It refuses when two residues would get the same identifiers.
9. **`mmcif_prune` (and `--prune`) after the first release.** Removing the rows of deleted molecules is well
   defined; keeping renames in step is not settled yet. A note for later: the edit tools change one scheme only
   (`mmcif_rplchain` changes `auth_asym_id` and leaves `label_asym_id`), so label ids could serve as the key for
   re-deriving author ids in the other tables; only residue and atom *names* are changed in both schemes
   (decision 4), and names are not keys.
10. **No gemmi path.** The core stays stdlib only. Checking a file against the PDBx/mmCIF dictionary is left to
    external validators (`mmcif_validate --help` mentions `gemmi validate`).
11. **Leftovers fixed in phase 3:** TER records in `mmcif_topdb`, quoted `'.'` / `'?'` read as text, and
    `mmcif_selaltloc` updating `_struct_conn`.

Decided 2026-10-24:

12. **gemmi is not in the final version.** No gemmi in the package, its help texts, the tests, the CI workflow
    or the `[test]` extra (done 2026-10-24, start of phase 4). The only gemmi check is
    `tools/check_with_gemmi.py`, a separate script run by hand: it runs the gemmi *program* (`gemmi validate`),
    which the user installs separately and can point to with `--gemmi PATH` (or `$GEMMI`).

Decided for phase 4 (2026-10-24):

13. **Name: mmcif-tools everywhere.** The PyPI name `mmcif-tools`, the import package `mmciftools`, commands
    `mmcif_*`, and `MMCIFTOOLS_SCHEME`. Two reasons the old names could not stay:
    - `ciftools` is already a PyPI package (a BinaryCIF library), with the same import name, and PyPI refuses
      new names that are too similar, so `cif-tools` would probably have been refused.
    - cod-tools (small-molecule crystallography, in Debian/Ubuntu) installs about 60 `cif_*` commands,
      including `cif_validate` and `cif_merge`.
14. **Independent project.** mmcif-tools is inspired by pdb-tools but not related to it, and is distributed
    on its own. pdb-tools is credited in README, NOTICE (as the Apache-2.0 licence asks, since the conventions
    and parts of the help texts follow pdb-tools) and CITATION.cff; nothing more. The help texts end with
    "This program is part of the mmcif-tools suite." (no longer "a port of pdb-tools"). The comparisons with
    pdb_X tools in the help texts stay, for users coming from pdb-tools.
15. **No search of existing pipelines.** The goal is that the tools work on any mmCIF file. Phase 4 adds an
    AlphaFold DB (ModelCIF) model and an entry with real symmetry operations (1STP) to the real-entry set.
16. **Docs:** the tool reference `docs/tools.md` is generated from the `--help` texts (`tools/make_docs.py`;
    a test and a CI step check that it is up to date). CHANGELOG.md and CITATION.cff are added. A user guide
    on GitHub Pages comes later. Done 2026-10-29: a Zensical site (`zensical.toml`; Zensical succeeds
    Material for MkDocs), built and published by `.github/workflows/docs.yml` (`zensical build --strict`,
    so broken links and anchors fail). The README stays the one place to edit the guide: `tools/make_docs.py`
    writes the other pages of `docs/` from its sections and from CHANGELOG.md (PAGES lists which section goes
    on which page), and `--check` covers them too. Then one website for everything (2026-10-29): a main page
    (`web/index.html`), the documentation under `/docs/` and mmcif-tools in the browser under `/app/`
    (`web/app/`: Pyodide, with the latest release from PyPI, which the workflow checks and pins; it also runs
    after each release), on its own domain, `mmcif-tools.dev` (GitHub forwards the github.io addresses). The
    old page addresses forward to `/docs/`, since the README of 0.1.1 on PyPI links to them. A separate repository for the in-browser page was tried first and dropped: two repositories
    and two sites were harder to follow than one.
17. **Version 0.1.0**, a test release: the first users test it from the repository; PyPI and conda-forge come
    after that. PyPI from 0.1.1 (2026-10-29, no tool changes): `.github/workflows/release.yml` builds the sdist
    and the wheel, checks them (`twine check --strict`, the smoke test on the installed wheel) and uploads them
    with Trusted Publishing, so no API token exists anywhere: a manual run goes to TestPyPI, a published
    GitHub release (tag `v<version>`, checked against `__version__`) to PyPI. The README's links are full URLs
    so that they work on PyPI too; `tools/make_docs.py` turns the ones to the website back into page links.
18. **Python 3.9 stays while it costs nothing** (the code is stdlib only and macOS's own python3 is 3.9); CI
    also covers the newest versions, up to 3.15.
19. **Four more mmCIF tools:** `mmcif_selentity` / `mmcif_delentity`, `mmcif_assembly`, `mmcif_info`.
20. **pdb-tools stays a test-only cross-check** (the `[test]` extra); it is not shipped.
21. **The test suite is not public** (2026-10-27, before the repository goes public for the pre-release).
    `tests/` stays in the developers' working copy (git-ignored) and in the history up to commit 36e75f3;
    the `[test]` extra is gone. CI runs `tools/smoke_test.py` (every command, a few pipelines) and the
    tool-reference check on Python 3.9-3.15, Linux, macOS and Windows.
22. **`mmcif_delmeta`** (2026-10-27): a command of its own for "drop the metadata that edits make stale",
    instead of an option of `mmcif_keepcoord` (which stays the counterpart of `pdb_keepcoord`). It keeps
    `_atom_site`, `_atom_site_anisotrop`, `_struct_conn`, `_entry`, `_cell`, `_symmetry`, `_exptl`,
    `_refine`, `_reflns`, `_em_3d_reconstruction`, `_pdbx_audit_revision_history`, `_entity`,
    `_entity_poly`, `_entity_poly_seq`, `_struct_asym`, `_chem_comp` and `_atom_type`; `-all` keeps the
    atoms, anisotropic B-factors and links only. 53 tools.
    - **AlphaFold/ColabFold templates** (same day, after a user report): ColabFold requires `_chem_comp.id`/`type`,
      `_struct_asym.id`/`entity_id` and `_entity_poly_seq.mon_id`; AlphaFold's parser also reads the method,
      resolution and release date, and fails (KeyError) on a `_struct_asym` chain of a protein entity with no
      atoms in model 1. So the kept set includes these categories, and `_struct_asym`, `_entity`,
      `_entity_poly` and `_entity_poly_seq` rows of chains and entities without atoms are removed;
      `_entity.pdbx_number_of_molecules` and `_entity_poly.pdbx_strand_id` are recomputed (unless unknown).
      `_chem_comp` is not pruned: AlphaFold needs the type of every residue of the sequence, resolved or not.
      It is therefore the *last* step of a pipeline. Checked by re-implementing those requirements (AlphaFold's
      own code was not run) on the stress-test entries.
    - **Other predictors and MD programs** (same day; their code read as text, not run): Boltz, Chai-1 and
      Protenix need nothing more. OpenFold3 requires `_pdbx_poly_seq_scheme` (asym_id -> entity_id), so it is
      rebuilt as AlphaFold 3 rebuilds a missing one: one row per `_entity_poly_seq` row and chain, with the author
      number of the first-model residue at that label_seq_id (the label one without author columns), or `.`
      (and `?` author columns) without atoms. A file without sequences or label ids keeps its table as it is.
      On the 300 stress entries, unedited, the rebuilt table equals the PDB's for residues with atoms, except
      where the PDB's records author names or numbers that `_atom_site` does not have (7NRF, 108D, 7Y7A).
      BioEmu2 needs `_pdbx_entity_nonpoly` for ligands (pruned by entity). OpenMM >= 8.3 and MDTraj >= 1.11 read
      `_chem_comp_bond` for ligand bonds (192 of 200 RCSB entries have it, with `_chem_comp_atom`); both are kept,
      minus the components whose heavy-atom names in `_atom_site` are not all in their definition (after
      `mmcif_uniqname`; hydrogens are left out, as NMR entries name terminal ones outside the CCD). PDBFixer and
      BioEmu2 read `_pdbx_struct_mod_residue`: rows are found by label ids and get the current author ids, or
      are matched on the shared ids, and removed when their residue is gone.

## 8. Phase 0 status (2026-10-23): done

**What exists:**
- `mmciftools/_cif.py`: the shared streaming layer, stdlib only.
  - **Events:** `DataBlock`, `Verbatim`, `Loop`, `Row`, `LoopEnd`.
  - **Functions:** `parse()`, `write()`/`Writer`, `quote()`, `get_scheme()`/`column()` (auth by default), and
    `open_input()` (handles gzip and CRLF).
  - **CLI:** `python -m mmciftools._cif` round-trips a file.
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
| Select / delete | `mmcif_selchain`, `mmcif_delchain`, `mmcif_selres`, `mmcif_delres`, `mmcif_selresname`, `mmcif_delresname`, `mmcif_selatom`, `mmcif_selelem`, `mmcif_selhetatm`, `mmcif_delhetatm`, `mmcif_selaltloc`, `mmcif_selmodel` |
| Edit | `mmcif_chain`, `mmcif_rplchain`, `mmcif_b`, `mmcif_occ`, `mmcif_reres`, `mmcif_shiftres` |
| Inspect | `mmcif_wc`, `mmcif_tofasta` |
| Fetch | `mmcif_fetch` |

All 21 are registered as console scripts in `pyproject.toml`. Each tool has pdb-tools' layout: the docstring is
its `--help` text, then `check_input`, `run` and `main`.

**Shared modules:**
- `mmciftools/_cli.py`: pdb-tools-style argument handling. It reads `[-option] [file]` or stdin, takes
  `--auth`/`--label` and `--help`, writes `ERROR!! ...` plus the usage with exit code 1, and exits quietly on a
  closed pipe.
- `mmciftools/_atoms.py`:
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
     made `mmcif_delhetatm` on the 80S ribosome (4V6X, 237k atoms) go from 2.5 to 1.3 s.
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
  - In the first pass, the keep-test was only run on link partners, which broke `mmcif_delres -::N`, a tool that
    counts residues (dangling links on 4V6X and 6VXX).
  - `mmcif_reres` mismatched the anisotropic rows of residues modelled as two residue types (1EJG PRO/SER and
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

The cross-check also showed that `mmcif_selmodel` should return a single-model file whole, as `pdb_selmodel`
does. That is fixed.

**Deliberate deviations from pdb-tools** (each is noted in the tool's `--help`):
- **Chain identifiers.** They are case-sensitive and may be longer than one character (`pdb_delchain` upper-cases
  them). mmCIF has no blank chain, so `mmcif_chain` needs an identifier.
- **`mmcif_selaltloc`.** It chooses one location per residue, keeps atoms without an altloc, and keeps the atom
  order.
  - `pdb_selaltloc` compares residue numbers as text and ignores insertion codes. On 1PPB it loses
    insertion-code residues; on 1EJG it drops backbone atoms with `-A`.
- **`mmcif_selhetatm`.** It keeps the metadata categories; `pdb_selhetatm` drops every non-HETATM record.
- **`mmcif_selelem`.** It matches elements ignoring case.
- **`mmcif_delres`.** It accepts a one-residue range (`-5:5`).
- **`mmcif_tofasta`.** It maps DNA (DA, DC, DG, DT) to bases, where `pdb_tofasta` gives X. It reads the first model
  only, never writes empty records, and adds `-seqres`.
- **`mmcif_wc`.** The `-a` / `-h` lists show atom names and HETATM residue names. (`pdb_wc -h` lists the first 3
  characters of atom names.)
- **`mmcif_fetch`.**
  - It exits 1 on a failed download, where `pdb_fetch` exits 0.
  - It accepts extended PDB ids and AlphaFold DB ids, and `-biounit` uses `<id>-assembly1.cif.gz`.
- **Residue numbers.** There are no limits (pdb-tools stops at 9999). `--label` renumbering warns that it breaks
  the link to `_entity_poly_seq`.
- **Help.** It is `--help` only, because `-h` is a real option (`mmcif_wc -h`, chain `h`).

**Known limitations:**
- These categories are not updated by any tool: `_struct_conf`, `_struct_sheet_range`, `_pdbx_poly_seq_scheme`,
  `_pdbx_struct_assembly_gen`, `_struct_site_gen`, `_pdbx_struct_conn_angle`, ...
- `mmcif_selaltloc` leaves `_struct_conn.pdbx_ptnr*_label_alt_id` as it was.
- Anisotropic rows placed *before* `_atom_site` are not pruned. No real file has been seen with that order.
- `mmcif_fetch` was checked live on one RCSB entry (byte-identical to the copy fetched earlier). The AlphaFold DB
  path is tested with a fake network only.
- With the real entries fetched and gemmi + pdb-tools installed, the full suite (93 tests) takes about 6.5 minutes.
  Most of that is the edit-every-value and gemmi checks on the two ribosomes and the pdb-tools cross-check.
  Without the real entries, it takes about 1 second.

## 10. Phase 2 status (2026-10-23): done

**The 21 remaining ports** (42 tools in total; `mmcif_chainbows` and `mmcif_swapids` moved to phase 3, decision 5):

| Group | Tools |
|---|---|
| Single file | `mmcif_delelem`, `mmcif_keepcoord`, `mmcif_head`, `mmcif_element`, `mmcif_rplresname`, `mmcif_fixinsert` (+ alias `mmcif_delinsertion`), `mmcif_reatom`, `mmcif_uniqname`, `mmcif_gap`, `mmcif_sort` |
| Models / several files | `mmcif_splitchain`, `mmcif_splitmodel`, `mmcif_mkensemble`, `mmcif_chkensemble`, `mmcif_merge`, `mmcif_intersect` |
| Format | `mmcif_tidy`, `mmcif_validate`, `mmcif_topdb`, `mmcif_frompdb` |

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
  - `mmcif_element -all` gives back the deposited `type_symbol` for **all 403,053 atoms** of the 10 entries (hemes,
    ions, glycans, lipids, nucleic acids);
  - `mmcif_topdb` → `mmcif_frompdb` gives back every author identifier on the 9 entries that fit PDB. It also gives
    back every label identifier on the 5 entries without missing residues;
  - 1D3Z split into its 10 models and rebuilt with `mmcif_mkensemble` passes `mmcif_chkensemble`;
  - `mmcif_validate` finds no errors in any entry;
  - `mmcif_tidy` keeps every value and is stable (running it twice changes nothing).
- **Cross-check against pdb-tools** (40 command lines × 14 structures = 630 comparisons): 606 identical, 15
  documented deviations, 9 not comparable. The 9 are `pdb_sort` refusing multi-model files, and the 9999-residue
  limit on 1JJ2.

**Deliberate deviations** (each in the tool's `--help`):
- **`mmcif_element`**:
  - It fills only missing elements by default (`-all` re-assigns), because deposited elements are authoritative.
  - mmCIF names have no column alignment, so it uses its own rules (ion = residue name; SE in MSE; two-letter
    elements in ligands; else the first letter).
- **`mmcif_fixinsert`.** The shift restarts at each chain instance (`label_asym_id`, the TER equivalent) and in each
  model. It accepts multi-character chains and negative numbers.
- **`mmcif_reatom`.** Ids continue across models (mmCIF ids are unique); ANISOU ids follow.
- **`mmcif_uniqname`.** It also updates the atom names in `_atom_site_anisotrop` and `_struct_conn`.
- **`mmcif_gap`.** `52 → 52A` is not a gap, and alternate-location Cα atoms are not compared.
- **`mmcif_sort`.** Multi-model files are sorted model by model, nothing is removed, and anisotropic rows follow
  their atoms.
- **Split tools.** Each output is a complete mmCIF file (metadata kept, block renamed). `mmcif_splitmodel` writes a
  single-model file too.
- **`mmcif_mkensemble` / `mmcif_merge`.**
  - The `_atom_site` loop has the union of the columns, and atom ids are renumbered to stay unique.
  - `mmcif_merge` merges into one data block and warns on chain clashes; `pdb_merge` concatenates.
- **`mmcif_chkensemble`.** The atom id is not part of an atom's identity, otherwise every ensemble would differ.
- **`mmcif_intersect`.** It keeps every model, and the first file's metadata.
- **`mmcif_head` / `mmcif_keepcoord`.** They keep a complete mmCIF file.
- **`mmcif_topdb`.** It writes the PDB format as specified: right-justified names and elements, MODEL columns
  11–14, charges as `1-`. It refuses structures that do not fit PDB, with no partial output.
- **`mmcif_frompdb`.** It assigns real label identifiers, plus `_entry`, `_cell` and `_symmetry` from CRYST1.
- **`mmcif_tidy` / `mmcif_validate`.** Both are mmCIF-specific: `pdb_tidy` and `pdb_validate` fix and check PDB
  columns, which mmCIF does not have.

**Known limitations:**
- **`mmcif_frompdb`:**
  - `label_seq_id` counts only the residues present in the file, because SEQRES is not used. It differs from the
    PDB entry after missing residues (1JJ2, 1KX5, 6VXX).
  - PDB entries order ligand instances and group glycans (branched entities) by their own conventions (2RH1,
    6VXX).
- **`mmcif_merge`.** It does not reconcile `_entity` / `_struct_asym`, and it drops the anisotropic and link rows of
  the files after the first.
- **`mmcif_topdb`.** It writes no TER, CONECT or header records.
- **`mmcif_validate`.** It does not check the file against the PDBx/mmCIF dictionary.
- **`mmcif_sort`.** It holds the atom rows in memory, as raw text.
- **Run time.** With the real entries fetched and gemmi + pdb-tools installed, the full suite (130 tests) takes
  about 7.5 minutes. Without the real entries, it takes a few seconds.

## 11. Phase 3 status (2026-10-24): done

**The 6 new tools** (48 in total):

| Group | Tools |
|---|---|
| Identifiers | `mmcif_label2auth` (`-chain`, `-res`), `mmcif_auth2label` |
| Categories / blocks | `mmcif_cat`, `mmcif_delcat`, `mmcif_keepcat`, `mmcif_blocks` |

- **`mmcif_label2auth`.** Copies `label_asym_id` → `auth_asym_id` and, for polymer residues, `label_seq_id` →
  `auth_seq_id` (insertion codes removed); ligands and waters keep their numbers. Links and anisotropic rows
  follow. It refuses when two residues would end up with the same chain and number.
- **`mmcif_auth2label`.** Copies author chains onto the label chains of polymers only. It refuses when an author
  chain holds two polymers, or when the new name is already the label chain of a ligand or water.
- **`mmcif_cat -_atom_site[.col,col]`.** One category as TSV: a header line, then one line per row; nulls as in the
  file; tabs, newlines and backslashes escaped. With several blocks, the header is repeated only when the
  columns change.
- **`mmcif_delcat` / `mmcif_keepcat`.** Comma-separated names, case-insensitive; a trailing `*` matches a prefix
  (`'-_pdbx_validate*'`). `mmcif_keepcat -_atom_site,_struct_conn` gives the same output as `mmcif_keepcoord`.
- **`mmcif_blocks`.** Lists the block names, or writes the named blocks; it stops reading once they are written
  (quick for one component of the whole CCD).

**Fixes:**
- **`mmcif_frompdb` and TER records** (found while planning). In a PDB file with blank chains whose numbering
  restarts after each TER (typical MD output), two molecules got the same label ids. Now each TER record ends a
  molecule, a blank chain gets its label chain as author chain (or the next free name), and later models without
  TER records are matched by chain and residue number.
- **`mmcif_topdb`** writes a TER record after the last residue of each polymer, with a serial number.
- **Quoted `'.'` and `'?'`.** The reader returns them as `_cif.Text`, a str that is not equal to the null, so
  `value in ('.', '?')` tests only match real nulls, and `quote()` writes them back in quotes (`mmcif_tidy`, rows
  formatted afresh).
- **`mmcif_selaltloc`** removes links (`_struct_conn`) to a removed conformer, and sets `pdbx_ptnr*_label_alt_id`
  to `.` on links to the kept one. None of the 10 real entries has such links, so this is tested on a fixture.
- **`parse()`** accepts category prefixes (`'_pdbx_*'`).

**Tested:**
- `tests/test_tools_phase3.py`: 31 unit and command-line tests; they also pass on Python 3.9.
- Real entries:
  - the new tools run on every entry with the phase-1 consistency checks;
  - `mmcif_label2auth` gives author ids equal to the label ids;
  - when `mmcif_auth2label` succeeds, it never merges two label chains;
  - `mmcif_cat` prints one line per atom, `mmcif_blocks` lists the entry, and `mmcif_keepcat` matches `mmcif_keepcoord`;
  - the PDB round trip still gives back the label ids, now through TER records.
- What the identifier tools do on the unedited entries:

  | Entry | `mmcif_auth2label` | `mmcif_label2auth -res` |
  |---|---|---|
  | 1D3Z, 1EHZ, 1EJG, 4HHB | no change | no change |
  | 1JJ2, 1KX5, 4V6X | changes the label chains | renumbers |
  | 1PPB | changes the label chains | refuses: `H 16` would become `H 1`, the number of the inhibitor 0G6 |
  | 2RH1 | no change | refuses: `A 267` would become `A 402`, the number of a sulfate |
  | 6VXX | no change | renumbers |

  Without `-res`, `mmcif_label2auth` (chains and numbers) works on all 10, because ligands then move to chains of
  their own.
- **Full suite:** 163 tests pass, in about 8.5 minutes, with the real entries, gemmi and pdb-tools installed. The
  pdb-tools cross-check is unchanged by the new TER records.

**Known limitations:**
- `mmcif_label2auth` and `mmcif_auth2label` leave the other categories that name chains or residues
  (`_struct_asym`, `_pdbx_poly_seq_scheme`, `_struct_conf`, ...) as they were; `mmcif_prune` is postponed.
- `mmcif_delcat` / `mmcif_keepcat` check nothing: removing `_atom_site` leaves `_atom_site_anisotrop` without atoms.
- `mmcif_cat` prints a quoted `'.'` and a null `.` the same way.

## 12. Phase 4 status (2026-10-24): test release 0.1.0

**Renamed to mmcif-tools** (§7, decisions 13-14): package `mmciftools`, 52 commands `mmcif_*`,
`MMCIFTOOLS_SCHEME`, NOTICE and help texts as an independent project inspired by pdb-tools. Files were moved
with `git mv`, so their history follows.

**gemmi out of the package** (§7, decision 12):
- gemmi is not imported or mentioned anywhere in the package, its help texts, the tests, the CI workflow or the
  `[test]` extra.
- The gemmi checks in the tests were replaced by standard-library ones:
  - every tool's output re-reads to the same text and has no syntax error for `mmcif_validate`;
  - a new test checks that the package imports only the Python standard library.
- `tools/check_with_gemmi.py` runs the gemmi program, which the user installs separately. It runs
  `gemmi validate` on each file, with options `--gemmi PATH`, `--dict`, `--parents` and `-` for stdin. It is not
  installed and not in the source distribution.

**Four new tools** (§7, decision 19):
- **`mmcif_selentity` / `mmcif_delentity`.** Select by entity number or type (polymer, non-polymer, branched,
  macrolide, water), with a fallback for files without `_entity`. On 1EHZ, `mmcif_delentity -non-polymer,water`
  keeps the whole tRNA (1652 atoms); `mmcif_delhetatm` cuts 323 atoms of modified nucleotides out of it.
- **`mmcif_assembly`.** Applies `_pdbx_struct_oper_list` to the chains of `_pdbx_struct_assembly_gen`,
  including products such as `(1-60)(61-88)` (applied right to left, as the dictionary defines). The naming
  follows the RCSB's assembly files: the identity copy keeps its names, the others become `A-2`, or `A-1x61` for
  a product.
  - Checked against the RCSB's own assembly files: 1STP (4 copies, 4004 atoms) and 2BUK (60 copies, 95,280
    atoms) have the same atoms, names and coordinates (largest difference 0.000 Å). Only the order differs: the
    RCSB lists ions copy by copy.
  - For 1M4X, `(1-60)(61-88)` gives 1680 operations. With our order, the 60 icosahedral copies made with each
    of the 28 other operations lie on one sphere (radius spread 3e-6 Å); the reverse order spreads them over
    229 Å.
  - gemmi 0.7.5 expands this expression to only 60 operations, so it could not serve as the reference.
- **`mmcif_info`.** Title, keywords, method, resolution, R-factors, dates, organisms, model quality (ModelCIF),
  entities and assemblies, as `Item:<tab>value` lines.

**Real-entry set: 12 entries.** 1STP (symmetry operations) and AF-P69905-F1 (an AlphaFold DB model, ModelCIF)
were added. `tests/fetch_real.py` fetches AlphaFold DB models through the AlphaFold DB API.

**Docs:**
- `docs/tools.md` is generated from the help texts by `tools/make_docs.py`. `tests/test_docs.py` and a CI step
  check that it is up to date and that the README lists every tool.
- README rewritten for users; CHANGELOG.md and CITATION.cff added (pdb-tools cited as a reference).

**Packaging:**
- **Metadata.** Version 0.1.0 (single source: `mmciftools.__version__`), SPDX licence `Apache-2.0` with
  LICENSE and NOTICE as licence files, authors, URLs, classifiers for Python 3.9 to 3.15. MANIFEST.in puts the
  fixtures, the docs, the changelog and the citation in the sdist.
- **Wheel.** Holds only `mmciftools`. Installed in clean virtual environments for Python 3.13 and macOS's 3.9
  (no network), all 52 commands answer `--help`, and a pipeline (`mmcif_assembly | mmcif_delentity -water |
  mmcif_wc`) works.
- **sdist.** Its own tests pass from the unpacked sdist: 177 tests, 18 skipped (no real entries, no `tools/`).
- **CI.** Python 3.9 (ubuntu-24.04), 3.10 to 3.15 (ubuntu-latest), and 3.15 on macOS and Windows. The package
  is installed with the test extras and tested again.

**Tests:** 177 tests pass in about 8 minutes, with the 12 real entries and the pdb-tools cross-check.

**Dictionary check with gemmi** (gemmi 0.7.5 program; a one-off run outside the repo, through
`tools/check_with_gemmi.py`):
- **What was checked:** every tool invocation of `tests/test_real.py`, plus `mmcif_assembly`,
  `mmcif_label2auth -res`, `mmcif_auth2label`, `mmcif_tidy`, `mmcif_keepcat`, `mmcif_merge` and the
  `mmcif_topdb | mmcif_frompdb` round trip, on the 12 entries: 42 outputs per entry.
- **Dictionaries:** the PDBx/mmCIF dictionary (v5.417) for PDB entries, and the ModelCIF dictionary
  (`mmcif_ma.dic`) for the AlphaFold model.
- **Result:** no output has a dictionary error that its input does not have (item types, allowed values,
  mandatory items, unique keys). The 11 PDB entries and all their outputs pass. The AlphaFold DB file has one
  violation of its own (`_pdbx_audit_revision_details.type` `Release`), which its outputs keep.
- **Refusals,** as expected: `mmcif_label2auth -res` on 1PPB and 2RH1 (a ligand already has one of the new
  numbers), and `mmcif_assembly` on the AlphaFold model (no assemblies).
- **Parent check (`--parents`):** compared with each input, it lists the categories that newly refer to
  something missing. These are the categories no tool keeps in step, and so the input list for `mmcif_prune`
  (phase 6):
  - **after removing atoms:** `_pdbx_validate_*`, `_struct_conf`, `_struct_sheet_range`,
    `_pdbx_struct_sheet_hbond`, `_struct_mon_prot_cis`, `_struct_site_gen`, `_pdbx_nonpoly_scheme`,
    `_pdbx_branch_scheme`, `_pdbx_nonpoly_atom_coordination*`, `_pdbx_modification_feature`,
    `_pdbx_struct_mod_residue`, `_pdbx_unobs_or_zero_occ_*`, `_pdbx_entity_instance_feature`,
    `_ndb_struct_na_base_pair(_step)`, and `_ma_qa_metric_local` (ModelCIF);
  - **after renaming:** the same categories, and `_atom_site` itself, whose ids no longer match
    `_pdbx_poly_seq_scheme` and `_struct_asym`;
  - **after `mmcif_keepcoord` / `mmcif_keepcat`:** `_atom_site` refers to the categories that were removed,
    as intended.

## 13. Stress test on 200 systems (2026-10-25)

The test lives outside the repository, in `../mmcif-tools-stress/`: the selection script, the runner, the downloaded
files (343 MB) and the reports (`reports/summary.md`, `systems.csv`, `all_runs.csv`, `tools.csv`, one CSV per
system).

**The systems.**
- **Selection:** 190 PDB entries drawn with a fixed seed from 44 RCSB Search API queries, plus 5 AlphaFold DB
  models and the PDBe "updated" files of 5 entries.
- **Methods:** X-ray (tiny, ultra-high resolution, typical, low resolution, large), cryo-EM (single particle,
  helical, subtomogram, 0.6-2.5M atoms), solution and solid-state NMR (up to 50+ models), neutron, electron
  crystallography, fiber, powder, solution scattering, EPR/theoretical models, multi-method entries.
- **Contents:** DNA, RNA, protein-DNA and protein-RNA complexes, DNA/RNA hybrids, PNA, D-amino acids, glycans,
  oligosaccharides only, peptides, membrane proteins, antibodies, FeS clusters, heme, selenomethionine, many
  ligands, UNK/UNL, BIRD molecules, more than 100 chains, icosahedral assemblies, chains of more than 3000
  residues, entries from before 1985 and from 2026, an integrative (IHM) model.
- **Size:** from 34 to 2,499,256 atoms.

**The runs.**
- **What ran:** every tool on every system, about 60 command lines each, through the installed wheel's
  commands.
- **Checks on each output:**
  - the exit status;
  - it re-reads to the same text;
  - no new unmatched anisotropic or link rows;
  - atom counts (tools that keep every atom, select/delete pairs, head, intersect, merge, assembly size, split +
    rebuild, PDB round trip);
  - no new gemmi dictionary error (PDBx, or ModelCIF for AlphaFold models; outputs up to 40 MB).
- `mmcif_fetch` downloaded every PDB entry and AlphaFold model.

**Found and fixed** (before the final run):
- **`mmcif_splitchain`:**
  - Chains whose ids differ only in case (`A`/`a`, `sA`/`sa`; 17 systems) overwrote each other's files on
    macOS's case-insensitive file system. File names are now unique ignoring case (`x_a~2.cif`), with a note.
  - It read the input once per chain (hours for 7Y7A's 1792 chains). It now writes up to 100 files per
    reading, with byte-identical output: checked on all fixtures and real entries, both schemes.
- **Missing columns (9A1M, an integrative model without auth_seq_id / pdbx_PDB_ins_code):**
  - `mmcif_fixinsert` / `mmcif_delinsertion` stopped with an error; such a file is now written unchanged.
  - `mmcif_label2auth` adds the missing author columns.
- **`mmcif_tofasta -seqres`:**
  - It needed `_pdbx_poly_seq_scheme`; it now falls back to `_entity_poly_seq`, with the same output on the
    12 real entries.
  - It wrote both residue types of a microheterogeneity (1EJG came out 48 residues instead of 46).
  - For a file with no polymer (3IRI, sugars only), it now says so.

**Final run** (fixed build):
- **Outcomes:** 12,424 runs: 12,152 ok, 70 refused, 202 skipped (one model: no `mmcif_mkensemble`; no PDB file:
  no `mmcif_frompdb`), 0 errors.
- **Checks:** all passed: atom counts 9,730, re-read 9,556, references 9,556, dictionary 8,812.
- **Refusals:** all documented ones:
  - `mmcif_topdb`: structures that do not fit the PDB format (26);
  - `mmcif_auth2label`: the new label chain is already a ligand's or water's (25);
  - `mmcif_label2auth -res`: two residues would get the same number (10);
  - `mmcif_assembly`: no assemblies (6);
  - `mmcif_tofasta -seqres`: no polymer (3).
- **Speed:** most tools take 7-11 µs per atom on large files, including Python start-up and gzip; `mmcif_tidy`
  takes 23. 7Y7A (2.5M atoms) takes 15-25 s per tool, and its 1792-chain split 83 s.
- **Memory:** the highest peaks are `mmcif_assembly` (2.0 GB on 7Y7A: it holds the asymmetric unit),
  `mmcif_validate` (1.4 GB) and `mmcif_frompdb` (0.7 GB). The other tools stream, at 15-300 MB.

**For later:**
- `mmcif_auth2label` refuses on 12% of entries. In PDB entries, author chain names often coincide with the
  label chains of ligands (e.g. `H` for an antibody heavy chain), so an option to rename those ligand label
  chains too may be worth adding.
- The memory of `mmcif_assembly` and `mmcif_validate` on the largest files could be reduced.

## 14. Second stress test: 100 systems aimed at the weak spots (2026-10-27)

Round 1 found most of its problems in chain identifiers, in files from outside the standard PDB pipeline, in
microheterogeneity, and in very large files. Round 2 aimed at those. It lives in
`../mmcif-tools-stress/round2/`: the selection and download scripts, `make_variants.py`, the data (148 MB) and
the reports.

**The systems:**
- **A, other sources (30):**
  - 10 integrative (PDB-IHM) entries;
  - 8 ModelArchive ModelCIF models (`ma-bak-cepc-*`; other ModelArchive projects could not be found by id);
  - 7 PDBe "updated" files, compared with the RCSB copy of the same entry;
  - 3 chemical components, 1 BIRD file, and 10 chemical components in one file (no `_atom_site`).
- **B, PDB entries for one weak spot each (40),** found from RCSB metadata:
  - 12 with chain ids differing only in case (they are common: 217 of 500 entries with 12-400 chains);
  - 3 with 500 or more chains (up to 716k atoms);
  - 5 with 100 or more models;
  - 5 with microheterogeneity;
  - 5 with residue numbers down to -397;
  - 5 antibodies with 15 or more insertion codes;
  - 5 with assemblies from several operation lists (up to 60) or a product, `(1-60)(61-88)` in 1M4X.
  - Residue numbers of 10,000 or more were not found in 600 entries; the modified copies cover them.
- **C, modified copies (30)** of 2KJK, 1X6Z, 4NHC, 3DIO, 9CVO and 9V1K (round-1 entries: NMR, altlocs,
  insertion codes, RNA, ligands, case-only chain ids):
  - **layout only** (Windows line endings, tabs, every value quoted, one value per line, comments in the loop,
    shuffled columns, upper-case tags, uncompressed without a final newline, rewritten by `mmcif_tidy`): every
    command must give the same atoms and values as on the original;
  - **two data blocks:** each half must equal the original;
  - **rewritten by the gemmi program:** the same atom counts;
  - **changed values or columns** (residue numbers +20000 / -1000, no author columns, no label columns, no
    model number, no element, no occupancy and B-factor, no `group_PDB`): the usual checks.

**The new checks.**
- *compare*: each command runs on the original file too, and the results are compared exactly (layout
  copies), doubled (two blocks), by counts (gemmi), or on the shared columns (PDBe vs RCSB).
- gemmi dictionary messages are compared without block names and quoted values.

**Found and fixed:**
- **4PTJ** (an ensemble whose models hold different residues, each with anisotropic rows): after `mmcif_reres`
  and `mmcif_fixinsert`, 1842 anisotropic rows of later models had model 1's residue numbers. Anisotropic rows
  now follow their atom, by atom id.
- **Upper-case tags:** `Row.moved()` matched column names exactly, so `mmcif_merge` left 15,690 duplicate ids
  (`_ATOM_SITE.ID`). Now case-insensitive.
- **Missing scheme columns:**
  - edits silently fell back to the other scheme (`--label` on a file without label columns renumbered the
    author ids, and the references then disagreed);
  - by default, an integrative model without author columns had its `label_seq_id` renumbered.
  - Edits now stop with a message that names the missing column and the way out (`mmcif_label2auth`).
- **gemmi-written links without residue names:** `mmcif_fixinsert` did not renumber them; its fallback key never
  matched.
- **`mmcif_head`** counted each data block separately, against its help.
- **Smaller ones:**
  - `mmcif_selaltloc` without an altloc column writes the file unchanged;
  - `mmcif_b` / `mmcif_occ` add their column;
  - missing-column errors name the mmCIF column;
  - the first-data-block behaviour of `mmcif_merge`, `mmcif_mkensemble` and `mmcif_topdb` is documented.

**Final round 2** (fixed build):
- **Outcomes:** 6,181 runs: 5,916 ok, 108 refused, 157 skipped, 0 errors.
- **Checks:** all passed: atom counts 4,733, re-read 4,701, references 4,701, dictionary 4,368, compare 1,427.
- **Refusals:** clear messages, e.g. `mmcif_topdb` on 5-character residue names (`A1A0K`, the new CCD
  identifiers), on multi-character chains and on more than 99,999 atoms; missing columns on the modified
  copies; files without atoms.
- **Unit tests:** 189 pass.

**Round 1 again, on the round-2 build** (a regression check):
- **Outcomes:** 12,424 runs, 0 errors, 0 failed checks.
- **The only change from round 1:** 3 more refusals, as intended. On the integrative entry 9A1M, which has no
  `auth_seq_id`, `mmcif_reres` and `mmcif_shiftres` (both forms) now refuse instead of renumbering
  `label_seq_id`.
