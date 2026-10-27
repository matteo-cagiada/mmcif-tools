# mmcif-tools

mmcif-tools is an independent project whose style and much of its tool set are inspired by
[pdb-tools](https://www.bonvinlab.org/pdb-tools/) (Rodrigues JPGLM, Teixeira JMC, Trellet M and Bonvin AMJJ,
*F1000Research* 2018, 7:1961). mmcif-tools ports and expands the pdb-tools functionality to the mmCIF (PDBx/mmCIF)
format.

As pdb-tools, it is a set of small, single-purpose command-line tools that read and write mmCIF files and chain
together with Unix pipes: select, delete, edit, convert and inspect structures without leaving the shell. The
package needs only the Python standard library (Python 3.9 or newer).

```bash
mmcif_fetch 1brs | mmcif_selchain -A,D | mmcif_delentity -non-polymer,water | mmcif_reres -1 > clean.cif
```

**CURRENT VERSION: Version 0.1.0, pre-release: 52 available tools**


## Tools

mmcif-tools reads the mmCIF text format, as written by the PDB (wwPDB, RCSB, PDBe), the AlphaFold DB (ModelCIF),
structure prediction and modelling programs and viewers, plain or gzip-compressed. BinaryCIF is currently not
supported.

The full `--help` text of every tool is in the [tool reference](docs/tools.md).

| Group | Tool | What it does |
|---|---|---|
| Fetch / convert | `mmcif_fetch` | Download an entry from the RCSB (or `-biounit`), or an AlphaFold DB model (`AF-P69905-F1`) |
| | `mmcif_topdb` / `mmcif_frompdb` | Convert to / from the PDB format (`mmcif_frompdb` assigns label identifiers) |
| Select / delete | `mmcif_selchain` / `mmcif_delchain` | Keep / remove chains (`-A,B`) |
| | `mmcif_selres` / `mmcif_delres` | Keep / remove residues by number (`-1:10,20:30`, `-::5`) |
| | `mmcif_selresname` / `mmcif_delresname` | Keep / remove residues by name (`-HOH`) |
| | `mmcif_selatom` | Keep atoms by name (`-CA,C,N,O`) |
| | `mmcif_selelem` / `mmcif_delelem` | Keep / remove atoms by element (`-C,N`, `-H`) |
| | `mmcif_selentity` / `mmcif_delentity` | Keep / remove entities by number or type (`-polymer`, `-non-polymer,water`) |
| | `mmcif_selhetatm` / `mmcif_delhetatm` | Keep / remove HETATM records (ligands, ions, waters, and modified residues) |
| | `mmcif_selaltloc` | Keep one alternate location (highest occupancy, or `-A`) |
| | `mmcif_selmodel` | Keep models (`-1,3`) |
| | `mmcif_head` | Keep the first N atoms (`-100`) |
| | `mmcif_keepcoord` | Keep only the coordinates (`_atom_site`, `_struct_conn`) |
| | `mmcif_intersect` | Keep the atoms of the first file found in all the others |
| Edit | `mmcif_chain` / `mmcif_rplchain` | Set all chains (`-C`) / rename one chain (`-A:B`) |
| | `mmcif_reres` / `mmcif_shiftres` | Renumber residues from a number / shift them by a number |
| | `mmcif_fixinsert` (`mmcif_delinsertion`) | Remove insertion codes, shifting the numbering |
| | `mmcif_reatom` | Renumber atom ids |
| | `mmcif_rplresname` | Rename a residue type (`-HIP:HIS`) |
| | `mmcif_uniqname` | Rename ligand atoms C1, C2, O1, … |
| | `mmcif_element` | Assign elements from atom names (`-all` to redo all) |
| | `mmcif_b` / `mmcif_occ` | Set B-factors / occupancies |
| Identifiers | `mmcif_label2auth` | Copy label chains and residue numbers onto the author ones (`-chain`, `-res`) |
| | `mmcif_auth2label` | Copy author chains onto the label chains of polymers |
| Models, files, assemblies | `mmcif_splitchain` / `mmcif_splitmodel` | Write one file per chain / model |
| | `mmcif_mkensemble` / `mmcif_chkensemble` | Build a multi-model file / check its models match |
| | `mmcif_merge` | Merge the atoms of several files into one |
| | `mmcif_assembly` | Build a biological assembly from the file's symmetry operations (`-1`) |
| Categories / blocks | `mmcif_cat` | Print one category as a table (TSV): `-_atom_site.id,Cartn_x` |
| | `mmcif_delcat` / `mmcif_keepcat` | Remove / keep whole categories (`-_struct_conf`, `'-_pdbx_validate*'`) |
| | `mmcif_blocks` | List the data blocks, or extract some (`-ATP`) |
| Inspect / format | `mmcif_info` | Title, method, resolution, organisms, entities, assemblies |
| | `mmcif_wc` | Counts: models, chains, residues, atoms, … |
| | `mmcif_tofasta` | Sequence in FASTA format (`-multi`, `-seqres`) |
| | `mmcif_gap` | Report chain breaks |
| | `mmcif_sort` | Sort atoms by chain and residue |
| | `mmcif_tidy` | Rewrite in a uniform, aligned layout |
| | `mmcif_validate` | Check syntax, coordinates and consistency |

## Usage examples and tips

```bash
# A clean protein for docking or simulations: one conformer, no ligands or waters, an aligned layout
mmcif_fetch 4hhb | mmcif_selaltloc | mmcif_delentity -non-polymer,water | mmcif_tidy > 4hhb_protein.cif

# Chain A with its heme, numbered from 1, as a PDB file for programs that need one
mmcif_fetch 4hhb | mmcif_selchain -A | mmcif_delentity -water | mmcif_reres -1 | mmcif_topdb > 4hhb_A.pdb

# The biological assembly, built from the symmetry operations in the file (chains A, A-2, A-3, A-4)
mmcif_fetch 1stp | mmcif_assembly > 1stp_tetramer.cif

# What an entry holds: method, resolution, organisms, entities, assemblies (and the pLDDT of AlphaFold models)
mmcif_fetch 4hhb | mmcif_info
mmcif_fetch AF-P69905-F1 | mmcif_info

# A table of atoms for awk, R or pandas (tab-separated, one header line)
mmcif_cat -_atom_site.auth_asym_id,auth_seq_id,label_comp_id,B_iso_or_equiv 4hhb.cif > atoms.tsv

# Residues numbered along the deposited sequence, and one author chain per molecule
mmcif_label2auth 4hhb.cif > 4hhb_by_sequence.cif

# A leaner file without the PDB's validation reports and secondary structure
mmcif_delcat '-_pdbx_validate*,_struct_conf,_struct_sheet_range' 4hhb.cif > 4hhb_lean.cif
```

- Each tool takes at most one option (starting with `-`) and one file, and reads stdin when no file is given, so
  tools chain with pipes. Negative numbers need a double dash (`mmcif_shiftres --5`).
- `--help` prints a tool's usage (not `-h`, which is a real option of some tools, e.g. `mmcif_wc -h`).
- Errors go to stderr as `ERROR!! ...`, with exit status 1.
- Rows a tool does not change are written back byte-for-byte. In changed rows, only the edited values are
  rewritten, and the other columns keep their alignment.
- To remove ligands and waters, select by entity (`mmcif_delentity -non-polymer,water`) rather than by HETATM
  record: `mmcif_delhetatm` also removes modified residues inside chains (selenomethionine, modified nucleotides).
- Quote options that contain `*`, so that the shell does not expand them: `'-_pdbx_validate*'`.
- Input files can be gzip-compressed (`.cif.gz`); outputs are plain text, so pipe them to `gzip` to compress.
- After editing, `mmcif_validate` checks the syntax, the coordinates and the links between tables.

## Coming from pdb-tools

Most tools have the pdb-tools name with `mmcif_` for `pdb_` and the same options. The differences:

- `pdb_tocif` / `pdb_fromcif` are `mmcif_frompdb` / `mmcif_topdb`.
- The segment tools (`pdb_seg`, `pdb_selseg`, `pdb_splitseg`) have no counterpart: mmCIF has no segment identifier,
  and `label_asym_id` names each molecule (`mmcif_selchain --label`, `mmcif_splitchain --label`).
- `pdb_chainbows`, `pdb_chainxseg` and `pdb_segxchain` have no counterpart: mmCIF has no TER records. `mmcif_frompdb`
  uses the TER records of a PDB file to separate molecules and names blank chains, and `mmcif_label2auth -chain`
  gives each molecule its own author chain.
- Where pdb-tools behaves differently (or has a bug), the tool's `--help` says so.

## Things to know about mmCIF

- **Two sets of identifiers.** mmCIF names chains and residues twice: author names (`auth_*`, what papers and
  viewers show) and the PDB's own names (`label_*`, one label chain per molecule, residues numbered along the
  sequence). Tools use author names by default; pass `--label` to use `label_*` instead, or set
  `MMCIFTOOLS_SCHEME=label` for every tool. Edits change one scheme only. Residue and atom names are the same in
  both, so `mmcif_rplresname` and `mmcif_uniqname` change both name columns.
- **Related tables stay in step.** Removing atoms also removes their rows in `_atom_site_anisotrop` and
  `_struct_conn`; renaming chains or renumbering residues updates the same two tables. Other categories that
  mention chains or residues (`_struct_conf`, `_pdbx_poly_seq_scheme`, ...) are left as they were.
- **Entities.** Every molecule belongs to an entity: polymer, non-polymer, branched (glycans) or water.
  `mmcif_selentity` / `mmcif_delentity` select by entity, which keeps modified residues (HETATM records inside a
  chain, such as MSE or the modified nucleotides of a tRNA) with their chain; `mmcif_delhetatm` removes them.
- **No fixed columns.** Chain identifiers can be longer than one character, and there are no limits on atom or
  residue numbers. `mmcif_topdb` refuses structures that do not fit the PDB format.

## Installation

```bash
git clone https://github.com/matteo-cagiada/mmcif-tools.git
cd mmcif-tools
python -m pip install .
```

There are no runtime dependencies. `python tools/smoke_test.py` then checks the installation: every command runs,
and a few pipelines give the expected results.

## Extra tools

`tools/` holds scripts that are not part of mmcif-tools and are not installed with it.

- **`tools/check_with_gemmi.py`** checks CIF files with [gemmi](https://github.com/project-gemmi/gemmi), an
  independent CIF reader: the syntax, and optionally a dictionary such as the PDBx/mmCIF one (`--dict`). It runs the
  gemmi program, which you install yourself (e.g. `pip install gemmi-program`); give its path with `--gemmi` if it
  is not on the PATH:

  ```bash
  python tools/check_with_gemmi.py --gemmi /path/to/gemmi out.cif
  mmcif_selchain -A 1abc.cif | python tools/check_with_gemmi.py --dict mmcif_pdbx_v50.dic -
  ```

- **`tools/make_docs.py`** writes the [tool reference](docs/tools.md) from the `--help` texts.
- **`tools/smoke_test.py`** checks an installation (see [Installation](#installation)).

## Citation and changes

How to cite mmcif-tools is in [CITATION.cff](CITATION.cff); the changes between versions are in
[CHANGELOG.md](CHANGELOG.md). [PLAN.md](PLAN.md) records the design decisions and the test results.

## Acknowledgements

Most of the code of mmcif-tools has been written with [Claude Code](https://claude.com/claude-code),
Anthropic's AI coding assistant, The author design the projects, tools and test to guarantee the quality of the tool.

## License

Apache-2.0. See [LICENSE](LICENSE).
