# mmcif-tools

### [mmcif-tools.dev](https://mmcif-tools.dev/)

Run the tools [in your browser](https://mmcif-tools.dev/app/), with nothing to install, or read
the [documentation](https://mmcif-tools.dev/docs/).

mmcif-tools is an independent project whose style and much of its tool set are inspired by
[pdb-tools](https://www.bonvinlab.org/pdb-tools/) (Rodrigues JPGLM, Teixeira JMC, Trellet M and Bonvin AMJJ,
*F1000Research* 2018, 7:1961). mmcif-tools ports and expands the pdb-tools functionality to the mmCIF (PDBx/mmCIF)
format.

As pdb-tools, it is a set of small, single-purpose command-line tools that read and write mmCIF files and chain
together with Unix pipes: select, delete, edit, convert and inspect structures without leaving the shell. The
package needs only the Python standard library (Python 3.9 or newer).

```bash
mmcif_fetch 1brs | mmcif_selchain -A,D | mmcif_delentity -non-polymer,water | mmcif_reres -1 | mmcif_delmeta > clean.cif
```

**CURRENT VERSION: Version 0.1.1, pre-release: 53 available tools**


## Tools

mmcif-tools reads the mmCIF text format, as written by the PDB (wwPDB, RCSB, PDBe), the AlphaFold DB (ModelCIF),
structure prediction and modelling programs and viewers, plain or gzip-compressed. BinaryCIF is currently not
supported.

The full `--help` text of every tool is in the [tool reference](https://mmcif-tools.dev/docs/tools/).

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
| | `mmcif_delmeta` | Delete the metadata that edits make stale; keep the rest (entry, method, entities, chains, ...) in step with the atoms (`-all`: atoms only) |
| | `mmcif_keepcoord` | Keep only the coordinates (`_atom_site`, `_struct_conn`), as `pdb_keepcoord` |
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

### Command examples

```bash
# A clean protein for docking or simulations: one conformer, no ligands or waters, an aligned layout
mmcif_fetch 4hhb | mmcif_selaltloc | mmcif_delentity -non-polymer,water | mmcif_delmeta | mmcif_tidy > 4hhb_protein.cif

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

# Delete only some categories: here the PDB's validation reports and secondary structure
mmcif_delcat '-_pdbx_validate*,_struct_conf,_struct_sheet_range' 4hhb.cif > 4hhb_lean.cif
```

### `mmcif_delmeta`

A PDB entry is mostly metadata: 4HHB has 72 categories, and only a few of them hold the atoms. Most of the rest
refers to residues by chain and number: validation reports, secondary structure, sequence schemes, binding sites,
assemblies, and more. The edit tools keep the atoms and their links consistent, but not these tables, so after
renumbering or deleting residues they describe a structure that is no longer there. Simulation and prediction
programs rarely read them anyway.

`mmcif_delmeta` deletes that metadata and keeps what is fundamental: the atoms, their links and anisotropic B-factors;
the entry id, method, resolution and release dates; the unit cell and space group; the entities, chains and
sequences; the ligands and modified residues; and the chemical components with their atoms and bonds. It also
brings these tables in step with the atoms: chains, entities and modified residues with no atoms left are removed,
the number of molecules and the chains of each entity are recomputed, modified residues get their current author
numbers, and the table that places each chain's residues on its sequence (`_pdbx_poly_seq_scheme`) is rebuilt.
For 4HHB that is 18 categories instead of 72, and the file is still valid mmCIF.

```bash
mmcif_fetch 4hhb | mmcif_delentity -non-polymer,water | mmcif_reres -10 | mmcif_reatom -4 | mmcif_delmeta > 4hhb_clean.cif
```

We suggest running `mmcif_delmeta` as the last step, after you have completed the structure manipulations.

The kept categories are also those that structure prediction programs read from template or input files (AlphaFold 2
and 3, ColabFold, OpenFold3, Boltz, Chai-1, Protenix, BioEmu), and MD preparation programs from structures (OpenMM,
PDBFixer).

`mmcif_delmeta -all` keeps only the atoms and links, and `mmcif_keepcat` keeps exactly the categories you list.

### Tips

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

#### Files for prediction and MD programs

- **AlphaFold 3** templates must hold a single protein chain: select it with `mmcif_selchain` (ligands may stay).
- **Boltz-2 and OpenFold3** name template chains by their label id (`_struct_asym.id`), not the author chain. To
  see which is which: `mmcif_cat -_atom_site.label_asym_id,auth_asym_id file.cif | sort -u`.
- **Renumber the author numbers only** (the default of `mmcif_reres` and `mmcif_shiftres`). These programs place
  residues on the sequence by `label_seq_id`, which `--label` renumbering changes.
- **Residues deleted inside a chain are still in the sequence**: prediction programs treat them as unresolved, and
  PDBFixer rebuilds them as missing residues. `mmcif_delmeta -all` drops the sequences as well.
- **OpenMM takes `_cell` as the periodic box.** Cryo-EM and NMR entries often carry a placeholder cell of 1 Å;
  remove it with `mmcif_delcat -_cell,_symmetry`, or set the box when solvating.

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

mmCIF (PDBx/mmCIF) is the format of the PDB archive: the wwPDB made it the standard in 2014, and structures too
large for the PDB format (more than 99,999 atoms or 62 chains) exist only as mmCIF. Two of its differences from PDB
files surprise most people: every atom has two sets of chain and residue identifiers, and the atoms are one table
among dozens. The examples below use 4HHB (`mmcif_fetch 4hhb > 4hhb.cif`).

### mmCIF and PDB files

The first atom of 4HHB in a PDB file, where each value has fixed columns:

```text
ATOM      1  N   VAL A   1      19.323  29.727  42.781  1.00 49.05           N
```

And in the mmCIF file, where the values are separated by spaces and named by the header lines above the rows:

```text
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   1    N  N   . VAL A 1 1   ? 19.323  29.727  42.781 1.00 49.05 ? 1   VAL A N   1
```

The *n*-th value of a row belongs to the *n*-th header line: `VAL` is the `label_comp_id`, the `A` after it the
`label_asym_id` (the chain), then `1` the `label_entity_id` and `1` the `label_seq_id` (the residue number). `.`
means "does not apply" and `?` "unknown". Files can have other columns, in another order, so programs (and these
tools) find a column by its name, never by its position.

The atoms are one table (a *category*) among many: 4HHB has 72, from `_entry` and `_cell` to the entities, the
sequences, the secondary structure and the validation reports. A table with a single row is written as name–value
pairs (`_entry.id 4HHB`) instead of a `loop_`.

### Two sets of identifiers: auth and label

Every atom names its chain and residue twice:

- **author identifiers** (`auth_asym_id`, `auth_seq_id`, near the end of the row) are the chains and numbers chosen
  by the authors. Papers, PDB files and most viewers (PyMOL and ChimeraX by default) use them.
- **label identifiers** (`label_asym_id`, `label_seq_id`, near the start of the row) are the PDB's own: one chain
  per molecule, and residues numbered 1, 2, 3, ... along the sequence. Ligands and waters have no sequence, so
  their `label_seq_id` is `.`. Structure prediction programs (AlphaFold, Boltz, ...) use the label numbers to place
  residues on the sequence.

The heme of the first chain of 4HHB is in chain `E`, with no residue number, by its label identifiers, and residue
142 of chain `A` by its author identifiers (as in the PDB file, `HEM A 142`):

```text
HETATM 4385 C  CHA . HEM E 3 .   ? 18.675  18.641  20.464 1.00 16.31 ? 142 HEM A CHA 1
```

`mmcif_cat` prints columns by name, so it shows both chains of every molecule:

```bash
mmcif_cat -_atom_site.label_asym_id,auth_asym_id,label_entity_id 4hhb.cif | uniq
```

```text
label_asym_id	auth_asym_id	label_entity_id
A	A	1
B	B	2
C	C	1
D	D	2
E	A	3
F	B	4
G	B	3
H	C	3
I	D	4
J	D	3
K	A	5
L	B	5
M	C	5
N	D	5
```

14 label chains, one per molecule (four proteins, four hemes, two phosphates and the waters of each chain), in 4
author chains, `A` to `D`, each holding a protein with its ligands and waters.

The tools use the author identifiers by default. Pass `--label` to use the label ones, or set
`MMCIFTOOLS_SCHEME=label` for every tool. **An edit changes one set only.** Residue and atom names are the same in
both, so `mmcif_rplresname` and `mmcif_uniqname` change both name columns.

### After renumbering: what changed

`mmcif_reres -101` numbers the residues from 101. The first atom of chain A, before and after:

```text
ATOM   1    N  N   . VAL A 1 1   ? 19.323  29.727  42.781 1.00 49.05 ? 1   VAL A N   1
ATOM   1    N  N   . VAL A 1 1   ? 19.323  29.727  42.781 1.00 49.05 ? 101 VAL A N   1
```

The `1` near the start of the row is the `label_seq_id`, and it stays 1: the new number is the `auth_seq_id`, near
the end. That is what you want: viewers show residue 101, and prediction programs still find it at position 1 of
the sequence. To see both numbers side by side:

```bash
mmcif_selchain -A 4hhb.cif | mmcif_reres -101 | mmcif_cat -_atom_site.label_seq_id,auth_seq_id,label_comp_id | uniq | head -5
```

```text
label_seq_id	auth_seq_id	label_comp_id
1	101	VAL
2	102	LEU
3	103	SER
4	104	PRO
```

The links in `_struct_conn` follow: the bond between the heme's iron and its histidine now names `HIS 187` and
`HEM 242` in its author columns. `mmcif_reres -101 --label` would renumber the `label_seq_id` instead, and warns
that this breaks the link between the atoms and the sequence.

### Reading your output

- **Find a column by its header line:** the *n*-th `_atom_site.` line names the *n*-th value of each atom row.
- **Print the columns you need** with `mmcif_cat`, as a table with one header line (for `less`, awk, R or pandas):
  `mmcif_cat -_atom_site.auth_asym_id,auth_seq_id,label_comp_id,label_atom_id out.cif`.
- **Summaries:** `mmcif_wc` counts models, chains, residues and atoms; `mmcif_info` shows the entry's metadata;
  `mmcif_validate` checks the syntax, the coordinates and the links between tables.
- **Compare with the input:** the tools write the rows they do not change byte-for-byte, and in changed rows only
  the edited values, so `diff 4hhb.cif out.cif` shows exactly what changed.
- **Look at it:** PyMOL, ChimeraX and Mol* open mmCIF files directly.

### Metadata goes stale: finish with mmcif_delmeta

The edit tools keep the atoms in step with their anisotropic B-factors (`_atom_site_anisotrop`) and their links
(`_struct_conn`). The other tables that name chains or residues are left as they were. After `mmcif_reres -101`,
the table that places each residue on the sequence still has the old numbers:

```bash
mmcif_selchain -A 4hhb.cif | mmcif_reres -101 | mmcif_cat -_pdbx_poly_seq_scheme.seq_id,mon_id,pdb_seq_num,auth_seq_num | head -3
```

```text
seq_id	mon_id	pdb_seq_num	auth_seq_num
1	VAL	1	1
2	LEU	2	2
```

and so do the secondary structure (`_struct_conf`), the validation reports and the binding sites. Ending the
pipeline with `mmcif_delmeta` deletes the tables that went stale and brings the ones it keeps in step with the atoms:
with `| mmcif_delmeta` before `mmcif_cat`, the same command prints `101` and `102`. Chain A of 4HHB then has 18
categories instead of 72. See [`mmcif_delmeta`](#mmcif_delmeta) for what it keeps.

### Entities

Every molecule belongs to an entity: polymer, non-polymer, branched (glycans) or water.

```bash
mmcif_cat -_entity.id,type,pdbx_description 4hhb.cif
```

```text
id	type	pdbx_description
1	polymer	Hemoglobin subunit alpha
2	polymer	Hemoglobin subunit beta
3	non-polymer	PROTOPORPHYRIN IX CONTAINING FE
4	non-polymer	PHOSPHATE ION
5	water	water
```

`mmcif_selentity` / `mmcif_delentity` select by entity, which keeps modified residues (HETATM records inside a
chain, such as MSE or the modified nucleotides of a tRNA) with their chain; `mmcif_delhetatm` removes them.

### No fixed columns

Chain identifiers can be longer than one character, and there are no limits on atom or residue numbers.
`mmcif_topdb` refuses structures that do not fit the PDB format.

## Installation

```bash
python -m pip install mmcif-tools
```

There are no runtime dependencies. For tab completion of the `mmcif_*` commands, open a new terminal after
installing (or run `rehash` in an open zsh terminal).

The extra tools below are in the repository, not in the package. To install from the repository instead,
`python tools/smoke_test.py` then checks the installation: every command runs, and a few pipelines give the
expected results.

```bash
git clone https://github.com/matteo-cagiada/mmcif-tools.git
cd mmcif-tools
python -m pip install .
python tools/smoke_test.py
```

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

- **`tools/make_docs.py`** writes the [tool reference](https://mmcif-tools.dev/docs/tools/) and
  the other pages of the [documentation](https://mmcif-tools.dev/docs/).
- **`tools/smoke_test.py`** checks an installation (see
  [Installation](https://mmcif-tools.dev/docs/installation/)).

## Citation and changes

How to cite mmcif-tools is in [CITATION.cff](https://github.com/matteo-cagiada/mmcif-tools/blob/main/CITATION.cff);
the changes between versions are in the [changelog](https://mmcif-tools.dev/docs/changelog/).
[PLAN.md](https://github.com/matteo-cagiada/mmcif-tools/blob/main/PLAN.md) records the design decisions and the test
results.

## Acknowledgements

Most of the code of mmcif-tools has been written with [Claude Code](https://claude.com/claude-code),
Anthropic's AI coding assistant. The author designed the project, tools and tests to guarantee the quality of the tool.
We thank the Bonvin Lab for developing and releasing pdb-tools, a fundamental tool for the authors' academic work.

## License

Apache-2.0. See [LICENSE](https://github.com/matteo-cagiada/mmcif-tools/blob/main/LICENSE).
