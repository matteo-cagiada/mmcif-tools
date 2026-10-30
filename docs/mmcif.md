<!-- Generated from README.md by tools/make_docs.py: edit README.md and run python tools/make_docs.py. -->

# Things to know about mmCIF

mmCIF (PDBx/mmCIF) is the format of the PDB archive: the wwPDB made it the standard in 2014, and structures too
large for the PDB format (more than 99,999 atoms or 62 chains) exist only as mmCIF. Two of its differences from PDB
files surprise most people: every atom has two sets of chain and residue identifiers, and the atoms are one table
among dozens. The examples below use 4HHB (`mmcif_fetch 4hhb > 4hhb.cif`).

## mmCIF and PDB files

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

## Two sets of identifiers: auth and label

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

## After renumbering: what changed

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

## Reading your output

- **Find a column by its header line:** the *n*-th `_atom_site.` line names the *n*-th value of each atom row.
- **Print the columns you need** with `mmcif_cat`, as a table with one header line (for `less`, awk, R or pandas):
  `mmcif_cat -_atom_site.auth_asym_id,auth_seq_id,label_comp_id,label_atom_id out.cif`.
- **Summaries:** `mmcif_wc` counts models, chains, residues and atoms; `mmcif_info` shows the entry's metadata;
  `mmcif_validate` checks the syntax, the coordinates and the links between tables.
- **Compare with the input:** the tools write the rows they do not change byte-for-byte, and in changed rows only
  the edited values, so `diff 4hhb.cif out.cif` shows exactly what changed.
- **Look at it:** PyMOL, ChimeraX and Mol* open mmCIF files directly.

## Metadata goes stale: finish with mmcif_delmeta

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
categories instead of 72. See [`mmcif_delmeta`](usage.md#mmcif_delmeta) for what it keeps.

## Entities

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

## No fixed columns

Chain identifiers can be longer than one character, and there are no limits on atom or residue numbers.
`mmcif_topdb` refuses structures that do not fit the PDB format.
