<!-- Generated from README.md by tools/make_docs.py: edit README.md and run python tools/make_docs.py. -->

# Usage examples and tips

## Command examples

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

## `mmcif_delmeta`

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

## Tips

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

### Files for prediction and MD programs

- **AlphaFold 3** templates must hold a single protein chain: select it with `mmcif_selchain` (ligands may stay).
- **Boltz-2 and OpenFold3** name template chains by their label id (`_struct_asym.id`), not the author chain. To
  see which is which: `mmcif_cat -_atom_site.label_asym_id,auth_asym_id file.cif | sort -u`.
- **Renumber the author numbers only** (the default of `mmcif_reres` and `mmcif_shiftres`). These programs place
  residues on the sequence by `label_seq_id`, which `--label` renumbering changes.
- **Residues deleted inside a chain are still in the sequence**: prediction programs treat them as unresolved, and
  PDBFixer rebuilds them as missing residues. `mmcif_delmeta -all` drops the sequences as well.
- **OpenMM takes `_cell` as the periodic box.** Cryo-EM and NMR entries often carry a placeholder cell of 1 Å;
  remove it with `mmcif_delcat -_cell,_symmetry`, or set the box when solvating.
