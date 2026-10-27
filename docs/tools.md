# mmcif-tools: tool reference

Version 0.1.0, 53 tools. Generated from the tools' `--help` texts by `python tools/make_docs.py`; do not edit by hand.

| Group | Tools |
|---|---|
| Fetch / convert | [`mmcif_fetch`](#mmcif_fetch), [`mmcif_topdb`](#mmcif_topdb), [`mmcif_frompdb`](#mmcif_frompdb) |
| Select / delete | [`mmcif_selchain`](#mmcif_selchain), [`mmcif_delchain`](#mmcif_delchain), [`mmcif_selres`](#mmcif_selres), [`mmcif_delres`](#mmcif_delres), [`mmcif_selresname`](#mmcif_selresname), [`mmcif_delresname`](#mmcif_delresname), [`mmcif_selatom`](#mmcif_selatom), [`mmcif_selelem`](#mmcif_selelem), [`mmcif_delelem`](#mmcif_delelem), [`mmcif_selhetatm`](#mmcif_selhetatm), [`mmcif_delhetatm`](#mmcif_delhetatm), [`mmcif_selentity`](#mmcif_selentity), [`mmcif_delentity`](#mmcif_delentity), [`mmcif_selaltloc`](#mmcif_selaltloc), [`mmcif_selmodel`](#mmcif_selmodel), [`mmcif_head`](#mmcif_head), [`mmcif_delmeta`](#mmcif_delmeta), [`mmcif_keepcoord`](#mmcif_keepcoord), [`mmcif_intersect`](#mmcif_intersect) |
| Edit | [`mmcif_chain`](#mmcif_chain), [`mmcif_rplchain`](#mmcif_rplchain), [`mmcif_reres`](#mmcif_reres), [`mmcif_shiftres`](#mmcif_shiftres), [`mmcif_fixinsert`](#mmcif_fixinsert), [`mmcif_delinsertion`](#mmcif_delinsertion), [`mmcif_reatom`](#mmcif_reatom), [`mmcif_rplresname`](#mmcif_rplresname), [`mmcif_uniqname`](#mmcif_uniqname), [`mmcif_element`](#mmcif_element), [`mmcif_b`](#mmcif_b), [`mmcif_occ`](#mmcif_occ) |
| Identifiers | [`mmcif_label2auth`](#mmcif_label2auth), [`mmcif_auth2label`](#mmcif_auth2label) |
| Models, files and assemblies | [`mmcif_splitchain`](#mmcif_splitchain), [`mmcif_splitmodel`](#mmcif_splitmodel), [`mmcif_mkensemble`](#mmcif_mkensemble), [`mmcif_chkensemble`](#mmcif_chkensemble), [`mmcif_merge`](#mmcif_merge), [`mmcif_assembly`](#mmcif_assembly) |
| Categories and data blocks | [`mmcif_cat`](#mmcif_cat), [`mmcif_delcat`](#mmcif_delcat), [`mmcif_keepcat`](#mmcif_keepcat), [`mmcif_blocks`](#mmcif_blocks) |
| Inspect / format | [`mmcif_info`](#mmcif_info), [`mmcif_wc`](#mmcif_wc), [`mmcif_tofasta`](#mmcif_tofasta), [`mmcif_gap`](#mmcif_gap), [`mmcif_sort`](#mmcif_sort), [`mmcif_tidy`](#mmcif_tidy), [`mmcif_validate`](#mmcif_validate) |

## Fetch / convert

### mmcif_fetch

Downloads a structure in mmCIF format and writes it to standard output.

```text
Downloads a structure in mmCIF format and writes it to standard output.

PDB entries (4-character codes, or extended ids such as pdb_00001abc) come
from the RCSB; -biounit downloads the first biological assembly instead of
the deposited coordinates. AlphaFold DB models (AF-<UniProt accession>, e.g.
AF-P69905-F1) are looked up through the AlphaFold DB API, so the latest
model version is used.

Usage:
    mmcif_fetch [-biounit] <pdb code | AlphaFold DB id>

Example:
    mmcif_fetch 1brs > 1brs.cif            # deposited coordinates, all 6 chains
    mmcif_fetch -biounit 1brs > 1brs.cif   # first biological assembly, 2 chains
    mmcif_fetch AF-P69905-F1 > hba.cif     # AlphaFold DB model

This program is part of the mmcif-tools suite.
```

### mmcif_topdb

Converts the coordinates of an mmCIF file to the PDB format.

```text
Converts the coordinates of an mmCIF file to the PDB format.

As pdb_fromcif, writes the ATOM/HETATM records of _atom_site (with MODEL /
ENDMDL records when there are several models) and END; other categories,
and data blocks after the first, are not converted. Chains, residue numbers and names are auth identifiers,
falling back to label ones when a file lacks them (with --label, label ones
are used). Atom serials restart at 1 in each model. Unlike pdb_fromcif, a
TER record follows the last residue of each polymer (atoms with a
label_seq_id, or ATOM records in files without one) and takes a serial
number, as the PDB format specifies.

Will not convert a structure that does not fit the PDB format: chains of
more than one character, more than 99,999 atoms in a model, residue numbers
outside -999..9999, residue names longer than 3 characters or atom names
longer than 4; nothing is written then. Unlike pdb_fromcif, residue names
and elements are right-justified and MODEL numbers placed as the PDB
format specifies, formal charges are written as 1-, 2+, ..., and missing
occupancies / B-factors are written as 1.00 / 0.00.

Usage:
    mmcif_topdb <cif file>

Example:
    mmcif_topdb 1abc.cif > 1abc.pdb

This program is part of the mmcif-tools suite.
```

### mmcif_frompdb

Converts a PDB file to mmCIF format.

```text
Converts a PDB file to mmCIF format.

Converts the coordinates (ATOM/HETATM records, models) into _atom_site,
plus _entry and, from CRYST1, _cell and _symmetry. Author identifiers
(auth_*) are those of the PDB file. Unlike pdb_tocif, which copies them into
the label_* columns, label identifiers are assigned as in PDB entries:

    label_asym_id    one per polymer chain, one per ligand, and one for the
                     waters of each chain (A, B, ..., Z, AA, BA, ...)
    label_entity_id  polymer chains with the same sequence share an entity;
                     each ligand type and water are entities of their own
    label_seq_id     1, 2, ... along each polymer chain; '.' elsewhere

A TER record ends a polymer: the residues of a chain up to a TER record
are one polymer, and the next TER record of the same chain ends another.
So a file whose chains are blank or repeated, such as MD output with one
TER record per molecule, gets one label_asym_id per molecule. In chains
without TER records, the polymer runs from the first to the last ATOM
record (so modified residues inside a chain belong to it).

Atoms with a blank chain get their label_asym_id as their author chain, or
the next free name when another chain of the file already has it (this
replaces pdb_chainbows | pdb_tocif).

For chains without missing residues this gives the label identifiers of the
PDB entry; but label_seq_id only counts the residues present in the file
(SEQRES is not used), so it differs from the entry after missing residues,
and PDB entries order ligand instances and group glycans (branched
entities) by their own conventions. Missing elements are assigned as by
mmcif_element. Atom ids run from 1 over the whole file; models keep their
MODEL numbers. ANISOU, CONECT and the other PDB records are not converted.

Usage:
    mmcif_frompdb <pdb file>

Example:
    mmcif_frompdb 1abc.pdb > 1abc.cif

This program is part of the mmcif-tools suite.
```

## Select / delete

### mmcif_selchain

Extracts one or more chains from an mmCIF file.

```text
Extracts one or more chains from an mmCIF file.

Chains are matched on auth_asym_id (the chain shown by viewers), or on
label_asym_id with --label. Chain identifiers are case-sensitive and may be
longer than one character. Rows of _atom_site_anisotrop and _struct_conn
that refer to removed atoms are removed too.

Usage:
    mmcif_selchain -<chain id>[,<chain id>...] [--label] <cif file>

Example:
    mmcif_selchain -C 1abc.cif     # selects chain C
    mmcif_selchain -A,C 1abc.cif   # selects chains A and C

This program is part of the mmcif-tools suite.
```

### mmcif_delchain

Deletes all atoms of one or more chains from an mmCIF file.

```text
Deletes all atoms of one or more chains from an mmCIF file.

Chains are matched on auth_asym_id (the chain shown by viewers), or on
label_asym_id with --label. Chain identifiers are case-sensitive (unlike
pdb_delchain, which upper-cases them) and may be longer than one character.
Rows of _atom_site_anisotrop and _struct_conn that refer to removed atoms
are removed too.

Usage:
    mmcif_delchain -<chain id>[,<chain id>...] [--label] <cif file>

Example:
    mmcif_delchain -A 1abc.cif     # removes chain A
    mmcif_delchain -A,B 1abc.cif   # removes chains A and B

This program is part of the mmcif-tools suite.
```

### mmcif_selres

Selects residues by their number, piecewise or in a range.

```text
Selects residues by their number, piecewise or in a range.

The range option has three components: start, end, and step. Start and end
are optional and if omitted the range will start at the first residue or
end at the last, respectively. As in pdb_selres, residues are matched by
number only (in every chain and model, with or without insertion code), and
the step counts residue occurrences in file order.

Residue numbers are auth_seq_id, or label_seq_id with --label. Rows of
_atom_site_anisotrop and _struct_conn that refer to removed atoms are
removed too.

Usage:
    mmcif_selres -[resid]:[resid]:[step] [--label] <cif file>

Example:
    mmcif_selres -1,2,4,6 1abc.cif     # extracts residues 1, 2, 4 and 6
    mmcif_selres -1:10 1abc.cif        # extracts residues 1 to 10
    mmcif_selres -1:10,20:30 1abc.cif  # extracts residues 1 to 10 and 20 to 30
    mmcif_selres -1: 1abc.cif          # extracts residues 1 to END
    mmcif_selres -:5 1abc.cif          # extracts residues from START to 5
    mmcif_selres -::5 1abc.cif         # extracts every 5th residue
    mmcif_selres -1:10:5 1abc.cif      # extracts every 5th residue from 1 to 10
    mmcif_selres --5:5 1abc.cif        # extracts residues -5 to 5

This program is part of the mmcif-tools suite.
```

### mmcif_delres

Deletes a range of residues from an mmCIF file.

```text
Deletes a range of residues from an mmCIF file.

The range option has three components: start, end, and step. Start and end
are optional and if omitted the range will start at the first residue or
end at the last, respectively. As in pdb_delres, the start and end values
are residue numbers (matched in every chain and model), while the step
counts residue occurrences in file order, regardless of their number.
Unlike pdb_delres, a range may start and end at the same residue (-5:5).

Residue numbers are auth_seq_id, or label_seq_id with --label. Rows of
_atom_site_anisotrop and _struct_conn that refer to removed atoms are
removed too.

Usage:
    mmcif_delres -[resid]:[resid]:[step] [--label] <cif file>

Example:
    mmcif_delres -1:10 1abc.cif     # deletes residues 1 to 10
    mmcif_delres -1: 1abc.cif       # deletes residues 1 to END
    mmcif_delres -:5 1abc.cif       # deletes residues from START to 5
    mmcif_delres -::5 1abc.cif      # deletes every 5th residue
    mmcif_delres -1:10:5 1abc.cif   # deletes every 5th residue from 1 to 10

This program is part of the mmcif-tools suite.
```

### mmcif_selresname

Selects all residues matching the given name(s) in an mmCIF file.

```text
Selects all residues matching the given name(s) in an mmCIF file.

Residue names are matched on auth_comp_id (label_comp_id with --label);
the names given are upper-cased, as in pdb_selresname. Rows of
_atom_site_anisotrop and _struct_conn that refer to removed atoms are
removed too.

Usage:
    mmcif_selresname -<name>[,<name>...] <cif file>

Example:
    mmcif_selresname -ALA 1abc.cif      # keeps only alanines
    mmcif_selresname -ASP,GLU 1abc.cif  # keeps (-) charged residues

This program is part of the mmcif-tools suite.
```

### mmcif_delresname

Removes all residues matching the given name(s) in an mmCIF file.

```text
Removes all residues matching the given name(s) in an mmCIF file.

Residue names are matched on auth_comp_id (label_comp_id with --label);
the names given are upper-cased, as in pdb_delresname. Rows of
_atom_site_anisotrop and _struct_conn that refer to removed atoms are
removed too.

Usage:
    mmcif_delresname -<name>[,<name>...] <cif file>

Example:
    mmcif_delresname -ALA 1abc.cif      # removes only alanines
    mmcif_delresname -HOH 1abc.cif      # removes waters
    mmcif_delresname -ASP,GLU 1abc.cif  # removes (-) charged residues

This program is part of the mmcif-tools suite.
```

### mmcif_selatom

Selects all atoms matching the given name(s) in an mmCIF file.

```text
Selects all atoms matching the given name(s) in an mmCIF file.

Atom names are matched on auth_atom_id (label_atom_id with --label); the
names given are upper-cased, as in pdb_selatom. mmCIF atom names carry no
padding, so -CA keeps alpha carbons and calcium ions alike, as pdb_selatom
does. Rows of _atom_site_anisotrop and _struct_conn that refer to removed
atoms are removed too.

Usage:
    mmcif_selatom -<name>[,<name>...] <cif file>

Example:
    mmcif_selatom -CA 1abc.cif        # keeps only alpha-carbon atoms
    mmcif_selatom -CA,C,N,O 1abc.cif  # keeps only backbone atoms

This program is part of the mmcif-tools suite.
```

### mmcif_selelem

Selects all atoms of the given element(s) in an mmCIF file.

```text
Selects all atoms of the given element(s) in an mmCIF file.

Elements are read from _atom_site.type_symbol and compared ignoring case
(pdb_selelem only matches upper-case elements). Rows of
_atom_site_anisotrop and _struct_conn that refer to removed atoms are
removed too.

Usage:
    mmcif_selelem -<element>[,<element>...] <cif file>

Example:
    mmcif_selelem -H 1abc.cif     # selects all hydrogens
    mmcif_selelem -N 1abc.cif     # selects all nitrogens
    mmcif_selelem -H,N 1abc.cif   # selects all hydrogens and nitrogens

This program is part of the mmcif-tools suite.
```

### mmcif_delelem

Deletes all atoms of the given element(s) from an mmCIF file.

```text
Deletes all atoms of the given element(s) from an mmCIF file.

Elements are read from _atom_site.type_symbol and compared ignoring case
(pdb_delelem only matches upper-case elements). Rows of
_atom_site_anisotrop and _struct_conn that refer to removed atoms are
removed too.

Usage:
    mmcif_delelem -<element>[,<element>...] <cif file>

Example:
    mmcif_delelem -H 1abc.cif     # deletes all hydrogens
    mmcif_delelem -H,D 1abc.cif   # deletes hydrogens and deuteriums

This program is part of the mmcif-tools suite.
```

### mmcif_selhetatm

Selects all HETATM records (ligands, ions, waters, ...) in an mmCIF file.

```text
Selects all HETATM records (ligands, ions, waters, ...) in an mmCIF file.

Atoms are selected on _atom_site.group_PDB. Unlike pdb_selhetatm, which
drops every other record, the rest of the file (metadata categories) is
kept; rows of _atom_site_anisotrop and _struct_conn that refer to removed
atoms are removed.

Usage:
    mmcif_selhetatm <cif file>

Example:
    mmcif_selhetatm 1abc.cif

This program is part of the mmcif-tools suite.
```

### mmcif_delhetatm

Removes all HETATM records (ligands, ions, waters, ...) from an mmCIF file.

```text
Removes all HETATM records (ligands, ions, waters, ...) from an mmCIF file.

Atoms are selected on _atom_site.group_PDB. Rows of _atom_site_anisotrop
and _struct_conn that refer to removed atoms are removed too.

Usage:
    mmcif_delhetatm <cif file>

Example:
    mmcif_delhetatm 1abc.cif

This program is part of the mmcif-tools suite.
```

### mmcif_selentity

Selects atoms by entity: by entity number (label_entity_id) or by type.

```text
Selects atoms by entity: by entity number (label_entity_id) or by type.

In mmCIF every molecule belongs to an entity (_entity), whose type is one
of:

    polymer       protein and nucleic-acid chains
    non-polymer   ligands and ions
    branched      oligosaccharides (glycans)
    macrolide     macrolides
    water         water

Unlike mmcif_selhetatm, this keeps modified residues with their chain: they
are HETATM records inside a polymer (MSE in proteins, PSU and 2MG in tRNA,
...). Numbers and types can be mixed (-polymer,5); types are
case-insensitive. Files without an _entity category (some non-PDB files) are
classified from their atoms: atoms with a label_seq_id are polymer, HOH,
WAT and DOD are water, and the rest non-polymer. Rows of
_atom_site_anisotrop and _struct_conn of removed atoms are removed too.

Usage:
    mmcif_selentity -<entity id or type>[,...] <cif file>

Example:
    mmcif_selentity -polymer 1abc.cif           # chains only, with modified residues
    mmcif_selentity -polymer,branched 1abc.cif  # chains and glycans
    mmcif_selentity -1,3 1abc.cif               # entities 1 and 3

This program is part of the mmcif-tools suite.
```

### mmcif_delentity

Deletes atoms by entity: by entity number (label_entity_id) or by type.

```text
Deletes atoms by entity: by entity number (label_entity_id) or by type.

Entity types are polymer, non-polymer (ligands and ions), branched
(glycans), macrolide and water; see mmcif_selentity. Unlike
mmcif_delhetatm, this never cuts modified residues (HETATM records inside
a polymer) out of their chain: -non-polymer,water removes the ligands, ions
and waters and nothing else. Rows of _atom_site_anisotrop and _struct_conn
of removed atoms are removed too.

Usage:
    mmcif_delentity -<entity id or type>[,...] <cif file>

Example:
    mmcif_delentity -water 1abc.cif
    mmcif_delentity -non-polymer,water 1abc.cif   # keeps chains and glycans
    mmcif_delentity -2 1abc.cif                   # removes entity 2

This program is part of the mmcif-tools suite.
```

### mmcif_selaltloc

Selects alternate locations (_atom_site.label_alt_id) for the entire file.

```text
Selects alternate locations (_atom_site.label_alt_id) for the entire file.

By default, for each residue the alternate location with the highest
occupancy is kept and the others are removed. With an option (e.g. -A),
location A is kept wherever a residue has it; atoms that have other
locations but no A (e.g. only B and C) are kept as they are. Kept atoms
lose their altloc label ('.'); occupancies are not changed.

The choice is made per residue, so a residue modelled as two different
residue types (e.g. SER/THR) is never mixed. Atoms keep their order, and
rows of _atom_site_anisotrop of removed atoms are removed too. Links
(_struct_conn) to a removed conformer are removed, and links to a kept one
lose their altloc label, as the atom does. A file without a label_alt_id
column has no alternate locations and is written unchanged.

Usage:
    mmcif_selaltloc [-<altloc>] <cif file>

Example:
    mmcif_selaltloc 1abc.cif      # picks locations with highest occupancy
    mmcif_selaltloc -A 1abc.cif   # picks alternate locations labelled 'A'

This program is part of the mmcif-tools suite.
```

### mmcif_selmodel

Extracts one or more models from an mmCIF file.

```text
Extracts one or more models from an mmCIF file.

Models are read from _atom_site.pdbx_PDB_model_num. As in pdb_selmodel, a
file with a single model (or without model numbers) is returned whole,
whatever model is asked for. Rows of _atom_site_anisotrop and _struct_conn
that refer to removed atoms are removed too.

Usage:
    mmcif_selmodel -<model id>[,<model id>...] <cif file>

Example:
    mmcif_selmodel -1 1abc.cif     # selects model 1
    mmcif_selmodel -1,3 1abc.cif   # selects models 1 and 3

This program is part of the mmcif-tools suite.
```

### mmcif_head

Keeps the first N atoms (_atom_site rows) of an mmCIF file.

```text
Keeps the first N atoms (_atom_site rows) of an mmCIF file.

Atoms are counted over the whole file, across models and data blocks.
Unlike pdb_head, which outputs nothing but the first N ATOM/HETATM lines,
the rest of the file (metadata categories) is kept, so the output is still
a complete mmCIF file;
rows of _atom_site_anisotrop and _struct_conn that refer to removed atoms
are removed.

Usage:
    mmcif_head -<num> <cif file>

Example:
    mmcif_head -100 1abc.cif   # first 100 atoms of the file

This program is part of the mmcif-tools suite.
```

### mmcif_delmeta

Deletes the metadata that edits make stale, keeping the coordinates and the metadata that stays true.

```text
Deletes the metadata that edits make stale, keeping the coordinates and the
metadata that stays true.

PDB entries carry dozens of categories that refer to residues by chain and
number: validation reports, secondary structure, sequence schemes, sites,
assemblies, ... The other tools do not update them, so after an edit they no
longer match the atoms, and most programs that read mmCIF ignore them anyway.
This tool keeps:

    _atom_site                 the atoms
    _atom_site_anisotrop       their anisotropic B-factors
    _struct_conn               the links (disulfides, metal sites, ...)
    _entry                     the entry id
    _cell, _symmetry           the unit cell and space group
    _entity, _entity_poly,     the molecules and their sequences
    _entity_poly_seq
    _chem_comp, _atom_type     the chemical components and elements

and deletes every other category. With -all, only the atoms, anisotropic
B-factors and links are kept (mmcif_keepcoord also drops the anisotropic
B-factors, as pdb_keepcoord drops ANISOU records).

Run it first in a pipeline, before the edits. Two kept items can still
disagree after edits: _entity keeps describing molecules whose atoms a later
tool removes, and _entity_poly.pdbx_strand_id names author chains, which
chain renames do not update.

Usage:
    mmcif_delmeta [-all] <cif file>

Example:
    mmcif_delmeta 1abc.cif > 1abc_clean.cif
    mmcif_delmeta 1abc.cif | mmcif_delentity -non-polymer,water > protein.cif
    mmcif_delmeta -all 1abc.cif > atoms_only.cif

This program is part of the mmcif-tools suite.
```

### mmcif_keepcoord

Removes all non-coordinate categories from an mmCIF file.

```text
Removes all non-coordinate categories from an mmCIF file.

Keeps the data_ line, _atom_site and _struct_conn (the mmCIF counterparts of
the ATOM/HETATM and CONECT records that pdb_keepcoord keeps); every other
category is removed, _atom_site_anisotrop included (pdb_keepcoord drops
ANISOU records too). To keep the metadata that edits do not make stale
(entry, unit cell, entities, sequences, chemical components), use
mmcif_delmeta instead.

Usage:
    mmcif_keepcoord <cif file>

Example:
    mmcif_keepcoord 1abc.cif

This program is part of the mmcif-tools suite.
```

### mmcif_intersect

Keeps only the atoms of the first mmCIF file that are present in all the others.

```text
Keeps only the atoms of the first mmCIF file that are present in all the
others.

As in pdb_intersect, atoms are judged equal if their name, altloc, residue
name, residue number, insertion code and chain are the same (auth
identifiers, or label ones with --label; models are not compared), and
coordinates are taken from the first file. Unlike pdb_intersect, the output
is the first file with the other atoms removed (metadata kept; rows of
_atom_site_anisotrop and _struct_conn of removed atoms removed).

Usage:
    mmcif_intersect <cif file> <cif file> [...]

Example:
    mmcif_intersect 1xyz.cif 1abc.cif

This program is part of the mmcif-tools suite.
```

## Edit

### mmcif_chain

Sets the chain identifier of all atoms in an mmCIF file.

```text
Sets the chain identifier of all atoms in an mmCIF file.

Changes auth_asym_id (label_asym_id with --label), in _atom_site and in the
matching columns of _atom_site_anisotrop and _struct_conn. mmCIF has no
blank chain, so unlike pdb_chain an identifier must be given; it may be
longer than one character.

Usage:
    mmcif_chain -<chain id> [--label] <cif file>

Example:
    mmcif_chain -C 1abc.cif

This program is part of the mmcif-tools suite.
```

### mmcif_rplchain

Replaces one chain identifier by another in an mmCIF file.

```text
Replaces one chain identifier by another in an mmCIF file.

Changes auth_asym_id (label_asym_id with --label), in _atom_site and in the
matching columns of _atom_site_anisotrop and _struct_conn. Identifiers are
case-sensitive and may be longer than one character; mmCIF has no blank
chain, so both must be given.

Usage:
    mmcif_rplchain -<from>:<to> [--label] <cif file>

Example:
    mmcif_rplchain -A:B 1abc.cif   # replaces chain A by chain B

This program is part of the mmcif-tools suite.
```

### mmcif_reres

Renumbers the residues of an mmCIF file starting from a given number (default 1).

```text
Renumbers the residues of an mmCIF file starting from a given number
(default 1).

As in pdb_reres, numbering continues across chains and restarts in each
model; a new residue starts whenever the residue name, chain, number or
insertion code changes from one atom to the next (so 52 and 52A become two
numbers), and insertion codes are kept. There is no upper limit.

Changes auth_seq_id by default, in _atom_site and in the matching columns of
_atom_site_anisotrop and _struct_conn. With --label it changes label_seq_id
of polymer residues instead; that number links atoms to the sequence in
_entity_poly_seq, so a warning is printed.

Usage:
    mmcif_reres -<number> [--label] <cif file>

Example:
    mmcif_reres -10 1abc.cif   # renumbers from 10
    mmcif_reres --1 1abc.cif   # renumbers from -1

This program is part of the mmcif-tools suite.
```

### mmcif_shiftres

Shifts the residue numbers of an mmCIF file by a constant value.

```text
Shifts the residue numbers of an mmCIF file by a constant value.

Changes auth_seq_id by default, in _atom_site and in the matching columns of
_atom_site_anisotrop and _struct_conn; insertion codes are kept. With
--label it shifts label_seq_id of polymer residues instead; that number
links atoms to the sequence in _entity_poly_seq, so a warning is printed.
There is no upper or lower limit.

Usage:
    mmcif_shiftres -<number> [--label] <cif file>

Example:
    mmcif_shiftres -10 1abc.cif   # adds 10 to the original numbering
    mmcif_shiftres --5 1abc.cif   # subtracts 5 from the original numbering

This program is part of the mmcif-tools suite.
```

### mmcif_fixinsert

Deletes insertion codes (_atom_site.pdbx_PDB_ins_code) and shifts the residue numbering (auth_seq_id) to make room.

```text
Deletes insertion codes (_atom_site.pdbx_PDB_ins_code) and shifts the
residue numbering (auth_seq_id) to make room.

As in pdb_fixinsert, the first residue carrying an insertion code just
loses it, and every later residue whose number was already seen moves the
numbering of the rest of the chain up by one: 51, 52A, 52B, 53 becomes
51, 52, 53, 54. With an option, only the listed residues lose their codes
(the numbering still shifts for the rest of the chain).

The shift restarts at each new chain instance (a change of chain or of
label_asym_id, e.g. between a protein and its ligands or waters: where
pdb_fixinsert meets a TER record) and in each model (pdb_fixinsert never
restarts between models). Chains may have several
characters and numbers may be negative: -A9,B12,AB-3. Insertion codes only
exist with author identifiers, so --label is not accepted. The matching
columns of _atom_site_anisotrop and _struct_conn are updated too. A file
without a pdbx_PDB_ins_code column (e.g. an integrative model) has no
insertion codes and is written unchanged.

Usage:
    mmcif_fixinsert [-<option>] <cif file>

Example:
    mmcif_fixinsert 1abc.cif            # deletes ALL insertion codes
    mmcif_fixinsert -A9,B12 1abc.cif    # deletes insertion codes of residue 9
                                        # of chain A and 12 of chain B

This program is part of the mmcif-tools suite.
```

### mmcif_delinsertion

Deletes insertion codes and shifts the residue numbering to make room.

```text
Deletes insertion codes and shifts the residue numbering to make room.

The same tool as mmcif_fixinsert, under the name of pdb_delinsertion (which
pdb-tools keeps as a deprecated alias of pdb_fixinsert). See mmcif_fixinsert
for the details.

Usage:
    mmcif_delinsertion [-<option>] <cif file>

Example:
    mmcif_delinsertion 1abc.cif            # deletes ALL insertion codes
    mmcif_delinsertion -A9,B12 1abc.cif    # deletes insertion codes of residue 9
                                           # of chain A and 12 of chain B

This program is part of the mmcif-tools suite.
```

### mmcif_reatom

Renumbers the atom ids (_atom_site.id) of an mmCIF file starting from a given number (default 1).

```text
Renumbers the atom ids (_atom_site.id) of an mmCIF file starting from a
given number (default 1).

The ids of _atom_site_anisotrop are changed to match. In mmCIF the atom id
must be unique in the whole file, so numbering continues across models
(pdb_reatom restarts at every MODEL), and there is no upper limit.

Usage:
    mmcif_reatom -<number> <cif file>

Example:
    mmcif_reatom -10 1abc.cif   # renumbers from 10
    mmcif_reatom --1 1abc.cif   # renumbers from -1

This program is part of the mmcif-tools suite.
```

### mmcif_rplresname

Replaces one residue name by another in an mmCIF file.

```text
Replaces one residue name by another in an mmCIF file.

Residue names are the same chemical identity in both identifier schemes,
and programs differ in which one they read, so both auth_comp_id and
label_comp_id are changed (whatever --auth/--label says), together with the
residue names in _atom_site_anisotrop and _struct_conn. Names are matched
exactly (case-sensitive) and may be longer than three characters. Other
categories that list residue names (_entity_poly_seq, _pdbx_poly_seq_scheme,
_chem_comp, ...) are not changed.

Usage:
    mmcif_rplresname -<from>:<to> <cif file>

Example:
    mmcif_rplresname -HIP:HIS 1abc.cif   # changes all HIP residues to HIS

This program is part of the mmcif-tools suite.
```

### mmcif_uniqname

Renames atoms sequentially (C1, C2, O1, ...) for each HETATM residue.

```text
Renames atoms sequentially (C1, C2, O1, ...) for each HETATM residue.

As in pdb_uniqname, only HETATM atoms are renamed (waters included), the
name is the element (_atom_site.type_symbol, see mmcif_element) followed by a
counter per element that restarts at each residue, and a missing element is
an error.

Atom names are the same chemical identity in both identifier schemes, so
both auth_atom_id and label_atom_id are changed, together with the atom
names of _atom_site_anisotrop and _struct_conn (pdb_uniqname leaves ANISOU
records alone). Residues are told apart by model, name, chain, number and
insertion code (auth identifiers, or label ones with --label).

Usage:
    mmcif_uniqname <cif file>

Example:
    mmcif_uniqname 1abc.cif

This program is part of the mmcif-tools suite.
```

### mmcif_element

Assigns the elements (_atom_site.type_symbol) of an mmCIF file from the atom names.

```text
Assigns the elements (_atom_site.type_symbol) of an mmCIF file from the
atom names.

By default only missing elements ('?', '.', or no type_symbol column) are
assigned, because deposited elements are authoritative; -all re-assigns
every atom (pdb_element always overwrites).

mmCIF atom names carry no column alignment (pdb_element tells ' CA ', an
alpha carbon, from 'CA  ', calcium, by it), so the rules are:
    1. an atom named like its residue is that element (ZN in ZN, CA in CA)
    2. SE in MSE / SEC is selenium
    3. in ligands (HETATM, not a standard residue) a name starting with a
       common two-letter element is that element (CL1, BR2, FE, ZN1, MG ...;
       not CA, NA, HG, CD, which in ligands are mostly carbon, nitrogen or
       hydrogen names, e.g. the heme nitrogens NA-ND)
    4. otherwise the first letter of the name, after leading digits
Atom names are auth_atom_id (label_atom_id with --label). Atoms for which no
element can be found get '?', and a warning is printed.

Usage:
    mmcif_element [-all] <cif file>

Example:
    mmcif_element 1abc.cif        # fills in missing elements
    mmcif_element -all 1abc.cif   # re-assigns every element

This program is part of the mmcif-tools suite.
```

### mmcif_b

Sets the temperature factor (_atom_site.B_iso_or_equiv) of all atoms (default 10.0).

```text
Sets the temperature factor (_atom_site.B_iso_or_equiv) of all atoms
(default 10.0).

Values are written with two decimals; a file without the column gets it.
Anisotropic factors (_atom_site_anisotrop) are not changed, as in pdb_b.

Usage:
    mmcif_b -<bfactor> <cif file>

Example:
    mmcif_b -10.0 1abc.cif

This program is part of the mmcif-tools suite.
```

### mmcif_occ

Sets the occupancy (_atom_site.occupancy) of all atoms (default 1.0).

```text
Sets the occupancy (_atom_site.occupancy) of all atoms (default 1.0).

Values are written with two decimals; a file without the column gets it.

Usage:
    mmcif_occ -<occupancy> <cif file>

Example:
    mmcif_occ -1.0 1abc.cif

This program is part of the mmcif-tools suite.
```

## Identifiers

### mmcif_label2auth

Copies the label identifiers of an mmCIF file onto the author ones.

```text
Copies the label identifiers of an mmCIF file onto the author ones.

By default copies both chains and residue numbers; an option copies only
one of them:

    -chain   auth_asym_id <- label_asym_id. Every molecule gets a chain of
             its own: in PDB entries, ligands and waters share the author
             chain of their protein, but each has its own label_asym_id.
    -res     auth_seq_id <- label_seq_id, for polymer residues. They are
             then numbered 1, 2, ... along the deposited sequence, so
             missing residues keep their place (unlike mmcif_reres), and
             insertion codes are removed. Ligands and waters have no
             label_seq_id and keep their numbers.

This is the mmCIF counterpart of pdb_chainxseg / pdb_segxchain and, with
-chain, of pdb_chainbows. Rows of _atom_site_anisotrop and _struct_conn
follow. Other categories that mention author chains or residue numbers
(_pdbx_poly_seq_scheme, _struct_conf, ...) are not updated. A file without
auth_asym_id or auth_seq_id columns (e.g. an integrative model) gets them;
residues without a label_seq_id are then numbered '?'. Nothing is
written if two residues would end up with the same chain and number, e.g.
with -res when a ligand of the chain already has one of the new numbers.

Usage:
    mmcif_label2auth [-chain | -res] <cif file>

Example:
    mmcif_label2auth 1abc.cif          # chains and residue numbers
    mmcif_label2auth -chain 1abc.cif   # one chain per molecule

This program is part of the mmcif-tools suite.
```

### mmcif_auth2label

Copies author chains onto the label chains of polymers (auth_asym_id -> label_asym_id) in an mmCIF file.

```text
Copies author chains onto the label chains of polymers (auth_asym_id ->
label_asym_id) in an mmCIF file.

For programs that read label_asym_id, e.g. after renaming chains with
mmcif_rplchain. Only chains are copied, and only for polymers (atoms with a
label_seq_id):

  - Ligands and waters keep their label chains. In mmCIF each of them is a
    molecule of its own, while in PDB entries they share the author chain
    of their protein, so copying would merge them.
  - Residue numbers are never copied: label_seq_id counts 1, 2, ... along
    the deposited sequence, which author numbers (with their offsets, gaps
    and insertion codes) do not. mmcif_label2auth -res goes the other way.

Nothing is written if an author chain holds more than one polymer, or if a
new label chain is already the label chain of another molecule. Rows of
_atom_site_anisotrop and _struct_conn follow; other categories that mention
label chains (_struct_asym, _pdbx_poly_seq_scheme, ...) are not updated, as
with mmcif_rplchain --label.

Usage:
    mmcif_auth2label <cif file>

Example:
    mmcif_rplchain -A:H 1abc.cif | mmcif_auth2label

This program is part of the mmcif-tools suite.
```

## Models, files and assemblies

### mmcif_splitchain

Splits an mmCIF file into several, each containing one chain.

```text
Splits an mmCIF file into several, each containing one chain.

As pdb_splitchain, writes <name>_<chain>.cif files in the current directory
(<name> is the input file name without its extension, or 'splitchains' when
reading stdin), overwriting existing ones, and prints nothing. Unlike
pdb_splitchain, each file is a complete mmCIF file: the metadata
categories are kept, the data block is renamed <block>_<chain>, and rows of
_atom_site_anisotrop and _struct_conn of other chains are removed.

Chains are auth_asym_id (label_asym_id with --label).

File names are unique ignoring case: on macOS and Windows, 1abc_A.cif and
1abc_a.cif would be the same file, so a chain whose file name differs from
an earlier one only in case gets a number (1abc_a~2.cif), with a note on
stderr. pdb_splitchain overwrites the file instead.

Usage:
    mmcif_splitchain <cif file>

Example:
    mmcif_splitchain 1abc.cif   # writes 1abc_A.cif, 1abc_B.cif, ...

This program is part of the mmcif-tools suite.
```

### mmcif_splitmodel

Splits an mmCIF file into several, each containing one model.

```text
Splits an mmCIF file into several, each containing one model.

As pdb_splitmodel, writes <name>_<model>.cif files in the current directory
(<name> is the input file name without its extension, or 'splitmodels'
when reading stdin), named by the model numbers as written, and prints
nothing. Unlike pdb_splitmodel, which writes nothing for a file without
MODEL records, a single-model file gives one file; and each file is a
complete mmCIF file: the metadata categories are kept, the data block is
renamed <block>_<model>, and rows of _atom_site_anisotrop and _struct_conn
of removed atoms are removed.

Usage:
    mmcif_splitmodel <cif file>

Example:
    mmcif_splitmodel 1abc.cif   # writes 1abc_1.cif, 1abc_2.cif, ...

This program is part of the mmcif-tools suite.
```

### mmcif_mkensemble

Merges several mmCIF files into one multi-model (ensemble) file.

```text
Merges several mmCIF files into one multi-model (ensemble) file.

As pdb_mkensemble, the n-th file becomes model n (all the atoms of its
first data block, whatever their model numbers), the provenance of each
model is written at the top (as comment lines), the other categories are
dropped except the _struct_conn of the first file, and the models are not
checked against each other (see mmcif_chkensemble).

The _atom_site loop has the columns of all files (missing values are
'?'), pdbx_PDB_model_num is added if needed, and atom ids are renumbered
from 1 over the whole ensemble, since mmCIF atom ids must be unique.

Usage:
    mmcif_mkensemble <cif file> <cif file> [...]

Example:
    mmcif_mkensemble 1abc.cif 1xyz.cif > ensemble.cif

This program is part of the mmcif-tools suite.
```

### mmcif_chkensemble

Checks that all models of a multi-model mmCIF file have the same composition.

```text
Checks that all models of a multi-model mmCIF file have the same
composition.

Composition is defined as same atoms/residues/chains: an atom is its chain,
residue name, number, insertion code, atom name and altloc (auth
identifiers, or label ones with --label). Unlike pdb_chkensemble, the atom
id is not part of it (mmCIF atom ids are unique across models, so every
ensemble would differ). Differences are written to stderr as by
pdb_chkensemble, one atom per line; the exit status is 1 if models differ.

Usage:
    mmcif_chkensemble <cif file>

Example:
    mmcif_chkensemble 1abc.cif

This program is part of the mmcif-tools suite.
```

### mmcif_merge

Merges the atoms of several mmCIF files into one.

```text
Merges the atoms of several mmCIF files into one.

pdb_merge concatenates files; concatenated mmCIF files would just be
several data blocks, so instead the output is the first file with the atoms
of all the others added to its _atom_site loop (which gets the columns of
all files; missing values are '?'). The atoms of the first file keep their
ids, those of the next files are renumbered after them. Only _atom_site is
taken from the other files (from their first data block, added to the
first data block of the first file): their metadata, _atom_site_anisotrop and
_struct_conn rows are not merged, and _entity / _struct_asym are not
updated.

Chains are not renamed: a warning is printed when two files use the same
chain (auth_asym_id, or label_asym_id with --label); see mmcif_rplchain.
Pipe the result through mmcif_tidy to align the columns.

Usage:
    mmcif_merge <cif file> <cif file> [...]

Example:
    mmcif_merge 1abc.cif 1xyz.cif > merged.cif

This program is part of the mmcif-tools suite.
```

### mmcif_assembly

Builds a biological assembly from the operations stored in an mmCIF file.

```text
Builds a biological assembly from the operations stored in an mmCIF file.

PDB entries describe their biological assemblies (_pdbx_struct_assembly) as
operations (_pdbx_struct_oper_list: a rotation matrix and a translation)
applied to label chains (_pdbx_struct_assembly_gen). This tool applies them
and writes the assembly's atoms in place of the asymmetric unit, offline and
for any file that carries these categories (mmcif_fetch -biounit downloads
the RCSB's prebuilt files instead). The option is the assembly id (default
1); mmcif_info lists the assemblies of a file.

Operation expressions follow the PDBx dictionary: 1, 1,2,5, (1-60), and
products such as (1-60)(61-88), whose operations are applied right to left.

Naming follows the RCSB's assembly files: the copy made by the identity
operation keeps its chain names, and every other copy gets the operation
added to its label and author chains (A-2; A-1x61 for a product). Atoms are
written chain by chain, each chain followed by its copies, in every model;
atom ids are renumbered from 1. Links (_struct_conn) are repeated in each
copy that holds both partners; links between copies or to symmetry mates
are not written. _atom_site_anisotrop is removed (its tensors would need
rotating). The other categories still describe the asymmetric unit.

Usage:
    mmcif_assembly [-<assembly id>] <cif file>

Example:
    mmcif_assembly 1abc.cif       # assembly 1
    mmcif_assembly -2 1abc.cif    # assembly 2

This program is part of the mmcif-tools suite.
```

## Categories and data blocks

### mmcif_cat

Prints one category of an mmCIF file as a table (tab-separated values).

```text
Prints one category of an mmCIF file as a table (tab-separated values).

The first line holds the column names (without the category prefix), then
each row of the category follows on a line of its own, so the output can go
to awk, sort, a spreadsheet or pandas. Columns can be chosen and ordered
after a dot (-_atom_site.id,Cartn_x). A category written as key-value pairs
(_cell.length_a 50.0, ...) is printed as one row.

Values are printed without their CIF quotes; nulls stay '.' and '?'. Tabs,
newlines and backslashes inside a value are written as \t, \n and \\.
In a file with several data blocks, the rows of every block are printed,
and the header is repeated only when the columns change. Category and
column names are case-insensitive, as in CIF.

Usage:
    mmcif_cat -<category>[.<column>,<column>...] <cif file>

Example:
    mmcif_cat -_atom_site 1abc.cif
    mmcif_cat -_atom_site.auth_asym_id,auth_seq_id,Cartn_x 1abc.cif
    mmcif_cat -_chem_comp components.cif

This program is part of the mmcif-tools suite.
```

### mmcif_delcat

Deletes whole categories from an mmCIF file.

```text
Deletes whole categories from an mmCIF file.

Categories are given by name, comma-separated; a name ending in '*' removes
every category that starts with it (quote it, so that the shell does not
expand it). Names are case-insensitive, and the leading '_' is optional.
Everything else is copied unchanged.

Nothing is checked or updated in the categories that remain: e.g. removing
_atom_site leaves _atom_site_anisotrop without its atoms.

Usage:
    mmcif_delcat -<category>[,<category>...] <cif file>

Example:
    mmcif_delcat -_atom_site_anisotrop 1abc.cif
    mmcif_delcat '-_pdbx_validate*' 1abc.cif      # PDB validation reports
    mmcif_delcat -_struct_conf,_struct_sheet_range 1abc.cif

This program is part of the mmcif-tools suite.
```

### mmcif_keepcat

Keeps only some categories of an mmCIF file.

```text
Keeps only some categories of an mmCIF file.

Categories are given by name, comma-separated; a name ending in '*' keeps
every category that starts with it (quote it, so that the shell does not
expand it). Names are case-insensitive, and the leading '_' is optional.
The data_ line of each block is kept; everything else (other categories,
comments) is removed. mmcif_keepcoord is the same as
mmcif_keepcat -_atom_site,_struct_conn.

Nothing is checked or added: e.g. keeping _atom_site_anisotrop without
_atom_site leaves it without its atoms.

Usage:
    mmcif_keepcat -<category>[,<category>...] <cif file>

Example:
    mmcif_keepcat -_atom_site 1abc.cif
    mmcif_keepcat '-_entity*,_struct_asym' 1abc.cif

This program is part of the mmcif-tools suite.
```

### mmcif_blocks

Lists or extracts the data blocks of a CIF file.

```text
Lists or extracts the data blocks of a CIF file.

Without an option, prints the name of each data block (data_<name>), one
per line. With block names, comma-separated, writes those blocks (in file
order) and nothing else. Names are case-insensitive, as in CIF; reading
stops once every requested block has been written, so extracting one
component from the whole Chemical Component Dictionary is quick.

Usage:
    mmcif_blocks [-<name>[,<name>...]] <cif file>

Example:
    mmcif_blocks components.cif              # lists the blocks
    mmcif_blocks -ATP,GTP components.cif     # writes the ATP and GTP blocks

This program is part of the mmcif-tools suite.
```

## Inspect / format

### mmcif_info

Summarizes the metadata of an mmCIF file: title, method, resolution, ...

```text
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
```

### mmcif_wc

Summarizes the contents of an mmCIF file, like the wc command in UNIX.

```text
Summarizes the contents of an mmCIF file, like the wc command in UNIX.

By default, this tool produces a general summary, but you can use several
options to produce focused but more detailed summaries:
    [m] - no. of models.
    [c] - no. of chains (plus per-model if multi-model file).
    [r] - no. of residues (plus per-model if multi-model file).
    [a] - no. of atoms (plus per-model if multi-model file).
    [h] - no. of HETATM (plus per-model if multi-model file).
    [o] - presence of disordered atoms (altloc).
    [i] - presence of insertion codes.

Counts follow pdb_wc: residues and atoms count ATOM records only, HETATM
records are counted separately, and the altloc and insertion-code flags
look at ATOM records. The lists printed with [a] and [h] give atom names
and HETATM residue names (pdb_wc prints padded atom names for both).

Usage:
    mmcif_wc [-<option>] <cif file>

Example:
    mmcif_wc 1abc.cif
    mmcif_wc -mc 1abc.cif

This program is part of the mmcif-tools suite.
```

### mmcif_tofasta

Extracts the residue sequence of an mmCIF file in FASTA format.

```text
Extracts the residue sequence of an mmCIF file in FASTA format.

Canonical amino acids and nucleotides (RNA and DNA) are represented by their
one-letter code while all others are represented by 'X'. As in pdb_tofasta,
the sequence is read from the atoms of the first model, ligands and waters
(HETATM) included, so pipe through mmcif_delhetatm to drop them.

The -multi option splits the different chains into different records in the
FASTA file. The -seqres option reads the full deposited sequence from
_pdbx_poly_seq_scheme instead, including residues without coordinates (from
_entity_poly_seq when a file has no _pdbx_poly_seq_scheme, with the chains
of _entity_poly.pdbx_strand_id, or of _struct_asym with --label). Where two
residue types share a position (microheterogeneity), the first is used. The
two can be combined: -multi,seqres.

Usage:
    mmcif_tofasta [-multi] [-seqres] <cif file>

Example:
    mmcif_tofasta 1abc.cif
    mmcif_tofasta -multi 1abc.cif
    mmcif_tofasta -multi,seqres 1abc.cif

This program is part of the mmcif-tools suite.
```

### mmcif_gap

Finds gaps between consecutive protein residues in an mmCIF file.

```text
Finds gaps between consecutive protein residues in an mmCIF file.

Detects gaps both by a distance criterion (consecutive alpha carbons more
than 4.0 A apart) and by discontinuous residue numbering, within each chain
and model. Only applies to protein residues (ATOM records), as pdb_gap. The
output lines are those of pdb_gap; unlike pdb_gap, 52 followed by 52A is
not a gap, and the second alpha carbon of a residue with alternate
locations is not compared with the first.

Chains and residue numbers are auth identifiers (label ones with --label).

Usage:
    mmcif_gap <cif file>

Example:
    mmcif_gap 1abc.cif

This program is part of the mmcif-tools suite.
```

### mmcif_sort

Sorts the atoms (_atom_site rows) of an mmCIF file.

```text
Sorts the atoms (_atom_site rows) of an mmCIF file.

As in pdb_sort, ATOM records come before HETATM records, and each group is
sorted by:
    chain, residue number, insertion code, altloc, atom id   (default, -CR)
    chain, altloc, atom id                                    (-C)
    chain order of appearance, residue number, insertion
    code, altloc, atom id                                     (-R)
so atoms keep their order within a residue. Chains sort as text. Unlike
pdb_sort, multi-model files are supported (each model is sorted on its
own, models stay in order), and nothing is removed. The rows of
_atom_site_anisotrop are put in the new order of their atoms.

Chains and residue numbers are auth identifiers (label ones with --label).

Usage:
    mmcif_sort [-<option>] <cif file>

Example:
    mmcif_sort 1abc.cif      # sorts by chain and residues
    mmcif_sort -C 1abc.cif   # sorts by chain only
    mmcif_sort -R 1abc.cif   # sorts by residue number/icode only

This program is part of the mmcif-tools suite.
```

### mmcif_tidy

Re-formats an mmCIF file in a uniform layout.

```text
Re-formats an mmCIF file in a uniform layout.

Every category is rewritten: a '#' line between categories, key-value
categories with their values aligned, loops with each column padded to its
widest value, values quoted only when needed (text fields on lines of
their own) and '\n' line endings. Comments between categories are kept
(the '#' separators are rewritten); comments inside a category are
dropped. The data itself is not changed. Use it after a pipeline of other
tools, which keep the original layout and only rewrite what they change.

Unlike pdb_tidy, which fixes PDB-format records (TER, END, line padding),
there is nothing to repair in mmCIF, so this tool only changes the layout.

Usage:
    mmcif_tidy <cif file>

Example:
    mmcif_tidy 1abc.cif

This program is part of the mmcif-tools suite.
```

### mmcif_validate

Validates the structure of an mmCIF file.

```text
Validates the structure of an mmCIF file.

Checks the CIF syntax and the coordinates (_atom_site); it does not check
the file against the PDBx/mmCIF dictionary (item types, allowed values,
mandatory items). Errors:

    - CIF syntax errors (with their line number)
    - no _atom_site category, or required _atom_site columns missing (id,
      type_symbol, Cartn_x/y/z, and an atom name, residue name, chain and
      residue number in at least one identifier scheme)
    - coordinates, occupancies or B-factors that are not numbers
    - group_PDB other than ATOM / HETATM
    - residue numbers (auth_seq_id, label_seq_id) or model numbers that are
      not integers ('.' and '?' are allowed)
    - the same atom id twice in a data block
    - rows of _atom_site_anisotrop or _struct_conn that match no atom

Warnings: only one identifier scheme (auth_* or label_*) present,
occupancies outside 0-1.

Prints one line per problem (the first 20 of each kind) and a summary.
The exit status is 1 if there are errors, 0 otherwise.

Usage:
    mmcif_validate <cif file>

Example:
    mmcif_validate 1abc.cif

This program is part of the mmcif-tools suite.
```
