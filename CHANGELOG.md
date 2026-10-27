# Changelog

All notable changes to mmcif-tools. Versions follow [semantic versioning](https://semver.org/); before 1.0,
command names and options may still change.

## 0.1.0 (2026-10-24): first test release

The first version for testers. It is not on PyPI or conda-forge yet; install it from the repository.

- **53 command-line tools** (`mmcif_*`), in seven groups: fetch and convert, select and delete, edit,
  identifiers, models/files/assemblies, categories and data blocks, inspect and format. See
  [docs/tools.md](docs/tools.md).
- **mmCIF-specific tools:**
  - `mmcif_selentity` / `mmcif_delentity`: select by entity number or type (polymer, non-polymer, branched,
    macrolide, water). Modified residues stay with their chain.
  - `mmcif_assembly`: build a biological assembly from `_pdbx_struct_oper_list` and
    `_pdbx_struct_assembly_gen`, including product expressions such as `(1-60)(61-88)`. The chain names follow
    the RCSB's assembly files.
  - `mmcif_info`: title, method, resolution, R-factors, dates, organisms, model quality, entities, assemblies.
  - `mmcif_label2auth` / `mmcif_auth2label`: copy one identifier scheme onto the other.
  - `mmcif_cat`, `mmcif_delcat`, `mmcif_keepcat`, `mmcif_blocks`: work on categories and data blocks.
- **Tools inspired by pdb-tools**, with the same options where they apply. Author identifiers by default,
  `--label` for the PDB's label identifiers.
- **A streaming mmCIF reader and writer** (standard library only). Rows a tool does not change are copied
  byte-for-byte. `_atom_site_anisotrop` and `_struct_conn` follow edits to the atoms.
- **Tested** on hand-written fixtures and 12 real entries (X-ray, NMR, cryo-EM, ribosomes, an AlphaFold DB model),
  with a cross-check against pdb-tools. Every tool's output on these entries was also checked against the
  PDBx/mmCIF (or ModelCIF) dictionary with gemmi, and has no error that its input does not have
  (outside the test suite).
- **Python 3.9 to 3.15**, on Linux, macOS and Windows.
- **Fixed after a 200-system stress test** (PDB entries of every method and size, AlphaFold DB models, PDBe
  files, an integrative model):
  - `mmcif_splitchain`: chains whose ids differ only in case (`A` and `a`) overwrote each other's files on macOS
    and Windows; such files now get a number (`1abc_a~2.cif`). It also writes up to 100 chains per reading of
    the input instead of one, with the same output: 1792 chains no longer take hours.
  - `mmcif_fixinsert` / `mmcif_delinsertion`: a file without an insertion code column is written unchanged
    (it was an error).
  - `mmcif_label2auth`: adds `auth_asym_id` / `auth_seq_id` when a file lacks them (integrative models).
  - `mmcif_tofasta -seqres`: reads `_entity_poly_seq` when a file has no `_pdbx_poly_seq_scheme`, writes one
    residue per position where two residue types share it, and explains when a file has no polymer.
- **Fixed after a second, targeted stress test** (100 systems: integrative models, ModelArchive and PDBe files,
  chemical component files, entries chosen for one weak spot each, and modified copies of entries):
  - `mmcif_reres` / `mmcif_fixinsert`: in files whose models hold different residues (4PTJ), the anisotropic
    rows of later models got model 1's numbers; they now follow their own atom.
  - `mmcif_merge` (and every tool that adds rows to a loop): with upper-case tags (`_ATOM_SITE.ID`), atom ids
    were not renumbered; tags are case-insensitive, as in CIF.
  - Edits with `--label` on a file without label columns (or without `--label` on one without author
    columns) changed the other scheme's columns; they now stop with a clear message.
  - `mmcif_fixinsert`: links whose rows have no residue name (as gemmi writes them) were not renumbered.
  - `mmcif_head` counted the atoms of each data block separately; it now counts over the whole file, as its
    help says.
  - `mmcif_selaltloc` writes a file without an altloc column unchanged; `mmcif_b` / `mmcif_occ` add their
    column when it is missing; "missing column" errors name the mmCIF column (`type_symbol`, `group_PDB`, ...).
  - `mmcif_merge`, `mmcif_mkensemble` and `mmcif_topdb` say in their help that they use the first data block.
- **`mmcif_delmeta`** deletes the metadata that edits make stale (validation reports, secondary structure, sequence
  schemes, sites, assemblies, ...) and keeps the atoms, links, entry, method and resolution, release dates, unit
  cell, entities, chains (`_struct_asym`), sequences and chemical components: what AlphaFold and ColabFold read
  from template files. Rows of chains and entities with no atoms left are removed, and the number of molecules and
  the chains of each entity are recomputed, so it is the last step of a pipeline. `-all` keeps only the atoms and
  links. (A first version deleted `_struct_asym`, and ColabFold rejected its files as templates.)
- **Checking an installation:** `python tools/smoke_test.py` runs every command and a few pipelines. The full test
  suite stays with the developers and is not part of the public repository.
- **Known limitation:** after an edit, categories other than `_atom_site`, `_atom_site_anisotrop` and
  `_struct_conn` that refer to atoms or residues (`_struct_conf`, `_pdbx_poly_seq_scheme`, `_pdbx_validate_*`,
  ...) are left as they were. A tool to bring them back in step is planned.
