# Changelog

All notable changes to mmcif-tools. Versions follow [semantic versioning](https://semver.org/); before 1.0,
command names and options may still change.

## 0.1.0 (2026-10-24): first test release

The first version for testers. It is not on PyPI or conda-forge yet; install it from the repository.

- **52 command-line tools** (`mmcif_*`), in seven groups: fetch and convert, select and delete, edit,
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
- **Known limitation:** after an edit, categories other than `_atom_site`, `_atom_site_anisotrop` and
  `_struct_conn` that refer to atoms or residues (`_struct_conf`, `_pdbx_poly_seq_scheme`, `_pdbx_validate_*`,
  ...) are left as they were. A tool to bring them back in step is planned.
