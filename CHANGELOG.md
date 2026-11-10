# Changelog

All notable changes to mmcif-tools. Versions follow [semantic versioning](https://semver.org/); before 1.0,
command names and options may still change.

## 0.2.0 (2026-11-10): metadata follows the edits

- **The metadata follows the edits.** Every tool keeps the tables that name chains, residues or atoms in step with
  the atoms, not only `_atom_site_anisotrop` and `_struct_conn`: the sequence schemes, the secondary structure,
  the sites, the validation reports, the assemblies, model quality (`_ma_qa_metric_local`), the SIFTS mappings of
  PDBe's updated files, the NCS domains, and more (50 tables).
  - **Renames** (`mmcif_chain`, `mmcif_rplchain`, `mmcif_reres`, `mmcif_shiftres`, `mmcif_fixinsert`,
    `mmcif_label2auth`, `mmcif_auth2label`, `mmcif_rplresname`, `mmcif_uniqname`) reach every table. After
    `mmcif_reres`, residues without atoms get `?` as their author number: their old number may now be another
    residue's. `mmcif_label2auth` gives them their label ids. `mmcif_rplresname` also renames the residue in the
    sequences (`_entity_poly_seq`) and its definition (`_chem_comp`, ...), or removes the old definition when the
    new name already has one.
  - **Removals** (the `sel` and `del` tools, `mmcif_head`, `mmcif_selaltloc`, `mmcif_intersect`, and the split
    tools): rows naming removed atoms, residues, chains or molecules are removed, the sequences included; helices,
    strands and alignments that lose their ends are shortened; lists of chains are pruned; sheets, sites,
    assemblies, sequence alignments and coordination spheres left without parts are removed; and
    `_entity.pdbx_number_of_molecules` and the one-letter sequences are updated (lines of up to 80 characters,
    broken between residues). Rows naming a removed model (`mmcif_selmodel`, `mmcif_splitmodel`) or conformer
    (`mmcif_selaltloc`) go too, and those naming a kept conformer lose its label, as the atoms do. Residues that had no atoms in the input stay, and so does a
    position of `_entity_poly_seq` that a kept chain lists without atoms. Atoms listed as never observed
    (`_pdbx_unobs_or_zero_occ_atoms`) go with their residue. Sequence alignments that name no chain
    (`_struct_ref_seq` in integrative models) follow the molecule's sequence. A molecule whose list of chains
    (`_entity_poly.pdbx_strand_id`) misses some of its copies keeps its row when only those are kept, and the
    list then names them. A residue deleted from the construct (`_struct_ref_seq_dif`, with no number) names
    only its chain, not the chain's waters. A row of a model (validation reports, cis peptides, base pairs)
    follows the residues of that model, as `mmcif_delres` and `mmcif_selres` with a step remove other residues
    in each model. Rows that name residues by their author ids only, with `?` label ids (TLS groups), follow
    and are shortened by those ids. A modification (`_pdbx_modification_feature`) goes with the atom of its bond
    (a disulfide's SG, a sugar's C1).
  - **Cost:** one more reading of the file. The tools that remove atoms take about twice as long (4HHB: 0.06 to
    0.15 s; a 238,000-atom ribosome, 4V6X: 1.6 to 2.3-3.2 s, with up to 50 MB more memory); renames take 10-25%
    longer.
  - **Checked** with the PDBx/mmCIF dictionary (version 5.418, with gemmi and its parent check): every tool's
    output on the 12 test entries has no error its input does not have, except when no polymer is left, where
    gemmi reports the empty `_entity_poly_seq`, and `mmcif_keepcoord`, which keeps only the coordinates. The
    stress test (300 entries and variants, about 60 commands each) was run again with the check; what it found
    is fixed. Files with two data blocks: only the first one's tables follow. Also checked: that no row is
    removed without a reason (each row of the input is tagged, and each one removed must name something with no
    atoms left), and tools run one after another: in either order, twice, a rename and its inverse, and
    through shell pipes.
- **`mmcif_validate`** warns about metadata out of step with the atoms, including rows of models with no atoms.
- **`mmcif_element`** adds a row to `_atom_type` for an element new to the file.
- **`mmcif_label2auth`**: an unobserved residue whose author number is also a water's or ligand's in the same
  author chain got that molecule's chain in `_pdbx_poly_seq_scheme` (8PJO); it now gets its own label ids.
- **`mmcif_uniqname`** names the conformers of an atom apart (S1, S2), as `pdb_uniqname` does; a row of another
  table naming one conformer now gets that conformer's new name.
- **`mmcif_delmeta`** is no longer needed to fix stale metadata: it deletes what prediction and simulation programs
  do not read. It now also keeps the tables the kept ones point to: `_struct_conn_type`, `_pdbx_nonpoly_scheme`
  and, for EM entries, `_em_image_processing` to `_em_experiment`. The README's "Metadata goes stale" is now
  "Metadata follows the edits".
- **DOI:** mmcif-tools is archived on Zenodo, [10.5281/zenodo.23966235](https://doi.org/10.5281/zenodo.23966235)
  (all versions; 0.1.2 is [10.5281/zenodo.23966236](https://doi.org/10.5281/zenodo.23966236)). The README has a
  DOI badge and explains how to cite the software, with a BibTeX entry, and CITATION.cff gives the DOI.

## 0.1.2 (2026-11-05): website, and archived on Zenodo

No changes to the tools.

- **One website**, <https://matteo-cagiada.github.io/mmcif-tools/>: a main page, the documentation under
  [docs/](https://matteo-cagiada.github.io/mmcif-tools/docs/) (the pages' old addresses forward there), and
  **mmcif-tools in the browser** under [app/](https://matteo-cagiada.github.io/mmcif-tools/app/): it runs the
  tools on your files with Pyodide (Python in WebAssembly), with nothing to install, and the files never leave
  your computer. It installs the latest release from PyPI. Its code is in `web/`.
- The README links to the website under its title.
- **Things to know about mmCIF**, rewritten with examples on 4HHB: mmCIF and PDB files side by side, the two sets
  of identifiers (`auth_*` and `label_*`), what renumbering changes, how to read an output, and why the metadata
  goes stale. On the website it comes right after Installation.
- **In the browser:** pipelines end with `mmcif_delmeta` unless you untick it, and the page links to the
  documentation where it helps (PDB and mmCIF files, the two sets of identifiers, reading the output).

## 0.1.1 (2026-10-29): on PyPI

No changes to the tools.

- **On PyPI:** `python -m pip install mmcif-tools`. Releases are built, checked and uploaded by
  `.github/workflows/release.yml`, with PyPI's Trusted Publishing.
- **Documentation website:** <https://matteo-cagiada.github.io/mmcif-tools/>, made from the README, the help texts
  and this changelog.
- The README installs from PyPI first, and its links are full URLs, so they also work on PyPI.

## 0.1.0 (2026-10-28): first test release

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
- **`mmcif_delmeta`** deletes the metadata that edits make stale (validation reports, secondary structure, sites,
  assemblies, ...) and keeps the atoms, links, entry, method and resolution, release dates, unit cell, entities,
  chains (`_struct_asym`), sequences, ligands (`_pdbx_entity_nonpoly`), modified residues and chemical components
  with their atoms and bonds: what structure prediction programs (AlphaFold 2 and 3, ColabFold, OpenFold3, Boltz,
  Chai-1, Protenix, BioEmu) and MD preparation programs (OpenMM, PDBFixer) read. The kept tables follow the atoms:
  chains, entities and modified residues with no atoms left are removed, modified residues get their current
  author numbers, the number of molecules and the chains of each entity are recomputed, and `_pdbx_poly_seq_scheme`
  is rebuilt (OpenFold3 needs it). So it is the last step of a pipeline. `-all` keeps only the atoms and links.
  (A first version deleted `_struct_asym`, and ColabFold rejected its files as templates.)
- **Checking an installation:** `python tools/smoke_test.py` runs every command and a few pipelines. The full test
  suite stays with the developers and is not part of the public repository.
- **Known limitation:** after an edit, categories other than `_atom_site`, `_atom_site_anisotrop` and
  `_struct_conn` that refer to atoms or residues (`_struct_conf`, `_pdbx_poly_seq_scheme`, `_pdbx_validate_*`,
  ...) are left as they were. A tool to bring them back in step is planned.
