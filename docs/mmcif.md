<!-- Generated from README.md by tools/make_docs.py: edit README.md and run python tools/make_docs.py. -->

# Things to know about mmCIF

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
