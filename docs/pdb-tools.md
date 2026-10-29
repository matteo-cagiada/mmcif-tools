<!-- Generated from README.md by tools/make_docs.py: edit README.md and run python tools/make_docs.py. -->

# Coming from pdb-tools

Most tools have the pdb-tools name with `mmcif_` for `pdb_` and the same options. The differences:

- `pdb_tocif` / `pdb_fromcif` are `mmcif_frompdb` / `mmcif_topdb`.
- The segment tools (`pdb_seg`, `pdb_selseg`, `pdb_splitseg`) have no counterpart: mmCIF has no segment identifier,
  and `label_asym_id` names each molecule (`mmcif_selchain --label`, `mmcif_splitchain --label`).
- `pdb_chainbows`, `pdb_chainxseg` and `pdb_segxchain` have no counterpart: mmCIF has no TER records. `mmcif_frompdb`
  uses the TER records of a PDB file to separate molecules and names blank chains, and `mmcif_label2auth -chain`
  gives each molecule its own author chain.
- Where pdb-tools behaves differently (or has a bug), the tool's `--help` says so.
