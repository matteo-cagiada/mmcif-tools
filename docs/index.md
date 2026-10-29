<!-- Generated from README.md by tools/make_docs.py: edit README.md and run python tools/make_docs.py. -->

# mmcif-tools

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

**CURRENT VERSION: Version 0.1.0, pre-release: 53 available tools**

## Tools

mmcif-tools reads the mmCIF text format, as written by the PDB (wwPDB, RCSB, PDBe), the AlphaFold DB (ModelCIF),
structure prediction and modelling programs and viewers, plain or gzip-compressed. BinaryCIF is currently not
supported.

The full `--help` text of every tool is in the [tool reference](tools.md).

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
