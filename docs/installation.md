<!-- Generated from README.md by tools/make_docs.py: edit README.md and run python tools/make_docs.py. -->

# Installation

```bash
git clone https://github.com/matteo-cagiada/mmcif-tools.git
cd mmcif-tools
python -m pip install .
```

There are no runtime dependencies. `python tools/smoke_test.py` then checks the installation: every command runs,
and a few pipelines give the expected results. For tab completion of the `mmcif_*` commands, open a new terminal
after installing (or run `rehash` in an open zsh terminal).

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

- **`tools/make_docs.py`** writes the [tool reference](tools.md) from the `--help` texts.
- **`tools/smoke_test.py`** checks an installation (see [Installation](installation.md)).
