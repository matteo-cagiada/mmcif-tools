<!-- Generated from README.md by tools/make_docs.py: edit README.md and run python tools/make_docs.py. -->

# Citation and license

## Citation and changes

If you use mmcif-tools in your research, please cite the software with its Zenodo DOI,
[10.5281/zenodo.23966235](https://doi.org/10.5281/zenodo.23966235). That DOI covers all versions and leads to the
latest one; each version also has a DOI of its own, listed on the Zenodo record. Please cite the version you used
whenever possible: it helps support the development of the software. You can get a formatted citation from the
"Cite this repository" button on [GitHub](https://github.com/matteo-cagiada/mmcif-tools) (made from
[CITATION.cff](https://github.com/matteo-cagiada/mmcif-tools/blob/main/CITATION.cff)) or from the
[Zenodo record](https://doi.org/10.5281/zenodo.23966235).

In BibTeX, for version 0.2.0:

```bibtex
@software{cagiada2026mmciftools,
  author    = {Cagiada, Matteo},
  title     = {{mmcif-tools: command-line tools for mmCIF files}},
  version   = {0.2.0},
  publisher = {Zenodo},
  year      = {2026},
  doi       = {10.5281/zenodo.24104355},
  url       = {https://doi.org/10.5281/zenodo.24104355}
}
```

For another version, use its number and the DOI on its Zenodo page (`python -m pip show mmcif-tools` prints the
version you have installed).

The changes between versions are in the [changelog](changelog.md).
[PLAN.md](https://github.com/matteo-cagiada/mmcif-tools/blob/main/PLAN.md) records the design decisions and the test
results.

## Acknowledgements

Most of the code of mmcif-tools has been written with [Claude Code](https://claude.com/claude-code),
Anthropic's AI coding assistant. The author designed the project, tools and tests to guarantee the quality of the tool.
We thank the Bonvin Lab for developing and releasing pdb-tools, a fundamental tool for the authors' academic work.

## License

Apache-2.0. See [LICENSE](https://github.com/matteo-cagiada/mmcif-tools/blob/main/LICENSE).
