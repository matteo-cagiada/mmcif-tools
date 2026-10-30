# The website: main page and mmcif-tools in the browser

The website <https://mmcif-tools.dev/> has three parts, which
[`.github/workflows/docs.yml`](../.github/workflows/docs.yml) builds and publishes together:

- `/`: the main page, `index.html` here;
- `/docs/`: the documentation, built by Zensical (`zensical.toml`, pages written by `tools/make_docs.py`);
- `/app/`: **mmcif-tools in the browser**, `app/` here. It runs mmcif-tools on mmCIF files without installing
  anything: open a `.cif` or `.cif.gz` file, chain tools into a pipeline
  (`mmcif_selchain -A | mmcif_delhetatm | mmcif_tidy`), look at the output and download it.

The rest of this file is about the in-browser page.

Everything runs in the browser. The page loads Python with [Pyodide](https://pyodide.org)
(CPython compiled to WebAssembly) and installs mmcif-tools from
[PyPI](https://pypi.org/project/mmcif-tools/) into it. The files you open are read on your
computer and never leave it: there is no server, and the only downloads are Pyodide (from
cdn.jsdelivr.net) and mmcif-tools (from PyPI).

What each tool does is in the [mmcif-tools documentation](https://mmcif-tools.dev/docs/).
For big files, many files or scripts, install mmcif-tools instead: `pip install mmcif-tools`.

## Run it locally

```bash
python -m http.server -d web
```

from the top of the repository, then open <http://localhost:8000/app/> (it installs the latest
release from PyPI). Opening `app/index.html` as a `file://` page does not work: the Web Worker and
the downloads of Python need http. `http.server` sends no cache headers, so after changing a file
reload the page with the cache bypassed (a change of the `#pipeline` part of the URL alone does not
reload it).

## A new mmcif-tools version

Nothing to change: after each release on PyPI (and on each push to `main`), the workflow takes the
latest release, runs a short pipeline with `runner.py` on it, and publishes the page with that
release pinned (`MMCIF_TOOLS_VERSION` in `worker.js`, which is `''`, the latest, in the
repository). GitHub Pages lets browsers cache the files for 10 minutes.

A new tool goes into `GROUPS` in `app/runner.py` as in `tools/make_docs.py` (a test checks that
they agree; a tool missing there is listed under "Other", and the workflow prints a warning). If it
is not "one mmCIF file in, one mmCIF file out", add it to `MULTI_INPUT`, `SPLIT` or `NOT_CIF` there
too. Tools that the installed release does not have yet are left out of the menu.

Pyodide is pinned as `PYODIDE_VERSION` in `app/worker.js`.

## Files

- `index.html`: the main page of the website.
- `app/index.html`: the in-browser page (layout and styles).
- `app/app.js`: the opened files, the pipeline editor (which ends with `mmcif_delmeta` unless it is
  unticked), the result and the download. The pipeline is
  kept in the URL (`#selchain -A | delhetatm | tidy`), so a link can open a ready pipeline.
- `app/worker.js`: loads Pyodide in a module Web Worker, installs mmcif-tools with micropip, and runs
  the pipelines. Running in a worker keeps the page responsive, and Stop simply ends the worker
  and starts a new one.
- `app/runner.py`: runs each tool through its own `main()`, as on the command line (`sys.argv` holds
  the arguments, stdout goes to a file, stderr is captured), and chains the steps through files.
  It also runs on CPython with mmcif-tools installed, which is how CI checks it and how the times
  below were compared: `python runner.py '<tool>[:<option>] | ...' <file>`.

## How it works

- The opened files are mounted read-only in Pyodide's file system (Emscripten's WORKERFS): Python
  reads them straight from the `File` objects, without loading them into memory or copying them.
  WORKERFS makes one call to the browser per read, and Python's text files read 8 kB at a time,
  which made reading a plain `.cif` file about 15 times slower than from memory; `worker.js`
  reads 4 MB blocks instead. `.cif.gz` files are read through `gzip`, which reads larger blocks
  anyway.
- Each step writes its output to Pyodide's in-memory file system and the next step reads it; the
  previous step's file is deleted as soon as the next one has run. The last output is handed to
  the page, which shows its first 200 kB and offers it as a download.
- The tool list, the option of each tool and the help texts come from the tools' docstrings, so
  the page follows the installed mmcif-tools version.

Not every tool fits "one mmCIF file in, one file out":

- `mmcif_fetch` is left out: it downloads files with `urllib`, which does not work in Pyodide.
  The page says so and points to the RCSB and the AlphaFold DB.
- `mmcif_splitchain` and `mmcif_splitmodel` write several files: they must be the last step, and
  the files are downloaded as a zip.
- `mmcif_merge`, `mmcif_mkensemble` and `mmcif_intersect` read every opened file (in the order of
  the list) and `mmcif_frompdb` reads a PDB file: they must be the first step.
- Tools that do not write mmCIF (`info`, `wc`, `gap`, `validate`, `chkensemble`, `cat`, `tofasta`,
  `topdb`) must be the last step. `validate` and `chkensemble` exit with status 1 when they find
  problems; their report is still shown.

## Speed and limits

Measured with mmcif-tools 0.1.0 on an Apple M3 Pro with 18 GB of memory: CPython 3.13 against
Pyodide 314.0.7 (Python 3.14) in Chromium, the tools alone (`wc`, `selchain -A`, `delhetatm`,
`b -10.0`, `tidy`, `sort` and `validate`, each on its own). The outputs were byte-identical.

| Input | Atoms | Uncompressed | CPython | Pyodide | Pyodide / CPython |
|---|---:|---:|---:|---:|---:|
| 4hhb.cif.gz | 4,779 | 0.8 MB | 0.26 s | 0.38 s | 1.46 |
| 6vxx.cif.gz | 23,694 | 2.6 MB | 0.93 s | 1.35 s | 1.45 |
| 1jj2.cif.gz | 98,543 | 11.9 MB | 5.03 s | 7.23 s | 1.44 |
| 4v6x.cif.gz | 237,685 | 28.8 MB | 11.7 s | 16.7 s | 1.43 |
| synthetic `.cif` (`wc`, `selchain -A`, `sort` only) | 1M | 100 MB | 7.6 s | 11.2 s | 1.48 |
| synthetic `.cif` (`wc`, `selchain -A`, `sort` only) | 2.5M | 252 MB | 19.6 s | 32.3 s | 1.65 |

- **Size.** Files of up to 2.5 million atoms (250 MB uncompressed) work. Beyond that, runs that
  keep state per atom slow down much more than on CPython (at 5 million atoms `mmcif_selchain -A`
  took 81 s against 9 s), and `mmcif_sort` and `mmcif_assembly`, which hold all the atoms in
  memory, can run out of the 4 GB a WebAssembly program can have: the page then says so and starts
  Python again.
- **Memory.** The input files are not loaded into memory, but each step's output is, in Pyodide's
  file system and then in the page (for the download), so a pipeline needs a few times the size of
  its output.
- **Loading.** The first visit downloads about 6 MB (compressed), which the browser then caches.
  Starting Python and installing mmcif-tools takes 5-8 s on each visit; after Stop, Python
  restarts in about 1.3 s.
- **Browsers.** The page needs module Web Workers (Chrome and Edge 80+, Safari 15+, Firefox 114+).
  It has been tried in Chromium. There a classic worker could not load Pyodide with
  `importScripts()`, which is why `worker.js` is a module worker and loads Pyodide with `import()`.
- **No interruption inside a tool.** Stop ends the worker (and Python with it): Pyodide's
  interrupt buffer would need a `SharedArrayBuffer`, which needs cross-origin isolation headers
  that GitHub Pages does not send.
