// Runs Python (Pyodide) and mmcif-tools in a Web Worker, so that the page
// stays responsive while a tool runs and a run can be stopped by ending the
// worker. The opened files are mounted read-only (WORKERFS), so they are read
// from disk as the tools stream through them and never copied or uploaded.
// It is a module worker (app.js starts it with type: 'module'): Pyodide is
// loaded with import(), which also works where importScripts() of a
// cross-origin script is blocked.

const PYODIDE_VERSION = '314.0.7';
const PYODIDE_URL = `https://cdn.jsdelivr.net/pyodide/v${PYODIDE_VERSION}/full/`;
// The mmcif-tools release installed from PyPI; '' installs the latest one. The
// website workflow (.github/workflows/docs.yml) sets the release it checked.
const MMCIF_TOOLS_VERSION = '';
const INPUT_DIR = '/work/in';
const RUN_DIR = '/work/run';
const READ_AHEAD = 4 << 20; // bytes read from an opened file at a time

let pyodide = null;
let runner = null;

async function init() {
  const start = performance.now();
  const { loadPyodide } = await import(PYODIDE_URL + 'pyodide.mjs');
  pyodide = await loadPyodide({ indexURL: PYODIDE_URL });
  await pyodide.loadPackage('micropip');
  const requirement = MMCIF_TOOLS_VERSION ? `mmcif-tools==${MMCIF_TOOLS_VERSION}` : 'mmcif-tools';
  try {
    await pyodide.pyimport('micropip').install(requirement);
  } catch (err) {
    // A Python error's message is its traceback: keep the last line.
    const reason = String(err && err.message ? err.message : err).trim().split('\n').pop();
    throw new Error(`Could not install ${requirement} from PyPI: ${reason}`);
  }

  const source = await fetch(new URL('runner.py', self.location.href));
  if (!source.ok) throw new Error('runner.py was not found next to the page.');
  pyodide.FS.mkdirTree('/work/lib');
  pyodide.FS.writeFile('/work/lib/runner.py', await source.text());
  pyodide.runPython("import sys; sys.path.insert(0, '/work/lib')");
  runner = pyodide.pyimport('runner');
  pyodide.FS.mkdirTree(INPUT_DIR);
  if (pyodide.FS.filesystems.WORKERFS) addReadAhead(pyodide.FS.filesystems.WORKERFS);

  self.postMessage({
    type: 'ready',
    tools: JSON.parse(runner.tools()),
    versions: { ...JSON.parse(runner.versions()), pyodide: pyodide.version },
    seconds: (performance.now() - start) / 1000,
  });
}

// WORKERFS makes one FileReaderSync call per read(), and Python's text files
// read 8 kB at a time (TextIOWrapper reads through read1(), which skips the
// buffer), so reading a file was ~15 times slower than from memory. Read
// READ_AHEAD bytes at a time instead, and serve the reads from the last block.
function addReadAhead(WORKERFS) {
  const reader = new FileReaderSync();
  WORKERFS.stream_ops.read = (stream, buffer, offset, length, position) => {
    const node = stream.node;
    if (position >= node.size) return 0;
    let block = node.readAhead;
    if (!block || position < block.start || position >= block.start + block.bytes.length) {
      const end = Math.min(node.size, position + Math.max(length, READ_AHEAD));
      const bytes = new Uint8Array(reader.readAsArrayBuffer(node.contents.slice(position, end)));
      block = node.readAhead = { start: position, bytes };
    }
    const start = position - block.start;
    const n = Math.min(length, block.bytes.length - start);
    buffer.set(block.bytes.subarray(start, start + n), offset);
    return n;
  };
}

// Make the opened files (a list of {name, file}) appear in INPUT_DIR.
function mountInputs(files) {
  const FS = pyodide.FS;
  try {
    FS.unmount(INPUT_DIR);
  } catch (err) {
    // nothing mounted yet
  }
  if (FS.filesystems.WORKERFS) {
    FS.mount(FS.filesystems.WORKERFS,
             { blobs: files.map((f) => ({ name: f.name, data: f.file })) }, INPUT_DIR);
    return;
  }
  // Without WORKERFS, copy the files into memory.
  for (const name of FS.readdir(INPUT_DIR)) {
    if (name !== '.' && name !== '..') FS.unlink(`${INPUT_DIR}/${name}`);
  }
  const reader = new FileReaderSync();
  for (const f of files) {
    FS.writeFile(`${INPUT_DIR}/${f.name}`, new Uint8Array(reader.readAsArrayBuffer(f.file)));
  }
}

function run(msg) {
  mountInputs(msg.files);
  const spec = {
    inputs: msg.inputs.map((name) => `${INPUT_DIR}/${name}`),
    steps: msg.steps,
    workdir: RUN_DIR,
  };
  const progress = (i) => self.postMessage({ type: 'progress', step: i });
  const result = JSON.parse(runner.run_pipeline(JSON.stringify(spec), progress));
  const output = result.output ? pyodide.FS.readFile(result.output) : null;
  // The page keeps the output; free the memory of the files written here.
  pyodide.runPython(`import shutil; shutil.rmtree(${JSON.stringify(RUN_DIR)}, ignore_errors=True)`);
  self.postMessage({ type: 'result', result, output }, output ? [output.buffer] : []);
}

self.onmessage = (event) => {
  if (event.data.type !== 'run') return;
  try {
    run(event.data);
  } catch (err) {
    // Python errors are caught by the runner, so this is the runtime failing (out of
    // memory, most likely); the page then starts a new worker.
    self.postMessage({ type: 'fatal', message: String(err && err.message ? err.message : err) });
  }
};

init().catch((err) => {
  self.postMessage({ type: 'fatal', message: String(err && err.message ? err.message : err), init: true });
});
