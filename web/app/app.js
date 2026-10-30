// The page: opened files, the pipeline of tools and the result. Python runs
// in worker.js; files are handed to it as File objects and never leave the
// browser.
'use strict';

const PREVIEW_BYTES = 200000; // of the output shown on the page
const DEFAULT_PIPELINE = 'info';

const $ = (id) => document.getElementById(id);

const state = {
  worker: null,
  ready: false,
  running: false,
  tools: [],
  byName: new Map(),
  files: [], // {name, file}
  input: 0, // index of the input file, for tools that read one file
  steps: [], // {tool, option, label, help}
  finish: { on: true, option: '', help: false }, // mmcif_delmeta after the steps
  sent: [], // the steps of the current run
  startedAt: 0,
  timer: null,
  downloadUrl: null,
};

// ---------------------------------------------------------------- formatting

function formatBytes(n) {
  if (n < 1000) return `${n} B`;
  if (n < 1e6) return `${(n / 1e3).toFixed(1)} kB`;
  if (n < 1e9) return `${(n / 1e6).toFixed(1)} MB`;
  return `${(n / 1e9).toFixed(2)} GB`;
}

function formatSeconds(s) {
  if (s < 1) return `${Math.round(s * 1000)} ms`;
  if (s < 60) return `${s.toFixed(2)} s`;
  return `${Math.floor(s / 60)} min ${Math.round(s % 60)} s`;
}

function shellQuote(word) {
  return /^[\w.,:+\-\/=@%^]+$/.test(word) ? word : `'${word.replace(/'/g, "'\\''")}'`;
}

function stem(name) {
  let base = name.replace(/\.gz$/i, '');
  const match = base.match(/\.(cif|mmcif|pdb|ent)$/i);
  if (match) base = base.slice(0, -match[0].length);
  return base;
}

function el(tag, props = {}, ...children) {
  const node = document.createElement(tag);
  Object.assign(node, props);
  node.append(...children);
  return node;
}

function setStatus(text, { busy = false, error = false } = {}) {
  $('status').classList.toggle('error', error);
  $('status').querySelector('.spinner').hidden = !busy;
  $('status-text').textContent = text;
}

// -------------------------------------------------------------------- worker

function startWorker() {
  state.ready = false;
  state.worker = new Worker('worker.js', { type: 'module' });
  state.worker.onmessage = (event) => onMessage(event.data);
  state.worker.onerror = (event) => {
    event.preventDefault();
    onMessage({ type: 'fatal', init: !state.ready, message: event.message || 'the worker failed' });
  };
  update();
}

function onMessage(msg) {
  if (msg.type === 'ready') {
    state.ready = true;
    if (!state.tools.length) {
      state.tools = msg.tools;
      state.byName = new Map(msg.tools.map((t) => [t.name, t]));
      try {
        loadPipeline(decodeURIComponent(location.hash.slice(1)));
      } catch (err) {
        state.steps = []; // a malformed #pipeline
      }
      if (!state.steps.length) loadPipeline(DEFAULT_PIPELINE);
      setupFinish();
      renderSteps();
    }
    const v = msg.versions;
    setStatus(`Ready (loaded in ${formatSeconds(msg.seconds)}).`);
    $('versions').textContent = `mmcif-tools ${v.mmciftools}, Python ${v.python}, Pyodide ${v.pyodide}.`;
    update();
  } else if (msg.type === 'progress') {
    state.progress = `Running step ${msg.step + 1} of ${state.sent.length}: mmcif_${state.sent[msg.step].tool}`;
    tick();
  } else if (msg.type === 'result') {
    finishRun();
    showResult(msg.result, msg.output);
    setStatus(msg.result.error ? 'Ready.' : `Done in ${formatSeconds(msg.result.seconds)}.`);
  } else if (msg.type === 'fatal') {
    const wasRunning = state.running;
    finishRun();
    state.worker.terminate();
    if (msg.init) {
      setStatus(`Could not load Python: ${msg.message}`, { error: true });
      return;
    }
    if (wasRunning) {
      showError(`Python stopped with an error, most likely because the browser ran out of memory ` +
                `for this file: ${msg.message}. Python is being restarted.`);
    }
    setStatus('Restarting Python…', { busy: true });
    startWorker();
  }
}

// --------------------------------------------------------------------- files

function addFiles(list) {
  for (const file of list) {
    let name = file.name;
    for (let n = 2; state.files.some((f) => f.name === name); n++) {
      name = file.name.replace(/^([^.]*)/, `$1~${n}`);
    }
    state.files.push({ name, file });
  }
  if (state.files.length && state.input >= state.files.length) state.input = 0;
  renderFiles();
  update();
}

function renderFiles() {
  const multi = firstTool()?.multi_input;
  $('files').replaceChildren(...state.files.map((f, i) => {
    const radio = el('input', {
      type: 'radio', name: 'input-file', checked: i === state.input,
      title: 'Run the pipeline on this file',
    });
    radio.hidden = multi || state.files.length < 2;
    radio.addEventListener('change', () => { state.input = i; update(); });
    const remove = el('button', { type: 'button', className: 'icon', title: 'Close this file', textContent: '×' });
    remove.setAttribute('aria-label', `Close ${f.name}`);
    remove.addEventListener('click', () => {
      state.files.splice(i, 1);
      if (state.input >= i && state.input > 0) state.input--;
      renderFiles();
      update();
    });
    return el('li', {}, radio, el('span', { className: 'name', textContent: f.name }),
              el('span', { className: 'size', textContent: formatBytes(f.file.size) }), remove);
  }));
  const note = $('files-note');
  note.hidden = state.files.length < 2;
  note.textContent = multi
    ? `mmcif_${firstTool().name} reads all ${state.files.length} files, in this order.`
    : 'The pipeline runs on the selected file.';
}

function setupDrop() {
  const drop = $('drop');
  $('file-input').addEventListener('change', (event) => {
    addFiles(event.target.files);
    event.target.value = '';
  });
  for (const type of ['dragenter', 'dragover']) {
    drop.addEventListener(type, (event) => { event.preventDefault(); drop.classList.add('over'); });
  }
  for (const type of ['dragleave', 'drop']) {
    drop.addEventListener(type, () => drop.classList.remove('over'));
  }
  drop.addEventListener('drop', (event) => {
    event.preventDefault();
    addFiles(event.dataTransfer.files);
  });
  // A file dropped next to the drop zone would replace the page.
  window.addEventListener('dragover', (event) => event.preventDefault());
  window.addEventListener('drop', (event) => event.preventDefault());
}

// ------------------------------------------------------------------ pipeline

function firstTool() {
  return state.steps.length ? state.byName.get(state.steps[0].tool) : null;
}

// The pipeline as command-line text without the files: "selchain -A | tidy".
function pipelineText() {
  return allSteps().map((s) => {
    const words = [s.tool];
    if (s.option) words.push('-' + s.option);
    if (s.label) words.push('--label');
    return words.join(' ');
  }).join(' | ');
}

function parsePipeline(text) {
  const steps = [];
  for (const part of text.split('|')) {
    const words = part.trim().split(/\s+/).filter(Boolean);
    if (!words.length) continue;
    const tool = words[0].replace(/^mmcif_/, '');
    if (!state.byName.has(tool)) continue;
    const step = { tool, option: '', label: false, help: false };
    for (const word of words.slice(1)) {
      if (word === '--label') step.label = true;
      else if (word.startsWith('-') && !step.option) step.option = word.slice(1);
    }
    steps.push(step);
  }
  return steps;
}

// The steps of a #pipeline. A last mmcif_delmeta is the finishing step; a
// pipeline that ends in mmCIF without it was made with that step unticked.
function loadPipeline(text) {
  state.steps = parsePipeline(text);
  const last = state.steps[state.steps.length - 1];
  if (state.steps.length > 1 && last.tool === 'delmeta') {
    state.steps.pop();
    state.finish.option = last.option;
  } else if (state.steps.length && !finishProblem()) {
    state.finish.on = false;
  }
}

// Why mmcif_delmeta cannot finish the pipeline (null if it can).
function finishProblem() {
  if (!state.byName.has('delmeta')) return 'Not available in this version of mmcif-tools.';
  const last = state.steps[state.steps.length - 1];
  if (!last) return 'Not used: there are no steps.';
  const tool = state.byName.get(last.tool);
  if (last.tool === 'delmeta') return 'Not used: the last step is already mmcif_delmeta.';
  if (tool.split) return `Not used: mmcif_${tool.name} writes several files.`;
  if (tool.last_only) return `Not used: mmcif_${tool.name} does not write mmCIF.`;
  return null;
}

// The steps that run: those of the list, then mmcif_delmeta if it is ticked and fits.
function allSteps() {
  if (!state.finish.on || finishProblem()) return state.steps;
  return [...state.steps, { tool: 'delmeta', option: state.finish.option, label: false }];
}

function setupFinish() {
  const tool = state.byName.get('delmeta');
  $('finish').hidden = !tool;
  if (!tool) return;
  $('finish-on').checked = state.finish.on;
  $('finish-option').value = state.finish.option;
  $('finish-option').placeholder = tool.placeholder + ' (optional)';
  $('finish-help').textContent = tool.help;
  $('finish-on').addEventListener('change', () => { state.finish.on = $('finish-on').checked; update(); });
  $('finish-option').addEventListener('input', () => { state.finish.option = $('finish-option').value.trim(); update(); });
  $('finish-option').addEventListener('keydown', (event) => { if (event.key === 'Enter') run(); });
  $('finish-help-button').addEventListener('click', () => {
    state.finish.help = !state.finish.help;
    $('finish-help').hidden = !state.finish.help;
    $('finish-help-button').setAttribute('aria-expanded', String(state.finish.help));
  });
}

// Why step i cannot run (null if it can); `bad` problems block the run.
function stepProblem(step, i) {
  const tool = state.byName.get(step.tool);
  const last = state.steps.length - 1;
  if (tool.unavailable) return { bad: true, text: `Not available here: this tool ${tool.unavailable}.` };
  if (i > 0 && tool.first_only) return { bad: true, text: 'This tool can only be the first step.' };
  if (i < last && tool.last_only) {
    return { bad: true, text: 'This tool does not write mmCIF, so it can only be the last step.' };
  }
  if (tool.option === 'required' && !step.option) return { bad: true, text: 'This tool needs an option.' };
  if (tool.multi_input && state.files.length < 2) {
    return { bad: false, text: 'This tool reads all the opened files: open two or more.' };
  }
  return null;
}

function toolSelect(step) {
  const select = el('select', { className: 'tool' });
  select.setAttribute('aria-label', 'Tool');
  let group = null;
  for (const tool of state.tools) {
    if (!group || group.label !== tool.group) {
      group = el('optgroup', { label: tool.group });
      select.append(group);
    }
    group.append(el('option', {
      value: tool.name,
      textContent: 'mmcif_' + tool.name + (tool.unavailable ? ' (not available)' : ''),
      disabled: Boolean(tool.unavailable),
    }));
  }
  select.value = step.tool;
  return select;
}

function stepItem(step, i) {
  const tool = state.byName.get(step.tool);
  const select = toolSelect(step);
  select.addEventListener('change', () => {
    Object.assign(step, { tool: select.value, option: '', label: false });
    renderSteps();
  });

  const option = el('input', {
    type: 'text', value: step.option, spellcheck: false, autocomplete: 'off',
    placeholder: tool.option === 'none' ? 'no option' : tool.placeholder + (tool.option === 'optional' ? ' (optional)' : ''),
    disabled: tool.option === 'none',
  });
  option.setAttribute('aria-label', `Option of mmcif_${tool.name}`);
  option.addEventListener('input', () => { step.option = option.value.trim(); update(); });
  option.addEventListener('keydown', (event) => { if (event.key === 'Enter') run(); });

  const label = el('label', { className: 'label-opt', title: 'Use label_* identifiers instead of auth_* ones' });
  const box = el('input', { type: 'checkbox', checked: step.label });
  box.addEventListener('change', () => { step.label = box.checked; update(); });
  label.append(box, ' --label');
  label.hidden = !tool.label;

  const helpButton = el('button', { type: 'button', className: 'icon', textContent: '?', title: 'Help' });
  helpButton.setAttribute('aria-expanded', String(step.help));
  const help = el('pre', { className: 'help', textContent: tool.help, hidden: !step.help });
  helpButton.addEventListener('click', () => {
    step.help = !step.help;
    help.hidden = !step.help;
    helpButton.setAttribute('aria-expanded', String(step.help));
  });

  const remove = el('button', { type: 'button', className: 'icon', textContent: '×', title: 'Remove this step' });
  remove.disabled = state.steps.length === 1;
  remove.addEventListener('click', () => { state.steps.splice(i, 1); renderSteps(); });

  const row = el('div', { className: 'step-row' },
    el('span', { className: 'step-num', textContent: String(i + 1) }), select,
    el('span', { className: 'option' }, el('span', { className: 'dash', textContent: '-' }), option),
    label, helpButton, remove);
  return el('li', { className: 'step' }, row,
    el('p', { className: 'summary', textContent: tool.summary }),
    el('p', { className: 'hint', hidden: true }), help);
}

function renderSteps() {
  $('steps').replaceChildren(...state.steps.map(stepItem));
  renderFiles();
  update();
}

function commandLine() {
  if (!state.steps.length) return '';
  const first = firstTool();
  let inputs;
  if (first.multi_input) inputs = state.files.map((f) => shellQuote(f.name));
  else inputs = state.files.length ? [shellQuote(state.files[state.input].name)] : ['<file>'];
  const steps = allSteps();
  const parts = steps.map((s, i) => {
    const words = ['mmcif_' + s.tool];
    if (s.option) words.push(shellQuote('-' + s.option));
    if (s.label) words.push('--label');
    if (i === 0) words.push(...inputs);
    return words.join(' ');
  });
  const last = state.byName.get(steps[steps.length - 1].tool);
  // the output is named after the last tool, or the one before a last mmcif_delmeta (as in runner.py)
  const named = steps[steps.length - (steps.length > 1 && last.name === 'delmeta' ? 2 : 1)].tool;
  const base = state.files.length ? stem(state.files[first.multi_input ? 0 : state.input].name) : 'out';
  let line = parts.join(' | ');
  if (!last.split) line += ' > ' + shellQuote(`${base}_${named}${last.ext}`);
  return line;
}

// Refresh everything that depends on the state: hints, buttons, command line.
function update() {
  let blocked = !state.ready || state.running || !state.files.length || !state.steps.length;
  document.querySelectorAll('#steps .step').forEach((item, i) => {
    const problem = stepProblem(state.steps[i], i);
    const hint = item.querySelector('.hint');
    hint.hidden = !problem;
    hint.textContent = problem ? problem.text : '';
    hint.classList.toggle('bad', Boolean(problem && problem.bad));
    if (problem && problem.bad) blocked = true;
  });
  const problem = finishProblem();
  const finishing = state.finish.on && !problem;
  $('finish').classList.toggle('unused', !finishing);
  $('finish-num').textContent = finishing ? String(state.steps.length + 1) : '–';
  $('finish-hint').hidden = !state.finish.on || !problem;
  $('finish-hint').textContent = problem || '';
  $('finish-option').disabled = !state.finish.on;
  $('run').disabled = blocked;
  $('run').hidden = state.running;
  $('stop').hidden = !state.running;
  $('add-step').disabled = !state.tools.length;
  $('command-box').hidden = !state.steps.length;
  $('command').textContent = commandLine();
  const hash = '#' + encodeURIComponent(pipelineText());
  if (state.steps.length && location.hash !== hash) history.replaceState(null, '', hash);
}

// ----------------------------------------------------------------------- run

function run() {
  if ($('run').disabled) return;
  const first = firstTool();
  const inputs = first.multi_input ? state.files.map((f) => f.name) : [state.files[state.input].name];
  state.sent = allSteps().map((s) => ({ tool: s.tool, option: s.option, label: s.label }));
  state.running = true;
  state.startedAt = performance.now();
  state.progress = 'Running';
  state.timer = setInterval(tick, 250);
  tick();
  update();
  state.worker.postMessage({
    type: 'run',
    files: state.files.map((f) => ({ name: f.name, file: f.file })),
    inputs,
    steps: state.sent,
  });
}

function tick() {
  if (!state.running) return;
  const seconds = (performance.now() - state.startedAt) / 1000;
  setStatus(`${state.progress}… ${seconds.toFixed(1)} s`, { busy: true });
}

function finishRun() {
  state.running = false;
  clearInterval(state.timer);
  update();
}

function stop() {
  finishRun();
  state.worker.terminate();
  setStatus('Stopped. Restarting Python…', { busy: true });
  startWorker();
}

// -------------------------------------------------------------------- result

function clearResult() {
  $('result').hidden = false;
  $('result-rows').replaceChildren();
  $('messages').replaceChildren();
  $('output-head').hidden = true;
  $('output-note').hidden = true;
  $('output').hidden = true;
  if (state.downloadUrl) URL.revokeObjectURL(state.downloadUrl);
  state.downloadUrl = null;
}

function showError(text) {
  clearResult();
  $('messages').append(el('div', { className: 'message bad', textContent: text }));
}

function showResult(result, output) {
  if (result.error) {
    showError(result.error);
    return;
  }
  clearResult();
  const rows = result.steps.map((s, i) => el('tr', {},
    el('td', { textContent: String(i + 1) }),
    el('td', {}, el('code', { textContent: s.args.map(shellQuote).join(' ') })),
    el('td', { className: s.status === 0 ? 'ok' : 'bad', textContent: s.status === 0 ? '0' : `${s.status} (failed)` }),
    el('td', { className: 'num', textContent: formatSeconds(s.seconds) }),
    el('td', { className: 'num', textContent: formatBytes(s.size) })));
  rows.push(el('tr', { className: 'total' }, el('td'), el('td', { textContent: 'Total' }), el('td'),
    el('td', { className: 'num', textContent: formatSeconds(result.seconds) }), el('td')));
  $('result-rows').replaceChildren(...rows);

  for (const s of result.steps) {
    if (!s.stderr) continue;
    $('messages').append(el('div', {
      className: 'message' + (s.status === 0 ? '' : ' bad'),
      textContent: `mmcif_${s.tool}:\n${s.stderr.replace(/\n$/, '')}`,
    }));
  }
  if (result.steps.length < state.sent.length) {
    const failed = result.steps.length;
    $('messages').append(el('div', {
      className: 'message bad',
      textContent: `Step ${failed} failed, so the steps after it did not run.`,
    }));
  }
  if (output) showOutput(result, output);
}

function showOutput(result, output) {
  const type = result.files ? 'application/zip' : 'text/plain';
  state.downloadUrl = URL.createObjectURL(new Blob([output], { type }));
  $('download').href = state.downloadUrl;
  $('download').download = result.name;
  $('download').textContent = `Download ${result.name}`;
  $('output-head').hidden = false;

  const pre = $('output');
  pre.hidden = false;
  if (result.files) {
    $('output-info').textContent = `${result.files.length} files, ${formatBytes(output.length)} zipped.`;
    pre.textContent = result.files.map(([name, size]) => `${name.padEnd(32)} ${formatBytes(size).padStart(10)}`).join('\n');
    return;
  }
  if (!output.length) {
    $('output-info').textContent = 'The output is empty.';
    pre.hidden = true;
    return;
  }
  let shown = output.subarray(0, PREVIEW_BYTES);
  if (shown.length < output.length) {
    const end = shown.lastIndexOf(10);
    if (end > 0) shown = shown.subarray(0, end + 1);
  }
  $('output-note').hidden = !/\.cif$/.test(result.name);
  const text = new TextDecoder().decode(shown);
  const lines = text.split('\n').length - (text.endsWith('\n') ? 1 : 0);
  $('output-info').textContent = shown.length < output.length
    ? `The first ${lines.toLocaleString()} lines (${formatBytes(shown.length)} of ${formatBytes(output.length)}).`
    : `${lines.toLocaleString()} lines, ${formatBytes(output.length)}.`;
  pre.textContent = text;
}

// ---------------------------------------------------------------------- init

// A new step goes before a last step that does not write mmCIF (a report).
$('add-step').addEventListener('click', () => {
  const last = state.steps[state.steps.length - 1];
  const at = last && state.byName.get(last.tool).last_only ? state.steps.length - 1 : state.steps.length;
  state.steps.splice(at, 0, { tool: 'selchain', option: '', label: false, help: false });
  renderSteps();
});
$('run').addEventListener('click', run);
$('stop').addEventListener('click', stop);
setupDrop();
startWorker();
