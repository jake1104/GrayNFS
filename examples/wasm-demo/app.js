const inputEl = document.getElementById("input");
const budgetEl = document.getElementById("budget");
const factorBtn = document.getElementById("factor");
const cancelBtn = document.getElementById("cancel");
const statusEl = document.getElementById("status");
const wallEl = document.getElementById("wall");
const factorsEl = document.getElementById("factors");
const methodEl = document.getElementById("method");
const degreeEl = document.getElementById("degree");
const polyEl = document.getElementById("poly");
const rootEl = document.getElementById("root");
const relationsEl = document.getElementById("relations");
const rankEl = document.getElementById("rank");
const elapsedEl = document.getElementById("elapsed");
const rawEl = document.getElementById("raw");

let worker = null;
let requestId = 0;
let runStarted = 0;

function spawnWorker() {
  if (worker) {
    worker.terminate();
  }
  worker = new Worker(new URL("./worker.js", import.meta.url), { type: "module" });
  worker.onmessage = onWorkerMessage;
  worker.onerror = (event) => {
    setStatus("error: " + (event.message || "worker failed"));
    setRunning(false);
  };
}

function setStatus(text) {
  statusEl.textContent = text;
}

function setRunning(running) {
  factorBtn.disabled = running;
  cancelBtn.disabled = !running;
}

function showResult(payload, wallMs) {
  let data;
  try {
    data = JSON.parse(payload);
  } catch {
    setStatus("error: invalid JSON from worker");
    return;
  }
  rawEl.textContent = JSON.stringify(data, null, 2);
  wallEl.textContent = wallMs + " ms";
  if (!data.ok) {
    setStatus("failed: " + data.error);
    return;
  }
  factorsEl.textContent = data.factors.join(" x ");
  methodEl.textContent = data.method;
  degreeEl.textContent = String(data.degree);
  polyEl.textContent = data.polynomial;
  rootEl.textContent = data.root;
  relationsEl.textContent = String(data.relations);
  rankEl.textContent = String(data.rank);
  elapsedEl.textContent = String(data.elapsedMs) + " ms";
  const product = data.factors.reduce((acc, f) => acc * BigInt(f), 1n);
  if (product === BigInt(data.input)) {
    setStatus("done: product verifies");
  } else {
    setStatus("done: PRODUCT MISMATCH");
  }
  setRunning(false);
}

function onWorkerMessage(event) {
  const { id, ok, body, error, wallMs } = event.data;
  if (id !== requestId) {
    return;
  }
  if (!ok) {
    setStatus("failed: " + error);
    setRunning(false);
    return;
  }
  showResult(body, wallMs);
}

spawnWorker();

factorBtn.addEventListener("click", () => {
  spawnWorker();
  const input = inputEl.value.trim();
  const budgetMs = Number(budgetEl.value);
  if (!input) {
    setStatus("error: empty input");
    return;
  }
  requestId += 1;
  runStarted = performance.now();
  setRunning(true);
  setStatus("running: no progress events, wait for done");
  worker.postMessage({ id: requestId, input, budgetMs });
});

cancelBtn.addEventListener("click", () => {
  spawnWorker();
  const waited = Math.round(performance.now() - runStarted);
  wallEl.textContent = waited + " ms (cancelled)";
  setStatus("cancelled by worker termination");
  setRunning(false);
});
