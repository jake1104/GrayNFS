import { readFileSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const pkgDir = resolve(here, "..", "crates", "graynfs-wasm", "pkg");
const gluePath = resolve(pkgDir, "graynfs_wasm.js");
const wasmPath = resolve(pkgDir, "graynfs_wasm_bg.wasm");

const SUCCESS_KEYS = [
  "ok",
  "input",
  "factors",
  "factorConfidence",
  "method",
  "degree",
  "polynomial",
  "root",
  "relations",
  "rank",
  "elapsedMs",
  "stages",
];
const STAGE_KEYS = ["stage", "degree", "ms", "detail"];
const FAILURE_KEYS = ["ok", "input", "error"];

const CASES = [
  {
    name: "small-semiprime",
    input: "10403",
    budgetMs: 3000,
    expectedFactors: ["101", "103"],
  },
  {
    name: "project-euler-3",
    input: "600851475143",
    budgetMs: 5000,
    expectedFactors: ["71", "839", "1471", "6857"],
  },
  {
    name: "hostile105-budget",
    input: "35658606617546370626847644936101",
    budgetMs: 1000,
    expectFailure: "resource limit exceeded",
  },
];

function sameKeys(object, expected) {
  const actual = Object.keys(object).sort();
  const wanted = [...expected].sort();
  return actual.length === wanted.length && actual.every((key, index) => key === wanted[index]);
}

function productOf(factors) {
  return factors.reduce((product, factor) => product * BigInt(factor), 1n);
}

const glue = await import(pathToFileURL(gluePath).href);
const init = glue.default;
const bytes = readFileSync(wasmPath);

let initForm = "options-object";
try {
  await init({ module_or_path: bytes });
} catch (error) {
  initForm = "positional";
  await init(bytes);
}

console.log(`version: ${glue.version()}`);
console.log(`init-form: ${initForm}`);
console.log(`wasm-bytes: ${bytes.length}`);

let passed = 0;
let observedSuccess = 0;
let failed = 0;

for (const testCase of CASES) {
  const started = Date.now();
  const raw = glue.factorize(testCase.input, testCase.budgetMs);
  const wallMs = Date.now() - started;
  let payload;
  try {
    payload = JSON.parse(raw);
  } catch (error) {
    failed += 1;
    console.log(`${testCase.name}: FAIL invalid-json ${String(error)}`);
    console.log(`raw: ${raw}`);
    continue;
  }
  console.log(`${testCase.name}: json ${JSON.stringify(payload)}`);
  console.log(`${testCase.name}: keys ${Object.keys(payload).sort().join(",")}`);
  console.log(`${testCase.name}: wall-ms ${wallMs}`);

  if (testCase.expectFailure) {
    if (payload.ok === true) {
      observedSuccess += 1;
      console.log(`${testCase.name}: OBSERVED-SUCCESS budget-ms ${testCase.budgetMs}`);
      continue;
    }
    if (!sameKeys(payload, FAILURE_KEYS)) {
      failed += 1;
      console.log(`${testCase.name}: FAIL failure-key-set ${Object.keys(payload).join(",")}`);
      continue;
    }
    if (payload.input !== testCase.input) {
      failed += 1;
      console.log(`${testCase.name}: FAIL failure-input ${String(payload.input)}`);
      continue;
    }
    if (typeof payload.error !== "string" || !payload.error.includes(testCase.expectFailure)) {
      failed += 1;
      console.log(`${testCase.name}: FAIL failure-error ${String(payload.error)}`);
      continue;
    }
    passed += 1;
    console.log(`${testCase.name}: PASS`);
    continue;
  }

  if (payload.ok !== true) {
    failed += 1;
    console.log(`${testCase.name}: FAIL ok-false ${String(payload.error)}`);
    continue;
  }
  if (!sameKeys(payload, SUCCESS_KEYS)) {
    failed += 1;
    console.log(`${testCase.name}: FAIL success-key-set ${Object.keys(payload).join(",")}`);
    continue;
  }
  const stages = Array.isArray(payload.stages) ? payload.stages : null;
  if (!stages || stages.some((stage) => !sameKeys(stage, STAGE_KEYS))) {
    failed += 1;
    console.log(`${testCase.name}: FAIL stage-key-set`);
    continue;
  }
  const factors = Array.isArray(payload.factors) ? payload.factors.map(String) : [];
  const expected = [...testCase.expectedFactors].sort();
  const actual = [...factors].sort();
  if (actual.length !== expected.length || actual.some((factor, index) => factor !== expected[index])) {
    failed += 1;
    console.log(`${testCase.name}: FAIL factors ${actual.join(",")}`);
    continue;
  }
  if (productOf(factors) !== BigInt(testCase.input)) {
    failed += 1;
    console.log(`${testCase.name}: FAIL product ${productOf(factors).toString()}`);
    continue;
  }
  if (typeof payload.elapsedMs !== "number" || payload.elapsedMs > testCase.budgetMs + 2000) {
    failed += 1;
    console.log(`${testCase.name}: FAIL elapsedMs ${String(payload.elapsedMs)}`);
    continue;
  }
  passed += 1;
  console.log(`${testCase.name}: PASS`);
}

console.log(`summary: ${passed}/${CASES.length} PASS observed-success ${observedSuccess} failed ${failed}`);
if (failed > 0 || observedSuccess > 0) {
  process.exitCode = 1;
}
