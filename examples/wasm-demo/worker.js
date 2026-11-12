import init, { factorize } from "../../crates/graynfs-wasm/pkg/graynfs_wasm.js";

const wasmUrl = new URL(
  "../../crates/graynfs-wasm/pkg/graynfs_wasm_bg.wasm",
  import.meta.url,
);
await init(wasmUrl);

self.onmessage = (event) => {
  const { id, input, budgetMs } = event.data;
  const start = Date.now();
  let payload;
  try {
    payload = factorize(input, budgetMs);
  } catch (error) {
    self.postMessage({ id, ok: false, error: String(error), wallMs: Date.now() - start });
    return;
  }
  self.postMessage({ id, ok: true, body: payload, wallMs: Date.now() - start });
};
