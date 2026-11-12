# SNFS special-form record

This directory records the SNFS special-form path of the GNFS engine.
`method` stays on the `gnfs` family; the SNFS parameters are carried in
`detail`/`polynomial` only.

- `corpus.json`: 6 special-form inputs covering degrees 3/4/5, both signs,
  `a == 1` and `a > 1`. Every entry proves the scaled identity
  `f(m) == a^(k-1) * N` with `f(m) mod N == 0` (and `f(m) == N` on the
  `a == 1` entries). Flags: `--skip-preprocessing --threads 1`, 120 s budget.
- `results.json`: end-to-end A/B runs of the same inputs through the real SNFS
  preset versus generic GNFS (`use_snfs = false`), reps of 3, threads 1.
  Every row is product-verified (`product_ok`).

Ineligible inputs fall back to generic GNFS with identical factors or the same
clean failure; `low_weight`-only detections never promote.
