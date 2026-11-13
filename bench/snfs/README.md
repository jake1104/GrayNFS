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

- `corpus-large.json`: 8 special-form inputs, 61-97 bit, degrees 3/4/5, both
  signs, `a == 1` and `a > 1`. Every kept entry passes the production gates
  (`is_snfs_eligible`, scaled identity `f(m) == a^(k-1) * N` with
  `f(m) mod N == 0`, scoped irreducibility); exclusions are recorded, never
  force-fitted. Flags: `--skip-preprocessing --threads 1`, 120 s budget.
- `results-large.json`: v0.2.5.7.2 same-budget showdown over the large corpus
  (same input, same 120 s budget, same default seed; SNFS preset vs generic
  `--no-snfs`; forced entry degree; 3 reps per arm; 48 runs total). Verdict:
  NO-GO (0 SNFS / 1 generic / 7 voided; bar was SNFS advantage on a majority
  of cells). Mechanism: usable-relation starvation at equal budget on the
  large forms (E1 SNFS usable 73/145/320 of target 3196 while raw stays high,
  `max_a` escalation to the deadline; no within-degree generic fallback),
  so the large-form scope does not promote.
- Scope: SNFS is a small-special-form path only. No new `bench/pins.json`
  cells were promoted from the large corpus, and no release-name claim
  follows from the toy corpus.
