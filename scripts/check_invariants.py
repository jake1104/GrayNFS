"""ECM/GNFS invariant guard (sibling to check_version.py and check_api_surface.py).

Pins the undocumented operating tables that decide every ECM/GNFS
operating point (probe/idea-docs-invariants.md M1/M2/M3/M8/M11):

M1  ecm_plan tiers in crates/graynfs-core/src/cascade.rs:
    digits <= 20 gives (B1 11000, B2 1900000, 74 curves);
    digits <= 25 gives (B1 50000, B2 13000000, 214 curves);
    above 25, one growth step per 5 digits multiplies B1 by 5,
    B2 by 7 and curves by 3 from the 25-digit base.

M2  scaled_anchor rows plus caps and growth steps in
    crates/graynfs-core/src/factor.rs:
    degree 3 uses anchor 27 bits with bounds (2000, 500, 1600, 700);
    degree 4 uses anchor 27 bits with bounds (2000, 2000, 6000, 1500);
    degree 5 uses anchor 43 bits with bounds (4000, 20000, 72000, 8000);
    caps rational 50000, algebraic 250000, sieve_a 500000,
    sieve_b 60000, large_prime 1000000;
    per-column growth steps rational 1, algebraic 5, sieve_a 3,
    sieve_b 2, large_prime 1 with shift = step * delta_bits / 32
    clamped to 20 (growth factor 1 << shift).

M3  sieve acceptance threshold buckets plus the fallback ladder in
    crates/graynfs-core/src/relations.rs:
    prime counts 0-200 give 4; 201-600 give 8; 601-1500 give 12;
    1501-5000 give 20; anything larger gives 28;
    the fallback ladder starts from the initial threshold and pushes
    4 (only when initial > 4) then 1 (only when initial > 1) then 0,
    sorts descending, deduplicates, and uses vec![0] when there are no
    sieve roots; threshold 0 is skipped when pairs exceed 2_000_000.

M8  cheap-slice plus Fermat constants and formula chain in
    crates/graynfs-core/src/cascade.rs:
    cheap_slice_ms returns remaining / 16;
    cheap_ecm_plan uses (B1 2000, B2 200000, 4-32 curves, full 32
    curves when no budget is known) with the per-tier wall-time rate
    table CHEAP_MS_PER_CURVE_BY_BITS ((max_bits 128, 8 ms), (u64 MAX,
    16 ms)) selected by input bits; curves fill the cheap slice at
    16 x the tier rate because the slice is remaining / 16, and the
    slice deadline stays the wall guard on overfill;
    cap_ecm_slice treats burn below live / 32 as negligible
    (NEGLIGIBLE_BURN_DIVISOR 32) and otherwise caps at live / 2;
    fermat_cap_for uses CEILING 200000, FLOOR 20000,
    REACH 65536.0, SCALE 85.0 with estimated = REACH * REACH /
    (8 * sqrt_n) from 2^(bits / 2), scaled by SCALE and clamped
    into [FLOOR, CEILING].

M11 primality constants in crates/graynfs-core/src/prime.rs:
    DETERMINISTIC_U64_BASES holds the 12 bases
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37;
    AKS_PROOF_MAX_BITS is 64.

The EXPECTED tables below record the current measured values. This guard
is change-detection, not truth: any edit to a pinned row or formula
constant fails with exit 1 and a diff, forcing the change to be reviewed
as an explicit invariant decision. Update the EXPECTED tables only inside
the same change that intentionally retunes the table, never silently.

Reads the source files, parses the rows with regular expressions,
and compares the parse output against EXPECTED. Usage:

    python scripts/check_invariants.py
    python scripts/check_invariants.py --root path/to/tree
    python scripts/check_invariants.py --cascade-path /tmp/cascade.rs --factor-path /tmp/factor.rs
    python scripts/check_invariants.py --relations-path /tmp/relations.rs --prime-path /tmp/prime.rs

Exit 0 when every parsed row matches, else exit 1. Stdlib only.
"""
import argparse
import re
import sys
from pathlib import Path


M1_EXPECTED_TIERS = (
    (20, 11000, 1900000, 74),
    (25, 50000, 13000000, 214),
)
M1_EXPECTED_BASE = (50000, 13000000, 214)
M1_EXPECTED_GROWTH = (5, 7, 3)
M1_EXPECTED_STEP_RULE = (25, 5)


M2_EXPECTED_ANCHORS = {
    3: (27, 2000, 500, 1600, 700),
    4: (27, 2000, 2000, 6000, 1500),
    5: (43, 4000, 20000, 72000, 8000),
}
M2_EXPECTED_CAPS = {
    "rational": 50000,
    "algebraic": 250000,
    "sieve_a": 500000,
    "sieve_b": 60000,
    "large_prime": 1000000,
}
M2_EXPECTED_STEPS = {
    "rational": 1,
    "algebraic": 5,
    "sieve_a": 3,
    "sieve_b": 2,
    "large_prime": 1,
}
M2_EXPECTED_SHIFT = (32, 20)


M3_EXPECTED_BUCKETS = (
    (0, 200, 4),
    (201, 600, 8),
    (601, 1500, 12),
    (1501, 5000, 20),
)
M3_EXPECTED_CATCHALL = 28
M3_EXPECTED_LADDER = {
    "empty": 0,
    "guards": (4, 1),
    "pushes": (4, 1, 0),
    "descending": True,
    "dedup": True,
}
M3_EXPECTED_PAIRS_CAP = 2_000_000


M8_EXPECTED_CHEAP_DIVISOR = 16
M8_EXPECTED_CHEAP_ECM = (2000, 200000, 4, 32)
M8_EXPECTED_CHEAP_TIER_TABLE = ((128, 8), (18446744073709551615, 16))
M8_EXPECTED_CLAMP_SHAPE = (True, True)
M8_EXPECTED_NEGLIGIBLE_DIVISOR = 32
M8_EXPECTED_FERMAT = (200000, 20000, 65536.0, 85.0)
M8_EXPECTED_FORMULA = (True, True, True, True)


M11_EXPECTED_BASE_COUNT = 12
M11_EXPECTED_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
M11_EXPECTED_AKS_MAX_BITS = 64


def read_text(path):
    return Path(path).read_text(encoding="utf-8")


def plain_number(token):
    return int(token.replace("_", ""))


def parse_ecm_plan(text):
    tiers = []
    tier_pattern = (
        r"if target_digits <= (\d+)\s*\{\s*"
        r"return \((\d[\d_]*),\s*(\d[\d_]*),\s*(\d[\d_]*)\);"
    )
    for match in re.finditer(tier_pattern, text):
        tiers.append(
            (
                int(match.group(1)),
                plain_number(match.group(2)),
                plain_number(match.group(3)),
                plain_number(match.group(4)),
            )
        )
    base_match = re.search(
        r"let mut bound_one = (\d[\d_]*)_u64;\s*"
        r"let mut bound_two = (\d[\d_]*)_u64;\s*"
        r"let mut curves = (\d[\d_]*)_u64;",
        text,
    )
    base = None
    if base_match:
        base = (
            plain_number(base_match.group(1)),
            plain_number(base_match.group(2)),
            plain_number(base_match.group(3)),
        )
    growth = []
    for name in ("bound_one", "bound_two", "curves"):
        match = re.search(
            re.escape(name)
            + r"\s*=\s*"
            + re.escape(name)
            + r"\.saturating_mul\((\d+)\)",
            text,
        )
        growth.append(int(match.group(1)) if match else None)
    step_match = re.search(r"saturating_sub\((\d+)\)\.div_ceil\((\d+)\)", text)
    step_rule = None
    if step_match:
        step_rule = (int(step_match.group(1)), int(step_match.group(2)))
    return {
        "tiers": tuple(tiers),
        "base": base,
        "growth": tuple(growth),
        "step_rule": step_rule,
    }


def parse_scaled_table(text):
    anchors = {}
    anchor_pattern = (
        r"(\d) => Some\(\((\d+),"
        r"\s*(\d[\d_]*),\s*(\d[\d_]*),\s*(\d[\d_]*),\s*(\d[\d_]*)\)\)"
    )
    for match in re.finditer(anchor_pattern, text):
        anchors[int(match.group(1))] = (
            int(match.group(2)),
            plain_number(match.group(3)),
            plain_number(match.group(4)),
            plain_number(match.group(5)),
            plain_number(match.group(6)),
        )
    caps = {}
    cap_names = (
        ("rational", "SCALED_RATIONAL_CAP"),
        ("algebraic", "SCALED_ALGEBRAIC_CAP"),
        ("sieve_a", "SCALED_SIEVE_A_CAP"),
        ("sieve_b", "SCALED_SIEVE_B_CAP"),
        ("large_prime", "SCALED_LARGE_PRIME_CAP"),
    )
    for key, const in cap_names:
        match = re.search(re.escape(const) + r": usize = (\d[\d_]*);", text)
        caps[key] = plain_number(match.group(1)) if match else None
    steps = {}
    step_pattern = r"let (\w+)_growth = scaled_growth\(anchor_bits, bits, (\d+)\);"
    for match in re.finditer(step_pattern, text):
        steps[match.group(1)] = int(match.group(2))
    shift_match = re.search(r"saturating_mul\(delta\) / (\d+)", text)
    min_match = re.search(r"1_usize << shift\.min\((\d+)\)", text)
    shift = None
    if shift_match and min_match:
        shift = (int(shift_match.group(1)), int(min_match.group(1)))
    return {
        "anchors": anchors,
        "caps": caps,
        "steps": steps,
        "shift": shift,
    }


def plain_float(token):
    return float(token.replace("_", ""))


def parse_threshold_buckets(text):
    block_match = re.search(
        r"let initial_threshold = match base\.algebraic_modulus_primes\.len\(\) \{(.*?)\};",
        text,
        re.DOTALL,
    )
    block = block_match.group(1) if block_match else text
    rows = []
    for match in re.finditer(r"(\d+)\.\.=(\d+)\s*=>\s*(\d+)", block):
        rows.append(
            (int(match.group(1)), int(match.group(2)), int(match.group(3)))
        )
    catch_match = re.search(r"_\s*=>\s*(\d+)", block)
    catchall = int(catch_match.group(1)) if catch_match else None
    return {"buckets": tuple(rows), "catchall": catchall}


def parse_threshold_ladder(text):
    ladder_match = re.search(
        r"let thresholds = if sieve_roots\.is_empty\(\) \{(.*?)\};",
        text,
        re.DOTALL,
    )
    ladder = ladder_match.group(0) if ladder_match else ""
    empty_match = re.search(r"vec!\[(\d+)\]", ladder)
    empty = int(empty_match.group(1)) if empty_match else None
    guards = tuple(
        int(match.group(1))
        for match in re.finditer(r"if initial_threshold > (\d+)", ladder)
    )
    pushes = tuple(
        int(match.group(1))
        for match in re.finditer(r"values\.push\((\d+)\)", ladder)
    )
    descending = (
        "sort_unstable_by(|left, right| right.cmp(left))" in ladder
    )
    dedup = "values.dedup()" in ladder
    pairs_match = re.search(r"threshold == 0 && pairs > ([\d_]+)", text)
    pairs_cap = plain_number(pairs_match.group(1)) if pairs_match else None
    return {
        "ladder": {
            "empty": empty,
            "guards": guards,
            "pushes": pushes,
            "descending": descending,
            "dedup": dedup,
        },
        "pairs_cap": pairs_cap,
    }


def parse_cheap_slice(text):
    match = re.search(
        r"remaining_ms\.map\(\|remaining\| remaining / (\d+)\)", text
    )
    return int(match.group(1)) if match else None


def parse_cheap_ecm(text):
    consts = []
    for name in ("BOUND_ONE", "BOUND_TWO", "MIN_CURVES", "MAX_CURVES"):
        match = re.search(
            r"const %s: u64 = " % name + r"([\d_]+);",
            text,
        )
        consts.append(plain_number(match.group(1)) if match else None)
    table = None
    table_match = re.search(
        r"CHEAP_MS_PER_CURVE_BY_BITS:\s*\[\(u64,\s*u128\);\s*\d+\]\s*=\s*\[(.*?)\];",
        text,
        re.DOTALL,
    )
    if table_match:
        rows = []
        for entry in re.finditer(
            r"\(\s*(u64::MAX|\d[\d_]*)\s*,\s*(\d[\d_]*)\s*\)",
            table_match.group(1),
        ):
            cap = entry.group(1)
            rows.append(
                (
                    18446744073709551615
                    if cap == "u64::MAX"
                    else plain_number(cap),
                    plain_number(entry.group(2)),
                )
            )
        table = tuple(rows)
    clamp_shape = (
        ".clamp(MIN_CURVES as u128, MAX_CURVES as u128)" in text,
        "None => MAX_CURVES" in text,
    )
    return {"consts": tuple(consts), "table": table,
            "clamp_shape": clamp_shape}


def parse_negligible_divisor(text):
    match = re.search(
        r"const NEGLIGIBLE_BURN_DIVISOR: u128 = ([\d_]+);", text
    )
    return plain_number(match.group(1)) if match else None


def parse_fermat_cap(text):
    int_consts = []
    for name in ("CEILING", "FLOOR"):
        match = re.search(
            re.escape("const %s: u64 = " % name) + r"([\d_]+);", text
        )
        int_consts.append(plain_number(match.group(1)) if match else None)
    float_consts = []
    for name in ("REACH", "SCALE"):
        match = re.search(
            re.escape("const %s: f64 = " % name) + r"([\d_\.]+);", text
        )
        float_consts.append(
            plain_float(match.group(1)) if match else None
        )
    formula = (
        "REACH * REACH / (8.0 * sqrt_n)" in text,
        "2_f64.powf(bits as f64 / 2.0)" in text,
        "(SCALE * estimated).clamp(0.0, CEILING as f64)" in text,
        "(FLOOR as f64).max(scaled).min(CEILING as f64)" in text,
    )
    return {
        "consts": (int_consts[0], int_consts[1], float_consts[0],
                   float_consts[1]),
        "formula": formula,
    }


def parse_primality(text):
    bases_match = re.search(
        r"pub const DETERMINISTIC_U64_BASES: \[u64; (\d+)\] = \[([^\]]+)\];",
        text,
    )
    count = None
    bases = None
    if bases_match:
        count = int(bases_match.group(1))
        bases = tuple(
            int(token.strip())
            for token in bases_match.group(2).split(",")
            if token.strip()
        )
    aks_match = re.search(
        r"pub const AKS_PROOF_MAX_BITS: u64 = (\d+);", text
    )
    aks = int(aks_match.group(1)) if aks_match else None
    return {"base_count": count, "bases": bases, "aks_max_bits": aks}


def check_invariants(cascade_path, factor_path, relations_path=None,
                      prime_path=None):
    results = []

    def record(name, found, expected):
        ok = found is not None and found == expected
        results.append(
            {"name": name, "found": found, "expected": expected, "ok": ok}
        )

    try:
        cascade_text = read_text(cascade_path)
    except OSError as exc:
        print("MISMATCH: cannot read %s: %s" % (cascade_path, exc))
        return 1
    try:
        factor_text = read_text(factor_path)
    except OSError as exc:
        print("MISMATCH: cannot read %s: %s" % (factor_path, exc))
        return 1
    if relations_path is None:
        relations_path = str(
            Path(cascade_path).parent / "relations.rs")
    if prime_path is None:
        prime_path = str(Path(cascade_path).parent / "prime.rs")
    try:
        relations_text = read_text(relations_path)
    except OSError as exc:
        print("MISMATCH: cannot read %s: %s" % (relations_path, exc))
        return 1
    try:
        prime_text = read_text(prime_path)
    except OSError as exc:
        print("MISMATCH: cannot read %s: %s" % (prime_path, exc))
        return 1

    ecm = parse_ecm_plan(cascade_text)
    record("M1 ecm_plan tier rows", ecm["tiers"], M1_EXPECTED_TIERS)
    record("M1 ecm_plan growth base", ecm["base"], M1_EXPECTED_BASE)
    record("M1 ecm_plan growth multipliers", ecm["growth"], M1_EXPECTED_GROWTH)
    record("M1 ecm_plan step rule", ecm["step_rule"], M1_EXPECTED_STEP_RULE)

    scaled = parse_scaled_table(factor_text)
    for degree in sorted(M2_EXPECTED_ANCHORS):
        record(
            "M2 scaled_anchor degree %d" % degree,
            scaled["anchors"].get(degree),
            M2_EXPECTED_ANCHORS[degree],
        )
    for key in M2_EXPECTED_CAPS:
        record(
            "M2 scaled cap %s" % key,
            scaled["caps"].get(key),
            M2_EXPECTED_CAPS[key],
        )
    for key in M2_EXPECTED_STEPS:
        record(
            "M2 scaled growth step %s" % key,
            scaled["steps"].get(key),
            M2_EXPECTED_STEPS[key],
        )
    record("M2 scaled shift rule", scaled["shift"], M2_EXPECTED_SHIFT)

    buckets = parse_threshold_buckets(relations_text)
    record("M3 threshold buckets", buckets["buckets"],
           M3_EXPECTED_BUCKETS)
    record("M3 threshold catchall", buckets["catchall"],
           M3_EXPECTED_CATCHALL)
    ladder = parse_threshold_ladder(relations_text)
    record("M3 fallback ladder", ladder["ladder"], M3_EXPECTED_LADDER)
    record("M3 zero-threshold pairs cap", ladder["pairs_cap"],
           M3_EXPECTED_PAIRS_CAP)

    record("M8 cheap slice divisor", parse_cheap_slice(cascade_text),
           M8_EXPECTED_CHEAP_DIVISOR)
    cheap = parse_cheap_ecm(cascade_text)
    record("M8 cheap ecm consts", cheap["consts"], M8_EXPECTED_CHEAP_ECM)
    record("M8 cheap tier table", cheap["table"],
           M8_EXPECTED_CHEAP_TIER_TABLE)
    record("M8 cheap ecm clamp shape", cheap["clamp_shape"],
           M8_EXPECTED_CLAMP_SHAPE)
    record("M8 negligible burn divisor",
           parse_negligible_divisor(cascade_text),
           M8_EXPECTED_NEGLIGIBLE_DIVISOR)
    fermat = parse_fermat_cap(cascade_text)
    record("M8 fermat consts", fermat["consts"], M8_EXPECTED_FERMAT)
    record("M8 fermat formula chain", fermat["formula"],
           M8_EXPECTED_FORMULA)

    primality = parse_primality(prime_text)
    record("M11 mr base count", primality["base_count"],
           M11_EXPECTED_BASE_COUNT)
    record("M11 mr bases", primality["bases"], M11_EXPECTED_BASES)
    record("M11 aks proof max bits", primality["aks_max_bits"],
           M11_EXPECTED_AKS_MAX_BITS)

    failures = [entry for entry in results if not entry["ok"]]
    for entry in results:
        status = "OK" if entry["ok"] else "MISMATCH"
        print(
            "[%s] %s: found=%s expected=%s"
            % (status, entry["name"], entry["found"], entry["expected"])
        )

    if failures:
        print("RESULT: MISMATCH %d of %d checks failed" % (len(failures), len(results)))
        return 1
    print("RESULT: all %d checks match" % len(results))
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Pin M1/M2/M3/M8/M11 invariant tables.")
    parser.add_argument(
        "--root",
        default=None,
        help="Repo root override (default: parent of scripts dir).",
    )
    parser.add_argument(
        "--cascade-path",
        default=None,
        help="Direct cascade.rs path override (tamper-copy testing).",
    )
    parser.add_argument(
        "--factor-path",
        default=None,
        help="Direct factor.rs path override (tamper-copy testing).",
    )
    parser.add_argument(
        "--relations-path",
        default=None,
        help="Direct relations.rs path override (tamper-copy testing).",
    )
    parser.add_argument(
        "--prime-path",
        default=None,
        help="Direct prime.rs path override (tamper-copy testing).",
    )
    args = parser.parse_args(argv)
    if args.root is not None:
        root = Path(args.root)
    else:
        root = Path(__file__).resolve().parent.parent
    cascade_path = args.cascade_path
    if cascade_path is None:
        cascade_path = root / "crates" / "graynfs-core" / "src" / "cascade.rs"
    factor_path = args.factor_path
    if factor_path is None:
        factor_path = root / "crates" / "graynfs-core" / "src" / "factor.rs"
    relations_path = args.relations_path
    if relations_path is None:
        relations_path = (
            root / "crates" / "graynfs-core" / "src" / "relations.rs")
    prime_path = args.prime_path
    if prime_path is None:
        prime_path = root / "crates" / "graynfs-core" / "src" / "prime.rs"
    return check_invariants(str(cascade_path), str(factor_path),
                            str(relations_path), str(prime_path))


if __name__ == "__main__":
    sys.exit(main())
