import argparse
import hashlib
import json
import os
import shlex
import statistics
import subprocess
import sys
import time


BIT100_N = "770817968652704265841785071059"
HOSTILE105_N = "35658606617546370626847644936101"


CELLS = [
    {
        "name": "bit100-deg2",
        "n": BIT100_N,
        "extra": ["--skip-preprocessing", "--degree", "2"],
        "method": "siqs",
        "relations": 1878,
        "rank": 945,
        "ref_ms": 292.5,
        "band": "subsecond",
    },
    {
        "name": "bit100-default",
        "n": BIT100_N,
        "extra": [],
        "method": "ecm",
        "relations": None,
        "rank": None,
        "ref_ms": 167.5,
        "band": "subsecond",
    },
    {
        "name": "hostile105-deg2",
        "n": HOSTILE105_N,
        "extra": ["--skip-preprocessing", "--degree", "2"],
        "method": "siqs",
        # 31.1 KEEP (2026-10-28): 1177/926 is standing identity, 3269/907 historical.
        "relations": 1177,
        "rank": 926,
        "ref_ms": 3432.0,
        "band": "multisecond",
    },
]


SUBSECOND_TOL_MS = 30.0
MULTISECOND_UPPER_RATIO = 1.25
WALL_SLACK_MS = 500.0
DEFAULT_REPS = 5
DEFAULT_THREADS = 16
DEFAULT_TIMEOUT_SECONDS = 60
SUBPROCESS_TIMEOUT_SECONDS = 300


def resolve_bin(explicit):
    if explicit:
        return shlex.split(explicit, posix=False)
    env = os.environ.get("GRAYNFS_BIN")
    if env:
        return shlex.split(env, posix=False)
    for candidate in ("target/release/graynfs.exe", "target/release/graynfs-cli.exe"):
        if os.path.exists(candidate):
            return [candidate]
    return ["target/release/graynfs.exe"]


def bin_identity(bin_cmd):
    argv0 = bin_cmd[0] if bin_cmd else ""
    try:
        mtime = os.path.getmtime(argv0)
    except OSError:
        mtime = "missing"
    try:
        digest = hashlib.sha256()
        with open(argv0, "rb") as handle:
            for chunk in iter(lambda: handle.read(65536), b""):
                digest.update(chunk)
        sha256 = digest.hexdigest()[:16]
    except OSError:
        sha256 = "missing"
    return {"argv0": argv0, "mtime": mtime, "sha256": sha256}


def band_limits(ref_ms, band):
    if band == "subsecond":
        hi = ref_ms + SUBSECOND_TOL_MS
        return 0.0, hi, "ref %.1f ms upper-only +%.1f ms (limit %.1f ms)" % (ref_ms, SUBSECOND_TOL_MS, hi)
    hi = ref_ms * MULTISECOND_UPPER_RATIO
    return 0.0, hi, "ref %.1f ms upper-only x%.2f (limit %.1f ms)" % (ref_ms, MULTISECOND_UPPER_RATIO, hi)


def run_rep(bin_cmd, cell, threads, timeout_seconds, rep_index):
    args = (
        bin_cmd
        + [cell["n"], "--threads", str(threads), "--timeout-seconds", str(timeout_seconds), "--json"]
        + cell["extra"]
    )
    target = int(cell["n"])
    started = time.perf_counter()
    try:
        proc = subprocess.run(args, capture_output=True, text=True, timeout=SUBPROCESS_TIMEOUT_SECONDS)
    except subprocess.TimeoutExpired:
        wall_ms = (time.perf_counter() - started) * 1000.0
        return {"ok": False, "error": "subprocess timeout", "wall_ms": wall_ms}
    wall_ms = (time.perf_counter() - started) * 1000.0
    if proc.returncode != 0:
        return {"ok": False, "error": "exit code %d" % proc.returncode, "wall_ms": wall_ms,
                "stdout_tail": proc.stdout[-500:], "stderr_tail": proc.stderr[-500:]}
    lines = [line for line in proc.stdout.splitlines() if line.strip()]
    if not lines:
        return {"ok": False, "error": "empty stdout", "wall_ms": wall_ms}
    try:
        payload = json.loads(lines[-1])
    except json.JSONDecodeError as exc:
        return {"ok": False, "error": "json decode: %s" % exc, "wall_ms": wall_ms}
    elapsed = payload.get("elapsed_ms")
    if not isinstance(elapsed, (int, float)):
        return {"ok": False, "error": "missing elapsed_ms", "wall_ms": wall_ms}
    raw_factors = payload.get("factors", [])
    try:
        factors = [int(x) for x in raw_factors]
    except (TypeError, ValueError):
        return {"ok": False, "error": "non-integer factors", "wall_ms": wall_ms,
                "elapsed_ms": float(elapsed)}
    product = 1
    for value in factors:
        product *= value
    verified = len(factors) >= 2 and all(value > 1 for value in factors) and product == target
    record = {
        "ok": True,
        "elapsed_ms": float(elapsed),
        "wall_ms": wall_ms,
        "method": payload.get("method"),
        "relations": payload.get("relations"),
        "rank": payload.get("rank"),
        "factors": factors,
        "verified": verified,
        "wall_ok": wall_ms <= float(elapsed) + WALL_SLACK_MS,
    }
    problems = []
    if not verified:
        problems.append("product mismatch")
    if cell["method"] is not None and record["method"] != cell["method"]:
        problems.append("method %s != %s" % (record["method"], cell["method"]))
    if cell["relations"] is not None and record["relations"] != cell["relations"]:
        problems.append("relations %s != %s" % (record["relations"], cell["relations"]))
    if cell["rank"] is not None and record["rank"] != cell["rank"]:
        problems.append("rank %s != %s" % (record["rank"], cell["rank"]))
    if not record["wall_ok"]:
        problems.append("wall %.1f ms exceeds elapsed %.1f ms + %.1f ms slack" % (
            wall_ms, float(elapsed), WALL_SLACK_MS))
    record["problems"] = problems
    record["rep"] = rep_index
    return record


def run_cell(bin_cmd, cell, reps, threads, timeout_seconds):
    times = []
    records = []
    for rep in range(1, reps + 1):
        record = run_rep(bin_cmd, cell, threads, timeout_seconds, rep)
        records.append(record)
        if record.get("ok"):
            times.append(record["elapsed_ms"])
            print("%s rep%d/%d elapsed_ms=%.1f wall_ms=%.1f method=%s relations=%s rank=%s verified=%s wall_ok=%s%s" % (
                cell["name"], rep, reps, record["elapsed_ms"], record["wall_ms"],
                record["method"], record["relations"], record["rank"],
                record["verified"], record["wall_ok"],
                "" if not record["problems"] else " PROBLEMS=" + ";".join(record["problems"])))
        else:
            print("%s rep%d/%d HARNESS_ERROR %s wall_ms=%.1f" % (
                cell["name"], rep, reps, record.get("error"), record.get("wall_ms", -1.0)))
        sys.stdout.flush()
    harness_errors = [r for r in records if not r.get("ok")]
    rep_problems = [p for r in records if r.get("ok") for p in r.get("problems", [])]
    median = statistics.median(times) if len(times) == reps else None
    lo, hi, desc = band_limits(cell["ref_ms"], cell["band"])
    band_ok = median is not None and lo <= median <= hi
    cell_pass = not harness_errors and not rep_problems and band_ok
    print("%s n=%d times=%s median=%s ref=%.1f band=[%s] band_ok=%s harness_errors=%d rep_problems=%d verdict=%s" % (
        cell["name"], reps,
        "[" + ",".join("%.1f" % t for t in times) + "]",
        "none" if median is None else "%.1f" % median,
        cell["ref_ms"], desc, band_ok, len(harness_errors), len(rep_problems),
        "PASS" if cell_pass else "FAIL"))
    sys.stdout.flush()
    return {
        "name": cell["name"],
        "times_ms": times,
        "median_ms": median,
        "ref_ms": cell["ref_ms"],
        "band": cell["band"],
        "band_lo_ms": None,
        "band_hi_ms": hi,
        "band_ok": band_ok,
        "harness_errors": [r.get("error") for r in harness_errors],
        "rep_problems": rep_problems,
        "verdict": "PASS" if cell_pass else "FAIL",
        "bin_identity": bin_identity(bin_cmd),
    }


def main():
    parser = argparse.ArgumentParser(description="Frozen 3-cell perf smoke gate for GrayNFS release CLI")
    parser.add_argument("--bin", default=None, help="Release binary path, or full command line when it contains spaces")
    parser.add_argument("--reps", type=int, default=DEFAULT_REPS, help="Reps per cell (default 5)")
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS, help="CLI worker threads (default 16)")
    parser.add_argument("--timeout-seconds", type=int, default=DEFAULT_TIMEOUT_SECONDS,
                        help="Per-run CLI budget in seconds (default 60)")
    parser.add_argument("--cells", default=None,
                        help="Comma-separated cell subset, default runs all frozen cells")
    parser.add_argument("--json-out", default=None, help="Path for machine-readable summary JSON")
    args = parser.parse_args()
    bin_cmd = resolve_bin(args.bin)
    wanted = None
    if args.cells:
        wanted = set(part.strip() for part in args.cells.split(",") if part.strip())
    selected = [cell for cell in CELLS if wanted is None or cell["name"] in wanted]
    if wanted:
        unknown = wanted - set(cell["name"] for cell in CELLS)
        if unknown:
            print("ERROR unknown cells: %s" % sorted(unknown))
            return 2
    print("perf_smoke cells=%s reps=%d bin=%s threads=%d timeout_seconds=%d wall_slack_ms=%.1f" % (
        [cell["name"] for cell in selected], args.reps, " ".join(bin_cmd),
        args.threads, args.timeout_seconds, WALL_SLACK_MS))
    sys.stdout.flush()
    results = [run_cell(bin_cmd, cell, args.reps, args.threads, args.timeout_seconds) for cell in selected]
    overall = all(item["verdict"] == "PASS" for item in results)
    print("OVERALL verdict=%s" % ("PASS" if overall else "FAIL"))
    if args.json_out:
        with open(args.json_out, "w") as handle:
            json.dump({"overall": "PASS" if overall else "FAIL",
                       "bin_identity": bin_identity(bin_cmd),
                       "cells": results}, handle, indent=2)
        print("summary written to %s" % args.json_out)
    if any(item["harness_errors"] for item in results):
        return 2
    return 0 if overall else 1


if __name__ == "__main__":
    sys.exit(main())
