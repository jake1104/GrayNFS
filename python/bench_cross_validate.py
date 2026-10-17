from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

script_directory = str(Path(__file__).resolve().parent)
sys.path = [entry for entry in sys.path if entry not in {"", script_directory}]

import graynfs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="logs/graynfs-benchmark.json")
    parser.add_argument("--output", default="logs/graynfs-benchmark-python.json")
    parser.add_argument("--timeout-seconds", type=int, default=30)
    return parser.parse_args()


def product(values: list[int]) -> int:
    result = 1
    for value in values:
        result *= value
    return result


def main() -> None:
    args = parse_args()
    source = json.loads(Path(args.input).read_text(encoding="utf-8"))
    results: list[dict[str, Any]] = []
    for observation in source:
        value = int(observation["input"])
        events: list[dict[str, Any]] = []
        started = time.perf_counter()
        status = "ok"
        error = None
        factors: list[int] = []
        try:
            factors = graynfs.factor_with_progress(
                value,
                events.append,
                rational_base_bound=5000,
                algebraic_base_bound=5000,
                sieve_a_bound=2500,
                sieve_b_bound=300,
                timeout_seconds=args.timeout_seconds,
                allow_fallback=True,
            )
        except Exception as exception:
            status = "error"
            error = str(exception)
        elapsed_ms = (time.perf_counter() - started) * 1000.0
        rust_factors = [int(value) for value in observation.get("factors", [])]
        expected = [int(value) for value in observation.get("expected", [])]
        verified = (
            status == "ok"
            and product(factors) == value
            and sorted(factors) == sorted(expected)
            and sorted(factors) == sorted(rust_factors)
        )
        results.append(
            {
                "name": observation["name"],
                "input": observation["input"],
                "bits": observation["bits"],
                "status": status,
                "error": error,
                "verified": verified,
                "total_ms": elapsed_ms,
                "python_factors": factors,
                "rust_factors": rust_factors,
                "expected": expected,
                "events": events,
            }
        )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print("case,status,verified,total_ms,event_count")
    for result in results:
        print(
            f"{result['name']},{result['status']},{result['verified']},"
            f"{result['total_ms']:.3f},{len(result['events'])}"
        )
    print(f"json: {output}")


if __name__ == "__main__":
    main()
