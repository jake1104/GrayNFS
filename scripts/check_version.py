"""Repo version-drift guard.

Compares the canonical workspace version against every user-visible copy.
Canonical source: root Cargo.toml [workspace.package] version.
Checked copies: pyproject.toml, README badge/wheel/example, test assertion.

Deliberately ignored: CHANGELOG history, Cargo.lock dependency versions,
task IDs like 2026.10.10.1.x, README v0.1.1/v0.1.0 history sections.
Only targeted patterns are matched, never a bare generic version scan.

Exit 0 when all copies equal the canonical version, else exit 1.
Stdlib only.
"""
import argparse
import re
import sys
from pathlib import Path

VERSION_RE = r"(\d+\.\d+\.\d+)(?!\.\d)"


def read_text(path):
    return Path(path).read_text(encoding="utf-8")


def section_body(text, section):
    pattern = re.compile(
        r"^\[" + re.escape(section) + r"\](.*?)(?=^\[|\Z)",
        re.MULTILINE | re.DOTALL,
    )
    match = pattern.search(text)
    return match.group(1) if match else ""


def extract_section_version(text, section):
    body = section_body(text, section)
    match = re.search(r'^\s*version\s*=\s*"([^"]+)"', body, re.MULTILINE)
    return match.group(1).strip() if match else None


def check_repo(root):
    root = Path(root)
    cargo_path = root / "Cargo.toml"
    pyproject_path = root / "pyproject.toml"
    readme_path = root / "README.md"
    test_path = root / "python" / "tests" / "test_basic.py"

    results = []

    def record(name, location, found, expected):
        ok = found is not None and found == expected
        results.append(
            {"name": name, "location": location, "found": found, "expected": expected, "ok": ok}
        )

    try:
        cargo_text = read_text(cargo_path)
    except OSError as exc:
        print("MISMATCH: cannot read %s: %s" % (cargo_path, exc))
        return 1
    cargo_version = extract_section_version(cargo_text, "workspace.package")
    if cargo_version is None:
        print("MISMATCH: no version under [workspace.package] in %s" % cargo_path)
        return 1
    print("canonical Cargo [workspace.package] version: %s" % cargo_version)

    try:
        pyproject_text = read_text(pyproject_path)
    except OSError as exc:
        print("MISMATCH: cannot read %s: %s" % (pyproject_path, exc))
        return 1
    pyproject_version = extract_section_version(pyproject_text, "project")
    record("pyproject.toml [project] version", str(pyproject_path), pyproject_version, cargo_version)

    try:
        readme_text = read_text(readme_path)
    except OSError as exc:
        print("MISMATCH: cannot read %s: %s" % (readme_path, exc))
        return 1

    badge = re.search(r"badge/version-v" + VERSION_RE, readme_text)
    record(
        "README badge version",
        "%s (badge/version-v)" % readme_path,
        badge.group(1) if badge else None,
        cargo_version,
    )

    wheels = re.findall(r"graynfs-" + VERSION_RE + r"-cp38-abi3", readme_text)
    if wheels:
        for index, value in enumerate(wheels):
            record(
                "README wheel version [%d]" % index,
                "%s (dist graynfs-<ver>-cp38-abi3)" % readme_path,
                value,
                cargo_version,
            )
    else:
        record(
            "README wheel version",
            "%s (dist graynfs-<ver>-cp38-abi3)" % readme_path,
            None,
            cargo_version,
        )

    example = re.search(
        r"graynfs\.version\(\)[^\n'\"]*['\"]" + VERSION_RE + r"['\"]", readme_text
    )
    record(
        "README version() example",
        "%s (graynfs.version() comment)" % readme_path,
        example.group(1) if example else None,
        cargo_version,
    )

    installed = re.search(r"installed\s+" + VERSION_RE + r"\s+wheel", readme_text)
    record(
        "README installed-wheel note",
        "%s (installed <ver> wheel)" % readme_path,
        installed.group(1) if installed else None,
        cargo_version,
    )

    notes = re.search(
        r"^##\s+Release notes\s+\(v" + VERSION_RE, readme_text, re.MULTILINE
    )
    record(
        "README latest release-notes heading",
        "%s (first ## Release notes (v<ver>))" % readme_path,
        notes.group(1) if notes else None,
        cargo_version,
    )

    try:
        test_text = read_text(test_path)
    except OSError as exc:
        print("MISMATCH: cannot read %s: %s" % (test_path, exc))
        return 1
    assertion = re.search(
        r"graynfs\.version\(\)\s*==\s*[\"']" + VERSION_RE + r"[\"']", test_text
    )
    record(
        "test_basic.py version assertion",
        "%s (graynfs.version() ==)" % test_path,
        assertion.group(1) if assertion else None,
        cargo_version,
    )

    failures = [entry for entry in results if not entry["ok"]]
    for entry in results:
        status = "OK" if entry["ok"] else "MISMATCH"
        found = entry["found"] if entry["found"] is not None else "<missing>"
        print(
            "[%s] %s: found=%s expected=%s (%s)"
            % (status, entry["name"], found, entry["expected"], entry["location"])
        )

    print("Ignored by design: CHANGELOG history, Cargo.lock dep versions,")
    print("task IDs like 2026.10.10.1.x, README v0.1.1/v0.1.0 history sections.")

    if failures:
        print("RESULT: MISMATCH %d of %d checks failed" % (len(failures), len(results)))
        return 1
    print("RESULT: all %d checks match %s" % (len(results), cargo_version))
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description="Check repo version drift.")
    parser.add_argument(
        "--root",
        default=None,
        help="Repo root override (default: parent of scripts dir).",
    )
    args = parser.parse_args(argv)
    if args.root is not None:
        root = Path(args.root)
    else:
        root = Path(__file__).resolve().parent.parent
    return check_repo(root)


if __name__ == "__main__":
    sys.exit(main())
