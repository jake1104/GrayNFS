"""API-surface drift guard (sibling to check_version.py).

Compares the live signatures against their checked-in copies:
- crates/graynfs-py/src/lib.rs #[pyo3(signature = ...)] vs
  python/graynfs/_graynfs.pyi (ordered parameter names, all four decimals).
- python/graynfs/__init__.py vs python/graynfs/__init__.pyi
  (factor, factor_with_progress, stage_timings, factor_report, factor_batch).
- factor_batch keyword tail vs factor keyword tail (ns + fail_fast head).
- README CLI flag table vs crates/graynfs-cli/src/main.rs long flags
  (one direction only: every documented flag must exist in the CLI;
  CLI-only flags such as --progress are reported, never failures).
- __all__ carries factor_batch.

Deliberately narrow: parameter names and flag presence only, never types,
defaults, or behavior. Exit 0 when all checks pass, else exit 1.
Stdlib only.
"""
import re
import sys
from pathlib import Path


def read_text(path):
    return Path(path).read_text(encoding="utf-8")


def split_top_level(body):
    parts = []
    depth = 0
    current = ""
    for char in body:
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        if char == "," and depth == 0:
            parts.append(current)
            current = ""
        else:
            current += char
    if current.strip():
        parts.append(current)
    return parts


def signature_names(body):
    names = []
    for part in split_top_level(body):
        part = part.strip()
        if not part or part == "*":
            continue
        names.append(part.split("=")[0].strip())
    return names


def _strip_string_literals(line):
    return re.sub(r"'(?:[^'\\]|\\.)*'|\"(?:[^\"\\]|\\.)*\"", '""', line)


def py_param_names(body):
    names = []
    for line in body.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("*") or stripped.startswith("#"):
            continue
        code = _strip_string_literals(stripped)
        depth = 0
        colon_index = -1
        for index, char in enumerate(code):
            if char in "([{":
                depth += 1
            elif char in ")]}":
                depth -= 1
            elif char == ":" and depth == 0:
                colon_index = index
                break
        if colon_index < 0:
            continue
        name = code[:colon_index].strip().rstrip(",")
        if not re.match(r"^\*{0,2}[A-Za-z_]\w*$", name):
            continue
        if name.startswith("*"):
            continue
        names.append(name)
    return [name for name in names if name]


def parse_lib_signatures(text):
    found = {}
    marker = re.compile(r"#\[pyo3\(signature\s*=\s*\(")
    name_pattern = re.compile(r"\)\]\s*\nfn\s+(\w+)\s*\(")
    for marker_match in marker.finditer(text):
        pos = marker_match.end()
        depth = 1
        index = pos
        in_single = False
        in_double = False
        escaped = False
        while index < len(text) and depth > 0:
            char = text[index]
            if escaped:
                escaped = False
            elif char == "\\" and (in_single or in_double):
                escaped = True
            elif char == "'" and not in_double:
                in_single = not in_single
            elif char == '"' and not in_single:
                in_double = not in_double
            elif not in_single and not in_double:
                if char == "(":
                    depth += 1
                elif char == ")":
                    depth -= 1
            index += 1
        body = text[pos:index - 1]
        name_match = name_pattern.match(text[index:index + 400])
        if name_match is None:
            continue
        found[name_match.group(1)] = signature_names(body)
    return found


def parse_stub_defs(text):
    found = {}
    pattern = re.compile(r"def\s+(\w+)\s*\((.*?)\)\s*->", re.DOTALL)
    for match in pattern.finditer(text):
        found[match.group(1)] = py_param_names(match.group(2))
    return found


def parse_cli_long_flags(text):
    match = re.search(r"struct Cli \{(.*?)\n\}", text, re.DOTALL)
    if match is None:
        return None
    flags = set()
    lines = match.group(1).splitlines()
    index = 0
    saw_long = False
    while index < len(lines):
        stripped = lines[index].strip()
        if stripped.startswith("#[arg("):
            block = stripped
            depth = block.count("(") - block.count(")")
            while depth > 0 and index + 1 < len(lines):
                index += 1
                block += "\n" + lines[index].strip()
                depth = block.count("(") - block.count(")")
            saw_long = re.search(r"\blong\b", _strip_string_literals(block)) is not None
            index += 1
            continue
        field = re.match(r"^([a-z][a-z0-9_]*)\s*:", stripped)
        if field and saw_long:
            flags.add("--" + field.group(1).replace("_", "-"))
        if field:
            saw_long = False
        index += 1
    return flags


def parse_readme_flags(text):
    return set(re.findall(r"^\|\s*`(--[a-z0-9-]+)`", text, re.MULTILINE))


def check_repo(root):
    root = Path(root)
    lib_path = root / "crates" / "graynfs-py" / "src" / "lib.rs"
    stub_path = root / "python" / "graynfs" / "_graynfs.pyi"
    init_path = root / "python" / "graynfs" / "__init__.py"
    init_stub_path = root / "python" / "graynfs" / "__init__.pyi"
    cli_path = root / "crates" / "graynfs-cli" / "src" / "main.rs"
    readme_path = root / "README.md"

    results = []

    def record(name, location, ok, detail=""):
        results.append({"name": name, "location": location, "ok": ok, "detail": detail})

    try:
        lib_text = read_text(lib_path)
        stub_text = read_text(stub_path)
        init_text = read_text(init_path)
        init_stub_text = read_text(init_stub_path)
        cli_text = read_text(cli_path)
        readme_text = read_text(readme_path)
    except OSError as exc:
        print("MISMATCH: cannot read input: %s" % exc)
        return 1

    lib_sigs = parse_lib_signatures(lib_text)
    stub_defs = parse_stub_defs(stub_text)
    for name in (
        "factor_decimal",
        "factor_with_progress_decimal",
        "stage_timings_decimal",
        "factor_report_decimal",
    ):
        live = lib_sigs.get(name)
        stub = stub_defs.get(name)
        if live is None:
            record("lib.rs %s" % name, str(lib_path), False, "signature not found")
        elif stub is None:
            record("_graynfs.pyi %s" % name, str(stub_path), False, "stub not found")
        else:
            record(
                "_graynfs.pyi %s matches lib.rs" % name,
                str(stub_path),
                stub == live,
                "stub=%s live=%s" % (stub, live) if stub != live else "",
            )

    init_defs = parse_stub_defs(init_text)
    init_stub_defs = parse_stub_defs(init_stub_text)
    for name in (
        "factor",
        "factor_with_progress",
        "stage_timings",
        "factor_report",
        "factor_batch",
    ):
        live = init_defs.get(name)
        stub = init_stub_defs.get(name)
        if live is None:
            record("__init__.py %s" % name, str(init_path), False, "def not found")
        elif stub is None:
            record("__init__.pyi %s" % name, str(init_stub_path), False, "stub not found")
        else:
            record(
                "__init__.pyi %s matches __init__.py" % name,
                str(init_stub_path),
                stub == live,
                "impl=%s stub=%s" % (live, stub) if stub != live else "",
            )

    factor_params = init_defs.get("factor", [])
    batch_params = init_defs.get("factor_batch", [])
    record(
        "factor_batch kwargs match factor kwargs",
        str(init_path),
        bool(factor_params)
        and bool(batch_params)
        and batch_params[:2] == ["ns", "fail_fast"]
        and batch_params[2:] == factor_params[1:],
        "batch=%s factor=%s" % (batch_params, factor_params)
        if not (
            bool(factor_params)
            and bool(batch_params)
            and batch_params[:2] == ["ns", "fail_fast"]
            and batch_params[2:] == factor_params[1:]
        )
        else "",
    )

    all_match = re.search(r"__all__\s*=\s*\[(.*?)\]", init_text, re.DOTALL)
    all_names = set(re.findall(r'"([^"]+)"', all_match.group(1))) if all_match else set()
    record(
        "__all__ carries factor_batch",
        str(init_path),
        "factor_batch" in all_names,
        "found=%s" % sorted(all_names) if "factor_batch" not in all_names else "",
    )

    cli_flags = parse_cli_long_flags(cli_text)
    readme_flags = parse_readme_flags(readme_text)
    if cli_flags is None:
        record("CLI long flags parse", str(cli_path), False, "struct Cli not found")
        cli_flags = set()
    stale = sorted(readme_flags - cli_flags)
    record(
        "README table flags exist in CLI",
        str(readme_path),
        not stale,
        "stale=%s" % stale if stale else "",
    )
    record(
        "README table carries --seed",
        str(readme_path),
        "--seed" in readme_flags,
        "--seed missing from flag table" if "--seed" not in readme_flags else "",
    )
    record(
        "CLI carries --seed",
        str(cli_path),
        "--seed" in cli_flags,
        "--seed missing from Cli struct" if "--seed" not in cli_flags else "",
    )
    only_cli = sorted(cli_flags - readme_flags)
    print("NOTE: CLI-only flags (prose-covered, table-optional): %s" % only_cli)

    failures = [entry for entry in results if not entry["ok"]]
    for entry in results:
        status = "OK" if entry["ok"] else "MISMATCH"
        extra = " (%s)" % entry["detail"] if entry["detail"] else ""
        print("[%s] %s (%s)%s" % (status, entry["name"], entry["location"], extra))

    if failures:
        print("RESULT: MISMATCH %d of %d checks failed" % (len(failures), len(results)))
        return 1
    print("RESULT: all %d checks pass" % len(results))
    return 0


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(description="Check API-surface drift.")
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
