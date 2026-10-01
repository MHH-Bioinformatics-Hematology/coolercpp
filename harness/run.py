#!/usr/bin/env python
"""Runs the coolercpp equivalence harness: every case through cooler 0.10.2
(oracle.py) and through coolercpp (the coolercpp-harness driver), compares the
outputs, measures peak RSS and CPU time of both processes, applies the gates
and writes a report.

    python harness/run.py --driver BUILD/harness/coolercpp-harness \\
        --hicx-data ~/src/HiCExplorer-v4/hicexplorer/test/test_data \\
        --out REPORT_DIR [--filter REGEX] [--cases FILE ...]

Path tokens in case files: @hicx/ (HiCExplorer test data), @cooler/ (this
repository's tests/data), @work/ (the side's own output directory),
@inputs/ (prepared input tables). Every @hicx/ and @cooler/ file must be listed
with its SHA-256 in harness/data_manifest.json.

Gates: coolercpp may not use more CPU time than Python, and its peak RSS must
stay within the case's rss_budget_mb, or within Python's peak RSS when the
case declares none.
"""

import argparse
import glob
import hashlib
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import compare  # noqa: E402

PROCESS_TIMEOUT_S = 3600


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


MEASURE = None  # path of coolercpp-measure, set from --measure


def run_measured(cmd, log_path):
    """Runs cmd through coolercpp-measure; returns exit status, CPU seconds and
    peak RSS of cmd itself (not of this runner, whose RSS a directly forked
    child would inherit)."""
    report = log_path + ".rusage.json"
    if os.path.exists(report):
        os.remove(report)
    with open(log_path, "wb") as log:
        proc = subprocess.Popen(
            [MEASURE, report, str(PROCESS_TIMEOUT_S), "--"] + cmd,
            stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
        )
        try:
            proc.wait()
        except BaseException:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
            raise
    if proc.returncode != 0 or not os.path.exists(report):
        return {"exit": -1, "cpu_s": float("nan"), "rss_mb": float("nan"), "wall_s": float("nan")}
    with open(report) as handle:
        return json.load(handle)


class Resolver:
    def __init__(self, hicx, manifest):
        self.roots = {"@hicx/": hicx, "@cooler/": os.path.join(REPO, "tests", "data")}
        self.manifest = manifest
        self.used = set()

    def substitute(self, value, work, inputs):
        if isinstance(value, dict):
            return {k: self.substitute(v, work, inputs) for k, v in value.items()}
        if isinstance(value, list):
            return [self.substitute(v, work, inputs) for v in value]
        if not isinstance(value, str):
            return value
        for token, root in self.roots.items():
            if value.startswith(token):
                self.used.add(value.split("::")[0])
                return os.path.join(root, value[len(token):])
        if value.startswith("@work/"):
            return os.path.join(work, value[len("@work/"):])
        if value.startswith("@inputs/"):
            return os.path.join(inputs, value[len("@inputs/"):])
        return value

    def verify(self, tokens):
        problems = []
        for token in sorted(tokens):
            path = self.substitute(token, "", "")
            entry = self.manifest.get(token)
            if entry is None:
                if token.endswith("does_not_exist.cool"):
                    continue
                problems.append(f"{token} is not listed in data_manifest.json")
            elif not os.path.exists(path):
                problems.append(f"{token} is missing ({path})")
            elif os.path.getsize(path) != entry["size"] or sha256(path) != entry["sha256"]:
                problems.append(f"{token} does not match its manifest checksum")
        return problems


def tokens_in(value, found):
    if isinstance(value, dict):
        for v in value.values():
            tokens_in(v, found)
    elif isinstance(value, list):
        for v in value:
            tokens_in(v, found)
    elif isinstance(value, str) and (value.startswith("@hicx/") or value.startswith("@cooler/")):
        found.add(value.split("::")[0])
    return found


def expand_chunks(case):
    for step in case.get("steps", []):
        chunks = step.get("pixel_chunks")
        if isinstance(chunks, dict):
            step["pixel_chunks"] = [
                f"@inputs/{chunks['input']}/pixels_{i}" for i in range(chunks["count"])
            ]


def normalize_dump(result):
    """Drops the normalized attributes from a dump's info item."""
    if result.get("kind") != "multi":
        return result
    info = result["items"].get("info")
    if info and info.get("kind") == "value" and info["value"].get("t") == "dict":
        info["value"]["v"] = [
            pair for pair in info["value"]["v"] if pair[0] not in compare.NORMALIZED_ATTRIBUTES
        ]
    return result


def load_result(directory):
    with open(os.path.join(directory, "result.json")) as handle:
        return json.load(handle)


def run_case(case, args, resolver, out_root):
    name = case["name"]
    work = os.path.join(out_root, "cases", re.sub(r"[^A-Za-z0-9_.-]", "_", name))
    shutil.rmtree(work, ignore_errors=True)
    inputs_dir = os.path.join(work, "inputs")
    os.makedirs(inputs_dir)
    record = {"name": name, "op": case["op"], "problems": [], "notes": []}

    missing = resolver.verify(tokens_in(case, set()))
    if missing:
        record["problems"].extend(missing)
        record["verdict"] = "FAIL"
        return record

    python = sys.executable
    oracle = os.path.join(HERE, "oracle.py")
    for input_name, spec in (case.get("inputs") or {}).items():
        spec_path = os.path.join(inputs_dir, f"{input_name}.json")
        with open(spec_path, "w") as handle:
            json.dump(resolver.substitute(spec, "", inputs_dir), handle)
        prep = run_measured(
            [python, oracle, "prep", spec_path, os.path.join(inputs_dir, input_name)],
            os.path.join(inputs_dir, f"{input_name}.log"),
        )
        if prep["exit"] != 0:
            record["problems"].append(f"input preparation {input_name} failed")
            record["verdict"] = "FAIL"
            return record

    expand_chunks(case)
    sides = {}
    for side in ("py", "cpp"):
        side_dir = os.path.join(work, side)
        files_dir = os.path.join(side_dir, "files")
        os.makedirs(files_dir)
        resolved = resolver.substitute(case, files_dir, inputs_dir)
        case_path = os.path.join(side_dir, "case.json")
        with open(case_path, "w") as handle:
            json.dump(resolved, handle)
        out_dir = os.path.join(side_dir, "out")
        os.makedirs(out_dir)
        cmd = [python, oracle, "run", case_path, out_dir] if side == "py" else [args.driver, case_path, out_dir]
        measure = run_measured(cmd, os.path.join(side_dir, "log.txt"))
        sides[side] = {"dir": out_dir, "files": files_dir, "case": resolved, "measure": measure}
        if measure["exit"] != 0 or not os.path.exists(os.path.join(out_dir, "result.json")):
            record["problems"].append(f"{side} process failed with exit status {measure['exit']}")

    record["py"] = sides["py"]["measure"]
    record["cpp"] = sides["cpp"]["measure"]
    if record["problems"]:
        record["verdict"] = "FAIL"
        return record

    py_doc = load_result(sides["py"]["dir"])
    cpp_doc = load_result(sides["cpp"]["dir"])
    record["py"]["op_s"] = py_doc["op_seconds"]
    record["cpp"]["op_s"] = cpp_doc["op_seconds"]
    cmp = compare.compare_results(
        py_doc["result"], cpp_doc["result"], sides["py"]["dir"], sides["cpp"]["dir"],
        messages=case.get("messages", "exact"),
    )
    if py_doc["result"]["kind"] == "error":
        record["notes"].append(f"both raise {py_doc['result']['type']}")

    # Any case that declares output files has them compared structurally and
    # read back by the other side: create cases, and the fileops cases that
    # copy, move and link groups.
    if case.get("outputs") and py_doc["result"]["kind"] != "error" and cmp.ok:
        outputs = case.get("outputs", [])
        for output in outputs:
            py_file = resolver.substitute(output, sides["py"]["files"], inputs_dir)
            cpp_file = resolver.substitute(output, sides["cpp"]["files"], inputs_dir)
            diffs = compare.compare_descriptions(
                compare.describe_file(py_file), compare.describe_file(cpp_file)
            )
            for diff in diffs:
                cmp.fail(f"structure {output}", diff)
        record["notes"].append(f"structure compared: {len(outputs)} file(s)")
        for uri in case.get("cross_read", outputs):
            label = re.sub(r"[^A-Za-z0-9_.-]", "_", uri)
            py_uri = resolver.substitute(uri, sides["py"]["files"], inputs_dir)
            cpp_uri = resolver.substitute(uri, sides["cpp"]["files"], inputs_dir)
            runs = {
                "py_reads_py": ([python, oracle, "run"], py_uri),
                "py_reads_cpp": ([python, oracle, "run"], cpp_uri),
                "cpp_reads_py": ([args.driver], py_uri),
            }
            dumps = {}
            for key, (prefix, target) in runs.items():
                directory = os.path.join(work, "cross", label, key)
                os.makedirs(directory)
                case_path = os.path.join(directory, "case.json")
                with open(case_path, "w") as handle:
                    json.dump({"op": "dump", "uri": target}, handle)
                result = run_measured(prefix + [case_path, directory], os.path.join(directory, "log.txt"))
                if result["exit"] != 0:
                    cmp.fail(f"cross read {uri}", f"{key} failed")
                    continue
                dumps[key] = (normalize_dump(load_result(directory)["result"]), directory)
            if len(dumps) == 3:
                base, base_dir = dumps["py_reads_py"]
                for key in ("py_reads_cpp", "cpp_reads_py"):
                    other, other_dir = dumps[key]
                    compare.compare_results(base, other, base_dir, other_dir, where=f"{key}:{uri}", cmp=cmp)
        record["notes"].append("cross read in both directions")

    record["problems"].extend(cmp.problems)
    record["float_class"] = cmp.float_class
    record["message_class"] = cmp.message_class

    budget = case.get("rss_budget_mb", record["py"]["rss_mb"])
    record["rss_budget_mb"] = budget
    record["cpu_gate"] = record["cpp"]["cpu_s"] <= record["py"]["cpu_s"]
    record["rss_gate"] = record["cpp"]["rss_mb"] <= budget
    if not record["cpu_gate"]:
        record["problems"].append(
            f"CPU gate: coolercpp {record['cpp']['cpu_s']:.3f} s > Python {record['py']['cpu_s']:.3f} s"
        )
    if not record["rss_gate"]:
        record["problems"].append(
            f"RSS gate: coolercpp {record['cpp']['rss_mb']:.1f} MB > budget {budget:.1f} MB"
        )
    record["verdict"] = "PASS" if not record["problems"] else "FAIL"
    return record


def write_report(records, out_root, args):
    passed = sum(1 for r in records if r["verdict"].lower() == "pass")
    failed = sum(1 for r in records if r["verdict"].lower() == "fail")
    lines = [
        "# coolercpp equivalence report",
        "",
        f"Date: {time.strftime('%Y-%m-%d')}",
        f"Cases: {len(records)}, pass: {passed}, fail: {failed}",
        "",
        "| case | verdict | floats | errors | Python CPU s | C++ CPU s | Python RSS MB | C++ RSS MB | RSS budget MB |",
        "|---|---|---|---|---:|---:|---:|---:|---:|",
    ]
    for r in records:
        py = r.get("py", {})
        cpp = r.get("cpp", {})
        lines.append(
            "| {name} | {verdict} | {floats} | {errors} | {pc} | {cc} | {pr} | {cr} | {b} |".format(
                name=r["name"],
                verdict=r["verdict"],
                floats=r.get("float_class") or "-",
                errors=r.get("message_class") or "-",
                pc=f"{py.get('cpu_s', float('nan')):.3f}",
                cc=f"{cpp.get('cpu_s', float('nan')):.3f}",
                pr=f"{py.get('rss_mb', float('nan')):.1f}",
                cr=f"{cpp.get('rss_mb', float('nan')):.1f}",
                b=f"{r.get('rss_budget_mb', float('nan')):.1f}",
            )
        )
    problems = [r for r in records if r["problems"]]
    if problems:
        lines += ["", "## Problems", ""]
        for r in problems:
            for p in r["problems"][:20]:
                lines.append(f"- {r['name']}: {p}")
    with open(os.path.join(out_root, "report.md"), "w") as handle:
        handle.write("\n".join(lines) + "\n")
    with open(os.path.join(out_root, "report.json"), "w") as handle:
        json.dump(records, handle, indent=1)
    return passed, failed


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--driver", required=True)
    parser.add_argument("--measure", default=None,
                        help="coolercpp-measure; defaults to the one next to the driver")
    parser.add_argument("--hicx-data", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--cases", nargs="*", default=None)
    parser.add_argument("--filter", default=None)
    parser.add_argument("--update-manifest", action="store_true",
                        help="rewrite data_manifest.json from the files the cases reference")
    args = parser.parse_args()
    global MEASURE
    MEASURE = args.measure or os.path.join(os.path.dirname(os.path.abspath(args.driver)), "coolercpp-measure")
    if not os.access(MEASURE, os.X_OK):
        sys.exit(f"coolercpp-measure not found at {MEASURE}")

    case_files = args.cases or sorted(glob.glob(os.path.join(HERE, "cases", "*.json")))
    cases = []
    for path in case_files:
        with open(path) as handle:
            cases.extend(json.load(handle)["cases"])
    names = [c["name"] for c in cases]
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        sys.exit(f"duplicate case names: {duplicates}")
    if args.filter:
        cases = [c for c in cases if re.search(args.filter, c["name"])]

    manifest_path = os.path.join(HERE, "data_manifest.json")
    manifest = {}
    if os.path.exists(manifest_path):
        with open(manifest_path) as handle:
            manifest = json.load(handle)
    resolver = Resolver(os.path.abspath(os.path.expanduser(args.hicx_data)), manifest)

    if args.update_manifest:
        tokens = set()
        for case in cases:
            tokens_in(case, tokens)
        for token in sorted(tokens):
            path = resolver.substitute(token, "", "")
            if os.path.exists(path):
                manifest[token] = {"sha256": sha256(path), "size": os.path.getsize(path)}
        with open(manifest_path, "w") as handle:
            json.dump(dict(sorted(manifest.items())), handle, indent=1)
            handle.write("\n")
        print(f"manifest: {len(manifest)} files")
        return 0

    os.makedirs(args.out, exist_ok=True)
    records = []
    for case in cases:
        record = run_case(json.loads(json.dumps(case)), args, resolver, args.out)
        records.append(record)
        py = record.get("py", {})
        cpp = record.get("cpp", {})
        print(
            f"{record['verdict']} {record['name']} floats={record.get('float_class')} "
            f"cpu py={py.get('cpu_s', 0):.3f}s cpp={cpp.get('cpu_s', 0):.3f}s "
            f"rss py={py.get('rss_mb', 0):.0f}MB cpp={cpp.get('rss_mb', 0):.0f}MB",
            flush=True,
        )
        for problem in record["problems"][:5]:
            print(f"    {problem}", flush=True)
        write_report(records, args.out, args)
    if not records:
        write_report(records, args.out, args)
        print("no cases selected")
        return 2
    passed, failed = write_report(records, args.out, args)
    print(f"cases {len(records)}: pass {passed}, fail {failed}; report {os.path.join(args.out, 'report.md')}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    signal.signal(signal.SIGTERM, lambda *_: sys.exit(143))
    sys.exit(main())
