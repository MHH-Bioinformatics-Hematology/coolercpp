#!/usr/bin/env python
"""Proves that harness cases can fail: applies each mutation of
harness/mutations.json to a scratch copy of the installed cooler package,
runs the cases the mutation targets with that copy first on PYTHONPATH, and
requires at least one of them to fail. The installed package is never touched.

    python harness/mutate.py --driver BUILD/harness/coolercpp-harness \\
        --hicx-data ~/src/HiCExplorer-v4/hicexplorer/test/test_data \\
        --scratch DIR
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--driver", required=True)
    parser.add_argument("--hicx-data", required=True)
    parser.add_argument("--scratch", required=True)
    args = parser.parse_args()

    import cooler

    installed = os.path.dirname(os.path.abspath(cooler.__file__))
    with open(os.path.join(HERE, "mutations.json")) as handle:
        mutations = json.load(handle)["mutations"]

    results = []
    for mutation in mutations:
        root = os.path.join(args.scratch, mutation["name"])
        shutil.rmtree(root, ignore_errors=True)
        package = os.path.join(root, "site", "cooler")
        shutil.copytree(installed, package, ignore=shutil.ignore_patterns("__pycache__"))
        target = os.path.join(package, mutation["file"])
        with open(target) as handle:
            source = handle.read()
        if source.count(mutation["old"]) != 1:
            results.append((mutation["name"], "INVALID (pattern not found exactly once)"))
            continue
        with open(target, "w") as handle:
            handle.write(source.replace(mutation["old"], mutation["new"]))
        env = dict(os.environ)
        env["PYTHONPATH"] = os.path.join(root, "site") + os.pathsep + env.get("PYTHONPATH", "")
        report = os.path.join(root, "report")
        proc = subprocess.run(
            [sys.executable, os.path.join(HERE, "run.py"), "--driver", args.driver,
             "--hicx-data", args.hicx_data, "--out", report, "--filter", mutation["cases"]],
            env=env, capture_output=True, text=True,
        )
        lines = proc.stdout.splitlines()
        verdicts = [line.split()[0] for line in lines if re.match(r"^(pass|fail)\b", line, re.I)]
        failed = sum(1 for v in verdicts if v.lower() == "fail")
        if not verdicts:
            status = "INVALID (no case ran)"
        elif failed > 0 and proc.returncode == 1:
            status = f"caught ({failed} of {len(verdicts)} cases fail)"
        else:
            status = f"MISSED ({len(verdicts)} cases pass)"
        results.append((mutation["name"], status))
        print(f"{mutation['name']}: {status}", flush=True)

    missed = [name for name, status in results if not status.startswith("caught")]
    print(f"mutations {len(results)}: caught {len(results) - len(missed)}, not caught {len(missed)}")
    return 0 if not missed else 1


if __name__ == "__main__":
    sys.exit(main())
