#!/usr/bin/env python3
"""Run the E03 spherical pulse refinement study, one case after the other.

The script writes each case's input, launches the executable and records the exit status,
the wall time and the peak memory of every rank. It writes measured data and run records
only: it computes no orders and draws no figures.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import time


LEVELS = (96, 144, 216, 324)


def git(*args):
    return subprocess.check_output(["git", "-C", Path(__file__).parent, *args], text=True)


def study_identity(args, experiment, binary_sha256):
    """The settings that every case of one study directory must share."""
    return {"experiment": experiment, "binary_sha256": binary_sha256, "ranks": args.ranks,
            "threads": args.threads, "profile_sharpness": args.profile_sharpness}


def prepare_study(args, experiment="E03"):
    """Pin the executable and record the working-source identity before the first run."""
    args.output.mkdir(parents=True, exist_ok=True)
    pinned = args.output / "experiment"
    with args.executable.open("rb") as executable:
        identity = study_identity(args, experiment,
                                  hashlib.file_digest(executable, "sha256").hexdigest())
    metadata = args.output / "study.json"
    if metadata.exists():
        previous = json.loads(metadata.read_text())
        previous.setdefault("experiment", "E03")
        previous.setdefault("profile_sharpness", 1.0)
        if any(previous[key] != value for key, value in identity.items()):
            raise ValueError("Continuation requires the same experiment, executable, ranks, "
                             "threads and profile sharpness")
        with pinned.open("rb") as executable:
            if hashlib.file_digest(executable, "sha256").hexdigest() != identity["binary_sha256"]:
                raise ValueError("The pinned executable changed")
        return pinned
    shutil.copy2(args.executable, pinned)
    (args.output / "source.patch").write_text(git("diff", "HEAD", "--binary"))
    (args.output / "source-status.txt").write_text(git("status", "--short"))
    metadata.write_text(json.dumps({**identity,
                                    "git_head_at_launch": git("rev-parse", "HEAD").strip(),
                                    "machine": platform.node(),
                                    "platform": platform.platform()}, indent=2) + "\n")
    return pinned


def case_input(args, n, experiment):
    """One ordinary YAML input; the sharpness is written only when it differs from one."""
    parameters = f"cells_per_axis: {n}"
    if experiment == "E03" and n == LEVELS[0]:
        parameters += ", write_plane: true"
    if args.profile_sharpness != 1:
        parameters += f", profile_sharpness: {args.profile_sharpness:g}"
    return (f"schema_version: 1\ncase: {experiment}\nunits: {{system: normalized}}\n"
            f"parameters: {{{parameters}}}\noutput: {{directory: data}}\n")


def run_case(args, executable, n, experiment="E03"):
    """Keep the input, combined log, rank memory peaks and any partial output."""
    run = args.output / f"n{n}"
    status = run / "status.json"
    if run.exists():
        if completed_run(run):
            return True
        print(f"Preserving incomplete or failed run {run}", flush=True)
        return False
    run.mkdir()
    inp = run / "input.yaml"
    inp.write_text(case_input(args, n, experiment))
    env = {**os.environ, "OMP_NUM_THREADS": str(args.threads), "OMP_PROC_BIND": "false"}
    command = ["mpiexec", "-n", str(args.ranks), "--bind-to", "none", "sh", "-c",
               'exec "$1" -f %M -o "rss-${OMPI_COMM_WORLD_RANK}.kib" "$2" --input "$3"',
               "pulse-rank", shutil.which("time"), str(executable), str(inp)]
    start = time.monotonic()
    print(f"Running {experiment} N={n} on {args.ranks} ranks", flush=True)
    with (run / "run.log").open("w") as log:
        code = subprocess.run(command, cwd=run, env=env, stdout=log,
                              stderr=subprocess.STDOUT).returncode
    peaks = [int(line) for p in sorted(run.glob("rss-*.kib"))
             for line in p.read_text().splitlines() if line.isdecimal()]
    record = {"exit_code": code, "wall_seconds": time.monotonic() - start,
              "rank_peak_rss_kib": peaks, "sum_rank_peaks_kib": sum(peaks)}
    status.write_text(json.dumps(record, indent=2) + "\n")
    print(f"{experiment} N={n}: {record}", flush=True)
    return code == 0


def completed_run(run):
    status = run / "status.json"
    return status.is_file() and json.loads(status.read_text())["exit_code"] == 0


def study_arguments(description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--levels", type=int, nargs="+", choices=LEVELS, default=LEVELS[:2])
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--profile-sharpness", type=float, default=1.0,
                        help="sharpness beta of the compact radial profile (default 1, the E01 profile)")
    args = parser.parse_args()
    if min(args.ranks, args.threads) < 1 or not shutil.which("time"):
        parser.error("Positive ranks and threads and GNU time on PATH are required")
    if not args.profile_sharpness > 0:
        parser.error("The profile sharpness must be positive")
    args.executable, args.output = args.executable.resolve(), args.output.resolve()
    return args


def main():
    args = study_arguments(__doc__)
    executable = prepare_study(args)
    success = True
    for n in sorted(set(args.levels)):
        success = run_case(args, executable, n) and success
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
