#!/usr/bin/env python3
"""Run and analyze the E03 spherical pulse refinement study."""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import time


LEVELS = (96, 144, 216, 324)
POLARIZATIONS = ("TEr", "TMr")
QUANTITIES = ("A_L2", "phi_L2", "E_L2", "E_Linf", "B_L2", "Q_L2", "Q_Linf", "G_L2", "G_Linf")


def git(*args):
    return subprocess.check_output(["git", "-C", Path(__file__).parent, *args], text=True)


def prepare_study(args, experiment="E03"):
    """Pin the executable and working-source identity before the first run."""
    args.output.mkdir(parents=True, exist_ok=True)
    pinned = args.output / "experiment"
    digest = hashlib.file_digest
    with args.executable.open("rb") as executable:
        identity = digest(executable, "sha256").hexdigest()
    metadata = args.output / "study.json"
    if metadata.exists():
        previous = json.loads(metadata.read_text())
        if previous.get("experiment", "E03") != experiment:
            raise ValueError("The output directory belongs to another experiment")
        if (previous["binary_sha256"], previous["ranks"], previous["threads"]) != (
            identity, args.ranks, args.threads
        ):
            raise ValueError("Continuation requires the same executable, ranks and threads")
        with pinned.open("rb") as executable:
            if digest(executable, "sha256").hexdigest() != identity:
                raise ValueError("The pinned executable changed")
        return pinned
    shutil.copy2(args.executable, pinned)
    shutil.copy2(__file__, args.output / "run_e03_study.py")
    (args.output / "source.patch").write_text(git("diff", "HEAD", "--binary"))
    (args.output / "source-status.txt").write_text(git("status", "--short"))
    metadata.write_text(json.dumps({"experiment": experiment,
                                    "git_head_at_launch": git("rev-parse", "HEAD").strip(),
                                    "binary_sha256": identity, "machine": platform.node(),
                                    "platform": platform.platform(), "ranks": args.ranks,
                                    "threads": args.threads}, indent=2) + "\n")
    return pinned


def run_case(args, executable, n, polarization, experiment="E03"):
    """Keep the input, combined log, rank memory peaks and any partial output."""
    run = args.output / f"n{n}-{polarization}"
    status = run / "status.json"
    if run.exists():
        if status.is_file() and json.loads(status.read_text())["exit_code"] == 0:
            return True
        print(f"Preserving incomplete or failed run {run}", flush=True)
        return False
    run.mkdir()
    inp = run / "input.yaml"
    plane = ", write_plane: true" if experiment == "E03" and n == LEVELS[0] else ""
    inp.write_text(f"schema_version: 1\ncase: {experiment}\nunits: {{system: normalized}}\n"
                   f"parameters: {{cells_per_axis: {n}, polarization: {polarization}{plane}}}\n"
                   "output: {directory: data}\n")
    env = {**os.environ, "OMP_NUM_THREADS": str(args.threads), "OMP_PROC_BIND": "false"}
    command = ["mpiexec", "-n", str(args.ranks), "--bind-to", "none", "sh", "-c",
               'exec "$1" -f %M -o "rss-${OMPI_COMM_WORLD_RANK}.kib" "$2" --input "$3"',
               "pulse-rank", shutil.which("time"), str(executable), str(inp)]
    start = time.monotonic()
    print(f"Running {experiment} {polarization} N={n} on {args.ranks} ranks", flush=True)
    with (run / "run.log").open("w") as log:
        code = subprocess.run(command, cwd=run, env=env, stdout=log,
                              stderr=subprocess.STDOUT).returncode
    peaks = [int(line) for p in sorted(run.glob("rss-*.kib"))
             for line in p.read_text().splitlines() if line.isdecimal()]
    record = {"exit_code": code, "wall_seconds": time.monotonic() - start,
              "rank_peak_rss_kib": peaks, "sum_rank_peaks_kib": sum(peaks)}
    status.write_text(json.dumps(record, indent=2) + "\n")
    print(f"{experiment} {polarization} N={n}: {record}", flush=True)
    return code == 0


def study_rows(root):
    """Derive the normalized values and paired orders from completed raw rows."""
    rows = []
    preceding = {}
    for n in LEVELS:
        for polarization in POLARIZATIONS:
            run = root / f"n{n}-{polarization}"
            if not completed_run(run):
                continue
            with (run / "data/measurements.csv").open() as stream:
                observed = list(csv.DictReader(stream))
            expected = 39 if polarization == "TEr" else 43
            if len(observed) != expected or not (run / "data/run.json").exists():
                raise ValueError(f"Incomplete E03 observations in {run}")
            for row in observed:
                scale = float(row["reference_scale"])
                value = float(row["value"]) / scale
                key = (polarization, row["quantity"], row["target_time"])
                earlier = preceding.get(key)
                order = ""
                if earlier and 2 * n == 3 * earlier[0] and min(earlier[1], value) > 1e-14:
                    order = math.log(earlier[1] / value) / math.log(n / earlier[0])
                preceding[key] = (n, value)
                rows.append({**row, "h": 12 / n, "normalized_value": value,
                             "observed_order": order})
    return rows


def completed_run(run):
    status = run / "status.json"
    return status.is_file() and json.loads(status.read_text())["exit_code"] == 0


def plot_errors(root, rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
    for col, polarization in enumerate(POLARIZATIONS):
        case = [r for r in rows if r["polarization"] == polarization]
        if not case:
            continue
        final = max(float(r["target_time"]) for r in case)
        for row_index, target in enumerate((0, final)):
            axis = axes[row_index, col]
            selected = [r for r in case if float(r["target_time"]) == target]
            for color, quantity in enumerate(QUANTITIES):
                values = [r for r in selected if r["quantity"] == quantity and r["normalized_value"] > 0]
                if values:
                    axis.loglog([r["h"] for r in values], [r["normalized_value"] for r in values],
                                "o-", color=f"C{color}", label=quantity)
            axis.set(title=f"{polarization}, {'startup' if target == 0 else 'final target'}",
                     xlabel="Cell width h", ylabel="Normalized norm")
            axis.grid(True, which="both", alpha=0.25)
            axis.legend(fontsize=8)
    figure.savefig(root / "e03-errors.png", dpi=160)
    plt.close(figure)


def plot_planes(root):
    import matplotlib.pyplot as plt
    import numpy as np

    n = LEVELS[0]
    for polarization in POLARIZATIONS:
        run = root / f"n{n}-{polarization}"
        path = run / "data/plane-z0.csv"
        if not completed_run(run) or not path.exists():
            continue
        data = np.genfromtxt(path, delimiter=",", names=True)
        steps = sorted(set(data["step"]))
        if len(steps) != 2 or len(data) != 2 * n * n:
            continue
        figure, axes = plt.subplots(2, 4, figsize=(12, 6), layout="constrained")
        figure.suptitle(f"E03 {polarization}, N={n}, z=0")
        for row_index, step in enumerate(steps):
            selected = data[data["step"] == step]
            for column, (prefix, field) in enumerate((("", "E"), ("", "B"),
                                                       ("error_", "E"), ("error_", "B"))):
                magnitude = np.sqrt(sum(selected[f"{prefix}{field}{axis}"] ** 2 for axis in "xyz"))
                axis = axes[row_index, column]
                plotted = axis.imshow(magnitude.reshape(n, n).T, origin="lower", extent=(-6, 6, -6, 6))
                label = f"error {field}" if prefix else field
                axis.set(title=f"|{label}|, t={selected['time'][0]:.3f}", xlabel="x", ylabel="y")
                figure.colorbar(plotted, ax=axis, shrink=0.75)
        figure.savefig(root / f"e03-plane-{polarization}.png", dpi=150)
        plt.close(figure)


def study_arguments(description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--levels", type=int, nargs="+", choices=LEVELS, default=LEVELS[:2])
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    if min(args.ranks, args.threads) < 1 or not shutil.which("time"):
        parser.error("Positive ranks and threads and GNU time on PATH are required")
    args.executable, args.output = args.executable.resolve(), args.output.resolve()
    return args


def main():
    args = study_arguments(__doc__)
    executable = prepare_study(args)
    success = True
    for n in sorted(set(args.levels)):
        for polarization in POLARIZATIONS:
            success = run_case(args, executable, n, polarization) and success
    rows = study_rows(args.output)
    if rows:
        analysis_source = Path(__file__).read_bytes()
        analysis_hash = hashlib.sha256(analysis_source).hexdigest()
        (args.output / f"analysis-{analysis_hash}.py").write_bytes(analysis_source)
        (args.output / "analysis-sha256.txt").write_text(analysis_hash + "\n")
        with (args.output / "convergence.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        plot_errors(args.output, rows)
        plot_planes(args.output)
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
