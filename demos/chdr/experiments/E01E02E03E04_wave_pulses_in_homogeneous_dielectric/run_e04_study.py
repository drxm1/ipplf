#!/usr/bin/env python3
"""Run E04 and integrate the exit and return diagnostics at the exact eta endpoints."""

import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import run_e03_study
from run_e03_study import LEVELS, POLARIZATIONS, completed_run, prepare_study, run_case, study_arguments


FACE_NAMES = ("xminus", "xplus", "yminus", "yplus", "zminus", "zplus")
LATE_TARGETS = (10, 12, 16)
PRECURSOR = (4.5, 5)
RETURN = (5, 16)


def load_trace(path):
    """Read a complete ordered trace before interpolating across its samples."""
    trace = np.genfromtxt(path, delimiter=",", names=True)
    if trace.ndim != 1 or len(trace) < 2 or not np.all(np.diff(trace["eta"]) > 0):
        raise ValueError(f"Invalid E04 time trace {path}")
    if any(not np.all(np.isfinite(trace[name])) for name in trace.dtype.names):
        raise ValueError(f"Nonfinite E04 time trace {path}")
    return trace


def window(eta, values, bounds, stride=1):
    """Clip the piecewise-linear trace at the exact endpoints in E04.typ."""
    begin, end = bounds
    if eta[0] > begin or eta[-1] < end:
        raise ValueError(f"E04 trace does not bracket {bounds}")
    if stride > 1:
        indices = np.unique(np.r_[np.arange(0, len(eta), stride), len(eta) - 1])
        eta, values = eta[indices], values[indices]
    inside = (eta > begin) & (eta < end)
    return (np.r_[begin, eta[inside], end],
            np.r_[np.interp(begin, eta, values), values[inside], np.interp(end, eta, values)])


def integral(eta, values, bounds, stride=1):
    coordinates, samples = window(eta, values, bounds, stride)
    return float(np.trapezoid(samples, coordinates))


def summarize_energy(energy):
    """E04.meas01 at the late targets and the maximum of the interpolated trace."""
    eta, values = energy["eta"], energy["remaining_fraction"]
    _, late = window(eta, values, (LATE_TARGETS[0], LATE_TARGETS[-1]))
    return {"startup_energy_fraction": values[0], "late_max": max(late),
            **{f"F_remain_{t}": float(np.interp(t, eta, values)) for t in LATE_TARGETS}}


def summarize_return(monitor):
    """E04.meas02 per-face contributions, return and precursor mean rates."""
    eta = monitor["eta"]
    inward = sum(monitor[f"{name}_in"] for name in FACE_NAMES)
    outward = sum(monitor[f"{name}_out"] for name in FACE_NAMES)
    by_face = {f"F_in_{name}": integral(eta, monitor[f"{name}_in"], RETURN)
               for name in FACE_NAMES}
    amount = integral(eta, inward, RETURN)
    _, precursor = window(eta, inward, PRECURSOR)
    return {"F_in": amount, "F_out": integral(eta, outward, RETURN),
            "F_in_stride2": integral(eta, inward, RETURN, stride=2),
            "precursor_mean_power": integral(eta, inward, PRECURSOR) / (PRECURSOR[1] - PRECURSOR[0]),
            "precursor_peak_power": max(precursor), "return_mean_power": amount / (RETURN[1] - RETURN[0]),
            **by_face}


def checkpoints(directory, n, polarization):
    """Normalize E04's recorded constraint and field rows without changing their clocks."""
    with (directory / "measurements.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 85:
        raise ValueError(f"Incomplete E04 observations in {directory}")
    for row in rows:
        value, scale = float(row["value"]), float(row["reference_scale"])
        if not np.isfinite(value) or not np.isfinite(scale) or scale <= 0:
            raise ValueError(f"Invalid E04 checkpoint scale or value in {directory}")
        yield {"cells_per_axis": n, "polarization": polarization, **row,
               "normalized_value": value / scale}


def save_csv(path, rows):
    if rows:
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def mesh_differences(results, observations):
    """Report adjacent changes without assuming decay to zero."""
    metrics = ("F_in", "late_max", "F_remain_10", "F_remain_12", "F_remain_16",
               "precursor_mean_power", "return_mean_power")
    series = [(r["cells_per_axis"], (r["polarization"], q, "", ""), r[q])
              for r in results for q in metrics]
    series.extend((r["cells_per_axis"], (r["polarization"], r["quantity"], r["region"], r["target_eta"]),
                   r["normalized_value"]) for r in observations)
    preceding, rows = {}, []
    for n, key, value in sorted(series):
        previous = preceding.get(key)
        if previous and 2 * n == 3 * previous[0]:
            pol, quantity, region, eta = key
            rows.append({"polarization": pol, "quantity": quantity, "region": region,
                         "target_eta": eta, "coarse_cells": previous[0], "fine_cells": n,
                         "coarse_value": previous[1], "fine_value": value,
                         "difference": value - previous[1]})
        preceding[key] = (n, value)
    return rows


def plot_meshes(root, results, observations):
    figure, axes = plt.subplots(3, 2, figsize=(10, 9), layout="constrained")
    for col, pol in enumerate(POLARIZATIONS):
        case = [r for r in results if r["polarization"] == pol]
        for key in ("F_remain_10", "F_remain_12", "F_remain_16"):
            axes[0, col].loglog([12 / r["cells_per_axis"] for r in case], [r[key] for r in case],
                                "o-", label=key)
        for key in ("precursor_mean_power", "return_mean_power"):
            axes[1, col].loglog([12 / r["cells_per_axis"] for r in case], [r[key] for r in case],
                                "o-", label=key)
        for region in ("interior", "boundary"):
            for quantity in ("Q_L2", "G_L2"):
                rows = [r for r in observations if r["polarization"] == pol and r["region"] == region
                        and r["quantity"] == quantity and float(r["target_eta"]) == 16]
                axes[2, col].loglog([12 / r["cells_per_axis"] for r in rows],
                                    [r["normalized_value"] for r in rows], "o-",
                                    label=f"{quantity} {region}")
        for axis, title in zip(axes[:, col], ("Late energy", "Mean inward power", "Constraints at eta=16")):
            axis.set(title=f"{pol}, {title}", xlabel="Cell width h", ylabel="Normalized value")
            axis.grid(True, which="both", alpha=0.25)
            axis.legend(fontsize=8)
    figure.savefig(root / "e04-meshes.png", dpi=160)
    plt.close(figure)


def plot_traces(root, traces):
    figure, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
    for n, pol, energy, monitor in traces:
        col = POLARIZATIONS.index(pol)
        axes[0, col].semilogy(energy["eta"], energy["remaining_fraction"], label=f"N={n}")
        inward = sum(monitor[f"{face}_in"] for face in FACE_NAMES)
        axes[1, col].semilogy(monitor["eta"], inward, label=f"N={n}")
    for col, pol in enumerate(POLARIZATIONS):
        for row, ylabel in enumerate(("Remaining energy fraction", "Normalized inward power")):
            axis = axes[row, col]
            axis.set(title=pol, xlabel="eta = vt/w", ylabel=ylabel)
            axis.grid(True, which="both", alpha=0.25)
            axis.legend(fontsize=8)
        for event in (3, 6 * np.sqrt(2) - 3, 6 * np.sqrt(3) - 3, 6 * np.sqrt(3) - 1):
            axes[0, col].axvline(event, color="0.65", linewidth=0.7)
        axes[1, col].axvspan(*PRECURSOR, color="0.8", alpha=0.4)
        axes[1, col].axvline(6, color="0.65", linewidth=0.7)
    figure.savefig(root / "e04-traces.png", dpi=160)
    plt.close(figure)


def analyze(root):
    results, observations, traces = [], [], []
    for n in LEVELS:
        for polarization in POLARIZATIONS:
            run = root / f"n{n}-{polarization}"
            if not completed_run(run):
                continue
            directory = run / "data"
            energy, monitor = (load_trace(directory / name) for name in ("energy.csv", "monitor.csv"))
            record = json.loads((directory / "run.json").read_text())
            steps = record["effective"]["protocol_steps"]
            if len(energy) != steps + 1 or int(energy["step"][-1]) != steps:
                raise ValueError(f"Incomplete E04 energy trace in {directory}")
            results.append({"cells_per_axis": n, "polarization": polarization,
                            "U0h": record["effective"]["initial_reference_energy"],
                            **summarize_energy(energy), **summarize_return(monitor)})
            observations.extend(checkpoints(directory, n, polarization))
            traces.append((n, polarization, energy, monitor))
    save_csv(root / "results.csv", results)
    save_csv(root / "checkpoints.csv", observations)
    save_csv(root / "mesh-differences.csv", mesh_differences(results, observations))
    if traces:
        plot_traces(root, traces)
        plot_meshes(root, results, observations)


def archive_analysis_sources(root):
    entries = {}
    for source in (Path(__file__), Path(run_e03_study.__file__)):
        content = source.read_bytes()
        digest = hashlib.sha256(content).hexdigest()
        archived = f"analysis-{digest}-{source.name}"
        (root / archived).write_bytes(content)
        entries[source.name] = {"sha256": digest, "file": archived}
    (root / "analysis-sources.json").write_text(json.dumps(entries, indent=2) + "\n")


def main():
    args = study_arguments(__doc__)
    executable = prepare_study(args, "E04")
    archive_analysis_sources(args.output)
    success = True
    for n in sorted(set(args.levels)):
        for polarization in POLARIZATIONS:
            success = run_case(args, executable, n, polarization, "E04") and success
    analyze(args.output)
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
