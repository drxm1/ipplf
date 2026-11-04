#!/usr/bin/env python3
"""Shared comparison tables and figures for recorded planar-interface experiments.

The C++ reference utility supplies coefficients, propagation roots and the flux
weight. Error and power definitions follow references.org, particularly its
continuum-power, phase-origin and reference-case sections. These routines do not
infer a physical angle or geometry from an input filename.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import importlib.util
import json
import math
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from read_run import PlaneRun, RunDataError, load_run, sha256_file, source_identity
from references import ReferenceTable, coefficient, evaluate


@dataclass(frozen=True)
class Comparison:
    """One measured run, its exact reference invocation and derived comparison rows."""

    identifier_m: str
    run_m: PlaneRun
    reference_m: ReferenceTable
    rows_m: tuple[dict, ...]
    directory_m: Path

    @property
    def central(self) -> dict:
        return self.rows_m[self.run_m.central_bin_m]


def argument_parser(description: str) -> argparse.ArgumentParser:
    """Use the same explicit input and fresh-output policy in each figure command."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("runs", nargs="*", type=Path, help="Completed run directories")
    parser.add_argument("--scan-manifest", action="append", type=Path, default=[],
                        help="Completed scan manifest; may be repeated for control scans")
    parser.add_argument("--reference-binary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path, help="New analysis directory")
    parser.add_argument("--phase-floor", type=float, default=1e-10,
                        help="Minimum coefficient magnitude for reporting phase")
    return parser


def write_csv(path: Path, rows: list[dict] | tuple[dict, ...], names=None) -> None:
    """Retain the numeric table used for each figure, including invalid-bin reasons."""
    names = names or list(dict.fromkeys(key for row in rows for key in row))
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=names)
        writer.writeheader()
        writer.writerows(rows)


def _parts(value: complex | None, name: str) -> dict:
    if value is None or not math.isfinite(value.real) or not math.isfinite(value.imag):
        return {name + "_real": None, name + "_imag": None}
    return {name + "_real": value.real, name + "_imag": value.imag}


def _phase(value: complex | None, floor: float) -> float | None:
    return None if value is None or abs(value) <= floor else math.atan2(value.imag, value.real)


def _coefficient_comparison(measured, expected, name: str, floor: float) -> dict:
    """Separate complex, magnitude and wrapped phase errors; avoid phase at zeros."""
    result = _parts(measured, name) | _parts(expected, name + "_continuum")
    result[name + "_phase"] = _phase(measured, floor)
    result[name + "_continuum_phase"] = _phase(expected, floor)
    usable = measured is not None and expected is not None
    result[name + "_complex_error"] = abs(measured - expected) if usable else None
    result[name + "_magnitude_error"] = abs(abs(measured) - abs(expected)) if usable else None
    phase_usable = usable and abs(measured) > floor and abs(expected) > floor
    result[name + "_phase_error"] = _phase(measured * expected.conjugate(), 0) if phase_usable else None
    return result


def _power(sample, reference: dict, usable: bool) -> dict:
    """Apply the C++ flux weight; this is a continuum flux estimate, not an invariant.

    See references.org, Power fractions (Ellingson §§5.6–5.7 and Eq. 3.98).
    """
    weight = reference["continuum_transmission_weight"]
    usable = usable and weight is not None
    reflected = abs(sample.reflection_m)**2 if usable else None
    transmitted = weight * abs(sample.transmission_m)**2 if usable else None
    return {"R": reflected, "T": transmitted,
            "R_plus_T_minus_one": reflected + transmitted - 1 if usable else None,
            "R_continuum": reference["continuum_R"],
            "T_continuum": reference["continuum_T"],
            "transmission_weight": weight}


def _sample_row(identifier: str, run: PlaneRun, index: int, reference: dict,
                phase_floor: float) -> dict:
    sample = run.samples_m[index]
    usable = sample.valid_m and reference["continuum_valid"]
    reason = ";".join(filter(None, (sample.reason_m, reference["continuum_reason"])))
    row = {"run_id": identifier, "bin": index, "central": index == run.central_bin_m,
           "omega_code": sample.omega_m, "angle_degrees": reference["continuum_angle_degrees"],
           "measured_valid": sample.valid_m, "comparison_valid": usable, "reason": reason,
           "lattice_valid": reference["lattice_valid"], "lattice_reason": reference["lattice_reason"]}
    for name, measured in (("r", sample.reflection_m), ("t", sample.transmission_m)):
        expected = coefficient(reference, "continuum", name)
        row.update(_coefficient_comparison(measured if usable else None, expected, name, phase_floor))
        lattice = coefficient(reference, "lattice", name)
        row.update(_parts(lattice, name + "_lattice"))
        row[name + "_lattice_error"] = abs(measured - lattice) if sample.valid_m and lattice is not None else None
    return row | _power(sample, reference, usable)


def _run_details(identifier: str, run: PlaneRun, directory: Path) -> dict:
    plane = run.plane
    return {"run_id": identifier, "run_directory": str(run.directory_m),
            "run_sha256": sha256_file(run.directory_m / "run.json"),
            "spectra_sha256": sha256_file(run.directory_m / "spectra.csv"),
            "build": run.record_m["build"], "effective_plane": plane,
            "spacing": run.spacing, "dt": run.dt, "output_factors": run.record_m["output_factors"],
            "reference_directory": str(directory / "reference"),
            "comparison_file": str(directory / "comparison.csv"),
            "valid_bin_count": sum(sample.valid_m for sample in run.samples_m),
            "bin_count": len(run.samples_m)}


def _prepare_run(path: Path, index: int, args, continuum_only: bool) -> Comparison:
    run = load_run(path)
    identifier = f"run-{index:03d}"
    directory = args.output / identifier
    directory.mkdir()
    table = evaluate(run, args.reference_binary, directory / "reference", continuum_only)
    rows = tuple(_sample_row(identifier, run, i, reference, args.phase_floor)
                 for i, reference in enumerate(table.rows_m))
    write_csv(directory / "comparison.csv", rows)
    result = Comparison(identifier, run, table, rows, directory)
    band_figure(result)
    return result


def _resolve_runs(args) -> list[Path]:
    paths = list(args.runs)
    for path in args.scan_manifest:
        manifest = json.loads(path.read_text())
        if manifest.get("scan_record_version") != 1 or manifest.get("status") != "completed":
            raise RunDataError(f"{path}: a completed scan is required; use explicit run paths for a deliberate subset")
        for attempt in manifest["attempts"]:
            if attempt.get("status") != "completed" or "output_directory" not in attempt:
                raise RunDataError(f"{path}: incomplete attempt in a completed scan")
            paths.append(Path(attempt["output_directory"]))
    if not paths:
        raise RunDataError("Supply run directories or a completed --scan-manifest")
    return list(dict.fromkeys(path.resolve() for path in paths))


def _prepare_output(args) -> list[Path]:
    if not math.isfinite(args.phase_floor) or args.phase_floor <= 0:
        raise RunDataError("phase-floor must be finite and positive")
    paths = _resolve_runs(args)
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    args.output_owned = True
    return paths


def _analysis_record(args, continuum_only: bool) -> dict:
    return {"schema_version": 1, "analysis_source": source_identity(),
            "phase_floor": args.phase_floor, "continuum_only": continuum_only, "runs": [],
            "scan_manifests": [{"path": str(path.resolve()), "sha256": sha256_file(path)}
                               for path in args.scan_manifest]}


def prepare_runs(args, continuum_only: bool = False) -> list[Comparison]:
    """Preserve partial results and a failed status if any input cannot be analysed."""
    paths = _prepare_output(args)
    manifest = _analysis_record(args, continuum_only)
    results = []
    try:
        for index, path in enumerate(paths):
            result = _prepare_run(path, index, args, continuum_only)
            results.append(result)
            manifest["runs"].append(_run_details(result.identifier_m, result.run_m, result.directory_m))
        manifest["status"] = "prepared"
        return results
    except Exception as error:
        manifest.update(status="failed", error=str(error))
        raise
    finally:
        (args.output / "analysis.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")


def plot_style() -> None:
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25,
                         "figure.dpi": 120, "savefig.dpi": 170, "axes.formatter.useoffset": False})


def save_figure(figure, directory: Path, name: str) -> None:
    """Save the shareable figure and vector version beside their numeric table."""
    figure.savefig(directory / (name + ".png"), bbox_inches="tight")
    figure.savefig(directory / (name + ".pdf"), bbox_inches="tight")
    plt.close(figure)


def run_label(run: PlaneRun) -> str:
    return (f"{run.plane['polarization']}, "
            f"$\\theta_0={run.plane['actual_angle_degrees']:.4g}^\\circ$, "
            f"$N_d={run.plane['cells_per_wavelength']:.5g}$")


def _column(rows, name: str) -> np.ndarray:
    return np.array([np.nan if row.get(name) is None else row[name] for row in rows], dtype=float)


def _band_coefficients(ax, rows, omega, quantity: str, phase: bool) -> None:
    if phase:
        expected = _column(rows, quantity + "_continuum_phase")
        measured = expected + _column(rows, quantity + "_phase_error")
        ylabel = f"phase({quantity}) [rad]"
    else:
        measured = np.hypot(_column(rows, quantity + "_real"), _column(rows, quantity + "_imag"))
        expected = np.hypot(_column(rows, quantity + "_continuum_real"), _column(rows, quantity + "_continuum_imag"))
        ylabel = f"|{quantity}|"
    valid = _column(rows, "comparison_valid").astype(bool)
    ax.plot(omega, np.where(valid, measured, np.nan), ".-", label="measured")
    ax.plot(omega, expected, "--", label="continuum reference")
    ax.set(xlabel="angular frequency [code units]", ylabel=ylabel)
    if phase:
        phase_limits(ax, measured, expected)
    ax.legend(fontsize=8)


def phase_limits(ax, *values) -> None:
    """Keep roundoff-scale phase scatter visible in context, without changing saved values."""
    finite = np.concatenate([np.asarray(value, dtype=float) for value in values])
    finite = finite[np.isfinite(finite)]
    if len(finite) and np.ptp(finite) < 0.1:
        centre = (finite.min() + finite.max()) / 2
        ax.set_ylim(centre - 0.05, centre + 0.05)


def resolution_ticks(ax, values) -> None:
    """Label the actual sampled resolutions directly on a logarithmic axis."""
    samples = sorted(set(values))
    ax.set_xticks(samples, labels=[f"{value:.5g}" for value in samples])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())


def band_figure(data: Comparison) -> None:
    """Show only valid measured bins; retained tables also expose every rejected bin."""
    rows = sorted(data.rows_m, key=lambda row: row["omega_code"])
    omega = _column(rows, "omega_code")
    figure, axes = plt.subplots(2, 2, figsize=(10, 6), constrained_layout=True)
    for column, quantity in enumerate(("r", "t")):
        _band_coefficients(axes[0, column], rows, omega, quantity, False)
        _band_coefficients(axes[1, column], rows, omega, quantity, True)
    omitted = sum(not row["comparison_valid"] for row in rows)
    figure.suptitle(run_label(data.run_m) + f"; {omitted}/{len(rows)} bins excluded")
    save_figure(figure, data.directory_m, "band")


def central_rows(data: list[Comparison]) -> list[dict]:
    """Collect actual central-bin quantities; an invalid central result stops aggregation."""
    rows = []
    for item in data:
        row, plane = dict(item.central), item.run_m.plane
        if not row["comparison_valid"]:
            raise RunDataError(f"{item.identifier_m}: unusable central bin: {row['reason']}")
        row.update(polarization=plane["polarization"], incidence=plane["incidence"],
                   epsilon=plane["epsilon"], requested_angle_degrees=plane["requested_angle_degrees"],
                   actual_angle_degrees=plane["actual_angle_degrees"],
                   cells_per_wavelength=plane["cells_per_wavelength"], h_x=item.run_m.spacing[0],
                   face_position=plane["face_position"], separation_factor=plane["separation_factor"],
                   timing_samples=plane["timing_samples"], dt=item.run_m.dt)
        rows.append(row)
    return rows


def observed_order(coarse_error: float, fine_error: float, coarse_h: float, fine_h: float):
    """Infer a two-level exponent from e=C*h**p; this is not an asymptotic-order claim.

    See references.org, Normal incidence and mesh refinement, and the qualifications
    discussed by Werner et al. (2013), sections 5–6.
    """
    if min(coarse_error, fine_error, coarse_h, fine_h) <= 0 or coarse_h == fine_h:
        return None
    return math.log(coarse_error / fine_error) / math.log(coarse_h / fine_h)


def matched_controls(data: list[Comparison], output: Path) -> list[dict]:
    """Compare rows with the same physics/spectrum but differing layout or gate controls."""
    comparisons = []
    for index, first in enumerate(data):
        for second in data[index + 1:]:
            if not _same_measurement(first.run_m, second.run_m):
                continue
            pairs = [(a, b) for a, b in zip(first.run_m.samples_m, second.run_m.samples_m)
                     if a.valid_m and b.valid_m]
            comparisons.append({"first": first.identifier_m, "second": second.identifier_m,
                "common_valid_bins": len(pairs),
                "max_r_difference": max((abs(a.reflection_m - b.reflection_m) for a, b in pairs), default=None),
                "max_t_difference": max((abs(a.transmission_m - b.transmission_m) for a, b in pairs), default=None)})
    write_csv(output / "matched-controls.csv", comparisons,
              ["first", "second", "common_valid_bins", "max_r_difference", "max_t_difference"])
    return comparisons


def _same_measurement(first: PlaneRun, second: PlaneRun) -> bool:
    keys = ("polarization", "incidence", "epsilon", "actual_angle_degrees", "kz", "omega0",
            "sigma_k", "tail_sigmas", "cutoff_fraction")
    if any(first.plane[key] != second.plane[key] for key in keys):
        return False
    if first.spacing != second.spacing or first.dt != second.dt:
        return False
    first_fraction = first.plane["face_position"] / first.spacing[0] % 1
    second_fraction = second.plane["face_position"] / second.spacing[0] % 1
    if not math.isclose(first_fraction, second_fraction, rel_tol=0, abs_tol=1e-12):
        return False
    return len(first.frequencies) == len(second.frequencies) and all(
        math.isclose(a, b, rel_tol=2e-13, abs_tol=0) for a, b in zip(first.frequencies, second.frequencies))


def baselines(data: list[Comparison], key) -> list[Comparison]:
    """Select the least enlarged sampling/layout row, retaining other rows as controls."""
    chosen = {}
    for item in data:
        plane = item.run_m.plane
        score = (plane["separation_factor"], plane["timing_samples"],
                 max(plane["gate_incident"], plane["gate_other"]))
        group = key(item)
        if group not in chosen or score < chosen[group][0]:
            chosen[group] = (score, item)
    return [item for _, item in chosen.values()]


def matching_cell_policy(first: PlaneRun, second: PlaneRun) -> bool:
    """Compare shape and dt/h_x, allowing floating representation roundoff only."""
    shape = zip((value / first.spacing[0] for value in first.spacing),
                (value / second.spacing[0] for value in second.spacing))
    return all(math.isclose(a, b, rel_tol=2e-12) for a, b in shape) and math.isclose(
        first.dt / first.spacing[0], second.dt / second.spacing[0], rel_tol=2e-12)


def update_status(output: Path, **values) -> None:
    path = output / "analysis.json"
    if path.exists():
        record = json.loads(path.read_text())
        record.update(values)
        path.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")


def execute_case(parser, analyse, continuum_only: bool = False) -> int:
    """Run a complete analysis, marking failed aggregation instead of claiming success."""
    args = parser.parse_args()
    args.output_owned = False
    try:
        plot_style()
        data = prepare_runs(args, continuum_only)
        analyse(args, data)
        script = Path(sys.argv[0]).resolve()
        hashes = {path.name: sha256_file(path) for path in (Path(__file__), script)}
        update_status(args.output, status="completed", analysis_files=hashes)
        print(args.output / "analysis.json")
        return 0
    except (OSError, ValueError, RuntimeError) as error:
        if args.output_owned:
            update_status(args.output, status="failed", error=str(error))
        print(f"analysis: {error}", file=sys.stderr)
        return 2


def frame_reader():
    """Reuse the existing VTK component/time-label reader."""
    path = Path(__file__).resolve().parents[2] / "apps/interface-movies/render_movie.py"
    spec = importlib.util.spec_from_file_location("interface_frame_reader", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def frame_sources(directory: Path, name: str) -> list[dict]:
    """Identify the exact wrapper and every VTK piece used for a retained profile."""
    wrapper = directory / name
    pieces = [directory / node.attrib["Source"] for node in ET.parse(wrapper).findall(".//Piece")]
    paths = [wrapper, *pieces]
    if any(not path.resolve().is_relative_to(directory.resolve()) for path in paths):
        raise RunDataError("Frame piece is outside its recorded run directory")
    return [{"path": str(path.relative_to(directory)), "sha256": sha256_file(path)} for path in paths]


def profile_data(data: Comparison) -> list[dict]:
    """Select three saved times before visible normal-edge contact and retain the criterion."""
    series = data.run_m.directory_m / "wave.pvti.series"
    if not series.exists():
        return []
    frames = json.loads(series.read_text())["files"]
    if len(frames) < 3:
        raise RunDataError(f"{data.identifier_m}: at least three saved frames are needed for profiles")
    reader = frame_reader()
    field = "E_raw" if data.run_m.plane["polarization"] == "TE" else "B_raw"
    profiles = [_one_profile(data, frame, reader, field) for frame in frames]
    selected, selection = _select_profile_times(profiles)
    (data.directory_m / "profile-selection.json").write_text(json.dumps(selection, indent=2) + "\n")
    return selected


def _select_profile_times(profiles: list[dict]) -> tuple[list[dict], dict]:
    peak = max(profile["peak"] for profile in profiles)
    if peak <= 0:
        raise RunDataError("The saved profile channel is zero in every frame")
    first_contact = next((i for i, profile in enumerate(profiles)
                          if profile["normal_edge_peak"] > 0.01 * peak), len(profiles))
    if first_contact < 3:
        raise RunDataError("Need three frames before visible normal-edge contact; use a finer field cadence")
    last = first_contact - 1
    midpoint = (profiles[0]["time"] + profiles[last]["time"]) / 2
    middle = min(range(1, last), key=lambda i: abs(profiles[i]["time"] - midpoint))
    selected = [profiles[i] for i in (0, middle, last)]
    samples = [{key: profile[key] for key in ("frame", "time", "peak", "normal_edge_peak", "sources")}
               for profile in profiles]
    record = {"purpose": "display selection, not a boundary-error bound", "edge_cells": 8,
              "threshold_fraction_of_global_peak": 0.01, "global_peak": peak, "samples": samples,
              "selected_frames": [item["frame"] for item in selected],
              "first_visible_contact": profiles[first_contact]["frame"] if first_contact < len(profiles) else None}
    return selected, record


def _one_profile(data: Comparison, frame: dict, reader, field: str) -> dict:
    values, grid = reader.read_cells(data.run_m.directory_m / frame["name"], field, "zx")
    if not np.isfinite(values[:, :, 1]).all():
        raise RunDataError(f"Nonfinite field profile in {frame['name']}")
    spacing, origin = grid.GetSpacing(), grid.GetOrigin()
    offsets = grid.GetFieldData().GetArray("component_offsets")
    offset = offsets.GetTuple(1 if field == "E_raw" else 4)[0]
    x = origin[0] + (np.arange(values.shape[0]) + 0.5 + offset) * spacing[0]
    factor = data.run_m.record_m["output_factors"]
    amplitude_scale = factor["electric_field" if field == "E_raw" else "magnetic_field"]
    field_values = values[:, :, 1] / amplitude_scale
    edge_peak = max(float(np.abs(field_values[:8]).max()), float(np.abs(field_values[-8:]).max()))
    return {"frame": frame["name"], "time": reader.field_time(grid, field) / factor["time"],
            "x": x / factor["length"], "field": field_values[:, 0],
            "quantity": "E_y" if field == "E_raw" else "B_y_raw",
            "sources": frame_sources(data.run_m.directory_m, frame["name"]),
            "peak": float(np.abs(field_values).max()), "normal_edge_peak": edge_peak}


def pulse_profiles(data: Comparison) -> bool:
    """Save actual samples, avoiding an analytic decomposition of a mixed field frame."""
    profiles = profile_data(data)
    if not profiles:
        return False
    figure, axes = plt.subplots(3, 1, figsize=(9, 7), sharex=True, constrained_layout=True)
    for ax, profile in zip(axes, profiles):
        ax.plot(profile["x"], profile["field"], linewidth=1)
        ax.axvline(data.run_m.plane["face_position"], color="black", linestyle=":", label="requested face")
        ax.set(ylabel=profile["quantity"], title=f"t={profile['time']:.7g} [code units]")
    axes[-1].set_xlabel("x [code units]")
    figure.suptitle(run_label(data.run_m))
    _save_profile_tables(data, profiles)
    save_figure(figure, data.directory_m, "profiles")
    return True


def _save_profile_tables(data: Comparison, profiles: list[dict]) -> None:
    rows = [{"frame": profile["frame"], "time": profile["time"], "x": x,
             "field": value, "quantity": profile["quantity"]}
            for profile in profiles for x, value in zip(profile["x"], profile["field"])]
    write_csv(data.directory_m / "profiles.csv", rows)
    sources = [{"frame": item["frame"], "time": item["time"], "quantity": item["quantity"],
                "files": item["sources"]} for item in profiles]
    (data.directory_m / "profile-sources.json").write_text(json.dumps(sources, indent=2) + "\n")
