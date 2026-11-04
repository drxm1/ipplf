#!/usr/bin/env python3
"""Compare tilted-plane spectra, raw-point fits, controls and matched aligned rows.

The continuum reference and spatial-fit audit come from references.py. This
script supplies no staircase oracle and repeats no Fresnel or projection formula.
See README.md and ../references.org for assumptions and interpretation limits.
"""
from copy import copy
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_common import (argument_parser, execute_case, observed_order, plt, prepare_runs,
                         resolution_ticks, run_label, save_figure, update_status, write_csv)
from read_run import RunDataError, sha256_file


def close(first, second):
    return math.isclose(first, second, rel_tol=2e-12, abs_tol=0)


def group_number(value):
    return float(format(value, ".12g"))


def finite(value):
    return value if value is not None and math.isfinite(value) else None


def cell_shape(run):
    return tuple(value / run.spacing[0] for value in run.spacing)


def require_tilted(data):
    if any(item.run_m.schema_m != 2 for item in data):
        raise RunDataError("Rotated runs require schema 2; pass aligned data with --aligned-run")


def central_rows(data):
    rows = []
    for item in data:
        row, run = dict(item.central), item.run_m
        row.update(polarization=run.plane["polarization"], incidence=run.plane["incidence"],
                   interface_tilt_degrees=run.plane["interface_tilt_degrees"],
                   actual_angle_degrees=run.plane["actual_angle_degrees"],
                   cells_per_wavelength=run.plane["cells_per_wavelength"], h_x=run.spacing[0],
                   requested_face_position=run.plane["requested_face_position"],
                   sigma_k=run.plane["sigma_k"], tail_sigmas=run.plane["tail_sigmas"])
        rows.append(row)
    return rows


def roots_and_fits(item):
    rows = []
    for index, sample in enumerate(item.run_m.samples_m):
        row = {"bin": index, "omega_code": sample.omega_m, "valid": sample.valid_m, "reason": sample.reason_m}
        row.update(incident_angles(item, index))
        for offset, name in enumerate(("incident", "reflected", "transmitted")):
            row[name + "_normal_k"] = finite(sample.normal_wavenumbers_m[offset])
            row[name + "_normal_group_velocity"] = finite(sample.normal_group_velocities_m[offset])
            row[name + "_fit_residual"] = finite(sample.fit_residuals_m[offset])
        rows.append(row)
    return rows


def incident_angles(item, index):
    """Distinguish continuum reference incidence from the lattice wavevector direction.

    The angle of the recorded normal/tangent components is a coordinate definition;
    their dispersive difference follows from Fallahi Eq. (3.17), cited in ../references.org.
    This is the phase direction, not the group-velocity direction.
    """
    sample, p = item.run_m.samples_m[index], item.run_m.plane
    reference = item.reference_m.rows_m[index]["continuum_angle_degrees"]
    reference = math.copysign(reference, p["kappa"]) if reference is not None else None
    normal = finite(sample.normal_wavenumbers_m[0])
    direction = -1 if p["incidence"] == "vacuum" else 1
    phase = math.degrees(math.atan2(p["kappa"], direction * normal)) if normal is not None else None
    return {"continuum_reference_angle_degrees": reference, "incident_lattice_phase_angle_degrees": phase,
            "lattice_minus_continuum_angle_degrees": phase - reference if None not in (phase, reference) else None}


def audit_projection(item, tolerance):
    """Check re-projection agreement only; this tolerance is not a physical error bound."""
    rows = item.reference_m.record_m["extraction_residuals"]
    write_csv(item.directory_m / "projection-audit.csv", rows,
              list(rows[0]) if rows else ["omega_code", "reflection_residual", "transmission_residual"])
    maximum = max((row[key] for row in rows for key in ("reflection_residual", "transmission_residual")), default=None)
    if maximum is not None and (not math.isfinite(maximum) or maximum > tolerance):
        raise RunDataError(f"{item.identifier_m}: raw-point re-projection exceeds the audit tolerance")
    return {"run_id": item.identifier_m, "audited_valid_bins": len(rows),
            "maximum_coefficient_discrepancy": maximum, "tolerance": tolerance,
            "basis": "recorded directional roots, actual lifted component coordinates"}


def fit_figure(item, rows):
    figure, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    ordered = sorted(rows, key=lambda row: row["omega_code"])
    for name in ("incident", "reflected", "transmitted"):
        axes[0].plot([row["omega_code"] for row in ordered], [row[name + "_fit_residual"] for row in ordered],
                     ".-", label=name)
    axes[0].axhline(item.run_m.plane["fit_residual_limit"], color="black", linestyle=":", label="bin validity limit")
    audit = item.reference_m.record_m["extraction_residuals"]
    for name in ("reflection", "transmission"):
        axes[1].plot([row["omega_code"] for row in audit], [row[name + "_residual"] for row in audit], ".", label=name)
    axes[0].set(xlabel="angular frequency [code units]", ylabel="spatial RMS residual / incident amplitude")
    axes[0].set_yscale("symlog", linthresh=1e-12)
    axes[1].set(xlabel="angular frequency [code units]", ylabel="raw-fit to published coefficient discrepancy")
    for ax in axes:
        ax.legend(fontsize=8)
    figure.suptitle(run_label(item.run_m) + "; residuals are not physical error bounds")
    save_figure(figure, item.directory_m, "projection")


def record_geometry(item):
    run, plane = item.run_m, item.run_m.plane
    names = ("interface_tilt_degrees", "requested_angle_degrees", "actual_angle_degrees", "omega0", "kappa",
             "cells_per_wavelength", "requested_face_position", "primitive_cells_x", "primitive_cells_z",
             "normal_period", "tangent_period", "tail_sigmas", "sigma_k", "cutoff_fraction", "band_sigmas",
             "probe_cells", "probe_half_width_cells", "timing_samples", "separation_factor",
             "floquet_margin", "retained_floquet_margin", "fit_residual_limit")
    row = {"run_id": item.identifier_m, "polarization": plane["polarization"], **{key: plane[key] for key in names}}
    row.update(h_x=run.spacing[0], h_y=run.spacing[1], h_z=run.spacing[2], dt=run.dt,
               actual_shape_y=cell_shape(run)[1], actual_shape_z=cell_shape(run)[2],
               discarded_field_weight_fraction=plane.get("discarded_field_weight_fraction"))
    return row


def record_raw_files(item):
    files = ("run.json", "spectra.csv", "probe-spectra.csv", "modes.csv")
    return {"run_id": item.identifier_m, "run_directory": str(item.run_m.directory_m),
            "files_sha256": {name: sha256_file(item.run_m.directory_m / name) for name in files}}


def physical_group(item):
    run, p = item.run_m, item.run_m.plane
    values = (p["epsilon"], p["interface_tilt_degrees"], p["actual_angle_degrees"], p["requested_face_position"],
              p["requested_cells_per_wavelength"], run.spacing[0], *cell_shape(run), run.dt)
    return (p["polarization"], p["incidence"], p["primitive_cells_x"], p["primitive_cells_z"],
            *(group_number(value) for value in values))


def baseline_score(item):
    p = item.run_m.plane
    return tuple(p[key] for key in ("tail_sigmas", "separation_factor", "probe_cells",
                                   "probe_half_width_cells", "timing_samples")) + (max(p["gate_incident"], p["gate_other"]),)


def choose_baselines(data):
    groups = {}
    for item in data:
        key = physical_group(item)
        if key not in groups or baseline_score(item) < baseline_score(groups[key]):
            groups[key] = item
    return list(groups.values())


def matching_frequencies(first, second):
    """Pair only equal recorded frequencies; there is no spectral interpolation."""
    return [(a, b) for a in first.samples_m for b in second.samples_m
            if a.valid_m and b.valid_m and close(a.omega_m, b.omega_m)]


def same_frequency_grid(first, second):
    return len(first.frequencies) == len(second.frequencies) and all(
        close(a, b) for a, b in zip(first.frequencies, second.frequencies))


def central_difference(first, second, attribute, comparable):
    if not comparable or not first.central["comparison_valid"] or not second.central["comparison_valid"]:
        return None
    a, b = first.run_m, second.run_m
    return abs(getattr(a.samples_m[a.central_bin_m], attribute) - getattr(b.samples_m[b.central_bin_m], attribute))


def control_difference(first, second):
    a, b = first.run_m, second.run_m
    comparable = close(a.plane["kappa"], b.plane["kappa"]) and close(a.plane["omega0"], b.plane["omega0"])
    pairs = matching_frequencies(a, b) if comparable else []
    changed = [key for key in ("sigma_k", "tail_sigmas", "band_sigmas", "separation_factor", "timing_samples",
                               "probe_cells", "probe_half_width_cells", "normal_period", "tangent_mode", "gate_incident", "gate_other")
               if a.plane[key] != b.plane[key]]
    row = {"baseline": first.identifier_m, "control": second.identifier_m,
           "central_reference_parameters_match": comparable, "changed_settings": ";".join(changed),
           "first_omega0": a.plane["omega0"], "second_omega0": b.plane["omega0"],
           "first_kappa": a.plane["kappa"], "second_kappa": b.plane["kappa"],
           "first_sigma_k": a.plane["sigma_k"], "second_sigma_k": b.plane["sigma_k"],
           "common_valid_bins": len(pairs), "all_frequencies_match": same_frequency_grid(a, b)}
    for name, attr in (("r", "reflection_m"), ("t", "transmission_m")):
        row["maximum_" + name + "_difference"] = max((abs(getattr(x, attr) - getattr(y, attr)) for x, y in pairs), default=None)
        row["central_" + name + "_difference"] = central_difference(first, second, attr, comparable)
    return row


def control_rows(data, baselines):
    chosen = {physical_group(item): item for item in baselines}
    return [control_difference(chosen[physical_group(item)], item) for item in data
            if item is not chosen[physical_group(item)]]


def ladder_group(item):
    p, run = item.run_m.plane, item.run_m
    values = (p["epsilon"], p["interface_tilt_degrees"], p["actual_angle_degrees"], p["omega0"],
              p["requested_face_position"], *cell_shape(run), run.dt / run.spacing[0],
              p["tail_sigmas"], p["cutoff_fraction"], p["band_sigmas"])
    return (p["polarization"], p["incidence"], p["primitive_cells_x"], p["primitive_cells_z"],
            *(group_number(value) for value in values))


def ladders(data):
    groups = {}
    for item in data:
        groups.setdefault(ladder_group(item), []).append(item)
    return [sorted(items, key=lambda item: item.run_m.spacing[0], reverse=True) for items in groups.values()]


def order_pair(coarse, fine):
    row = {"coarse_run": coarse.identifier_m, "fine_run": fine.identifier_m,
           "coarse_h": coarse.run_m.spacing[0], "fine_h": fine.run_m.spacing[0]}
    for quantity in ("r", "t"):
        for kind in ("complex", "magnitude"):
            key = quantity + "_" + kind + "_error"
            errors = (coarse.central[key], fine.central[key])
            order = observed_order(*errors, row["coarse_h"], row["fine_h"]) if None not in errors else None
            row[quantity + "_" + kind + "_order"] = order
    return row


def draw_errors(ax, items, quantity, geometry):
    for kind, colour in (("complex", "tab:blue"), ("magnitude", "tab:orange")):
        style = "o-" if geometry == "tilted" else "s--"
        ax.loglog([item.run_m.plane["cells_per_wavelength"] for item in items],
                  [item.central[quantity + "_" + kind + "_error"] for item in items], style,
                  color=colour, markerfacecolor="none" if geometry == "aligned" else colour,
                  label=geometry + " " + kind)


def convergence_figure(items, aligned, output, index):
    if len(items) < 2:
        return
    figure, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    matches = [match for item in items if (match := aligned_match(item, aligned)) is not None]
    for ax, quantity in zip(axes, ("r", "t")):
        draw_errors(ax, items, quantity, "tilted")
        if matches:
            draw_errors(ax, matches, quantity, "aligned")
        ax.set(xlabel="actual dielectric cells per wavelength", ylabel=f"absolute {quantity} error")
        resolution_ticks(ax, [item.run_m.plane["cells_per_wavelength"] for item in items])
        ax.legend()
    p = items[0].run_m.plane
    figure.suptitle(f"{p['polarization']}: interface tilt {p['interface_tilt_degrees']:g}°, incidence {p['actual_angle_degrees']:g}°")
    save_figure(figure, output, f"convergence-{index:02d}")


def prepare_aligned(args):
    if not args.aligned_run:
        return []
    settings = copy(args)
    settings.runs, settings.scan_manifest = args.aligned_run, []
    settings.output = args.output / "aligned-controls"
    data = prepare_runs(settings, continuum_only=True)
    if any(item.run_m.schema_m != 1 for item in data):
        raise RunDataError("Aligned controls require schema 1")
    update_status(settings.output, status="completed", purpose="aligned control comparisons")
    return data


def aligned_candidate(first, second):
    a, b = first.run_m, second.run_m
    same = all(a.plane[key] == b.plane[key] for key in ("polarization", "incidence", "epsilon"))
    numbers = (a.plane["actual_angle_degrees"], *a.spacing, a.dt)
    other = (b.plane["actual_angle_degrees"], *b.spacing, b.dt)
    return same and all(close(x, y) for x, y in zip(numbers, other))


def aligned_match(item, data):
    candidates = [other for other in data if aligned_candidate(item, other)
                  and close(item.run_m.plane["omega0"], other.run_m.plane["omega0"])
                  and close(item.run_m.plane["kappa"], other.run_m.plane["kz"])]
    score = lambda other: (other.run_m.plane["separation_factor"], other.run_m.plane["timing_samples"])
    return min(candidates, key=score) if candidates else None


def aligned_pair(first, second):
    a, b = first.run_m, second.run_m
    match = close(a.plane["omega0"], b.plane["omega0"]) and close(a.plane["kappa"], b.plane["kz"])
    row = {"tilted_run": first.identifier_m, "aligned_run": second.identifier_m,
           "central_reference_parameters_match": match, "tilted_omega0": a.plane["omega0"],
           "tilted_central_valid": first.central["comparison_valid"], "aligned_central_valid": second.central["comparison_valid"],
           "aligned_omega0": b.plane["omega0"], "tilted_kappa": a.plane["kappa"], "aligned_kz": b.plane["kz"],
           "tilted_sigma_k": a.plane["sigma_k"], "aligned_sigma_k": b.plane["sigma_k"],
           "same_gaussian_normal_width": close(a.plane["sigma_k"], b.plane["sigma_k"]),
           "tilted_tail_sigmas": a.plane["tail_sigmas"], "aligned_tail_sigmas": b.plane["tail_sigmas"]}
    row.update(aligned_spectrum_settings(a, b))
    for quantity in ("r", "t"):
        row[quantity + "_tilted_complex_error"] = first.central[quantity + "_complex_error"]
        row[quantity + "_aligned_complex_error"] = second.central[quantity + "_complex_error"]
    return row


def aligned_spectrum_settings(tilted, aligned):
    """Expose envelope/truncation differences even when the central reference agrees."""
    keys = ("mode_cutoff_sigmas", "cutoff_fraction", "band_sigmas", "envelope_half_width")
    return {prefix + key: run.plane.get(key)
            for prefix, run in (("tilted_", tilted), ("aligned_", aligned)) for key in keys}


def final_records(args, data, selected, aligned, audits):
    controls = control_rows(data, selected)
    matched = [aligned_pair(a, b) for a in selected for b in aligned if aligned_candidate(a, b)]
    write_csv(args.output / "controls.csv", controls, list(controls[0]) if controls else ["baseline", "control"])
    write_csv(args.output / "aligned-comparisons.csv", matched, list(matched[0]) if matched else ["tilted_run", "aligned_run"])
    update_status(args.output, baseline_runs=[item.identifier_m for item in selected],
                  projection_audits=audits, control_comparisons=controls, aligned_comparisons=matched,
                  raw_sources=[record_raw_files(item) for item in data],
                  qualification="Measured comparisons; controls and residuals require interpretation, not automatic validation",
                  limitations=["No scalar represented staircase face", "No aligned lattice-interface oracle",
                               "Spectral width can differ even at matched central frequency", "Fit residual is not an error bound"])


def analyse(args, data):
    require_tilted(data)
    if not math.isfinite(args.audit_tolerance) or args.audit_tolerance <= 0:
        raise RunDataError("audit-tolerance must be finite and positive")
    audits = []
    for item in data:
        rows = roots_and_fits(item)
        write_csv(item.directory_m / "roots-and-fits.csv", rows)
        audits.append(audit_projection(item, args.audit_tolerance))
        fit_figure(item, rows)
    selected, aligned = choose_baselines(data), prepare_aligned(args)
    write_csv(args.output / "central.csv", central_rows(data))
    write_csv(args.output / "geometry.csv", [record_geometry(item) for item in data])
    groups = ladders(selected)
    orders = [order_pair(a, b) for items in groups for a, b in zip(items, items[1:])]
    write_csv(args.output / "orders.csv", orders, list(orders[0]) if orders else ["coarse_run", "fine_run"])
    for index, items in enumerate(groups):
        convergence_figure(items, aligned, args.output, index)
    final_records(args, data, selected, aligned, audits)


if __name__ == "__main__":
    parser = argument_parser(__doc__)
    parser.add_argument("--aligned-run", action="append", type=Path, default=[], help="Matching aligned control; repeat as needed")
    parser.add_argument("--audit-tolerance", type=float, default=1e-10, help="Maximum raw-fit to published coefficient discrepancy")
    raise SystemExit(execute_case(parser, analyse, continuum_only=True))
