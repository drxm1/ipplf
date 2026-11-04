#!/usr/bin/env python3
"""Compare central plane coefficients and continuum flux estimates over incidence angle.

The C++ utility supplies every continuum reference, including the dense angle
curves. See ../references.org for fixed-H TM signs and Brewster conventions.
"""
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_common import (argument_parser, baselines, central_rows, execute_case,
                         matched_controls, np, phase_limits, plt, save_figure, update_status, write_csv)
from read_run import RunDataError
from references import angle_curve, coefficient


def require_same_scan(data):
    first = data[0].run_m
    for item in data:
        run = item.run_m
        for key in ("epsilon", "incidence", "requested_cells_per_wavelength"):
            if run.plane[key] != first.plane[key]:
                raise RunDataError("Angle rows must match media, incidence side and requested resolution")
        if run.spacing != first.spacing or not math.isclose(run.dt, first.dt, rel_tol=2e-12):
            raise RunDataError("Angle rows must retain one cell shape, cell size and time step")


def groups(rows):
    result = []
    for polarization in ("TE", "TM"):
        selected = sorted((row for row in rows if row["polarization"] == polarization),
                          key=lambda row: row["actual_angle_degrees"])
        if selected:
            result.append((polarization, selected))
    return result


def angle_samples(data):
    """Cover the measured signs, retaining the existing positive-only curve range."""
    measured = {item.run_m.plane["actual_angle_degrees"] for item in data}
    start = min(-80, min(measured)) if min(measured) < 0 else 0
    stop = max(80, max(measured)) if min(measured) >= 0 or max(measured) > 0 else 0
    return sorted(set(np.linspace(start, stop, 321)) | measured)


def signed_curve_rows(table, angles):
    """Retain signed plot coordinates beside the nonnegative reference angle.

    The Fresnel roots and fixed-component coefficients depend on the square of
    the tangential wavenumber. Their invariance under its sign reversal follows
    from the TE/TM derivations in ../references.org (Ellingson, sections 5.6–5.7).
    """
    return tuple(row | {"plot_angle_degrees": angle} for row, angle in zip(table.rows_m, angles))


def curve_provenance(item, table, angles):
    return {"polarization": item.run_m.plane["polarization"], "anchor_run": item.identifier_m,
            "omega_code": item.run_m.plane["omega0"], "directory": str(table.directory_m),
            "reference_record": str(table.directory_m / "reference.json"),
            "plot_angles_degrees": angles, "reference_angles_degrees": [abs(angle) for angle in angles],
            "angle_mapping": "signed plot angle maps to absolute continuum reference angle"}


def dense_curves(data, args):
    """Preserve one explicit continuum calculation per polarization and its anchor run."""
    anchors = {}
    for item in data:
        anchors.setdefault(item.run_m.plane["polarization"], item)
    angles = angle_samples(data)
    curves, provenance = {}, []
    for polarization, item in anchors.items():
        table = angle_curve(item.run_m, args.reference_binary, [abs(angle) for angle in angles],
                            args.output / ("continuum-" + polarization))
        curves[polarization] = signed_curve_rows(table, angles)
        provenance.append(curve_provenance(item, table, angles))
    update_status(args.output, reference_curve="dense C++ continuum angle calculation", reference_curves=provenance)
    return curves


def curve_values(rows, quantity, phase_floor=None):
    """Use the recorded reference validity and omit undefined phase near zeros."""
    angles, values = [], []
    for row in rows:
        value = coefficient(row, "continuum", quantity)
        angles.append(row["plot_angle_degrees"])
        if not row["continuum_valid"] or value is None or (phase_floor is not None and abs(value) <= phase_floor):
            values.append(np.nan)
        else:
            values.append(abs(value) if phase_floor is None else math.atan2(value.imag, value.real))
    if phase_floor is not None:
        for index in range(len(values) - 1, 0, -1):
            if abs(values[index] - values[index - 1]) > math.pi / 2:
                values[index] = np.nan
    return angles, values


def amplitude_figure(rows, curves, output):
    figure, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for polarization, selected in groups(rows):
        angle = [row["actual_angle_degrees"] for row in selected]
        for ax, quantity in zip(axes, ("r", "t")):
            values = [math.hypot(row[quantity + "_real"], row[quantity + "_imag"]) for row in selected]
            reference_angle, expected = curve_values(curves[polarization], quantity)
            line, = ax.plot(angle, values, "o", label=polarization + " measured")
            ax.plot(reference_angle, expected, "--", color=line.get_color(), label=polarization + " Fresnel")
            ax.set(xlabel="actual central incidence angle [degrees]", ylabel=f"|{quantity}|")
            ax.legend(fontsize=8)
    figure.suptitle("Vacuum–dielectric plane interface: central amplitudes")
    save_figure(figure, output, "angle-amplitudes")


def power_figure(rows, curves, output):
    figure, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for polarization, selected in groups(rows):
        angle = [row["actual_angle_degrees"] for row in selected]
        for quantity, marker in (("R", "o"), ("T", "s")):
            line, = axes[0].plot(angle, [row[quantity] for row in selected], marker,
                                label=polarization + " " + quantity + " measured")
            reference = curves[polarization]
            axes[0].plot([row["plot_angle_degrees"] for row in reference],
                         [row["continuum_" + quantity] for row in reference], "--",
                         color=line.get_color(), label=polarization + " " + quantity + " Fresnel")
        axes[1].plot(angle, [row["R_plus_T_minus_one"] for row in selected], "o-", label=polarization)
    axes[0].set(xlabel="actual central incidence angle [degrees]", ylabel="continuum flux estimate")
    axes[1].set(xlabel="actual central incidence angle [degrees]", ylabel="R + T − 1")
    axes[1].axhline(0, color="black", linestyle=":", linewidth=0.8)
    for ax in axes:
        ax.legend(fontsize=8)
    save_figure(figure, output, "angle-power")


def phase_figure(rows, curves, output, phase_floor):
    figure, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for polarization, selected in groups(rows):
        angle = [row["actual_angle_degrees"] for row in selected]
        for ax, quantity in zip(axes, ("r", "t")):
            reference = [row[quantity + "_continuum_phase"] for row in selected]
            error = [row[quantity + "_phase_error"] for row in selected]
            phase = [ref + err if ref is not None and err is not None else np.nan
                     for ref, err in zip(reference, error)]
            reference_angle, expected = curve_values(curves[polarization], quantity, phase_floor)
            line, = ax.plot(angle, phase, "o", label=polarization + " measured")
            ax.plot(reference_angle, expected, "--", color=line.get_color(), label=polarization + " Fresnel")
            ax.set(xlabel="actual central incidence angle [degrees]", ylabel=f"phase({quantity}) [rad]")
            ax.legend()
    for ax in axes:
        phase_limits(ax, *(line.get_ydata() for line in ax.lines))
    figure.suptitle("Phase relative to the continuum branch; undefined near reference zeros")
    save_figure(figure, output, "angle-phase")


def analyse(args, data):
    require_same_scan(data)
    selected = baselines(data, lambda item: (item.run_m.plane["polarization"],
                                            item.run_m.plane["actual_angle_degrees"]))
    rows = central_rows(selected)
    curves = dense_curves(selected, args)
    write_csv(args.output / "central.csv", central_rows(data))
    write_csv(args.output / "angle-baselines.csv", rows)
    amplitude_figure(rows, curves, args.output)
    power_figure(rows, curves, args.output)
    phase_figure(rows, curves, args.output, args.phase_floor)
    controls = matched_controls(data, args.output)
    update_status(args.output, baseline_runs=[item.identifier_m for item in selected],
                  isolation_comparisons=controls,
                  actual_resolution_range=[min(row["cells_per_wavelength"] for row in rows),
                                           max(row["cells_per_wavelength"] for row in rows)])


if __name__ == "__main__":
    raise SystemExit(execute_case(argument_parser(__doc__), analyse))
