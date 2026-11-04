#!/usr/bin/env python3
"""Compare central plane coefficients and continuum flux estimates over incidence angle.

The C++ utility supplies every continuum reference. A curve connecting central
reference values is labelled as sampled; it is not a separately evaluated dense
Fresnel curve. See ../references.org for fixed-H TM signs and Brewster conventions.
"""
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_common import (argument_parser, baselines, central_rows, execute_case,
                         matched_controls, np, plt, save_figure, update_status, write_csv)
from read_run import RunDataError


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


def amplitude_figure(rows, output):
    figure, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for polarization, selected in groups(rows):
        angle = [row["actual_angle_degrees"] for row in selected]
        for ax, quantity in zip(axes, ("r", "t")):
            values = [math.hypot(row[quantity + "_real"], row[quantity + "_imag"]) for row in selected]
            expected = [math.hypot(row[quantity + "_continuum_real"], row[quantity + "_continuum_imag"])
                        for row in selected]
            line, = ax.plot(angle, values, "o", label=polarization + " measured")
            ax.plot(angle, expected, "x--", color=line.get_color(), label=polarization + " Fresnel samples")
            ax.set(xlabel="actual central incidence angle [degrees]", ylabel=f"|{quantity}|")
            ax.legend(fontsize=8)
    figure.suptitle("Vacuum–dielectric plane interface: central amplitudes")
    save_figure(figure, output, "angle-amplitudes")


def power_figure(rows, output):
    figure, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for polarization, selected in groups(rows):
        angle = [row["actual_angle_degrees"] for row in selected]
        for quantity, marker in (("R", "o"), ("T", "s")):
            axes[0].plot(angle, [row[quantity] for row in selected], marker + "-",
                         label=polarization + " " + quantity)
        axes[1].plot(angle, [row["R_plus_T_minus_one"] for row in selected], "o-", label=polarization)
    axes[0].set(xlabel="actual central incidence angle [degrees]", ylabel="continuum flux estimate")
    axes[1].set(xlabel="actual central incidence angle [degrees]", ylabel="R + T − 1")
    axes[1].axhline(0, color="black", linestyle=":", linewidth=0.8)
    for ax in axes:
        ax.legend()
    save_figure(figure, output, "angle-power")


def phase_figure(rows, output):
    figure, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for polarization, selected in groups(rows):
        angle = [row["actual_angle_degrees"] for row in selected]
        for ax, quantity in zip(axes, ("r", "t")):
            reference = [row[quantity + "_continuum_phase"] for row in selected]
            error = [row[quantity + "_phase_error"] for row in selected]
            phase = [ref + err if ref is not None and err is not None else np.nan
                     for ref, err in zip(reference, error)]
            line, = ax.plot(angle, phase, "o-", label=polarization + " measured")
            ax.plot(angle, reference, "x--", color=line.get_color(), label=polarization + " Fresnel samples")
            ax.set(xlabel="actual central incidence angle [degrees]", ylabel=f"phase({quantity}) [rad]")
            ax.legend()
    figure.suptitle("Phase relative to the continuum branch; undefined near reference zeros")
    save_figure(figure, output, "angle-phase")


def analyse(args, data):
    require_same_scan(data)
    selected = baselines(data, lambda item: (item.run_m.plane["polarization"],
                                            item.run_m.plane["actual_angle_degrees"]))
    rows = central_rows(selected)
    write_csv(args.output / "central.csv", central_rows(data))
    write_csv(args.output / "angle-baselines.csv", rows)
    amplitude_figure(rows, args.output)
    power_figure(rows, args.output)
    phase_figure(rows, args.output)
    controls = matched_controls(data, args.output)
    update_status(args.output, baseline_runs=[item.identifier_m for item in selected],
                  isolation_comparisons=controls, reference_curve="central samples connected",
                  actual_resolution_range=[min(row["cells_per_wavelength"] for row in rows),
                                           max(row["cells_per_wavelength"] for row in rows)])


if __name__ == "__main__":
    raise SystemExit(execute_case(argument_parser(__doc__), analyse))
