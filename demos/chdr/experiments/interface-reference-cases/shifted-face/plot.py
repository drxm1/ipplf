#!/usr/bin/env python3
"""Compare shifted-face complex and magnitude errors at the requested physical plane.

The requested phase origin is retained. The unshifted aligned-lattice oracle is
not substituted for a shifted requested face. The phase-registration deduction
and its limits are derived in ../references.org, Shifted interface and phase origin.
"""
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_common import (argument_parser, baselines, central_rows, execute_case,
                         matched_controls, matching_cell_policy, observed_order, plt,
                         resolution_ticks, save_figure, update_status, write_csv)
from read_run import RunDataError


def shift_fraction(item, base_face):
    return (item.run_m.plane["face_position"] - base_face) / item.run_m.spacing[0]


def require_shifted(data, base_face):
    first = data[0].run_m
    for item in data:
        run, plane = item.run_m, item.run_m.plane
        if plane["polarization"] != "TE" or abs(plane["actual_angle_degrees"]) > 1e-10:
            raise RunDataError("The shifted-face ladder requires normal-incidence TE")
        for key in ("epsilon", "incidence", "omega0"):
            if isinstance(plane[key], str):
                same = plane[key] == first.plane[key]
            else:
                same = math.isclose(plane[key], first.plane[key], rel_tol=2e-12)
            if not same:
                raise RunDataError("Shifted rows must match media, incidence side and physical frequency")
        if abs(base_face / run.spacing[0] - round(base_face / run.spacing[0])) > 1e-10:
            raise RunDataError("base-face must be aligned on every supplied mesh")
        if "represented_face_position" not in plane:
            raise RunDataError("Shifted-face data require the represented face position")
        if not matching_cell_policy(run, first):
            raise RunDataError("Shifted rows must retain one cell shape and time-step policy")


def grouped_rows(data, base_face):
    result = central_rows(data)
    for row, item in zip(result, data):
        row["base_face"] = base_face
        row["shift_fraction"] = shift_fraction(item, base_face)
        row["represented_face"] = item.run_m.plane.get("represented_face_position")
    return result


def shifted_orders(rows):
    result = []
    shifts = sorted(set(round(row["shift_fraction"], 10) for row in rows))
    for shift in shifts:
        ladder = sorted((row for row in rows if math.isclose(row["shift_fraction"], shift, abs_tol=1e-10)),
                        key=lambda row: row["h_x"], reverse=True)
        for coarse, fine in zip(ladder, ladder[1:]):
            row = {"shift_fraction": shift, "coarse_run": coarse["run_id"], "fine_run": fine["run_id"]}
            for quantity in ("r", "t"):
                for kind in ("complex", "magnitude"):
                    key = quantity + "_" + kind + "_error"
                    row[quantity + "_" + kind + "_order"] = observed_order(
                        coarse[key], fine[key], coarse["h_x"], fine["h_x"])
            result.append(row)
    return result


def convergence_figure(rows, output):
    figure, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    shifts = sorted(set(round(row["shift_fraction"], 10) for row in rows))
    for shift in shifts:
        ladder = sorted((row for row in rows if math.isclose(row["shift_fraction"], shift, abs_tol=1e-10)),
                        key=lambda row: row["cells_per_wavelength"])
        for column, quantity in enumerate(("r", "t")):
            for row_index, kind in enumerate(("complex", "magnitude")):
                ax = axes[row_index, column]
                ax.loglog([row["cells_per_wavelength"] for row in ladder],
                          [row[quantity + "_" + kind + "_error"] for row in ladder],
                          "o-", label=f"shift/h = {shift:g}")
                ax.set(xlabel="actual dielectric cells per wavelength", ylabel=f"{quantity}: {kind} error")
                resolution_ticks(ax, [row["cells_per_wavelength"] for row in rows])
                ax.legend(fontsize=8)
    figure.suptitle("Shifted face: continuum comparison at the requested phase origin")
    save_figure(figure, output, "shifted-convergence")


def analyse(args, data):
    if not math.isfinite(args.base_face):
        raise RunDataError("base-face must be finite")
    require_shifted(data, args.base_face)
    selected = baselines(data, lambda item: (item.run_m.spacing[0], round(shift_fraction(item, args.base_face), 10)))
    rows = grouped_rows(selected, args.base_face)
    write_csv(args.output / "central.csv", grouped_rows(data, args.base_face))
    write_csv(args.output / "shift-baselines.csv", rows)
    orders = shifted_orders(rows)
    write_csv(args.output / "orders.csv", orders, ["shift_fraction", "coarse_run", "fine_run",
              "r_complex_order", "r_magnitude_order", "t_complex_order", "t_magnitude_order"])
    convergence_figure(rows, args.output)
    controls = matched_controls(data, args.output)
    update_status(args.output, base_face=args.base_face, observed_orders=orders,
                  baseline_runs=[item.identifier_m for item in selected], isolation_comparisons=controls,
                  reference_scope="continuum at requested face; no unshifted lattice substitution")


if __name__ == "__main__":
    parser = argument_parser(__doc__)
    parser.add_argument("--base-face", type=float, default=500.0,
                        help="Common aligned physical anchor before adding fractional-cell shifts")
    raise SystemExit(execute_case(parser, analyse, continuum_only=True))
