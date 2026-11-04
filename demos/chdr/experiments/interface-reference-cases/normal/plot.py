#!/usr/bin/env python3
"""Compare normal TE runs with continuum references and a mesh-refinement ladder.

The reference and convergence qualifications are in ../references.org. Numerical
tables retain complex and magnitude errors separately; a two-level exponent is
reported only as an observed value, not an asserted asymptotic order.
"""
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_common import (argument_parser, baselines, central_rows, execute_case,
                         matched_controls, matching_cell_policy, observed_order, plt, pulse_profiles,
                         save_figure, update_status, write_csv)
from read_run import RunDataError


def require_normal(data):
    first = data[0].run_m
    for item in data:
        run, plane = item.run_m, item.run_m.plane
        if plane["polarization"] != "TE" or abs(plane["actual_angle_degrees"]) > 1e-10:
            raise RunDataError("The normal ladder requires TE at normal incidence")
        for key in ("epsilon", "incidence"):
            if plane[key] != first.plane[key]:
                raise RunDataError("Normal ladder media and incidence side must match")
        if not math.isclose(plane["omega0"], first.plane["omega0"], rel_tol=2e-12):
            raise RunDataError("Normal refinement must keep physical central frequency fixed; use cell_size")
        if not matching_cell_policy(run, first):
            raise RunDataError("Normal ladder cell shape and time-step policy must match")


def order_rows(rows):
    result = []
    for coarse, fine in zip(rows, rows[1:]):
        item = {"coarse_run": coarse["run_id"], "fine_run": fine["run_id"],
                "coarse_h": coarse["h_x"], "fine_h": fine["h_x"]}
        for quantity in ("r", "t"):
            for kind in ("complex", "magnitude"):
                name = quantity + "_" + kind + "_error"
                item[quantity + "_" + kind + "_order"] = observed_order(
                    coarse[name], fine[name], coarse["h_x"], fine["h_x"])
        result.append(item)
    return result


def error_figure(rows, output):
    figure, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for ax, quantity in zip(axes, ("r", "t")):
        for kind, marker in (("complex", "o"), ("magnitude", "s")):
            ax.loglog([row["cells_per_wavelength"] for row in rows],
                      [row[quantity + "_" + kind + "_error"] for row in rows],
                      marker + "-", label=kind + " error")
        ax.set(xlabel="actual dielectric cells per wavelength", ylabel=f"absolute {quantity} error")
        ax.legend()
    figure.suptitle("Normal TE: continuum error at the recorded central frequency")
    save_figure(figure, output, "normal-convergence")


def write_ladder(output, data, rows):
    write_csv(output / "central.csv", central_rows(data))
    write_csv(output / "ladder.csv", rows)
    orders = order_rows(rows)
    write_csv(output / "orders.csv", orders,
              ["coarse_run", "fine_run", "coarse_h", "fine_h", "r_complex_order",
               "r_magnitude_order", "t_complex_order", "t_magnitude_order"])
    return orders


def analyse(args, data):
    require_normal(data)
    selected = baselines(data, lambda item: item.run_m.spacing[0])
    selected.sort(key=lambda item: item.run_m.spacing[0], reverse=True)
    rows = central_rows(selected)
    orders = write_ladder(args.output, data, rows)
    error_figure(rows, args.output)
    controls = matched_controls(data, args.output)
    profiles = [] if args.skip_profiles else [item.identifier_m for item in selected if pulse_profiles(item)]
    if not args.skip_profiles and len(profiles) != len(selected):
        raise RunDataError("Normal baseline profiles require saved frames; use --skip-profiles only for spectrum-only analysis")
    update_status(args.output, baseline_runs=[item.identifier_m for item in selected],
                  observed_orders=orders, isolation_comparisons=controls, profile_runs=profiles,
                  profiles_explicitly_skipped=args.skip_profiles,
                  qualification="Observed errors; boundary controls must be assessed separately")


if __name__ == "__main__":
    parser = argument_parser(__doc__)
    parser.add_argument("--skip-profiles", action="store_true",
                        help="Explicit spectrum-only analysis when no field frames are available")
    raise SystemExit(execute_case(parser, analyse))
