#!/usr/bin/env python3
"""Compare independent propagating dielectric-side rows using existing references.

Coefficient and flux definitions are in ../references.org, from Ellingson
sections 5.6-5.7 and the documented aligned-lattice recurrence derivation.
No coefficient formula is implemented here. The two supplied rows have different
materials, cell shapes and resolutions; no common angle curve is inferred.
"""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_common import argument_parser, central_rows, execute_case, update_status, write_csv
from read_run import RunDataError


def require_dielectric_rows(data):
    """Keep evanescent or unsupported geometry out of the selected primary comparison."""
    for item in data:
        if (item.run_m.plane["geometry"] != "aligned_plane"
                or item.run_m.plane["incidence"] != "dielectric"):
            raise RunDataError("This comparison requires aligned dielectric-side incidence")
        for sample, reference in zip(item.run_m.samples_m, item.reference_m.rows_m):
            if sample.valid_m and any(reference[f"{kind}_{side}_root"] != "propagating"
                                     for kind in ("continuum", "lattice")
                                     for side in ("incident", "transmitted")):
                raise RunDataError("A selected dielectric-side bin must propagate in both media")


def case_rows(data):
    """Retain each row's distinct material, requested resolution and actual cell shape."""
    rows = central_rows(data)
    for row, item in zip(rows, data):
        row["requested_cells_per_wavelength"] = item.run_m.plane["requested_cells_per_wavelength"]
        row["h_y"], row["h_z"] = item.run_m.spacing[1:]
        row["input_preset"] = item.run_m.record_m["input"].get("parameters", {}).get("preset", "")
    return rows


def analyse(args, data):
    require_dielectric_rows(data)
    write_csv(args.output / "central.csv", case_rows(data))
    update_status(args.output, comparison_scope="independent dielectric-side rows",
                  reference_scope="propagating continuum and aligned-lattice coefficients",
                  shared_angle_curve=False, observed_orders=[],
                  qualification="Per-bin comparisons; assess each run's isolation separately")


if __name__ == "__main__":
    raise SystemExit(execute_case(argument_parser(__doc__), analyse))
