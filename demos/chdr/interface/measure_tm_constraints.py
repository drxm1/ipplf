#!/usr/bin/env python3
"""Record the charge-free TM constraint from saved x-z movie frames.

Compute the backward divergence of harmonic-face epsilon times E, using the
solver's staggered positions [chew2014generalizedgauge, section 2, Eqs. (4), (11);
ryu2016potentialfdtd, Appendix A]. This is a diagnostic of div(epsilon E), not
a correction or an assertion that componentwise Mur preserves Gauss's law.
The input must be y-independent with periodic y and open x/z faces.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from render_movie import field_time, read_cells, read_metadata


def face_permittivity(epsilon, axis):
    """Match positive-face harmonic coefficients and constant material extension at open faces."""
    neighbour = np.roll(epsilon, -1, axis=axis)
    end = ([-1, slice(None)] if axis == 0 else [slice(None), -1])
    neighbour[tuple(end)] = epsilon[tuple(end)]
    return 2.0 / (1.0 / epsilon + 1.0 / neighbour)


def residual(electric, faces, spacing):
    """Backward divergence where the two contributing face samples were saved."""
    dx, dz = faces[0] * electric[:, :, 0], faces[1] * electric[:, :, 2]
    divergence = ((dx - np.roll(dx, 1, axis=0)) / spacing[0]
                  + (dz - np.roll(dz, 1, axis=1)) / spacing[2])
    valid = np.ones(divergence.shape, dtype=bool)
    valid[0, :], valid[:, 0] = False, False
    interior = np.zeros_like(valid)
    interior[2:-2, 2:-2] = True
    values = [float(np.abs(divergence[mask]).max())
              for mask in (valid, interior, valid & ~interior)]
    return values, max(float(np.abs(dx).max()), float(np.abs(dz).max()))


def measure(directory, output):
    """Write all saved-time residuals, separating the interior from the boundary band.

    The lowest x/z cell layers are omitted: their backward stencils need unsaved
    exterior E. The interior excludes two cells on all open faces. Other valid
    samples form the boundary band. Scale by initial max|epsilon E| / min(h_x,h_z),
    so the dimensionless column does not grow merely because the late pulse decays.
    """
    metadata = read_metadata(directory)
    if metadata.get("periodic") != "0,1,0" or metadata.get("polarization") != "TM":
        raise ValueError("Expected a y-independent TM movie with open x/z and periodic y")
    epsilon, grid = read_cells(directory / "material.pvti", "epsilon_r", "zx")
    faces = [face_permittivity(epsilon[:, :, 0], axis) for axis in (0, 1)]
    frames = json.loads((directory / "wave.pvti.series").read_text())["files"]
    with output.open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(["frame", "electric_time", "div_valid_max", "div_interior_max",
                         "div_boundary_max", "displacement_peak", "div_scaled_initial"])
        for index, frame in enumerate(frames):
            electric, grid = read_cells(directory / frame["name"], "E_raw", "zx")
            if not np.isfinite(electric).all() or np.any(electric[:, :, 1] != 0):
                raise ValueError("Expected finite in-plane electric fields")
            maxima, peak = residual(electric, faces, grid.GetSpacing())
            if index == 0:
                scale = peak / min(grid.GetSpacing()[0], grid.GetSpacing()[2])
                if scale <= 0:
                    raise ValueError("Initial electric field is zero")
            writer.writerow([index, field_time(grid, "E_raw"), *maxima, peak, maxima[0] / scale])
    print(f"Wrote {len(frames)} constraint samples to {output}")


def main():
    """Select saved output; no simulation or boundary state is modified."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    measure(args.directory, args.output or args.directory / "constraints.csv")


if __name__ == "__main__":
    main()
