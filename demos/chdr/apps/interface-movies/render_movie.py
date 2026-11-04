#!/usr/bin/env python3
"""Render saved interface fields and energy diagnostics without a solver or display."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter
import numpy as np
from vtkmodules.vtkIOXML import vtkXMLPImageDataReader
from vtkmodules.util.numpy_support import vtk_to_numpy


def plane_info(plane):
    """Return horizontal axis, vertical axis and out-of-plane vector component."""
    return (2, 0, 1) if plane == "zx" else (0, 1, 2)


def read_cells(path, name, plane="xy"):
    """Read logical VTK cells and return the requested display slice plus metadata."""
    reader = vtkXMLPImageDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    array = grid.GetCellData().GetArray(name)
    if array is None:
        raise ValueError(f"{path}: missing cell array {name}")
    shape = tuple(n - 1 for n in grid.GetDimensions())
    values = vtk_to_numpy(array).reshape((*shape[::-1], -1))
    if plane == "zx":
        values = values[:, 0, :, :].transpose(1, 0, 2)
    else:
        values = values[0]
    return values, grid


def field_time(grid, field):
    """Use the writer's explicit time label for the selected field."""
    name = "E_time" if field == "E_raw" else "B_time"
    return grid.GetFieldData().GetArray(name).GetTuple1(0)


def scan_frames(directory, frames, field, plane, normal_edges_only=False):
    """Find a global colour limit and the first boundary signal above 1% of that peak."""
    component = plane_info(plane)[2]
    peak, boundary = 0.0, []
    for index, frame in enumerate(frames):
        values, grid = read_cells(directory / frame["name"], field, plane)
        sample = values[:, :, component]
        if not np.isfinite(sample).all():
            raise ValueError(f"Nonfinite field in {frame['name']}")
        peak = max(peak, float(np.abs(sample).max()))
        edge = max(np.abs(sample[0]).max(), np.abs(sample[-1]).max())
        if not normal_edges_only:
            edge = max(edge, np.abs(sample[:, 0]).max(), np.abs(sample[:, -1]).max())
        boundary.append((field_time(grid, field), edge))
        if index % 60 == 0 or index + 1 == len(frames):
            print(f"Scanned {index + 1}/{len(frames)} frames", flush=True)
    if peak <= 0:
        raise ValueError("The selected out-of-plane field is zero in every frame")
    contact = next((time for time, edge in boundary if edge >= 0.01 * peak), None)
    return peak, contact


def display_geometry(grid, field, plane):
    """Return cell centres and the component-offset-corrected image extent."""
    horizontal, vertical, component = plane_info(plane)
    spacing, origin = grid.GetSpacing(), grid.GetOrigin()
    cells = tuple(n - 1 for n in grid.GetDimensions())
    axes = (horizontal, vertical)
    centres = [origin[d] + (np.arange(cells[d]) + 0.5) * spacing[d] for d in axes]
    extent = [value for d in axes
              for value in (origin[d], origin[d] + cells[d] * spacing[d])]
    offsets = vtk_to_numpy(grid.GetFieldData().GetArray("component_offsets"))
    offset = offsets[component + (3 if field == "B_raw" else 0)]
    shifted = [value + offset[axes[i // 2]] * spacing[axes[i // 2]]
               for i, value in enumerate(extent)]
    return centres, extent, shifted


def repeated_geometry(grid, field, plane, repeats):
    """Extend the horizontal coordinate by whole periods, preserving cell offsets."""
    centres, extent, shifted = display_geometry(grid, field, plane)
    period = extent[1] - extent[0]
    centres[0] = np.concatenate([centres[0] + index * period for index in range(repeats)])
    extent[1] += (repeats - 1) * period
    shifted[1] += (repeats - 1) * period
    return centres, extent, shifted


def repeated_component(values, plane, repeats):
    """Repeat already computed samples without interpolation or a new simulation."""
    return np.tile(values[:, :, plane_info(plane)[2]], (1, repeats))


def draw_field(ax, values, grid, epsilon, field, plane, limit, length_unit="", repeats=1):
    """Draw one field slice and material outline in physical coordinate order."""
    horizontal, vertical, component = plane_info(plane)
    centres, extent, shifted = repeated_geometry(grid, field, plane, repeats)
    artist = ax.imshow(repeated_component(values, plane, repeats), origin="lower", extent=shifted,
                       cmap="RdBu_r", vmin=-limit, vmax=limit, aspect="equal")
    if epsilon.min() != epsilon.max():
        ax.contour(*centres, np.tile(epsilon[:, :, 0], (1, repeats)),
                   levels=[(epsilon.min() + epsilon.max()) / 2],
                   colors="black", linewidths=0.8)
    ax.set(xlabel="xyz"[horizontal] + length_unit, ylabel="xyz"[vertical] + length_unit,
           xlim=extent[:2], ylim=extent[2:])
    return artist


def energy_scale(energy):
    """Prefer the first saved discrete total; a zero first sample uses the saved peak."""
    first = abs(float(energy["total"][0]))
    if first > 0:
        return first, "first saved discrete total"
    return max(float(np.abs(energy["total"]).max()), 1e-30), "peak discrete total"


def draw_energy(ax, energy, contact, time_unit="c = 1", contact_label="boundary signal"):
    """Distinguish the cross-level discrete diagnostic from the aligned field-energy sum."""
    if energy.size == 0:
        ax.set(xlabel=f"time ({time_unit})", ylabel="energy")
        ax.text(0.5, 0.5, "No energy samples saved", transform=ax.transAxes, ha="center")
        return
    scale, label = energy_scale(energy)
    curves = (("vacuum", "discrete: vacuum", "-", 0.5),
              ("dielectric", "discrete: dielectric", "-", 0.5),
              ("total", "discrete: total", "-", 1.0),
              ("field_energy", "aligned field energy", "--", 1.0))
    for name, title, style, alpha in curves:
        if name in energy.dtype.names:
            ax.plot(energy["time"], energy[name] / scale, style, label=title, alpha=alpha)
    if contact is not None:
        ax.axvline(contact, color="gray", linestyle=":", label=contact_label + ": 1% of peak")
    ax.set(xlabel=f"time ({time_unit})", ylabel=f"energy / {label}")
    ax.legend(loc="best", ncols=2, fontsize=7)


def read_metadata(directory):
    """Read the CSV case metadata; legacy outputs have no plane declaration."""
    with (directory / "energy.csv").open() as stream:
        line = stream.readline()
    return dict(token.split("=", 1) for token in line.split() if "=" in token)


def output_units(directory, field):
    """Use the common record for SI labels; older saved runs use normalized units."""
    path = directory / "run.json"
    system = json.loads(path.read_text())["output_factors"]["system"] if path.exists() else "normalized"
    if system == "si":
        return "s", "V/m" if field == "E_raw" else "T", " (m)"
    return "c = 1", "normalized units", ""


def validate_repeats(directory, plane, repeats):
    """Allow displayed z copies only for a recorded aligned, z-periodic plane run.

    The option changes the view, never the computed box or measured coefficients.
    Equal image aspect and actual component coordinates preserve the displayed
    wavefront orientation. The supplied periodic record is the authority.
    """
    if repeats < 1:
        raise ValueError("--transverse-repeats must be positive")
    if repeats == 1:
        return
    record = json.loads((directory / "run.json").read_text())
    effective = record.get("effective", {})
    geometry = effective.get("plane", {})
    periodic = effective.get("box", {}).get("periodic", [])
    if (plane != "zx" or geometry.get("geometry") != "aligned_plane"
            or geometry.get("incidence_plane") != "xz" or len(periodic) != 3 or periodic[2] != 1):
        raise ValueError("Transverse repeats require a recorded aligned plane with periodic z")


class MovieRenderer:
    """Read frames on demand and share plotting conventions across animation and previews."""

    def __init__(self, args):
        self.set_view(args)
        self.read_series()
        self.make_figure()

    def set_view(self, args):
        """Resolve metadata, units and the explicitly requested display repetition."""
        self.directory, self.field = args.directory, args.field
        metadata = read_metadata(args.directory)
        if self.field is None:
            self.field = "B_raw" if metadata.get("polarization") == "TM" else "E_raw"
        self.time_unit, self.field_unit, self.length_unit = output_units(self.directory, self.field)
        self.plane = args.plane or metadata.get("plane", "xy")
        self.repeats = args.transverse_repeats
        validate_repeats(self.directory, self.plane, self.repeats)
        self.set_caveat(metadata)
        self.label = self.field[0] + "xyz"[plane_info(self.plane)[2]]

    def set_caveat(self, metadata):
        """Identify preliminary legacy cases and distinguish a repeated view from the box."""
        if metadata.get("case") == "plane":
            self.caveat = "Periodic plane pulse"
        elif metadata.get("polarization") == "TM":
            self.caveat = "Preliminary: nonzero Gauss residual under Mur"
            if metadata.get("case", "").startswith("diagonal"):
                self.caveat += "\nDiagonal interface and boundary treatment not validated"
        elif metadata.get("case") == "diagonal":
            self.caveat = "Preliminary: diagonal interface and\nboundary treatment not validated"
        else:
            self.caveat = "Preliminary"
        self.contact_label = "normal-edge signal" if self.repeats > 1 else "boundary signal"
        if self.repeats > 1:
            self.caveat += f"\nDisplay only: {self.repeats} periodic z copies; equal axis scales"

    def read_series(self):
        """Load frame metadata and diagnostics, then determine one fixed colour scale."""
        self.frames = json.loads((self.directory / "wave.pvti.series").read_text())["files"]
        if not self.frames:
            raise ValueError("The series contains no frames")
        self.energy = np.atleast_1d(np.genfromtxt(
            self.directory / "energy.csv", delimiter=",", names=True, skip_header=1))
        self.epsilon, _ = read_cells(self.directory / "material.pvti", "epsilon_r", self.plane)
        self.limit, self.contact = scan_frames(
            self.directory, self.frames, self.field, self.plane, self.repeats > 1)

    def frame(self, index):
        """Read one saved frame without retaining the complete time series in memory."""
        return read_cells(self.directory / self.frames[index]["name"], self.field, self.plane)

    def make_figure(self):
        """Set up the field and energy panels, with only the field and time marker animated."""
        values, grid = self.frame(0)
        self.fig, (self.field_ax, energy_ax) = plt.subplots(
            2, 1, figsize=(7, 8), height_ratios=(3, 1))
        self.artist = draw_field(self.field_ax, values, grid, self.epsilon,
                                 self.field, self.plane, self.limit, self.length_unit, self.repeats)
        self.fig.colorbar(self.artist, ax=self.field_ax, label=f"{self.label} ({self.field_unit})")
        draw_energy(energy_ax, self.energy, self.contact, self.time_unit, self.contact_label)
        if self.energy.size == 0 and len(self.frames) > 1:
            energy_ax.set_xlim(self.frames[0]["time"], self.frames[-1]["time"])
        self.marker = energy_ax.axvline(field_time(grid, self.field), color="black", linestyle="-",
                                        linewidth=1, label="current frame")
        energy_ax.legend(loc="best", ncols=2, fontsize=7)
        self.fig.suptitle(self.caveat)
        self.fig.tight_layout(rect=(0, 0, 1, 0.94))

    def update(self, index):
        """Update the selected field component and its actual sample-time label."""
        values, grid = self.frame(index)
        self.artist.set_data(repeated_component(values, self.plane, self.repeats))
        time = field_time(grid, self.field)
        self.field_ax.set_title(
            f"{self.label} at t = {time:.3g} ({self.time_unit}); black line: dielectric face")
        self.marker.set_xdata([time, time])
        return self.artist, self.marker

    def save(self, output):
        """Write MP4 or GIF plus a representative field frame, overview and energy plot."""
        animation = FuncAnimation(self.fig, self.update, frames=len(self.frames), interval=100)
        writer = PillowWriter(fps=12) if output.suffix.lower() == ".gif" else FFMpegWriter(
            fps=12, bitrate=2400, extra_args=["-pix_fmt", "yuv420p"])
        animation.save(output, writer=writer, dpi=120)
        self.update(len(self.frames) // 3)
        self.fig.savefig(output.with_suffix(".png"), dpi=160)
        plt.close(self.fig)
        self.save_overview(output)
        fig, ax = plt.subplots(figsize=(7, 3))
        draw_energy(ax, self.energy, self.contact, self.time_unit, self.contact_label)
        fig.tight_layout()
        fig.savefig(output.with_name(output.stem + "-energy.png"), dpi=160)
        plt.close(fig)
        self.save_settings(output)

    def save_settings(self, output):
        """Retain the exact view choice and source record beside the derived movie."""
        source = self.directory / "run.json"
        record = {"run_directory": str(self.directory.resolve()), "field": self.field,
                  "plane": self.plane, "transverse_repeats": self.repeats,
                  "display_only": True, "axis_aspect": "equal", "colour_limit": self.limit,
                  "edge_display_criterion": self.contact_label + ": 1% of global sampled peak",
                  "edge_display_criterion_is_error_bound": False, "first_edge_time": self.contact,
                  "frames": self.frames, "movie": output.name,
                  "renderer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        if source.exists():
            record["run_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
        output.with_name(output.stem + "-render.json").write_text(json.dumps(record, indent=2) + "\n")

    def save_overview(self, output):
        """Show initial, interaction and late-time fields using the movie's fixed colour scale."""
        indices = (0, len(self.frames) // 3, 2 * len(self.frames) // 3, len(self.frames) - 1)
        fig, axes = plt.subplots(1, 4, figsize=(14, 3.7), layout="constrained")
        fig.suptitle(self.caveat)
        for ax, index in zip(axes, indices):
            values, grid = self.frame(index)
            artist = draw_field(ax, values, grid, self.epsilon,
                                self.field, self.plane, self.limit, self.length_unit, self.repeats)
            ax.set_title(f"t = {field_time(grid, self.field):.3g} ({self.time_unit})")
        fig.colorbar(artist, ax=list(axes), label=f"{self.label} ({self.field_unit})", shrink=0.8)
        fig.savefig(output.with_name(output.stem + "-overview.png"), dpi=140)
        plt.close(fig)


def main():
    """Select a saved run and render it; plane metadata preserves older x-y datasets."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--field", choices=("E_raw", "B_raw"), default=None,
                        help="Override the field selected from polarization metadata")
    parser.add_argument("--plane", choices=("zx", "xy"), default=None)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--transverse-repeats", type=int, default=1,
                        help="Display whole z-period copies of an aligned periodic plane (default: 1)")
    args = parser.parse_args()
    output = args.output or args.directory / "movie.mp4"
    if output.suffix.lower() not in (".mp4", ".gif"):
        parser.error("--output must end in .mp4 or .gif")
    MovieRenderer(args).save(output)


if __name__ == "__main__":
    main()
