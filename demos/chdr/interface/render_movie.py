#!/usr/bin/env python3
"""Render saved interface fields and energy diagnostics without a solver or display."""
import argparse
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


def scan_frames(directory, frames, field, plane):
    """Find a global colour limit and the first boundary signal above 1% of that peak."""
    component = plane_info(plane)[2]
    peak, boundary = 0.0, []
    for index, frame in enumerate(frames):
        values, grid = read_cells(directory / frame["name"], field, plane)
        sample = values[:, :, component]
        if not np.isfinite(sample).all():
            raise ValueError(f"Nonfinite field in {frame['name']}")
        peak = max(peak, float(np.abs(sample).max()))
        edge = max(np.abs(sample[0]).max(), np.abs(sample[-1]).max(),
                   np.abs(sample[:, 0]).max(), np.abs(sample[:, -1]).max())
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


def draw_field(ax, values, grid, epsilon, field, plane, limit):
    """Draw one field slice and material outline in physical coordinate order."""
    horizontal, vertical, component = plane_info(plane)
    centres, extent, shifted = display_geometry(grid, field, plane)
    artist = ax.imshow(values[:, :, component], origin="lower", extent=shifted,
                       cmap="RdBu_r", vmin=-limit, vmax=limit)
    if epsilon.min() != epsilon.max():
        ax.contour(*centres, epsilon[:, :, 0],
                   levels=[(epsilon.min() + epsilon.max()) / 2],
                   colors="black", linewidths=0.8)
    ax.set(xlabel="xyz"[horizontal], ylabel="xyz"[vertical],
           xlim=extent[:2], ylim=extent[2:])
    return artist


def energy_scale(energy):
    """Prefer the initial discrete total; zero-start datasets use their peak instead."""
    initial = abs(float(energy["total"][0]))
    if initial > 0:
        return initial, "initial discrete total"
    return max(float(np.abs(energy["total"]).max()), 1e-30), "peak discrete total"


def draw_energy(ax, energy, contact):
    """Distinguish the cross-level discrete diagnostic from the aligned field-energy sum."""
    scale, label = energy_scale(energy)
    curves = (("vacuum", "discrete: vacuum", "-", 0.5),
              ("dielectric", "discrete: dielectric", "-", 0.5),
              ("total", "discrete: total", "-", 1.0),
              ("field_energy", "aligned field energy", "--", 1.0))
    for name, title, style, alpha in curves:
        if name in energy.dtype.names:
            ax.plot(energy["time"], energy[name] / scale, style, label=title, alpha=alpha)
    if contact is not None:
        ax.axvline(contact, color="gray", linestyle=":", label="boundary signal: 1% of peak")
    ax.set(xlabel="time (c = 1)", ylabel=f"energy / {label}")
    ax.legend(loc="best", ncols=2, fontsize=7)


def read_metadata(directory):
    """Read the simple key=value run record; legacy outputs have no plane declaration."""
    with (directory / "energy.csv").open() as stream:
        line = stream.readline()
    return dict(token.split("=", 1) for token in line.split() if "=" in token)


class MovieRenderer:
    """Read frames on demand and share plotting conventions across animation and previews."""

    def __init__(self, args):
        self.directory, self.field = args.directory, args.field
        metadata = read_metadata(args.directory)
        self.plane = args.plane or metadata.get("plane", "xy")
        self.caveat = ("Preliminary: diagonal interface and\nboundary treatment not validated"
                       if metadata.get("case") == "diagonal" else "Preliminary")
        self.frames = json.loads((args.directory / "wave.pvti.series").read_text())["files"]
        if not self.frames:
            raise ValueError("The series contains no frames")
        self.energy = np.atleast_1d(np.genfromtxt(
            args.directory / "energy.csv", delimiter=",", names=True, skip_header=1))
        self.epsilon, _ = read_cells(args.directory / "material.pvti", "epsilon_r", self.plane)
        self.limit, self.contact = scan_frames(
            self.directory, self.frames, self.field, self.plane)
        self.label = self.field[0] + "xyz"[plane_info(self.plane)[2]]
        self.make_figure()

    def frame(self, index):
        """Read one saved frame without retaining the complete time series in memory."""
        return read_cells(self.directory / self.frames[index]["name"], self.field, self.plane)

    def make_figure(self):
        """Set up the field and energy panels, with only the field and time marker animated."""
        values, grid = self.frame(0)
        self.fig, (self.field_ax, energy_ax) = plt.subplots(
            2, 1, figsize=(7, 8), height_ratios=(3, 1))
        self.artist = draw_field(self.field_ax, values, grid, self.epsilon,
                                 self.field, self.plane, self.limit)
        self.fig.colorbar(self.artist, ax=self.field_ax, label=f"{self.label} (normalized units)")
        draw_energy(energy_ax, self.energy, self.contact)
        self.marker = energy_ax.axvline(self.energy["time"][0], color="gray", linestyle=":")
        self.fig.suptitle(self.caveat)
        self.fig.tight_layout(rect=(0, 0, 1, 0.94))

    def update(self, index):
        """Update the selected field component and its actual sample-time label."""
        values, grid = self.frame(index)
        self.artist.set_data(values[:, :, plane_info(self.plane)[2]])
        time = field_time(grid, self.field)
        self.field_ax.set_title(f"{self.label} at t = {time:.3f}; black line: dielectric face")
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
        draw_energy(ax, self.energy, self.contact)
        fig.tight_layout()
        fig.savefig(output.with_name(output.stem + "-energy.png"), dpi=160)
        plt.close(fig)

    def save_overview(self, output):
        """Show initial, interaction and late-time fields using the movie's fixed colour scale."""
        indices = (0, len(self.frames) // 3, 2 * len(self.frames) // 3, len(self.frames) - 1)
        fig, axes = plt.subplots(1, 4, figsize=(14, 3.7), layout="constrained")
        fig.suptitle(self.caveat)
        for ax, index in zip(axes, indices):
            values, grid = self.frame(index)
            artist = draw_field(ax, values, grid, self.epsilon,
                                self.field, self.plane, self.limit)
            ax.set_title(f"t = {field_time(grid, self.field):.2f}")
        fig.colorbar(artist, ax=list(axes), label=self.label, shrink=0.8)
        fig.savefig(output.with_name(output.stem + "-overview.png"), dpi=140)
        plt.close(fig)


def main():
    """Select a saved run and render it; plane metadata preserves older x-y datasets."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--field", choices=("E_raw", "B_raw"), default="E_raw")
    parser.add_argument("--plane", choices=("zx", "xy"), default=None)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    output = args.output or args.directory / "movie.mp4"
    if output.suffix.lower() not in (".mp4", ".gif"):
        parser.error("--output must end in .mp4 or .gif")
    MovieRenderer(args).save(output)


if __name__ == "__main__":
    main()
