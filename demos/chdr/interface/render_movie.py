#!/usr/bin/env python3
"""Render saved interface fields and the conserved-energy series without a solver or display."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
from vtkmodules.vtkIOXML import vtkXMLPImageDataReader
from vtkmodules.util.numpy_support import vtk_to_numpy


def read_cells(path, name):
    """Read logical VTK cells, x fastest, and return their z=0 slice plus metadata."""
    reader = vtkXMLPImageDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    shape = tuple(n - 1 for n in grid.GetDimensions())
    values = vtk_to_numpy(grid.GetCellData().GetArray(name))
    values = values.reshape((*shape[::-1], -1))[0]
    return values, grid


def main():
    """Use the writer's explicit times and one fixed colour scale for the complete movie."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--field", choices=("E_raw", "B_raw"), default="E_raw")
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    frames = json.loads((args.directory / "wave.pvti.series").read_text())["files"]
    energy = np.genfromtxt(args.directory / "energy.csv", delimiter=",", names=True, skip_header=1)
    epsilon, grid = read_cells(args.directory / "material.pvti", "epsilon_r")
    spacing = grid.GetSpacing()
    nx, ny = epsilon.shape[1], epsilon.shape[0]
    extent = (0, nx * spacing[0], 0, ny * spacing[1])
    initial, first_grid = read_cells(args.directory / frames[0]["name"], args.field)
    offsets = vtk_to_numpy(first_grid.GetFieldData().GetArray("component_offsets"))
    offset = offsets[2 if args.field == "E_raw" else 5]
    field_extent = tuple(value + offset[axis // 2] * spacing[axis // 2]
                         for axis, value in enumerate(extent))
    limit = np.max(np.abs(initial[:, :, 2]))
    fig, (field_ax, energy_ax) = plt.subplots(2, 1, figsize=(7, 8), height_ratios=(3, 1))
    artist = field_ax.imshow(initial[:, :, 2], origin="lower", extent=field_extent,
                            cmap="RdBu_r", vmin=-limit, vmax=limit)
    field_ax.contour((np.arange(nx) + .5) * spacing[0],
                     (np.arange(ny) + .5) * spacing[1], epsilon[:, :, 0],
                     levels=[(epsilon.min() + epsilon.max()) / 2], colors="black", linewidths=.8)
    field_ax.set(xlabel="x", ylabel="y")
    fig.colorbar(artist, ax=field_ax, label=args.field[0] + "z (normalized units)")
    for name in ("vacuum", "dielectric", "total"):
        energy_ax.plot(energy["time"], energy[name] / energy["total"][0], label=name)
    energy_ax.set(xlabel="time (c = 1)", ylabel="vector invariant / initial")
    energy_ax.legend(loc="best", ncols=3)
    marker = energy_ax.axvline(energy["time"][0], color="gray", linestyle=":")
    fig.suptitle("Preliminary")
    fig.tight_layout(rect=(0, 0, 1, .94))

    def update(index):
        values, grid = read_cells(args.directory / frames[index]["name"], args.field)
        artist.set_data(values[:, :, 2])
        time_name = "E_time" if args.field == "E_raw" else "B_time"
        time = grid.GetFieldData().GetArray(time_name).GetTuple1(0)
        field_ax.set_title(f"{args.field[0]}z at t = {time:.3f}; black line: dielectric face")
        marker.set_xdata([time, time])
        return artist, marker

    output = args.output or args.directory / "movie.gif"
    animation = FuncAnimation(fig, update, frames=len(frames), interval=100, blit=False)
    animation.save(output, writer=PillowWriter(fps=10))
    update(len(frames) - 1)
    fig.savefig(output.with_suffix(".png"), dpi=160)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 3))
    for name in ("vacuum", "dielectric", "total"):
        ax.plot(energy["time"], energy[name] / energy["total"][0], label=name)
    ax.set(xlabel="time (c = 1)", ylabel="vector invariant / initial",
           title="Preliminary")
    ax.legend(ncols=3)
    fig.tight_layout()
    fig.savefig(output.with_name(output.stem + "-energy.png"), dpi=160)


if __name__ == "__main__":
    main()
