"""Render one saved PVTI-series frame with fixed view/color settings, without a GUI.

Run with pvpython --force-offscreen-rendering --disable-registry.
API: https://www.paraview.org/paraview-docs/latest/python/paraview.simple.html
The PNG is a viewing aid; numerical acceptance uses the separate reader tests.
"""
import argparse
import json
from pathlib import Path

from paraview import servermanager
from paraview.simple import (ColorBy, CreateView, GetColorTransferFunction,
                            GetParaViewVersion, GetScalarBar, OpenDataFile,
                            Render, ResetCamera, SaveScreenshot, Show, Text)
from vtkmodules.vtkCommonCore import vtkVersion


def arguments():
    """Select one exact saved time and a signed vector component."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("series", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--time", type=float, default=0.25)
    parser.add_argument("--array", choices=("E_raw", "B_raw"), default="E_raw")
    parser.add_argument("--component", choices=("X", "Y", "Z"), default="X")
    parser.add_argument("--range", type=float, nargs=2, default=(-1.0, 1.0), dest="limits")
    options = parser.parse_args()
    if options.output.exists() or not options.limits[0] < options.limits[1]:
        parser.error("choose a new PNG path and an increasing color range")
    return options


def read_frame(options):
    """Use the explicit saved time; never substitute an inferred filename ordinal."""
    reader = OpenDataFile(str(options.series.resolve()))
    if reader is None or options.time not in list(reader.TimestepValues):
        raise ValueError("Requested time is not present in the series")
    reader.UpdatePipeline(options.time)
    image = servermanager.Fetch(reader)
    array = image.GetCellData().GetArray(options.array)
    if array is None or array.GetNumberOfComponents() != 3:
        raise ValueError("Selected three-component CellData vector is missing")
    return reader, image


def metadata(image):
    """Read labels directly from the dataset, including the separate E/B timestamps."""
    names = ("potential_time", "E_time", "B_time", "dt", "step", "normalization")
    result = {}
    for name in names:
        array = image.GetFieldData().GetArray(name)
        if array is None:
            raise ValueError(f"Missing metadata: {name}")
        result[name] = array.GetTuple1(0)
    return result


def configure_view(view, image):
    """Set an orthographic +y viewpoint: +x upward and +z rightward."""
    bounds = image.GetBounds()
    center = [(bounds[2 * d] + bounds[2 * d + 1]) / 2 for d in range(3)]
    lengths = [bounds[2 * d + 1] - bounds[2 * d] for d in range(3)]
    view.ViewSize = [960, 360]
    view.UseColorPaletteForBackground = 0
    view.Background = [1.0, 1.0, 1.0]
    view.OrientationAxesVisibility = 0
    view.CenterAxesVisibility = 0
    ResetCamera(view)
    view.CameraParallelProjection = 1
    view.CameraPosition = [center[0], center[1] + 4 * max(lengths), center[2]]
    view.CameraFocalPoint = center
    view.CameraViewUp = [1.0, 0.0, 0.0]
    view.CameraParallelScale = 0.7 * max(lengths[0], lengths[2] / (960 / 360))


def color_surface(reader, view, options):
    """Fix component, range and blue/white/red colors instead of automatic rescaling."""
    display = Show(reader, view)
    display.Representation = "Surface"
    display.Ambient, display.Diffuse, display.Specular = 1.0, 0.0, 0.0
    ColorBy(display, ("CELLS", options.array, options.component))
    lookup = GetColorTransferFunction(options.array)
    low, high = options.limits
    lookup.RGBPoints = [low, 0.23, 0.30, 0.75, (low + high) / 2, 0.87, 0.87, 0.87,
                        high, 0.71, 0.02, 0.15]
    lookup.ColorSpace = "RGB"
    lookup.RescaleTransferFunction(low, high)
    display.SetScalarBarVisibility(view, True)
    configure_legend(lookup, view, options)


def configure_legend(lookup, view, options):
    """Keep the fixed signed-color legend below the field without clipping its title."""
    legend = GetScalarBar(lookup, view)
    legend.Title = f"{options.array}.{options.component} (normalized)"
    legend.ComponentTitle = ""
    legend.TitleColor, legend.LabelColor = [0.1, 0.1, 0.1], [0.1, 0.1, 0.1]
    legend.Orientation = "Horizontal"
    legend.WindowLocation = "Any Location"
    legend.Position = [0.27, 0.06]
    legend.ScalarBarLength = 0.45
    legend.ScalarBarThickness = 12
    legend.TitleFontSize, legend.LabelFontSize = 12, 12


def annotate(view, values):
    """Show actual array-time labels and the orientation convention in the image."""
    label = Text()
    label.Text = (f"potential time {values['potential_time']:.8g}; "
                  f"E time {values['E_time']:.8g}; B time {values['B_time']:.8g}\n"
                  "+x upward; +z rightward; raw normalized fields")
    display = Show(label, view)
    display.WindowLocation = "Upper Left Corner"
    display.FontSize = 12
    display.Color = [0.1, 0.1, 0.1]


def render_provenance(options, values, view):
    """Check camera orientation explicitly; PNG pixels still require visual inspection."""
    client = view.GetClientSideObject()
    window, camera = client.GetRenderWindow(), client.GetRenderer().GetActiveCamera()
    matrix = camera.GetViewTransformMatrix()
    if not (matrix.GetElement(0, 2) > 0.999 and matrix.GetElement(1, 0) > 0.999):
        raise RuntimeError("Camera must have +x upward and +z rightward")
    record = dict(values, array=options.array, component=options.component, limits=options.limits,
                  vtk=vtkVersion.GetVTKVersion(), paraview=str(GetParaViewVersion()),
                  series=str(options.series.resolve()), gui_walkthrough=False,
                  render_window=window.GetClassName(), camera_position=list(view.CameraPosition),
                  camera_up=list(view.CameraViewUp), camera_orientation_checked=True,
                  camera_clipping_range=list(camera.GetClippingRange()),
                  camera_parallel_scale=camera.GetParallelScale(), image_pixels_checked=False)
    record["opengl"] = [line for line in window.ReportCapabilities().splitlines()
                        if line.lower().startswith(("opengl vendor", "opengl renderer", "opengl version"))]
    return record


def main():
    """Save PNG plus provenance, explicitly recording that this is a headless view."""
    options = arguments()
    reader, image = read_frame(options)
    values = metadata(image)
    view = CreateView("RenderView")
    view.ViewTime = options.time
    color_surface(reader, view, options)
    configure_view(view, image)
    annotate(view, values)
    Render(view)
    options.output.parent.mkdir(parents=True, exist_ok=True)
    SaveScreenshot(str(options.output), view, ImageResolution=[960, 360])
    record = render_provenance(options, values, view)
    options.output.with_suffix(".json").write_text(json.dumps(record, indent=2) + "\n")
    print(f"Saved headless image: {options.output}")


if __name__ == "__main__":
    main()
