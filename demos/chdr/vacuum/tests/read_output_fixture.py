"""Read the independent 72-cell sentinel using installed VTK and ParaView readers.

Expected physical coordinates/data come from the declared fixture, not the writer.
See ../ViewingGuide.org, output-check section; VTK XML and ParaView
temporal-file-series official specifications are cited in ../VacuumOutput.h.
This verifies reader interoperability, not a GUI viewing walkthrough.
"""
import json
import math
import struct
import sys
from pathlib import Path

from vtkmodules.vtkCommonCore import vtkVersion
from vtkmodules.vtkIOXML import vtkXMLImageDataReader, vtkXMLPImageDataReader
from paraview import servermanager
from paraview.simple import Delete, GetParaViewVersion, OpenDataFile

ORIGIN = (-2.0, 3.0, 5.0)
SPACING = (0.5, 1.25, 2.0)
TIMES = (0.0, 0.125, 0.375)
STEPS = (0, 1, 3)


def close(actual, expected):
    """Require finite numeric agreement independent of serialization formatting."""
    assert math.isfinite(actual) and math.isclose(actual, expected, rel_tol=1e-14, abs_tol=1e-14), (actual, expected)


def exact(actual, expected):
    """Require equal binary64 representations after the actual VTK reader conversion."""
    assert math.isfinite(actual) and struct.pack(">d", actual) == struct.pack(">d", expected), (actual.hex(), float(expected).hex())


def expected_vectors(i, j, k, time):
    """Independent Python nextafter/pi construction of the declared precision sentinel."""
    electric = (math.nextafter(1 + i + 100 * time, math.inf),
                math.nextafter(math.pi * (10 + j), math.inf),
                math.nextafter((100 + k) / math.pi, 0.0))
    magnetic = (math.nextafter(-10 - 2 * i, -math.inf),
                math.nextafter(-math.pi * (100 + 3 * j), -math.inf),
                math.nextafter((-1000 - 5 * k) / math.pi, 0.0))
    return electric, magnetic


def check_metadata(image, step, time):
    """Verify dataset/array timestamps, step and the declared normalization code."""
    values = {"TimeValue": time, "potential_time": time,
              "E_time": math.nextafter(time - 0.0625, -math.inf),
              "B_time": math.nextafter(time + 0.03125, math.inf),
              "dt": 0.125, "step": step, "normalization": 1}
    for name, expected in values.items():
        array = image.GetFieldData().GetArray(name)
        assert array is not None, name
        assert array.GetNumberOfTuples() == 1 and array.GetNumberOfComponents() == 1, name
        exact(array.GetTuple1(0), expected)


def physical_index(image, cell):
    """Infer a global logical index from the reader's actual physical cell bounds."""
    bounds = [0.0] * 6
    image.GetCellBounds(cell, bounds)
    centers = [(bounds[2 * d] + bounds[2 * d + 1]) / 2 for d in range(3)]
    indices = tuple(round((centers[d] - ORIGIN[d]) / SPACING[d] - 0.5) for d in range(3))
    for d in range(3):
        close(centers[d], ORIGIN[d] + (indices[d] + 0.5) * SPACING[d])
        close(bounds[2 * d + 1] - bounds[2 * d], SPACING[d])
    return indices


def check_cells(image, first_z, end_z, time):
    """Check every signed vector and require exactly one cell at each expected center."""
    assert image.GetNumberOfCells() == 12 * (end_z - first_z)
    arrays = [image.GetCellData().GetArray(name) for name in ("E_raw", "B_raw")]
    for array in arrays:
        assert array is not None and array.GetNumberOfComponents() == 3
        assert array.GetDataTypeAsString() == "double"
        assert array.GetNumberOfTuples() == image.GetNumberOfCells()
    found = set()
    for cell in range(image.GetNumberOfCells()):
        i, j, k = physical_index(image, cell)
        assert (i, j, k) not in found
        found.add((i, j, k))
        expected = expected_vectors(i, j, k, time)
        for array, vector in zip(arrays, expected):
            for actual, target in zip(array.GetTuple3(cell), vector):
                exact(actual, target)
    assert found == {(i, j, k) for k in range(first_z, end_z) for j in range(3) for i in range(4)}


def check_image(image, first_z, end_z, step, time):
    """Validate extent, global origin/spacing, fields and metadata after actual reading."""
    assert image.GetExtent() == (0, 4, 0, 3, first_z, end_z), image.GetExtent()
    for actual, expected in zip(image.GetOrigin(), ORIGIN):
        close(actual, expected)
    for actual, expected in zip(image.GetSpacing(), SPACING):
        close(actual, expected)
    assert image.GetPointData().GetNumberOfArrays() == 0
    check_metadata(image, step, time)
    check_cells(image, first_z, end_z, time)


def read_xml(path, parallel):
    """Use VTK's serial/parallel image reader and turn reader errors into failures."""
    reader = vtkXMLPImageDataReader() if parallel else vtkXMLImageDataReader()
    errors = []
    reader.AddObserver("ErrorEvent", lambda *args: errors.append("VTK ErrorEvent"))
    reader.SetFileName(str(path))
    reader.Update()
    assert not errors and reader.GetErrorCode() == 0, (path, errors, reader.GetErrorCode())
    return reader.GetOutput()


def check_xml_frames(directory, ranks):
    """Read each standalone piece and each PVTI, including cells on both interface sides."""
    for step, time in zip(STEPS, TIMES):
        image = read_xml(directory / f"wave-{step}.pvti", True)
        check_image(image, 0, 6, step, time)
        for rank in range(ranks):
            image = read_xml(directory / f"wave-{step}-rank-{rank}.vti", False)
            check_image(image, rank * (6 // ranks), (rank + 1) * (6 // ranks), step, time)


def check_series(directory):
    """Ask ParaView to expose exact unequal series times and read every selected state."""
    path = directory / "wave.pvti.series"
    content = json.loads(path.read_text())
    assert content["file-series-version"] == "1.0"
    assert content["files"] == [{"name": f"wave-{s}.pvti", "time": t} for s, t in zip(STEPS, TIMES)]
    reader = OpenDataFile(str(path))
    assert reader is not None
    assert list(reader.TimestepValues) == list(TIMES), reader.TimestepValues
    for step, time in zip(STEPS, TIMES):
        reader.UpdatePipeline(time)
        check_image(servermanager.Fetch(reader), 0, 6, step, time)
    Delete(reader)


def check_failure(directory, mode):
    """Previously complete series remains readable; a failed frame is never listed."""
    content = json.loads((directory / "wave.pvti.series").read_text())
    assert content["files"] == [{"name": "wave-0.pvti", "time": 0.0}]
    assert (directory / "wave-0.pvti").is_file()
    if mode != "series":
        assert not (directory / "wave-1.pvti").exists()
    check_image(read_xml(directory / "wave-0.pvti", True), 0, 6, 0, 0)


def main():
    """Run a successful-readback case or inspect an intentionally failed publication."""
    directory, argument = Path(sys.argv[1]), sys.argv[2]
    print(f"Reader versions: VTK {vtkVersion.GetVTKVersion()}, ParaView {GetParaViewVersion()}")
    if argument in ("piece", "wrapper", "series", "nonfinite"):
        check_failure(directory, argument)
        print(f"PASS retained complete frame after {argument} publication failure")
    else:
        check_xml_frames(directory, int(argument))
        check_series(directory)
        print(f"PASS all 72 cells: bit-exact nextafter/pi vectors and explicit times, extents, halo exclusion: {argument} piece(s)")


if __name__ == "__main__":
    main()
