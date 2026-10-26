#!/usr/bin/env pvpython
"""Read every owned wave cell with VTK and compare to independent discrete/continuum fields.

References: Fallahi, MITHRA 2.0, arXiv:2009.13645v1, (3.8)-(3.17).
VTK reader correctness is checked separately by the distinguishable output fixture.
"""
import argparse
import cmath
import itertools
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path

import vtk
from wave_checks import K, load_rows, check_rows, mode_recurrence, check_phase


def field_metadata(grid):
    names = ['TimeValue', 'potential_time', 'E_time', 'B_time', 'dt', 'step', 'normalization']
    return {name: grid.GetFieldData().GetArray(name).GetTuple1(0) for name in names}

def cell_keys(grid):
    centers = vtk.vtkCellCenters()
    centers.SetInputData(grid)
    centers.Update()
    points = centers.GetOutput().GetPoints()
    keys = []
    for index in range(points.GetNumberOfPoints()):
        xyz = points.GetPoint(index)
        logical = [(xyz[d] - grid.GetOrigin()[d]) / grid.GetSpacing()[d] - 0.5 for d in range(3)]
        key = tuple(round(value) for value in logical)
        assert max(abs(value - integer) for value, integer in zip(logical, key)) < 1e-12
        keys.append(key)
    return keys

def load_frame(path):
    reader = vtk.vtkXMLPImageDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    assert reader.GetErrorCode() == 0 and grid.GetNumberOfCells() > 0
    arrays = {name: [grid.GetCellData().GetArray(name).GetTuple3(i) for i in range(grid.GetNumberOfCells())]
              for name in ('E_raw', 'B_raw')}
    assert all(math.isfinite(value) for array in arrays.values() for vector in array for value in vector)
    return {'origin': grid.GetOrigin(), 'spacing': grid.GetSpacing(),
            'extent': grid.GetExtent(), 'keys': cell_keys(grid),
            'metadata': field_metadata(grid), **arrays}

def expected_fields(frame):
    meta = frame['metadata']
    z = [frame['origin'][2] + (key[2] + 0.5) * frame['spacing'][2] for key in frame['keys']]
    wave = [cmath.exp(1j * K * value) for value in z]
    current, previous = mode_recurrence(int(meta['step']), meta['dt'], frame['spacing'][2])
    electric = [(-(value * (current - previous)).imag / (K * meta['dt']), 0., 0.) for value in wave]
    factor = math.sin(K * frame['spacing'][2]) / (K * frame['spacing'][2])
    magnetic = [(0., factor * (value * current).real, 0.) for value in wave]
    continuum_e = [math.cos(K * (value - meta['E_time'])) for value in z]
    continuum_b = [math.cos(K * (value - meta['B_time'])) for value in z]
    return electric, magnetic, continuum_e, continuum_b

def relative_error(actual, expected):
    denominator = math.fsum(value * value for value in expected)
    assert math.isfinite(denominator) and denominator > 0
    return math.sqrt(math.fsum((left - right) ** 2 for left, right in zip(actual, expected)) / denominator)

def maximum_difference(left, right):
    return max(abs(a - b) for va, vb in zip(left, right) for a, b in zip(va, vb))

def check_frame(frame, nz, time, row):
    meta = frame['metadata']
    assert tuple(frame['extent']) == (0, nz // 8, 0, nz // 8, 0, nz)
    assert tuple(frame['origin']) == (0, 0, 0) and tuple(frame['spacing']) == (1 / nz,) * 3
    assert meta['potential_time'] == meta['B_time'] == meta['TimeValue'] == time
    assert meta['dt'] == 1 / (2 * nz) and meta['E_time'] == time - meta['dt'] / 2
    assert meta['step'] == row['step'] and meta['normalization'] == 1
    assert meta['step'] * meta['dt'] == time
    keys = frame['keys']
    assert sorted(keys) == list(itertools.product(range(nz // 8), range(nz // 8), range(nz)))
    electric, magnetic, continuum_e, continuum_b = expected_fields(frame)
    discrete = max(maximum_difference(frame['E_raw'], electric), maximum_difference(frame['B_raw'], magnetic))
    errors = [relative_error([v[0] for v in frame['E_raw']], continuum_e),
              relative_error([v[1] for v in frame['B_raw']], continuum_b)]
    assert discrete <= 1e-10 and max(errors) <= (0.003 if nz == 128 else 0.01)
    assert abs(errors[0] - row['relativeE']) < 1e-12 and abs(errors[1] - row['relativeB']) < 1e-12
    return {'discrete_max': float(discrete), 'relativeE': errors[0], 'relativeB': errors[1]}

def check_piece_extents(path, nz, ranks):
    pieces = ET.parse(path).getroot().find('PImageData').findall('Piece')
    expected = [(0, nz // 8, 0, nz // 8, rank * nz // ranks, (rank + 1) * nz // ranks)
                for rank in range(ranks)]
    actual = [tuple(map(int, piece.attrib['Extent'].split())) for piece in pieces]
    assert actual == expected

def check_case(run):
    nz = 128 if run['case'] == 'fine-serial' else 64
    rows = load_rows(run['csv'])
    check_rows(rows, nz)
    directory = Path(run['output'])
    series = json.loads((directory / 'wave.pvti.series').read_text())
    assert [item['time'] for item in series['files']] == [0, .25, .5, .75, 1]
    frames, errors = [], []
    for entry in series['files']:
        check_piece_extents(directory / entry['name'], nz, 2 if run['case'] == 'coarse-mpi2' else 1)
        frame = load_frame(directory / entry['name'])
        errors.append(check_frame(frame, nz, entry['time'], rows[int(frame['metadata']['step'])]))
        frames.append(frame)
    return rows, frames, errors, check_phase(rows, nz)

def compare_mpi(serial, parallel):
    assert len(serial[0]) == len(parallel[0])
    diagnostics = max(abs(left[key] - right[key]) for left, right in zip(serial[0], parallel[0])
                      for key in left if left[key] is not None)
    fields = 0.0
    for left, right in zip(serial[1], parallel[1]):
        assert sorted(left['keys']) == sorted(right['keys'])
        for name in ('E_raw', 'B_raw'):
            ordered_l = [value for _, value in sorted(zip(left['keys'], left[name]))]
            ordered_r = [value for _, value in sorted(zip(right['keys'], right[name]))]
            fields = max(fields, maximum_difference(ordered_l, ordered_r))
    assert diagnostics <= 1e-11 and fields <= 1e-11
    return {'diagnostics_max_difference': diagnostics, 'fields_max_difference': fields}


def summarize(cases):
    coarse, fine = cases['coarse-serial'], cases['fine-serial']
    ratios = {name: max(e[name] for e in coarse[2]) / max(e[name] for e in fine[2])
              for name in ('relativeE', 'relativeB')}
    assert all(3.0 <= ratio <= 5.0 for ratio in ratios.values())
    return {'vtk_version': vtk.vtkVersion.GetVTKVersion(), 'refinement_ratios': ratios,
            'mpi_comparison': compare_mpi(coarse, cases['coarse-mpi2']),
            'cases': {name: {'frames': case[2], 'maximum_discrete_error_all_steps':
                            max(row['discreteError'] for row in case[0]),
                            'final_phase_speed': case[0][-1]['phaseSpeed'],
                            'phase_coefficient_max_error': case[3]}
                      for name, case in cases.items()}, 'passed': True}


def main():
    if not __debug__:
        raise RuntimeError('Run verification without Python optimization.')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--summary', type=Path, required=True)
    args = parser.parse_args()
    record = json.loads(args.manifest.read_text())
    assert record['source_unchanged'] and len(record['runs']) == 3
    assert all(run['exit_code'] == 0 for run in record['runs'])
    summary = summarize({run['case']: check_case(run) for run in record['runs']})
    summary['source_hashes'] = record['source_before']
    args.summary.write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
