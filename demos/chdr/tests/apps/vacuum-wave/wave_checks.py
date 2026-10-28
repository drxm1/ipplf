"""Independent benchmark checks; no imports from the C++ diagnostic implementation.

The scalar recurrence and signed phase follow our axial substitution into
Fallahi, MITHRA 2.0, arXiv:2009.13645v1, equations (3.8)-(3.17), printed pp.12-14.
The benchmark uses c=1, wavelength1, uniform grids and the declared raw E/B times.
"""
import cmath
import csv
import math
from pathlib import Path

K = 2.0 * math.pi


def load_rows(path):
    with Path(path).open() as stream:
        rows = list(csv.DictReader(stream))
    return [{key: float(value) if value else None for key, value in row.items()} for row in rows]

def check_rows(rows, nz):
    assert len(rows) == 2 * nz + 1
    for step, row in enumerate(rows):
        assert row['step'] == step and row['count'] == (nz // 8) ** 2 * nz
        assert row['dt'] == 1 / (2 * nz)
        assert row['timeA'] == step / (2 * nz) and row['timeB'] == row['timeA']
        assert row['timeE'] == row['timeA'] - 1 / (4 * nz)
        assert all(value is None or math.isfinite(value) for value in row.values())
        assert row['discreteError'] <= 1e-10 and row['zeroError'] <= 1e-12
        assert row['sourceError'] <= 1e-12 and row['phaseAmplitude'] > 0.5
        if step:
            assert -math.pi < row['phaseIncrement'] < 0
            assert abs(row['phaseSpeed'] - 1.0) <= 0.002
        else:
            assert row['phaseSpeed'] is None

def mode_recurrence(step, dt, spacing):
    current, previous = 1.0 + 0j, cmath.exp(1j * K * dt)
    coefficient = 2.0 - 4.0 * (dt / spacing) ** 2 * math.sin(K * spacing / 2.0) ** 2
    for _ in range(step):
        current, previous = coefficient * current - previous, current
    return current, previous

def check_phase(rows, nz):
    previous, accumulated, maximum = None, 0.0, 0.0
    for step, row in enumerate(rows):
        expected = -1j * mode_recurrence(step, 1 / (2 * nz), 1 / nz)[0]
        difference = abs(complex(row['phaseReal'], row['phaseImag']) - expected)
        assert difference <= 1e-10
        maximum = max(maximum, difference)
        if step:
            increment = cmath.phase(expected / previous)
            accumulated += increment
            assert abs(increment - row['phaseIncrement']) <= 1e-11
            assert abs(row['phaseSpeed'] + accumulated / (K * row['timeA'])) <= 1e-11
        previous = expected
    return maximum


def checkpoint_errors(rows, nz):
    selected = [row for row in rows if row['step'] % (nz // 2) == 0]
    assert len(selected) == 5
    maxima = {name: max(row[name] for row in selected)
              for name in ('relativeE', 'relativeB')}
    assert all(0.0 <= value <= (0.003 if nz == 128 else 0.01)
               for value in maxima.values())
    return maxima


def compare_csv_runs(runs):
    rows = {run['case']: load_rows(run['csv']) for run in runs}
    for name, values in rows.items():
        nz = 128 if name == 'fine-serial' else 64
        check_rows(values, nz)
        check_phase(values, nz)
    coarse = checkpoint_errors(rows['coarse-serial'], 64)
    fine = checkpoint_errors(rows['fine-serial'], 128)
    checkpoint_errors(rows['coarse-mpi2'], 64)
    ratios = {name: coarse[name] / fine[name] for name in coarse}
    assert all(3.0 <= value <= 5.0 for value in ratios.values())
    differences = [abs(a[key] - b[key])
                   for a, b in zip(rows['coarse-serial'], rows['coarse-mpi2'])
                   for key in a if a[key] is not None]
    assert max(differences) <= 1e-11
    return {'refinement_ratios': ratios, 'mpi_max_difference': max(differences)}
