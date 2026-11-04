#!/usr/bin/env python3
"""Check YAML/SI field and spectrum publication, without a new accuracy claim.

Each of chdr-vacuum-wave, interface-movies and planar-interface runs one case twice on one
rank: once in code units (4 steps, or the default run of the plane case) and once for the
equal SI duration, with a length scale of 2 mm and a charge-density scale of 3e-6 C/m^3.
Every SI number (record factors, VTK origin, spacing and time labels, E_raw, B_raw, CSV
columns, plane spectra) must equal the code-unit number times its factor within a relative
tolerance of 2e-12. File names, step lists, extents and record fields must be equal, and the
time step of the run record must equal the VTK dt exactly. CTest entry: chdr.input.runs.si,
which exists only when CMake finds Python 3.

Scale substitutions: fitzpatrick2014ClassicalElectromagnetism, Eqs. (1.1), (1.3),
(1.11)-(1.16), (4.70), (6.85); constants: mohr2025Codata2022, Table XXXIII.
The expectations below are independent of the C++ conversion helpers.
"""
import csv
import json
import math
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

LENGTH, DENSITY, LIGHT_SPEED, EPSILON = 0.002, 0.000003, 299792458.0, 8.8541878188e-12
TIME = LENGTH / LIGHT_SPEED
ELECTRIC = DENSITY * LENGTH / EPSILON
MAGNETIC = ELECTRIC / LIGHT_SPEED
ENERGY = EPSILON * ELECTRIC ** 2 * LENGTH ** 3


def close(actual, expected, label):
    assert math.isfinite(float(actual)) and math.isclose(float(actual), float(expected),
        rel_tol=2e-12, abs_tol=1e-280), (label, actual, expected)


def run(binary, directory, movie, duration=None, plane=False):
    directory.mkdir()
    output = directory / 'output'
    units = ('system: normalized' if duration is None else
             f'system: si, length_m: {LENGTH}, charge_density_C_m3: {DENSITY}')
    timing = 'steps: 4' if duration is None else f'duration: {duration:.17g}'
    case = 'plane' if plane else 'brick' if movie else 'coarse'
    document = (f'schema_version: 1\ncase: {case}\nunits: {{{units}}}\n'
                f'output: {{directory: {json.dumps(str(output))}}}\n')
    if not plane or duration is not None:
        document += f'run: {{{timing}}}\n'
    if not plane:
        document += ('observations:\n  fields: {every: 3, include_final: true}\n'
                     f'  diagnostics: {{every: {2 if movie else 1}, include_final: true}}\n')
    path = directory / 'input.yaml'
    path.write_text(document + ('parameters: {smoke: true}\n' if movie else ''))
    result = subprocess.run([binary, '--input', str(path)], capture_output=True, text=True,
        env=dict(os.environ, OMP_NUM_THREADS='1', OMP_PROC_BIND='false'), timeout=90)
    assert result.returncode == 0, result.stdout + result.stderr
    return output, json.loads((output / 'run.json').read_text()), result.stdout


def vtk(path, binary):
    data = path.read_bytes()
    marker = b'<AppendedData encoding="raw">_'
    header, payload = data.split(marker, 1) if marker in data else (data, b'')
    tree = ET.fromstring(header + b'</VTKFile>' if payload else header)
    assert bool(payload) == binary and tree.attrib['byte_order'] == 'LittleEndian'
    arrays = {}
    for node in tree.findall('.//CellData/DataArray'):
        if payload:
            offset = int(node.attrib['offset'])
            size, = struct.unpack_from('<Q', payload, offset)
            arrays[node.attrib['Name']] = struct.unpack_from(f'<{size // 8}d', payload, offset + 8)
        else:
            arrays[node.attrib['Name']] = list(map(float, node.text.split()))
    metadata = {node.attrib['Name']: list(map(float, node.text.split()))
                for node in tree.findall('.//FieldData/DataArray')}
    return tree.find('ImageData').attrib, metadata, arrays


def check_frames(normal, physical, record, movie):
    series = [json.loads((path / 'wave.pvti.series').read_text())['files']
              for path in (normal, physical)]
    assert [row['name'] for row in series[0]] == [f'wave-{k}.pvti' for k in (0, 3, 4)]
    assert [row['name'] for row in series[1]] == [row['name'] for row in series[0]]
    dt, cells = record['effective']['solver']['time_step'], record['effective']['box']['cells']
    for left, right in zip(*series):
        close(right['time'], left['time'] * TIME, 'series time')
        frames = []
        for directory, row in ((normal, left), (physical, right)):
            wrapper = ET.parse(directory / row['name'])
            frames.append(vtk(directory / wrapper.find('.//Piece').attrib['Source'], movie))
        (geometry, metadata, arrays), (scaled_geometry, scaled_metadata, scaled_arrays) = frames
        assert geometry['WholeExtent'] == scaled_geometry['WholeExtent']
        extent = list(map(int, geometry['WholeExtent'].split()))
        assert [extent[2*d+1] - extent[2*d] for d in range(3)] == cells
        for key in ('Origin', 'Spacing'):
            for a, b in zip(map(float, geometry[key].split()), map(float, scaled_geometry[key].split())):
                close(b, a * LENGTH, key)
        step = metadata['step'][0]
        assert metadata['dt'] == [dt], 'run record must preserve the round-trip VTK timestep'
        for key, expected in {'dt': dt, 'TimeValue': step*dt, 'potential_time': step*dt,
                              'E_time': (step-0.5)*dt, 'B_time': step*dt}.items():
            close(metadata[key][0], expected, key)
            close(scaled_metadata[key][0], expected * TIME, key)
        assert metadata['normalization'] == [1] and scaled_metadata['normalization'] == [0]
        assert metadata.get('component_offsets') == scaled_metadata.get('component_offsets')
        for key, factor in (('E_raw', ELECTRIC), ('B_raw', MAGNETIC)):
            assert len(arrays[key]) == len(scaled_arrays[key]) == 3 * math.prod(cells)
            for a, b in zip(arrays[key], scaled_arrays[key]):
                close(b, a * factor, key)


def check_csv(left, right, movie):
    rows = [list(csv.DictReader(line for line in text.splitlines() if not line.startswith('#')))
            for text in (left, right)]
    assert [int(row['step']) for row in rows[0]] == ([0, 2, 4] if movie else list(range(5)))
    assert len(rows[0]) == len(rows[1])
    for normal, physical in zip(*rows):
        for key, value in normal.items():
            if value == '':
                assert physical[key] == ''
                continue
            factor = TIME if key in ('time', 'dt', 'timeA', 'timeE', 'timeB') else 1
            if not movie and key in ('energyE', 'energyB', 'energy'):
                factor = ENERGY
            if key == 'phaseSpeed':
                factor = LIGHT_SPEED
            close(physical[key], float(value) * factor, key)
        if movie:
            close(physical['field_energy_J'], float(normal['field_energy']) * ENERGY, 'energy')


def check_records(record, si_record, steps):
    for key in ('solver', 'box', 'observations', 'completed_solves', 'requested_steps'):
        assert record['effective'][key] == si_record['effective'][key], key
    assert record['effective']['completed_solves'] == record['effective']['requested_steps'] == steps
    for item, system in ((record, 'normalized'), (si_record, 'si')):
        assert item['status'] == 'completed' and item['execution']['ranks'] == 1
        assert item['input']['units']['system'] == item['output_factors']['system'] == system
        assert item['effective']['value_system'] == 'normalized'
    for key, factor in {'length': LENGTH, 'time': TIME, 'electric_field': ELECTRIC,
                        'magnetic_field': MAGNETIC, 'physical_field_energy': ENERGY}.items():
        close(si_record['output_factors'][key], factor, key)


def check_csv_units(stdout, si_stdout, record, si_record):
    assert not stdout.startswith('#') and si_stdout.startswith('# ')
    csv_dt = float(next(csv.DictReader(stdout.splitlines()))['dt'])
    assert record['effective']['solver']['time_step'] == csv_dt
    labels = dict(item.split('=', 1) for item in si_stdout.splitlines()[0][2:].split())
    expected = {'time_units': 's', 'phase_speed_units': 'm/s',
                'physical_energy_columns': 'energyE,energyB,energy',
                'physical_energy_units': 'J', 'energy_invariant_units': 'normalized'}
    for key, value in expected.items():
        assert labels[key] == si_record['effective']['output']['stdout_csv'][key] == value
    for key, value in (('length_m', LENGTH), ('charge_density_C_m3', DENSITY)):
        close(labels[key], value, key)
        close(si_record['effective']['output']['stdout_csv'][key], value, key)


def check_case(binary, directory, movie):
    normal, record, stdout = run(binary, directory / 'normalized', movie)
    physical, si_record, si_stdout = run(binary, directory / 'si', movie,
                                        4 * record['effective']['solver']['time_step'] * TIME)
    check_records(record, si_record, 4)
    assert record['effective']['observations']['fields']['every'] == 3
    assert record['effective']['observations']['diagnostics']['every'] == (2 if movie else 1)
    assert record['effective']['observations']['diagnostics']['include_final'] == 1
    check_frames(normal, physical, record, movie)
    check_csv((normal / 'energy.csv').read_text() if movie else stdout,
              (physical / 'energy.csv').read_text() if movie else si_stdout, movie)
    if not movie:
        check_csv_units(stdout, si_stdout, record, si_record)


def check_plane(binary, directory):
    normal, record, _ = run(binary, directory / 'normalized', False, plane=True)
    duration = record['effective']['completed_solves'] * record['effective']['solver']['time_step']
    physical, si_record, _ = run(binary, directory / 'si', False, duration * TIME, plane=True)
    rows = [list(csv.DictReader((path / 'spectra.csv').read_text().splitlines()))
            for path in (normal, physical)]
    assert len(rows[0]) == len(rows[1]) == 3
    for left, right in zip(*rows):
        close(right['angular_frequency_rad_per_s'], float(left['angular_frequency']) / TIME,
              'angular frequency')
        for key in ('r_real', 'r_imag', 't_real', 't_imag', 'omega_code',
                    'incident_real', 'incident_imag', 'total_real', 'total_imag',
                    'transmitted_real', 'transmitted_imag'):
            close(right[key], left[key], key)
        assert right['valid'] == left['valid'] == '1'
        assert right['reason'] == left['reason'] == ''
    check_records(record, si_record, record['effective']['requested_steps'])
    for item, system in ((record, 'normalized'), (si_record, 'si')):
        assert item['input']['case'] == 'plane'
        spectra = item['effective']['output']['spectra']
        assert spectra['amplitudes_units'] == 'dimensionless'
        assert spectra['angular_frequency_units'] == ('rad/s' if system == 'si' else 'normalized')
        assert spectra['raw_units'] == 'normalized'
        assert spectra['raw_phase_origin'] == 'probe_coordinates'
        assert spectra['coefficient_phase_origin'] == 'requested_face'
        assert item['effective']['plane']['schema_version'] == 1
        assert item['effective']['plane']['central_bin'] == 0
        assert item['effective']['completed_solves'] > 0


if __name__ == '__main__':
    assert not sys.flags.optimize, 'Run checks without Python optimization'
    assert len(sys.argv) == 4, 'Pass vacuum-wave, interface-movies and planar-interface paths'
    with tempfile.TemporaryDirectory(prefix='chdr-si-') as temporary:
        for index, binary in enumerate(sys.argv[1:3]):
            directory = Path(temporary) / str(index)
            directory.mkdir()
            check_case(binary, directory, bool(index))
        directory = Path(temporary) / 'plane'
        directory.mkdir()
        check_plane(sys.argv[3], directory)
    print('YAML/SI vacuum ASCII, movie binary and plane spectrum comparisons passed')
