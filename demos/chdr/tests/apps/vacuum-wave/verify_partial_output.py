#!/usr/bin/env pvpython
"""Read the CLI's output-enabled zero/one-step cases with the actual VTK reader.

The full-period checker remains unchanged. We reuse its independent frame and
phase checks, derived there from Fallahi, MITHRA 2.0, arXiv:2009.13645v1,
Eqs. (3.8)-(3.17), for the fixed coarse grid and the declared raw E/B times.
"""
import argparse
import json
from pathlib import Path

import vtk
from verify_wave_output import check_frame, check_piece_extents, load_frame
from wave_checks import check_phase, load_rows


def check_partial(run, last_step):
    """Require exactly the initial/final saved states and inspect every owned cell."""
    assert run['exit_code'] == 0 and run['passed'] and run['output']
    rows = load_rows(run['stdout'])
    assert len(rows) == last_step + 1
    assert [row['step'] for row in rows] == list(range(last_step + 1))
    assert [row['timeA'] for row in rows] == [step / 128 for step in range(last_step + 1)]
    directory = Path(run['output'])
    series = json.loads((directory / 'wave.pvti.series').read_text())
    assert series['file-series-version'] == '1.0'
    expected = [{'name': f'wave-{step}.pvti', 'time': step / 128}
                for step in range(last_step + 1)]
    assert series['files'] == expected
    errors = []
    for step, entry in enumerate(series['files']):
        path = directory / entry['name']
        check_piece_extents(path, 64, 1)
        errors.append(check_frame(load_frame(path), 64, entry['time'], rows[step]))
    return {'steps': last_step, 'frames': errors,
            'phase_coefficient_max_error': check_phase(rows, 64)}


def main():
    """Require successful CLI evidence and persist actual-reader results separately."""
    if not __debug__:
        raise RuntimeError('Run verification without Python optimization.')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--summary', type=Path, required=True)
    args = parser.parse_args()
    record = json.loads(args.manifest.read_text())
    assert record['source_unchanged'] and record['passed']
    runs = {run['case']: run for run in record['runs']}
    cases = {'output-initial': 0, 'output-first-step': 1}
    summary = {'vtk_version': vtk.vtkVersion.GetVTKVersion(), 'passed': True,
               'cases': {name: check_partial(runs[name], step) for name, step in cases.items()},
               'source_hashes': record['source_before']}
    args.summary.write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
