#!/usr/bin/env python3
"""Exercise the actual driver CLI, including a coordinated MPI rejection.

This checks parsing and process exit status. Mid-run numerical gates are tested
in GoogleTest. The output-enabled zero/one-step cases supply a separate actual
reader test through cli-checks.json; their process success alone is not readback.
"""
import argparse
import json
import os
from pathlib import Path
import tempfile

from run_wave_checks import fingerprint, invoke


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--mpiexec', required=True)
    parser.add_argument('--numproc-flag', default='-n')
    parser.add_argument('--mpi-before', action='append', default=[])
    parser.add_argument('--mpi-after', action='append', default=[])
    return parser.parse_args()


def cases(args):
    mpi = [args.mpiexec, args.numproc_flag, '2'] + args.mpi_before
    return [('initial', [], ['--steps', '0'], True),
            ('fine-initial', [], ['--fine', '--steps', '0'], True),
            ('timer-option', [], ['--timer-fences', 'on', '--steps', '0'], True),
            ('debug-option', [], ['--debug', '--steps', '0'], True),
            ('unknown', [], ['--unrecognized-vacuum-option'], False),
            ('missing-output', [], ['--output'], False),
            ('output-initial', [], ['--steps', '0', '--output', 'output-initial'], True),
            ('output-first-step', [], ['--steps', '1', '--output', 'output-first-step'], True),
            ('missing-steps', [], ['--steps'], False),
            ('negative-steps', [], ['--steps', '-1'], False),
            ('suffix-steps', [], ['--steps', '1x'], False),
            ('overflow-steps', [], ['--steps', '999999999999999999999'], False),
            ('period-limit-mpi2', mpi, args.mpi_after + ['--steps', '999'], False)]


def output_path(directory, options):
    """Resolve a supplied field directory against the recorded driver working directory."""
    if '--output' not in options:
        return None
    index = options.index('--output') + 1
    return str(directory / options[index]) if index < len(options) else None


def run_case(args, directory, specification):
    name, prefix, options, success = specification
    argv = prefix + [str(args.binary)] + options
    out_path, err_path = directory / (name + '.stdout'), directory / (name + '.stderr')
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS='1', OMP_PROC_BIND='false')
    with out_path.open('w') as out, err_path.open('w') as err:
        with (directory / 'stdin.txt').open('r') as stdin:
            status = invoke(argv, directory, environment, out, err, stdin)
    reason = '--steps exceeds the one-period benchmark' if prefix else 'vacuum rank '
    expected = status == 0 if success else status not in (0, 124) and reason in err_path.read_text()
    return {'case': name, 'argv': argv, 'cwd': str(directory), 'exit_code': status,
            'expected_success': success, 'passed': expected,
            'stdout': str(out_path), 'stderr': str(err_path),
            'environment': {'OMP_NUM_THREADS': '1', 'OMP_PROC_BIND': 'false'},
            'stdin': 'x\n', 'output': output_path(directory, options)}


def main():
    args = arguments()
    args.binary, args.results = args.binary.resolve(), args.results.resolve()
    args.results.mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix='cli-', dir=args.results))
    (directory / 'stdin.txt').write_text('x\n')
    record = {'source_before': fingerprint(args.binary)}
    record['runs'] = [run_case(args, directory, case) for case in cases(args)]
    record['source_after'] = fingerprint(args.binary)
    record['source_unchanged'] = record['source_before'] == record['source_after']
    record['passed'] = record['source_unchanged'] and all(r['passed'] for r in record['runs'])
    (directory / 'cli-checks.json').write_text(json.dumps(record, indent=2) + '\n')
    (args.results / 'cli-checks.json').write_text(json.dumps(record, indent=2) + '\n')
    for run in record['runs']:
        print(run['case'], 'PASS' if run['passed'] else 'FAIL', run['exit_code'])
    return 0 if record['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
