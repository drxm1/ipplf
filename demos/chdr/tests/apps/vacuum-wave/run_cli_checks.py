#!/usr/bin/env python3
## @file run_cli_checks.py
# @brief CTest entry chdr.vacuum.cli: 13 YAML/command-line runs of chdr-vacuum-wave.
# @details Six runs must exit with status 0: coarse/fine initialization, initialization
# with IPPL's timer/debug options, and zero/one-step field output. Six invalid selectors
# or YAML inputs on one rank and run.steps: 999 on two ranks must exit with a status other
# than 0 and 124 and print "rank " (one rank) or "exceeds the one-period benchmark"
# (two ranks) to standard error. Each run is limited to 120 s. The
# source hashes before and after the runs must be equal. The record cli-checks.json is the
# manifest that chdr.vacuum.readback.partial reads.
"""Exercise the chdr-vacuum-wave command line, including a coordinated MPI rejection.

This checks parsing and process exit status. Mid-run numerical gates are tested
in GoogleTest. The output-enabled zero/one-step cases supply a separate actual
reader test through cli-checks.json; their process success alone is not readback.
"""
import argparse
import json
import os
from pathlib import Path
import tempfile

from run_wave_checks import fingerprint, invoke, write_input


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--mpiexec', required=True)
    parser.add_argument('--numproc-flag', default='-n')
    parser.add_argument('--mpi-before', action='append', default=[])
    parser.add_argument('--mpi-after', action='append', default=[])
    return parser.parse_args()


def yaml_inputs(directory):
    """Prepare exact step values, including the invalid values exercised by the parser."""
    settings = [('initial', 'coarse', 0), ('fine-initial', 'fine', 0),
                ('output-initial', 'coarse', 0), ('output-first-step', 'coarse', 1),
                ('negative-steps', 'coarse', -1), ('suffix-steps', 'coarse', '1x'),
                ('overflow-steps', 'coarse', 999999999999999999999),
                ('period-limit-mpi2', 'coarse', 999)]
    inputs = {}
    for name, case, steps in settings:
        output = directory / name if name in ('output-initial', 'output-first-step') else None
        options = write_input(directory / (name + '.yaml'), case=case, steps=steps, output=output)
        inputs[name] = (options, str(output) if output is not None else None)
    return inputs


def cases(args, directory):
    mpi = [args.mpiexec, args.numproc_flag, '2'] + args.mpi_before
    inputs = yaml_inputs(directory)
    initial = inputs['initial'][0]
    return [('initial', [], initial, True, None),
            ('fine-initial', [], inputs['fine-initial'][0], True, None),
            ('timer-option', [], ['--timer-fences', 'on'] + initial, True, None),
            ('debug-option', [], ['--debug'] + initial, True, None),
            ('removed-fine-option', [], ['--fine'], False, None),
            ('missing-input-file', [], ['--input'], False, None),
            ('output-initial', [], inputs['output-initial'][0], True, inputs['output-initial'][1]),
            ('output-first-step', [], inputs['output-first-step'][0], True,
             inputs['output-first-step'][1]),
            ('missing-input', [], [], False, None),
            ('negative-steps', [], inputs['negative-steps'][0], False, None),
            ('suffix-steps', [], inputs['suffix-steps'][0], False, None),
            ('overflow-steps', [], inputs['overflow-steps'][0], False, None),
            ('period-limit-mpi2', mpi, args.mpi_after + inputs['period-limit-mpi2'][0], False, None)]


def run_case(args, directory, specification):
    name, prefix, options, success, output = specification
    argv = prefix + [str(args.binary)] + options
    out_path, err_path = directory / (name + '.stdout'), directory / (name + '.stderr')
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS='1', OMP_PROC_BIND='false')
    with out_path.open('w') as out, err_path.open('w') as err:
        with (directory / 'stdin.txt').open('r') as stdin:
            status = invoke(argv, directory, environment, out, err, stdin)
    reason = 'exceeds the one-period benchmark' if prefix else 'rank '
    expected = status == 0 if success else status not in (0, 124) and reason in err_path.read_text()
    return {'case': name, 'argv': argv, 'cwd': str(directory), 'exit_code': status,
            'expected_success': success, 'passed': expected,
            'stdout': str(out_path), 'stderr': str(err_path),
            'environment': {'OMP_NUM_THREADS': '1', 'OMP_PROC_BIND': 'false'},
            'stdin': 'x\n', 'output': output}


def main():
    args = arguments()
    args.binary, args.results = args.binary.resolve(), args.results.resolve()
    args.results.mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix='cli-', dir=args.results))
    (directory / 'stdin.txt').write_text('x\n')
    record = {'source_before': fingerprint(args.binary)}
    record['runs'] = [run_case(args, directory, case) for case in cases(args, directory)]
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
