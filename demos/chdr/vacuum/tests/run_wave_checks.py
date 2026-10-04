#!/usr/bin/env python3
"""Run the fixed coarse, fine and two-rank benchmark; verify all CSV rows.

Expected quantities and tolerances follow README.md and wave_checks.py's cited
reference. This unit compares global diagnostics across ranks; field-by-field
MPI agreement and VTK/readback acceptance remain deferred to the output unit.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import tempfile

from wave_checks import compare_csv_runs


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--mpiexec', required=True)
    parser.add_argument('--numproc-flag', default='-n')
    parser.add_argument('--mpi-before', action='append', default=[])
    parser.add_argument('--mpi-after', action='append', default=[])
    return parser.parse_args()


def fingerprint(binary):
    source = Path(__file__).resolve().parent.parent
    paths = [p for p in source.rglob('*') if p.is_file()
             and (p.suffix in ('.h', '.cpp', '.py', '.cmake') or p.name == 'CMakeLists.txt')]
    paths.append(binary)
    return {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}


def command(args, case):
    extra = ['--fine'] if case == 'fine-serial' else []
    base = [str(args.binary)] + extra
    if case != 'coarse-mpi2':
        return base
    return ([args.mpiexec, args.numproc_flag, '2'] + args.mpi_before
            + base[:1] + args.mpi_after + base[1:])


def stop_group(process):
    os.killpg(process.pid, signal.SIGTERM)
    try:
        process.wait(timeout=2)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.wait()


def invoke(argv, cwd, environment, out, err, stdin=None):
    with subprocess.Popen(argv, cwd=cwd, env=environment, stdout=out, stderr=err, stdin=stdin,
                          start_new_session=True) as process:
        try:
            return process.wait(timeout=120)
        except subprocess.TimeoutExpired:
            stop_group(process)
            return 124


def run_case(args, directory, case):
    csv = directory / (case + '.csv')
    stderr = directory / (case + '.stderr.log')
    argv = command(args, case)
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS='1', OMP_PROC_BIND='false')
    with csv.open('w') as out, stderr.open('w') as err:
        status = invoke(argv, directory, environment, out, err)
    return {'case': case, 'argv': argv, 'cwd': str(directory), 'exit_code': status,
            'environment': {'OMP_NUM_THREADS': '1', 'OMP_PROC_BIND': 'false'},
            'csv': str(csv), 'stderr': str(stderr)}


def save_record(args, directory, record):
    record['source_after'] = fingerprint(args.binary)
    record['source_unchanged'] = record['source_before'] == record['source_after']
    text = json.dumps(record, indent=2) + '\n'
    (directory / 'wave-runs.json').write_text(text)
    temporary = args.results / 'wave-runs.json.tmp'
    temporary.write_text(text)
    temporary.replace(args.results / 'wave-runs.json')


def execute(args, directory, record):
    for case in ('coarse-serial', 'fine-serial', 'coarse-mpi2'):
        run = run_case(args, directory, case)
        record['runs'].append(run)
        if run['exit_code']:
            raise RuntimeError(f"{case} exited {run['exit_code']}: {run['stderr']}")
    record['csv_checks'] = compare_csv_runs(record['runs'])


def main():
    if not __debug__:
        raise RuntimeError('Run verification without Python optimization.')
    args = arguments()
    args.binary, args.results = args.binary.resolve(), args.results.resolve()
    args.results.mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix='wave-', dir=args.results))
    record = {'utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'source_before': fingerprint(args.binary), 'runs': [],
              'mpi_comparison': 'global diagnostics only',
              'deferred': ['field-by-field MPI comparison', 'VTK output and real-reader checks']}
    try:
        execute(args, directory, record)
    finally:
        save_record(args, directory, record)
    assert record['source_unchanged'], 'Sources changed during the benchmark.'
    print(json.dumps(record['csv_checks'], indent=2))


if __name__ == '__main__':
    main()
