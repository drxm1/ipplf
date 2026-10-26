#!/usr/bin/env python3
"""Run the C++ finite-window probe test at both rank counts and compare its output.

The test FourierProbeSpectra.FrozenGatesAndNaNGhosts writes the record. It records a
scale R from the time step, window sum and tone amplitudes, defined by checkProbeValue().
Every complex output is compared using that scale with relative tolerance 1e-14.
Field evolution and the closed-form Fourier reference remain in the C++ test.
"""
import argparse
import math
import os
from pathlib import Path
import struct
import subprocess

RECORD_TEST = "FourierProbeSpectra.FrozenGatesAndNaNGhosts"


def expected_keys():
    """The declared fixture has two planes and six points, eleven channels and six frequencies."""
    return {(probe, cell, channel, frequency)
            for probe in range(8) for cell in range(128 if probe < 2 else 1)
            for channel in range(11) for frequency in range(6)}


def arguments():
    """Retain CMake's MPI launch flags, including flags with leading hyphens."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", required=True)
    parser.add_argument("--mpiexec", required=True)
    parser.add_argument("--numproc-flag", default="-n")
    parser.add_argument("--mpi-preflag", action="append", default=[])
    parser.add_argument("--mpi-postflag", action="append", default=[])
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_record(path, ranks):
    """Read every cell/channel/frequency once; duplicate or missing data is a failure."""
    result = {}
    lines = path.read_text().splitlines()
    if not lines or lines[0] != f"# communicator_size {ranks}":
        raise RuntimeError(f"Missing or incorrect communicator size in {path}")
    for line in lines[1:]:
        probe, cell, channel, frequency, scale, real, imaginary = line.split()
        key = tuple(map(int, (probe, cell, channel, frequency)))
        if key in result:
            raise RuntimeError(f"Duplicate probe datum {key}")
        result[key] = tuple(map(float, (scale, real, imaginary)))
    if result.keys() != expected_keys():
        raise RuntimeError(f"Incomplete or unexpected probe coverage in {path}")
    return result


def test_command(args, ranks):
    test_filter = f"--gtest_filter={RECORD_TEST}"
    if ranks == 1:
        return [args.executable, "--gtest_color=no", test_filter]
    return [args.mpiexec, args.numproc_flag, "2", *args.mpi_preflag,
            args.executable, *args.mpi_postflag, "--gtest_color=no", test_filter]


def run_case(args, ranks):
    """Execute the real accumulator, preserve its log, and bound each run to 60 seconds."""
    path = args.output_dir / f"probe-ranks{ranks}.txt"
    path.unlink(missing_ok=True)
    environment = dict(os.environ, CHDR_PROBE_OUTPUT=str(path.resolve()),
                       OMP_NUM_THREADS="1", OMP_PROC_BIND="false")
    completed = subprocess.run(test_command(args, ranks), env=environment, text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               timeout=60, check=False)
    (args.output_dir / f"probe-ranks{ranks}.log").write_text(completed.stdout)
    print(completed.stdout, end="")
    if completed.returncode:
        raise RuntimeError(f"{RECORD_TEST} failed with {ranks} ranks: exit {completed.returncode}")
    if f"[       OK ] {RECORD_TEST}" not in completed.stdout:
        raise RuntimeError(f"{RECORD_TEST} did not pass with {ranks} ranks")
    return read_record(path, ranks)


def compare_value(key, serial, parallel):
    """Use the C++ record's R, including when the transform itself is nearly zero."""
    scale, real, imaginary = serial
    mpi_scale, mpi_real, mpi_imaginary = parallel
    if not math.isfinite(scale) or not math.isfinite(mpi_scale) or scale != mpi_scale or not scale > 0:
        raise RuntimeError(f"Invalid or inconsistent normalization at {key}")
    if not all(map(math.isfinite, (real, imaginary, mpi_real, mpi_imaginary))):
        raise RuntimeError(f"Nonfinite probe output at {key}")
    relative = math.hypot(real - mpi_real, imaginary - mpi_imaginary) / scale
    if relative > 1e-14:
        raise RuntimeError(f"Rank mismatch {relative}/R at {key}")
    identical = struct.pack("dd", real, imaginary) == struct.pack("dd", mpi_real, mpi_imaginary)
    return relative, identical


def compare_records(serial, parallel):
    if serial.keys() != parallel.keys():
        raise RuntimeError("Serial and MPI probe outputs contain different cells/channels/frequencies")
    maximum, bitwise = 0.0, True
    for key, value in serial.items():
        relative, identical = compare_value(key, value, parallel[key])
        maximum = max(maximum, relative)
        bitwise &= identical
    print(f"Probe spectra rank comparison: {len(serial)} complex values, "
          f"max difference/R={maximum:.17g}, "
          f"bitwise_equal={bitwise}")


def main():
    args = arguments()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    serial = run_case(args, 1)
    parallel = run_case(args, 2)
    compare_records(serial, parallel)


if __name__ == "__main__":
    main()
