#!/usr/bin/env python3
"""Run one solver configuration of E11, E12, E13 or E14, one mesh after the other.

The first argument names the experiment; the other arguments are those of the E03 launcher
(run_e03_study.py), which writes the inputs, launches the executable and records the exit
status, wall time and peak memory of every rank. --solver-configuration is required.
For the two staggered configurations, E11 and E13 always save the z=0 sections at the initial
state and the final field target on every level (write_plane), which E11.meas01 compares.
The script writes measured data and run records only.

Usage: run_dipole_study.py E11|E12|E13|E14 --executable FILE --output DIRECTORY
       --solver-configuration O1VS|Homogeneous|Standard|StandardR2Prime
       [run_e03_study.py arguments]
"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent
                       / "E01E02E03E04_wave_pulses_in_homogeneous_dielectric"))
from run_e03_study import prepare_study, run_case, study_arguments  # noqa: E402

EXPERIMENTS = ("E11", "E12", "E13", "E14")


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in EXPERIMENTS:
        raise SystemExit(__doc__)
    experiment = sys.argv.pop(1)
    args = study_arguments(__doc__)
    if not args.solver_configuration:
        raise SystemExit("--solver-configuration is required for " + experiment)
    staggered = args.solver_configuration in ("O1VS", "Homogeneous")
    args.write_plane = args.write_plane or (experiment in ("E11", "E13") and staggered)
    executable = prepare_study(args, experiment)
    success = True
    for n in sorted(set(args.levels)):
        success = run_case(args, executable, n, experiment) and success
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
