#!/usr/bin/env python3
"""Run the E04 spherical pulse exit study, one mesh after the other.

The script uses the launcher of run_e03_study.py with `case: E04`. It writes measured data
and run records only: it integrates no traces and draws no figures.
"""

from run_e03_study import prepare_study, run_case, study_arguments


def main():
    args = study_arguments(__doc__)
    executable = prepare_study(args, "E04")
    success = True
    for n in sorted(set(args.levels)):
        success = run_case(args, executable, n, "E04") and success
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
