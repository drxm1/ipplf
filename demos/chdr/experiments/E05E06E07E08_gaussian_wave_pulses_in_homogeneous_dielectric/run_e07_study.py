#!/usr/bin/env python3
"""Run the E07 study (E03 with the Gaussian profile), one mesh after the other.

The script uses the launcher of E03 (run_e03_study.py) with `case: E07`. Do not pass
--profile-sharpness: E07 rejects that parameter. --plane-every writes z=0 sections of E and
R2 B for two-dimensional videos. The script writes measured data and run records only.
"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent
                       / "E01E02E03E04_wave_pulses_in_homogeneous_dielectric"))
from run_e03_study import prepare_study, run_case, study_arguments  # noqa: E402


def main():
    args = study_arguments(__doc__)
    executable = prepare_study(args, "E07")
    success = True
    for n in sorted(set(args.levels)):
        success = run_case(args, executable, n, "E07") and success
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
