"""Run a real pvpython readback or report its absent dependency as CTest skip 77."""
import argparse
from pathlib import Path
import subprocess


def main():
    """Preserve reader failure status; an unavailable reader is never a passing check."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reader")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if not args.reader or not Path(args.reader).is_file():
        print("SKIP: pvpython unavailable; actual wave-reader acceptance is unfinished")
        return 77
    return subprocess.call([args.reader, "--disable-registry"] + args.command)


if __name__ == "__main__":
    raise SystemExit(main())
