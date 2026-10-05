"""Run one writer/real-reader case in a fresh directory; retain commands and evidence."""
import argparse
import hashlib
import json
import os
import subprocess
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent


def arguments():
    """The command after -- is an executable or the configured MPI launcher plus executable."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reader", required=True)
    parser.add_argument("--output-parent", type=Path, required=True)
    parser.add_argument("--ranks", type=int, choices=(1, 2), required=True)
    parser.add_argument("--case", choices=("valid", "piece", "wrapper", "series", "nonfinite"), required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    options = parser.parse_args()
    if options.command[:1] == ["--"]:
        options.command = options.command[1:]
    if not options.command:
        parser.error("a fixture command is required after --")
    return options


def test_environment(directory):
    """Confine settings/cache writes and keep the CPU/thread footprint bounded."""
    environment = dict(os.environ, OMP_NUM_THREADS="1", OMP_PROC_BIND="false",
                       OPENBLAS_NUM_THREADS="1", VTK_SMP_MAX_THREADS="1")
    for key, name in (("XDG_CONFIG_HOME", "config"), ("XDG_CACHE_HOME", "cache")):
        path = directory / name
        path.mkdir()
        environment[key] = str(path)
    return environment


def execute(command, directory, environment, label):
    """Save unabridged child output and a truthful exit status for one bounded command."""
    result = subprocess.run(command, cwd=directory, env=environment, capture_output=True,
                            text=True, timeout=90)
    log = directory / f"{label}.log"
    log.write_text(result.stdout + result.stderr)
    print(log.read_text(), end="", flush=True)
    return {"command": command, "exit_code": result.returncode, "log": str(log)}


def source_hashes():
    """Identify the local fixture, reader, orchestration and writer source revision."""
    files = (HERE / "output_fixture.cpp", HERE / "read_output_fixture.py",
             HERE / "run_output_fixture.py", HERE.parent / "VacuumOutput.h")
    return {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in files}


def run_case(options, directory):
    """Run the real writer first, then require VTK/ParaView readback to succeed."""
    environment = test_environment(directory)
    data = directory / "data"
    command = options.command + [str(data)]
    if options.case != "valid":
        command.append(options.case)
    results = [execute(command, directory, environment, "writer")]
    if results[-1]["exit_code"] == 0:
        mode = str(options.ranks) if options.case == "valid" else options.case
        reader = [options.reader, "--disable-registry", str(HERE / "read_output_fixture.py")]
        results.append(execute(reader + [str(data), mode], directory, environment, "reader"))
    record = {"case": options.case, "ranks": options.ranks, "sources": source_hashes(),
              "commands": results, "omp_threads": 1, "gui_walkthrough": False}
    (directory / "record.json").write_text(json.dumps(record, indent=2) + "\n")
    print(f"Evidence: {directory / 'record.json'}", flush=True)
    return 0 if len(results) == 2 and all(item["exit_code"] == 0 for item in results) else 1


def main():
    """Missing pvpython is an explicit CTest skip, never successful reader acceptance."""
    options = arguments()
    if not options.reader or not Path(options.reader).is_file():
        print("SKIP: pvpython unavailable; actual output-reader acceptance is unfinished")
        return 77
    options.output_parent.mkdir(parents=True, exist_ok=True)
    prefix = f"{options.case}-{options.ranks}rank-"
    directory = Path(tempfile.mkdtemp(prefix=prefix, dir=options.output_parent)).resolve()
    return run_case(options, directory)


if __name__ == "__main__":
    raise SystemExit(main())
