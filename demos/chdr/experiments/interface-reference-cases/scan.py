#!/usr/bin/env python3
"""Run explicit plane YAML inputs with bounded local resources and retained failures.

Only output/directory is changed. Physical input text is preserved byte for byte;
the application remains responsible for resolving numerical settings. PyYAML is
used for node locations, not to deserialize and rewrite the physical parameters.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import threading
import time
import uuid

from read_run import load_run, sha256_file, source_identity


@dataclass(frozen=True)
class ScanSettings:
    binary_m: Path
    output_root_m: Path
    ranks_m: int = 1
    threads_m: int = 1
    jobs_m: int = 1
    max_cores_m: int | None = None
    mpi_m: str = "mpiexec"
    mpi_arguments_m: tuple[str, ...] = ()
    timeout_m: float | None = None


class ScanCancelled(Exception):
    """An attempt was interrupted or prevented from starting by scan cancellation."""


class AttemptProcesses:
    """Serialize launch registration against cancellation, including queued attempts."""
    def __init__(self):
        self.lock_m = threading.Lock()
        self.cancelled_m = False
        self.processes_m = set()
        self.interrupted_m = set()

    def launch(self, command: list[str], **options) -> subprocess.Popen:
        with self.lock_m:
            if self.cancelled_m:
                raise ScanCancelled("Scan cancelled before this attempt started")
            process = subprocess.Popen(command, start_new_session=True, **options)
            self.processes_m.add(process)
            return process

    def interrupted(self, process: subprocess.Popen) -> bool:
        with self.lock_m:
            return process in self.interrupted_m

    def remove(self, process: subprocess.Popen) -> None:
        with self.lock_m:
            self.processes_m.discard(process)

    def cancel(self) -> list[dict]:
        with self.lock_m:
            self.cancelled_m = True
            processes = tuple(self.processes_m)
            self.interrupted_m.update(processes)
        return _stop_groups(processes)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def _members(node, label: str) -> dict:
    if node.id != "mapping":
        raise ValueError(label + " must be a YAML mapping")
    result = {}
    for key, value in node.value:
        if key.id != "scalar" or key.value in result:
            raise ValueError(label + " has a nonscalar or duplicate key")
        result[key.value] = value
    return result


def output_input(text: str, output: Path) -> str:
    """Replace one scalar by source positions, preserving all physical YAML spelling."""
    import yaml
    if any(getattr(event, "anchor", None) is not None for event in yaml.parse(text)):
        raise ValueError("Scan inputs must be explicit YAML without anchors or aliases")
    root = yaml.compose(text)
    if root is None:
        raise ValueError("Empty YAML input")
    members = _members(root, "input")
    if "case" not in members or members["case"].value != "plane":
        raise ValueError("This runner requires an explicit case: plane")
    if "output" in members:
        return _replace_directory(text, _members(members["output"], "output"), output)
    return _append_output(text, root, output)


def _replace_directory(text: str, output: dict, directory: Path) -> str:
    node = output.get("directory")
    if node is None or node.id != "scalar":
        raise ValueError("An existing output mapping must contain a scalar directory")
    start, end = node.start_mark.index, node.end_mark.index
    replacement = json.dumps(str(directory))
    if text[start:end].endswith("\n"):
        replacement += "\n"
    return text[:start] + replacement + text[end:]


def _append_output(text: str, root, directory: Path) -> str:
    value = json.dumps({"directory": str(directory)})
    if root.flow_style:
        import yaml
        index = root.end_mark.index - 1
        end = root.value[-1][1].end_mark.index
        trailing_comma = any(isinstance(token, yaml.tokens.FlowEntryToken)
                             and end <= token.start_mark.index < index for token in yaml.scan(text))
        separator = " " if trailing_comma else ", "
        return text[:index] + separator + "output: " + value + text[index:]
    index = root.end_mark.index
    newline = "" if index == 0 or text[index - 1] == "\n" else "\n"
    return text[:index] + newline + "output: " + value + "\n" + text[index:]


def _available_cores() -> int:
    return len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)


def _validate(settings: ScanSettings) -> None:
    if not hasattr(os, "killpg"):
        raise ValueError("The scan runner requires POSIX process-group cleanup")
    counts = (settings.ranks_m, settings.threads_m, settings.jobs_m)
    if any(value < 1 for value in counts):
        raise ValueError("ranks, threads and jobs must be positive")
    budget = _available_cores() if settings.max_cores_m is None else settings.max_cores_m
    if budget < 1 or budget > _available_cores() or settings.jobs_m * settings.ranks_m * settings.threads_m > budget:
        raise ValueError("jobs * ranks * threads exceeds the allowed local CPU budget")
    if settings.timeout_m is not None and (not math.isfinite(settings.timeout_m) or settings.timeout_m <= 0):
        raise ValueError("timeout must be positive")
    if not settings.binary_m.is_file() or not os.access(settings.binary_m, os.X_OK):
        raise ValueError("--binary must name an executable file")
    if settings.ranks_m > 1 and shutil.which(settings.mpi_m) is None:
        raise ValueError("MPI launcher not found: " + settings.mpi_m)
    if settings.ranks_m == 1 and settings.mpi_arguments_m:
        raise ValueError("MPI arguments require ranks > 1")


def _command(settings: ScanSettings, input_path: Path) -> list[str]:
    command = [str(settings.binary_m.resolve()), "--input", str(input_path)]
    if settings.ranks_m > 1:
        command = [shutil.which(settings.mpi_m), *settings.mpi_arguments_m, "-n", str(settings.ranks_m)] + command
    return command


def _attempt(root: Path, source: Path, index: int, scan_id: str) -> tuple[Path, dict]:
    stem = re.sub(r"[^A-Za-z0-9_-]+", "_", source.stem).strip("_") or "row"
    row_id = f"{index:04d}-{stem}"
    directory = (root / row_id / scan_id).resolve()
    if not directory.is_relative_to(root):
        raise ValueError("Attempt directory escapes the selected output root")
    directory.mkdir(parents=True, exist_ok=False)
    record = {"attempt_record_version": 1, "row_id": row_id, "attempt_id": scan_id,
              "source_input": str(source.resolve()), "directory": str(directory),
              "status": "preparing", "physics_validation": "not_performed"}
    _atomic_json(directory / "attempt.json", record)
    return directory, record


def _prepare(source: Path, directory: Path, record: dict, settings: ScanSettings) -> None:
    original = source.read_bytes()
    (directory / "input.original.yaml").write_bytes(original)
    record["original_input_sha256"] = hashlib.sha256(original).hexdigest()
    text = output_input(original.decode("utf-8"), directory / "output")
    input_path = directory / "input.yaml"
    input_path.write_text(text)
    record["executed_input_sha256"] = sha256_file(input_path)
    record["command"] = _command(settings, input_path)
    record["working_directory"] = str(Path.cwd())
    record["status"] = "prepared"


def _signal_group(process: subprocess.Popen, action: int) -> None:
    try:
        os.killpg(process.pid, action)
    except ProcessLookupError:
        pass


def _wait_group(process: subprocess.Popen, seconds: float) -> bool:
    deadline = time.monotonic() + seconds
    while time.monotonic() < deadline:
        process.poll()
        try:
            os.killpg(process.pid, 0)
        except ProcessLookupError:
            return True
        time.sleep(.02)
    return False


def _stop_group(process: subprocess.Popen) -> bool:
    """Bound cleanup of the local group, including ranks whose launcher already exited."""
    _signal_group(process, signal.SIGTERM)
    if not _wait_group(process, 2):
        _signal_group(process, signal.SIGKILL)
    process.wait(timeout=2)
    return _wait_group(process, 2)


def _wait_groups(processes: tuple, seconds: float) -> tuple:
    deadline = time.monotonic() + seconds
    remaining = processes
    while remaining and time.monotonic() < deadline:
        remaining = tuple(process for process in remaining if _group_present(process))
        if remaining:
            time.sleep(.02)
    return remaining


def _group_present(process: subprocess.Popen) -> bool:
    process.poll()
    try:
        os.killpg(process.pid, 0)
        return True
    except ProcessLookupError:
        return False


def _stop_groups(processes: tuple) -> list[dict]:
    """Signal every active group before waiting, so cleanup does not scale with job count."""
    for process in processes:
        _signal_group(process, signal.SIGTERM)
    for process in _wait_groups(processes, 2):
        _signal_group(process, signal.SIGKILL)
    remaining = _wait_groups(processes, 2)
    result = []
    for process in processes:
        try:
            process.wait(timeout=0)
            result.append({"pid": process.pid, "returncode": process.returncode,
                           "process_group_still_present": process in remaining})
        except subprocess.TimeoutExpired:
            result.append({"pid": process.pid, "cleanup_error": "launcher_did_not_exit_after_SIGKILL"})
    return result


def _command_returncode(process: subprocess.Popen, timeout: float | None,
                        active: AttemptProcesses) -> int:
    deadline = None if timeout is None else time.monotonic() + timeout
    while True:
        if active.interrupted(process):
            raise ScanCancelled("Scan cancelled while this attempt was running")
        remaining = .1 if deadline is None else max(0, min(.1, deadline - time.monotonic()))
        try:
            result = process.wait(timeout=remaining)
            if active.interrupted(process):
                raise ScanCancelled("Scan cancelled while this attempt was running")
            return result
        except subprocess.TimeoutExpired:
            if deadline is not None and time.monotonic() >= deadline:
                raise subprocess.TimeoutExpired(process.args, timeout)


def _wait_process(process: subprocess.Popen, record: dict, settings: ScanSettings,
                  active: AttemptProcesses) -> int:
    try:
        record["returncode"] = _command_returncode(process, settings.timeout_m, active)
        return record["returncode"]
    except subprocess.TimeoutExpired:
        if active.interrupted(process):
            raise ScanCancelled("Scan cancelled while this attempt was running")
        disappeared = _stop_group(process)
        record.update(returncode=process.returncode, timeout_cleanup="group_signalled_launcher_reaped",
                      process_group_still_present=not disappeared)
        raise
    except ScanCancelled:
        raise
    except BaseException:
        _stop_group(process)
        raise


def _wait_command(record: dict, settings: ScanSettings, output, errors, environment: dict,
                  active: AttemptProcesses) -> int:
    process = active.launch(record["command"], stdout=output, stderr=errors, env=environment,
                            cwd=record["working_directory"])
    try:
        return _wait_process(process, record, settings, active)
    except ScanCancelled:
        record["returncode"] = process.poll()
        raise
    finally:
        active.remove(process)


def _execute(directory: Path, record: dict, settings: ScanSettings, binary_hash: str,
             active: AttemptProcesses) -> None:
    if sha256_file(settings.binary_m) != binary_hash:
        raise ValueError("Executable changed before this row; use an immutable binary copy")
    environment = dict(os.environ, OMP_NUM_THREADS=str(settings.threads_m), OMP_PROC_BIND="false")
    record.update(status="running", started_utc=_now(), binary_sha256=binary_hash,
                  ranks=settings.ranks_m, threads=settings.threads_m, timeout_seconds=settings.timeout_m)
    _atomic_json(directory / "attempt.json", record)
    with (directory / "stdout.txt").open("x") as output, (directory / "stderr.txt").open("x") as errors:
        record["returncode"] = _wait_command(record, settings, output, errors, environment, active)
    if record["returncode"]:
        raise ValueError(f"Simulation exited with code {record['returncode']}")
    if sha256_file(settings.binary_m) != binary_hash:
        raise ValueError("Executable changed during this row; result identity is ambiguous")


def _completed(directory: Path, record: dict) -> None:
    run = load_run(directory / "output")
    record.update(status="completed", output_validation="record_and_spectrum_schema",
                  output_directory=str(run.directory_m), build=run.record_m.get("build", {}),
                  valid_bins=sum(sample.valid_m for sample in run.samples_m),
                  bins=len(run.samples_m))
    record["output_sha256"] = {name: sha256_file(run.directory_m / name)
                               for name in ("run.json", "spectra.csv")}


def _run_attempt(directory: Path, record: dict, settings: ScanSettings, binary_hash: str,
                 active: AttemptProcesses | None = None) -> dict:
    started = time.monotonic()
    try:
        _execute(directory, record, settings, binary_hash, active or AttemptProcesses())
        _completed(directory, record)
    except ScanCancelled as error:
        record.update(status="cancelled", error_type=type(error).__name__, error=str(error))
    except Exception as error:
        record.update(status="failed", error_type=type(error).__name__, error=str(error))
    finally:
        record.update(finished_utc=_now(), elapsed_seconds=time.monotonic() - started)
        _atomic_json(directory / "attempt.json", record)
    return record


def _prepare_all(inputs: list[Path], root: Path, scan_id: str, settings: ScanSettings) -> list:
    prepared = []
    for index, source in enumerate(inputs):
        directory, record = _attempt(root, source, index, scan_id)
        try:
            _prepare(source, directory, record, settings)
        except Exception as error:
            record.update(status="failed", error_type=type(error).__name__, error=str(error))
        _atomic_json(directory / "attempt.json", record)
        prepared.append((directory, record))
    return prepared


def run_scan(inputs: list[Path], settings: ScanSettings) -> Path:
    """Allocate new attempts, preserve inputs first, then run within the explicit CPU budget."""
    _validate(settings)
    if not inputs:
        raise ValueError("At least one explicit input is required")
    root = settings.output_root_m.resolve()
    root.mkdir(parents=True, exist_ok=True)
    scan_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ-") + uuid.uuid4().hex[:8]
    manifest = root / ("scan-" + scan_id + ".json")
    prepared = _prepare_all(inputs, root, scan_id, settings)
    record = {"scan_record_version": 1, "scan_id": scan_id, "status": "running",
              "binary": str(settings.binary_m.resolve()), "binary_sha256": sha256_file(settings.binary_m),
              "source_at_launch": source_identity(), "started_utc": _now(),
              "ranks": settings.ranks_m, "threads": settings.threads_m, "jobs": settings.jobs_m,
              "max_cores": _available_cores() if settings.max_cores_m is None else settings.max_cores_m,
              "timeout_seconds": settings.timeout_m,
              "physics_validation": "not_performed", "attempts": [item[1] for item in prepared]}
    _atomic_json(manifest, record)
    _finish_scan(prepared, settings, record, manifest)
    return manifest


def _pool_results(pool, prepared: list, settings: ScanSettings, record: dict,
                  manifest: Path, active: AttemptProcesses) -> None:
    futures = {pool.submit(_run_attempt, directory, dict(item), settings, record["binary_sha256"], active): index
               for index, (directory, item) in enumerate(prepared) if item["status"] == "prepared"}
    for future in as_completed(futures):
        record["attempts"][futures[future]] = future.result()
        _atomic_json(manifest, record)


def _execute_pool(prepared: list, settings: ScanSettings, record: dict, manifest: Path) -> bool:
    active = AttemptProcesses()
    with ThreadPoolExecutor(max_workers=settings.jobs_m) as pool:
        try:
            _pool_results(pool, prepared, settings, record, manifest, active)
        except KeyboardInterrupt:
            signal.signal(signal.SIGINT, signal.SIG_IGN)
            record.update(status="cancelling", cancellation="SIGINT")
            record["cancellation_cleanup"] = active.cancel()
            _atomic_json(manifest, record)
        except BaseException:
            active.cancel()
            raise
    return active.cancelled_m


def _refresh_cancelled(prepared: list, record: dict) -> None:
    for index, (directory, _) in enumerate(prepared):
        item = json.loads((directory / "attempt.json").read_text())
        if item["status"] in ("preparing", "prepared"):
            item.update(status="cancelled", error_type="ScanCancelled",
                        error="Scan cancelled before this attempt was submitted", finished_utc=_now())
            _atomic_json(directory / "attempt.json", item)
        record["attempts"][index] = item


def _finish_scan(prepared: list, settings: ScanSettings, record: dict, manifest: Path) -> None:
    previous = signal.getsignal(signal.SIGINT)
    try:
        cancelled = _execute_pool(prepared, settings, record, manifest)
        if cancelled:
            _refresh_cancelled(prepared, record)
        completed = all(item["status"] == "completed" for item in record["attempts"])
        record["status"] = "cancelled" if cancelled else ("completed" if completed else "failed")
        record["finished_utc"] = _now()
        _atomic_json(manifest, record)
    finally:
        if threading.current_thread() is threading.main_thread():
            signal.signal(signal.SIGINT, previous)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    for name in ("ranks", "threads", "jobs"):
        parser.add_argument("--" + name, type=int, default=1)
    parser.add_argument("--max-cores", type=int)
    parser.add_argument("--mpi", default="mpiexec")
    parser.add_argument("--mpi-arg", action="append", default=[])
    parser.add_argument("--timeout", type=float)
    return parser


def main() -> int:
    parser = _parser()
    args = parser.parse_args()
    settings = ScanSettings(args.binary.resolve(), args.output_root, args.ranks, args.threads,
                            args.jobs, args.max_cores, args.mpi, tuple(args.mpi_arg), args.timeout)
    try:
        manifest = run_scan(args.inputs, settings)
        print(manifest)
        status = json.loads(manifest.read_text())["status"]
        return 130 if status == "cancelled" else (0 if status == "completed" else 1)
    except (OSError, ValueError) as error:
        parser.exit(2, f"scan: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
