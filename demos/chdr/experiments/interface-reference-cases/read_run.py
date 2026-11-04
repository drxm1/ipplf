#!/usr/bin/env python3
"""Read recorded plane measurements without reconstructing requested setup choices.

Effective quantities are authoritative normalized values. Schema zero is the
explicit first aligned-row record (no raw spectra); schema one adds raw spectra
and their conventions. See references.org for physical definitions.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import subprocess


class RunDataError(ValueError):
    """A completed record is missing, unsupported or internally inconsistent."""


@dataclass(frozen=True)
class SpectrumSample:
    """One normalized frequency, measured face ratios and optional raw transforms."""

    omega_m: float
    reflection_m: complex | None
    transmission_m: complex | None
    valid_m: bool
    reason_m: str
    incident_m: complex | None = None
    total_m: complex | None = None
    transmitted_m: complex | None = None


@dataclass(frozen=True)
class PlaneRun:
    """An accepted record with spectra in their original bin order."""

    directory_m: Path
    record_m: dict
    samples_m: tuple[SpectrumSample, ...]
    schema_m: int
    central_bin_m: int

    @property
    def plane(self) -> dict:
        return self.record_m["effective"]["plane"]

    @property
    def spacing(self) -> tuple[float, float, float]:
        return tuple(self.record_m["effective"]["box"]["spacing"])

    @property
    def dt(self) -> float:
        return self.record_m["effective"]["solver"]["time_step"]

    @property
    def frequencies(self) -> tuple[float, ...]:
        return tuple(sample.omega_m for sample in self.samples_m)


def sha256_file(path: Path) -> str:
    """Hash the exact bytes used, without loading a large executable into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_identity(source_root: Path | None = None) -> dict:
    """Record the analysis checkout separately from a binary's compiled identity."""
    root = source_root or Path(__file__).resolve().parents[4]
    try:
        prefix = ["git", "-C", str(root)]
        revision = subprocess.check_output(prefix + ["rev-parse", "HEAD"], text=True,
                                           stderr=subprocess.DEVNULL).strip()
        status = subprocess.check_output(prefix + ["status", "--porcelain"], text=True,
                                         stderr=subprocess.DEVNULL)
        return {"revision": revision, "dirty": bool(status), "status": status}
    except (OSError, subprocess.CalledProcessError):
        return {"revision": "unrecorded", "dirty": "unrecorded"}


def _at(record: dict, *path: str):
    current = record
    for part in path:
        if not isinstance(current, dict) or part not in current:
            raise RunDataError("Missing metadata: " + "/".join(path))
        current = current[part]
    return current


def _finite(value, label: str, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise RunDataError(f"{label} must be numeric")
    result = float(value)
    if not math.isfinite(result) or (positive and result <= 0):
        raise RunDataError(f"{label} must be finite" + (" and positive" if positive else ""))
    return result


def _integer(value, label: str, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise RunDataError(f"{label} must be an integer >= {minimum}")
    return value


def _object(record: dict, *path: str) -> dict:
    value = _at(record, *path)
    if not isinstance(value, dict):
        raise RunDataError("Expected an object at " + "/".join(path))
    return value


def _equal(value, expected, label: str) -> None:
    if value != expected:
        raise RunDataError(f"Unsupported {label}: {value!r}; expected {expected!r}")


def _vector(record: dict, key: str, integer: bool = False) -> tuple:
    values = _at(record, "effective", "box", key)
    if not isinstance(values, list) or len(values) != 3:
        raise RunDataError(f"effective/box/{key} must contain three values")
    check = _integer if integer else _finite
    return tuple(check(value, f"box/{key}", 1 if integer else True) for value in values)


def _plane_metadata(record: dict) -> tuple[dict, int]:
    plane = _object(record, "effective", "plane")
    schema = _integer(plane.get("schema_version", 0), "plane/schema_version")
    if schema not in (0, 1):
        raise RunDataError(f"Unsupported plane schema version: {schema}")
    _equal(plane.get("geometry"), "aligned_plane", "plane/geometry")
    _equal(plane.get("incidence_plane"), "xz", "plane/incidence_plane")
    _equal(plane.get("frequency_units"), "normalized", "plane/frequency_units")
    _equal(plane.get("coordinate_units"), "normalized", "plane/coordinate_units")
    if plane.get("polarization") not in ("TE", "TM"):
        raise RunDataError("plane/polarization must be TE or TM")
    _equal(plane.get("component"), "E_y" if plane["polarization"] == "TE" else "H_y",
           "plane/component")
    if plane.get("incidence") not in ("vacuum", "dielectric"):
        raise RunDataError("plane/incidence must identify vacuum or dielectric incidence")
    return plane, schema


def _plane_numbers(record: dict, plane: dict) -> None:
    positive = ("epsilon", "omega0", "kx0", "sigma_k", "tail_sigmas",
                "requested_cells_per_wavelength", "cells_per_wavelength", "separation_factor")
    other = ("kz", "face_position", "incident_probe_x", "other_probe_x", "probe_z",
             "requested_angle_degrees", "actual_angle_degrees", "cutoff_fraction")
    for name in positive + other:
        _finite(_at(plane, name), "plane/" + name, name in positive)
    for name in ("gate_incident", "gate_other"):
        _integer(_at(plane, name), "plane/" + name, 1)
    _integer(_at(plane, "timing_samples"), "plane/timing_samples")
    solves = _integer(_at(record, "effective", "completed_solves"), "completed_solves")
    if solves < max(plane["gate_incident"], plane["gate_other"]) - 1:
        raise RunDataError("Completed solves do not cover the recorded probe gates")


def _spectra_conventions(record: dict, schema: int) -> float:
    factors = _object(record, "output_factors")
    if factors.get("system") not in ("normalized", "si"):
        raise RunDataError("Unsupported output unit system")
    factor = _finite(_at(factors, "frequency"), "output_factors/frequency", True)
    output = _object(record, "effective", "output", "spectra")
    _equal(output.get("amplitudes_units"), "dimensionless", "spectra/amplitudes_units")
    _equal(output.get("angular_frequency_units"),
           "rad/s" if factors["system"] == "si" else "normalized", "spectra/frequency units")
    if schema == 1:
        required = {"raw_units": "normalized", "transform": "dt sum w f exp(+i omega t)",
                    "raw_phase_origin": "probe_coordinates",
                    "coefficient_phase_origin": "requested_face"}
        for key, value in required.items():
            _equal(output.get(key), value, "spectra/" + key)
    if factors["system"] == "normalized" and factor != 1:
        raise RunDataError("Normalized output must have frequency factor one")
    return factor


def _metadata(record: dict) -> tuple[dict, int, float]:
    _equal(_integer(_at(record, "record_version"), "record_version", 1), 1, "record_version")
    _equal(_at(record, "status"), "completed", "run status")
    _equal(_at(record, "application"), "planar-interface", "application")
    _equal(_at(record, "input", "case"), "plane", "case")
    _equal(_at(record, "effective", "value_system"), "normalized", "effective/value_system")
    _vector(record, "cells", True)
    _vector(record, "spacing")
    _finite(_at(record, "effective", "solver", "time_step"), "solver/time_step", True)
    plane, schema = _plane_metadata(record)
    _plane_numbers(record, plane)
    return plane, schema, _spectra_conventions(record, schema)


def _csv_number(row: dict, key: str, finite: bool = True) -> float:
    try:
        value = float(row[key])
    except (KeyError, TypeError, ValueError) as error:
        raise RunDataError(f"Missing or invalid spectra column {key}") from error
    if finite and not math.isfinite(value):
        raise RunDataError(f"Non-finite spectra column {key}")
    return value


def _complex(row: dict, prefix: str, required_finite: bool = True) -> complex | None:
    if not required_finite and all(row.get(prefix + suffix) == "" for suffix in ("_real", "_imag")):
        return None
    parts = [_csv_number(row, prefix + suffix, required_finite) for suffix in ("_real", "_imag")]
    return complex(*parts) if all(map(math.isfinite, parts)) else None


def _validity(row: dict, schema: int) -> tuple[bool, str]:
    if schema == 0:
        return True, ""
    value, reason = row["valid"].strip().lower(), row["reason"].strip()
    if value not in ("0", "1", "true", "false"):
        raise RunDataError("valid must be 0/1 or true/false")
    valid = value in ("1", "true")
    if (not valid and not reason) or (valid and reason not in ("", "ok")):
        raise RunDataError("Inconsistent valid/reason columns")
    return valid, reason


def _raw_transform(row: dict, prefix: str, valid: bool) -> complex:
    """Preserve non-finite raw values on explicitly invalid bins for diagnosis."""
    return complex(*(_csv_number(row, prefix + suffix, valid) for suffix in ("_real", "_imag")))


def _sample(row: dict, schema: int, frequency: str, factor: float) -> SpectrumSample:
    scaled = _csv_number(row, frequency)
    omega = _csv_number(row, "omega_code") if schema else scaled / factor
    if omega <= 0 or not math.isclose(scaled, omega * factor, rel_tol=2e-12, abs_tol=0):
        raise RunDataError("Frequency column disagrees with normalized omega/output factor")
    valid, reason = _validity(row, schema)
    reflection, transmission = (_complex(row, key, valid) for key in ("r", "t"))
    raw = tuple(_raw_transform(row, key, valid) for key in ("incident", "total", "transmitted")) if schema else (None,) * 3
    return SpectrumSample(omega, reflection, transmission, valid, reason, *raw)


def _read_spectra(path: Path, schema: int, system: str, factor: float) -> tuple:
    frequency = "angular_frequency_rad_per_s" if system == "si" else "angular_frequency"
    required = {frequency, "r_real", "r_imag", "t_real", "t_imag"}
    if schema:
        required |= {"omega_code", "valid", "reason"}
        required |= {key + suffix for key in ("incident", "total", "transmitted") for suffix in ("_real", "_imag")}
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        names = reader.fieldnames or []
        if len(names) != len(set(names)) or not required.issubset(names):
            raise RunDataError("Missing or duplicate spectrum columns")
        rows = list(reader)
    if not rows or any(None in row or any(value is None for value in row.values()) for row in rows):
        raise RunDataError("Empty or malformed spectrum rows")
    samples = tuple(_sample(row, schema, frequency, factor) for row in rows)
    if len({sample.omega_m for sample in samples}) != len(samples):
        raise RunDataError("Duplicate spectrum frequencies")
    return samples


def load_run(directory: Path | str) -> PlaneRun:
    """Validate one completed aligned run; never infer its physics from its directory name."""
    directory = Path(directory).resolve()
    try:
        record = json.loads((directory / "run.json").read_text())
        plane, schema, factor = _metadata(record)
        samples = _read_spectra(directory / "spectra.csv", schema,
                                record["output_factors"]["system"], factor)
        central = _integer(_at(plane, "central_bin"), "plane/central_bin") if schema else 0
        if central >= len(samples) or not math.isclose(samples[central].omega_m,
                                                       plane["omega0"], rel_tol=2e-12):
            raise RunDataError("Central bin does not agree with plane/omega0")
        return PlaneRun(directory, record, samples, schema, central)
    except (OSError, json.JSONDecodeError) as error:
        raise RunDataError(f"Cannot read completed run {directory}: {error}") from error


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    args = parser.parse_args()
    try:
        run = load_run(args.run)
        print(json.dumps({"run": str(run.directory_m), "plane_schema": run.schema_m,
                          "bins": len(run.samples_m), "valid_bins": sum(s.valid_m for s in run.samples_m),
                          "central_bin": run.central_bin_m, "physics_validation": "not_performed"}, indent=2))
        return 0
    except RunDataError as error:
        parser.exit(2, f"read_run: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
