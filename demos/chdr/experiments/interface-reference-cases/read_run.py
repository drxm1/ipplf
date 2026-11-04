#!/usr/bin/env python3
"""Read recorded plane measurements without reconstructing requested setup choices.

Effective quantities are authoritative normalized values. Schema zero is the
explicit first aligned-row record (no raw spectra); schema one adds raw spectra.
Schema two records tilted plane-wave fits and lifted probe spectra separately.
See references.org for physical definitions.
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
    face_amplitudes_m: tuple[complex, complex, complex] | None = None
    normal_wavenumbers_m: tuple[float, float, float] | None = None
    normal_group_velocities_m: tuple[float, float, float] | None = None
    fit_residuals_m: tuple[float, float, float] | None = None


@dataclass(frozen=True)
class ProbeSample:
    """One actual component sample, its periodic lift and raw temporal transforms."""

    bin_m: int
    side_m: str
    point_m: int
    position_m: tuple[float, float, float]
    winding_m: tuple[int, int]
    normal_position_m: float
    tangent_position_m: float
    incident_m: complex
    total_m: complex


@dataclass(frozen=True)
class PlaneRun:
    """An accepted record with spectra in their original bin order."""

    directory_m: Path
    record_m: dict
    samples_m: tuple[SpectrumSample, ...]
    schema_m: int
    central_bin_m: int
    probes_m: tuple[ProbeSample, ...] = ()

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
    if schema not in (0, 1, 2):
        raise RunDataError(f"Unsupported plane schema version: {schema}")
    _equal(plane.get("geometry"), "periodic_tilted_plane" if schema == 2 else "aligned_plane", "plane/geometry")
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
    positive = ("epsilon", "omega0", "sigma_k", "tail_sigmas",
                "requested_cells_per_wavelength", "cells_per_wavelength", "separation_factor")
    other = ("requested_angle_degrees", "actual_angle_degrees", "cutoff_fraction")
    if plane.get("schema_version") != 2:
        positive += ("kx0",)
        other += ("kz", "face_position", "incident_probe_x", "other_probe_x", "probe_z")
    for name in positive + other:
        _finite(_at(plane, name), "plane/" + name, name in positive)
    for name in ("gate_incident", "gate_other"):
        _integer(_at(plane, name), "plane/" + name, 1)
    _integer(_at(plane, "timing_samples"), "plane/timing_samples")
    solves = _integer(_at(record, "effective", "completed_solves"), "completed_solves")
    if solves < max(plane["gate_incident"], plane["gate_other"]) - 1:
        raise RunDataError("Completed solves do not cover the recorded probe gates")


def _plane_vector(plane: dict, name: str, positive: bool = False) -> tuple:
    values = _at(plane, name)
    if not isinstance(values, list) or len(values) != 3:
        raise RunDataError("plane/" + name + " must contain three values")
    return tuple(_finite(value, "plane/" + name, positive) for value in values)


def _tilted_basis(plane: dict) -> None:
    normal, tangent = (_plane_vector(plane, name) for name in ("normal", "tangent"))
    _check_face_anchor(plane, normal)
    _plane_vector(plane, "requested_cell_shape", True)
    angle = _finite(_at(plane, "interface_tilt_degrees"), "plane/interface_tilt_degrees")
    expected = (math.cos(math.radians(angle)), 0, math.sin(math.radians(angle)))
    if any(not math.isclose(a, b, rel_tol=2e-12, abs_tol=2e-12) for a, b in zip(normal, expected)):
        raise RunDataError("Tilt angle and recorded normal disagree")
    if any(not math.isclose(a, b, rel_tol=2e-12, abs_tol=2e-12)
           for a, b in zip(tangent, (-normal[2], 0, normal[0]))):
        raise RunDataError("Tilted tangent and normal basis disagree")
    _equal(_at(plane, "projection"), "unweighted_complex_plane_wave_fit", "plane/projection")
    _equal(_at(plane, "projection_basis"), "exp(i*(normal_k*u+kappa*v))", "plane/projection_basis")
    _equal(_at(plane, "projection_normalization"), "mean_over_points", "plane/projection_normalization")


def _check_face_anchor(plane: dict, normal: tuple) -> None:
    """The supplied phase anchor must lie on the requested signed normal-offset plane."""
    anchor = _plane_vector(plane, "anchor")
    offset = _finite(_at(plane, "requested_face_position"), "plane/requested_face_position")
    distance = sum(a * b for a, b in zip(normal, anchor))
    tolerance = 2e-12 * max(1, *(abs(value) for value in anchor))
    if not math.isclose(distance, offset, rel_tol=2e-12, abs_tol=tolerance):
        raise RunDataError("Tilted phase anchor does not lie on the requested face")


def _tilted_numbers(plane: dict) -> None:
    _tilted_basis(plane)
    positive = ("normal_period", "tangent_period", "primitive_tangent_period", "normal_dielectric_width",
                "normal_incident_width", "mode_cutoff_sigmas", "fit_residual_limit", "incident_floor")
    for name in positive + ("kappa", "floquet_margin", "retained_floquet_margin", "band_sigmas"):
        _finite(_at(plane, name), "plane/" + name, name in positive)
    for name in ("primitive_cells_x", "primitive_cells_z", "period_cells_multiplier", "retained_mode_count"):
        _integer(_at(plane, name), "plane/" + name, 1)
    mode = _at(plane, "tangent_mode")
    if isinstance(mode, bool) or not isinstance(mode, int):
        raise RunDataError("plane/tangent_mode must be an integer")
    if not math.isclose(plane["normal_dielectric_width"] + plane["normal_incident_width"],
                         plane["normal_period"], rel_tol=2e-12):
        raise RunDataError("Tilted material widths do not fill the normal period")
    _finite(_at(plane, "probe_cells"), "plane/probe_cells", True)
    _integer(_at(plane, "probe_half_width_cells"), "plane/probe_half_width_cells")
    for side in ("incident", "transmitted"):
        _integer(_at(plane, "probe_point_count_" + side), "plane/probe_point_count_" + side, 1)


def _raw_conventions(output: dict, schema: int) -> None:
    if schema in (1, 2):
        required = {"raw_units": "normalized", "transform": "dt sum w f exp(+i omega t)",
                    "raw_phase_origin": "probe_coordinates",
                    "coefficient_phase_origin": "requested_face"}
        if schema == 2:
            required.update(raw_phase_origin="lifted_probe_coordinates", face_amplitudes_units="normalized",
                            coefficient_phase_origin="requested_face_anchor")
        for key, value in required.items():
            _equal(output.get(key), value, "spectra/" + key)


def _spectra_conventions(record: dict, schema: int) -> float:
    factors = _object(record, "output_factors")
    if factors.get("system") not in ("normalized", "si"):
        raise RunDataError("Unsupported output unit system")
    factor = _finite(_at(factors, "frequency"), "output_factors/frequency", True)
    output = _object(record, "effective", "output", "spectra")
    _equal(output.get("amplitudes_units"), "dimensionless", "spectra/amplitudes_units")
    _equal(output.get("angular_frequency_units"),
           "rad/s" if factors["system"] == "si" else "normalized", "spectra/frequency units")
    _raw_conventions(output, schema)
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
    if schema == 2:
        _tilted_numbers(plane)
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


def _tilted_sample(row: dict, omega: float, reflection, transmission, valid: bool, reason: str) -> SpectrumSample:
    faces = tuple(_raw_transform(row, name + "_face", valid) for name in ("incident", "reflected", "transmitted"))
    roots = tuple(_csv_number(row, name + "_normal_k", valid) for name in ("incident", "reflected", "transmitted"))
    speeds = tuple(_csv_number(row, name + "_normal_group_velocity", valid)
                   for name in ("incident", "reflected", "transmitted"))
    residuals = tuple(_csv_number(row, "fit_residual_" + name, valid) for name in ("incident", "r", "t"))
    return SpectrumSample(omega, reflection, transmission, valid, reason,
                          face_amplitudes_m=faces, normal_wavenumbers_m=roots,
                          normal_group_velocities_m=speeds, fit_residuals_m=residuals)


def _sample(row: dict, schema: int, frequency: str, factor: float) -> SpectrumSample:
    scaled = _csv_number(row, frequency)
    omega = _csv_number(row, "omega_code") if schema else scaled / factor
    if omega <= 0 or not math.isclose(scaled, omega * factor, rel_tol=2e-12, abs_tol=0):
        raise RunDataError("Frequency column disagrees with normalized omega/output factor")
    valid, reason = _validity(row, schema)
    reflection, transmission = (_complex(row, key, valid) for key in ("r", "t"))
    if schema == 2:
        return _tilted_sample(row, omega, reflection, transmission, valid, reason)
    raw = tuple(_raw_transform(row, key, valid) for key in ("incident", "total", "transmitted")) if schema else (None,) * 3
    return SpectrumSample(omega, reflection, transmission, valid, reason, *raw)


def _spectrum_columns(schema: int, frequency: str) -> set:
    required = {frequency, "r_real", "r_imag", "t_real", "t_imag"}
    if schema:
        required |= {"omega_code", "valid", "reason"}
    if schema == 1:
        required |= {key + suffix for key in ("incident", "total", "transmitted") for suffix in ("_real", "_imag")}
    if schema == 2:
        for name in ("incident", "reflected", "transmitted"):
            required |= {name + suffix for suffix in ("_face_real", "_face_imag", "_normal_k", "_normal_group_velocity")}
        required |= {"fit_residual_incident", "fit_residual_r", "fit_residual_t"}
    return required


def _read_csv(path: Path, required: set) -> list:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        names = reader.fieldnames or []
        if len(names) != len(set(names)) or not required.issubset(names):
            raise RunDataError("Missing or duplicate spectrum columns")
        rows = list(reader)
    if not rows or any(None in row or any(value is None for value in row.values()) for row in rows):
        raise RunDataError("Empty or malformed spectrum rows")
    return rows


def _read_spectra(path: Path, schema: int, system: str, factor: float) -> tuple:
    frequency = "angular_frequency_rad_per_s" if system == "si" else "angular_frequency"
    rows = _read_csv(path, _spectrum_columns(schema, frequency))
    samples = tuple(_sample(row, schema, frequency, factor) for row in rows)
    if len({sample.omega_m for sample in samples}) != len(samples):
        raise RunDataError("Duplicate spectrum frequencies")
    return samples


def _csv_integer(row: dict, name: str) -> int:
    text = row[name]
    try:
        result = int(text)
    except ValueError as error:
        raise RunDataError("Invalid integer column: " + name) from error
    return result


def _probe_sample(row: dict, samples: tuple) -> ProbeSample:
    index, point = (_csv_integer(row, name) for name in ("bin", "point"))
    side = row["side"]
    if not 0 <= index < len(samples) or point < 0 or side not in ("incident", "transmitted"):
        raise RunDataError("Invalid tilted probe bin, point or side")
    if not math.isclose(_csv_number(row, "omega_code"), samples[index].omega_m, rel_tol=2e-12):
        raise RunDataError("Probe frequency does not match its spectrum bin")
    position = tuple(_csv_number(row, name) for name in ("x", "y", "z"))
    winding = tuple(_csv_integer(row, name) for name in ("lift_x", "lift_z"))
    valid = samples[index].valid_m
    return ProbeSample(index, side, point, position, winding, _csv_number(row, "u"), _csv_number(row, "v"),
                       _raw_transform(row, "incident", valid), _raw_transform(row, "total", valid))


def _check_probe_basis(probe: ProbeSample, record: dict) -> None:
    plane, box = record["effective"]["plane"], record["effective"]["box"]
    winding = (probe.winding_m[0], 0, probe.winding_m[1])
    lifted = tuple(x + shift * count * h - origin for x, shift, count, h, origin
                   in zip(probe.position_m, winding, box["cells"], box["spacing"], plane["anchor"]))
    for basis, value in ((plane["normal"], probe.normal_position_m), (plane["tangent"], probe.tangent_position_m)):
        expected = sum(a * b for a, b in zip(basis, lifted))
        tolerance = 2e-12 * max(1, *(abs(value) for value in lifted))
        if not math.isclose(expected, value, rel_tol=2e-12, abs_tol=tolerance):
            raise RunDataError("Tilted probe coordinates disagree with their lift and phase basis")


def _read_probes(directory: Path, record: dict, samples: tuple) -> tuple:
    required = {"bin", "omega_code", "side", "point", "x", "y", "z", "lift_x", "lift_z", "u", "v",
                "incident_real", "incident_imag", "total_real", "total_imag"}
    probes = tuple(_probe_sample(row, samples) for row in _read_csv(directory / "probe-spectra.csv", required))
    unique, coverage = set(), set()
    for probe in probes:
        key = (probe.bin_m, probe.side_m, probe.point_m)
        if key in unique:
            raise RunDataError("Duplicate tilted probe sample")
        unique.add(key)
        coverage.add((probe.bin_m, probe.side_m))
        _check_probe_basis(probe, record)
    expected = {(index, side) for index in range(len(samples)) for side in ("incident", "transmitted")}
    if coverage != expected:
        raise RunDataError("Every tilted bin needs both incident-side and transmitted-side probes")
    _check_probe_counts(probes, record["effective"]["plane"], len(samples))
    return probes


def _check_probe_counts(probes: tuple, plane: dict, bins: int) -> None:
    groups = {(index, side): set() for index in range(bins) for side in ("incident", "transmitted")}
    for probe in probes:
        groups[probe.bin_m, probe.side_m].add((probe.point_m, probe.position_m, probe.winding_m,
                                            probe.normal_position_m, probe.tangent_position_m))
    for side in ("incident", "transmitted"):
        expected = groups[0, side]
        if len(expected) != plane["probe_point_count_" + side]:
            raise RunDataError("Tilted probe count differs from the recorded geometry")
        if any(groups[index, side] != expected for index in range(bins)):
            raise RunDataError("Tilted probe geometry changes between frequency bins")


def _check_tilted_frequencies(plane: dict, samples: tuple) -> None:
    frequencies = _at(plane, "angular_frequencies")
    if not isinstance(frequencies, list) or len(frequencies) != len(samples):
        raise RunDataError("Tilted frequency metadata must cover every recorded bin")
    if any(not math.isclose(_finite(value, "plane/angular_frequencies", True), sample.omega_m, rel_tol=2e-12)
           for value, sample in zip(frequencies, samples)):
        raise RunDataError("Tilted frequency metadata differs from the spectrum")


def _check_tilted_bins(plane: dict, samples: tuple) -> None:
    _check_tilted_frequencies(plane, samples)
    direction = -1 if plane["incidence"] == "vacuum" else 1
    for sample in samples:
        if not sample.valid_m:
            continue
        if min(plane["floquet_margin"], plane["retained_floquet_margin"]) <= 0:
            raise RunDataError("A valid tilted bin needs positive recorded Floquet exclusion margins")
        if any(value < 0 or value > plane["fit_residual_limit"] for value in sample.fit_residuals_m):
            raise RunDataError("A valid tilted bin exceeds the recorded fit-residual limit")
        if any(value * sign <= 0 for value, sign in zip(sample.normal_group_velocities_m,
                                                       (direction, -direction, direction))):
            raise RunDataError("Tilted root group velocities have inconsistent propagation directions")
        if sample.face_amplitudes_m[0] == 0:
            raise RunDataError("A valid tilted bin has zero fitted incident amplitude")


def _assemble_run(directory: Path, record: dict) -> PlaneRun:
    plane, schema, factor = _metadata(record)
    samples = _read_spectra(directory / "spectra.csv", schema, record["output_factors"]["system"], factor)
    central = _integer(_at(plane, "central_bin"), "plane/central_bin") if schema else 0
    if central >= len(samples) or not math.isclose(samples[central].omega_m, plane["omega0"], rel_tol=2e-12):
        raise RunDataError("Central bin does not agree with plane/omega0")
    probes = ()
    if schema == 2:
        _check_tilted_bins(plane, samples)
        probes = _read_probes(directory, record, samples)
    return PlaneRun(directory, record, samples, schema, central, probes)


def load_run(directory: Path | str) -> PlaneRun:
    """Validate one completed run; never infer its physics from its directory name."""
    directory = Path(directory).resolve()
    try:
        record = json.loads((directory / "run.json").read_text())
        return _assemble_run(directory, record)
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
