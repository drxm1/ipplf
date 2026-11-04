#!/usr/bin/env python3
"""Call the C++ reference calculator and preserve its exact invocation and results.

No Fresnel or lattice-interface coefficient is implemented in Python. The optional
raw-spectrum audit only repeats the documented propagation-phase conversion.
"""
from __future__ import annotations

import argparse
import cmath
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import subprocess

from read_run import PlaneRun, RunDataError, load_run, sha256_file, source_identity


@dataclass(frozen=True)
class ReferenceTable:
    directory_m: Path
    record_m: dict
    rows_m: tuple[dict, ...]


def coefficient(row: dict, reference: str, quantity: str) -> complex | None:
    """Read a declared complex reference column without changing its basis."""
    real, imag = row[f"{reference}_{quantity}_real"], row[f"{reference}_{quantity}_imag"]
    return None if real is None or imag is None else complex(real, imag)


def _lattice_applicability(run: PlaneRun, continuum_only: bool) -> None:
    if run.schema_m == 2 and not continuum_only:
        raise RunDataError("A tilted staircase has no aligned lattice oracle; request continuum-only")
    if continuum_only:
        return
    solver = run.record_m["effective"]["solver"]
    expected = {"stencil": "Standard", "spatial": "staggered", "formulation": "O1VS"}
    if any(solver.get(key) != value for key, value in expected.items()):
        raise RunDataError("The exact aligned lattice reference requires Standard/staggered/O1VS")
    face = run.plane["face_position"] / run.spacing[0]
    if abs(face - round(face)) > 8 * math.ulp(face):
        raise RunDataError("A shifted requested face needs continuum-only or a separately phase-shifted lattice reference")


def _arguments(run: PlaneRun, binary: Path, sampling: dict, continuum_only: bool) -> list[str]:
    plane = run.plane
    incident, transmitted = ((1.0, plane["epsilon"]) if plane["incidence"] == "vacuum"
                             else (plane["epsilon"], 1.0))
    values = {"polarization": plane["polarization"], "epsilon-incident": repr(incident),
              "epsilon-transmitted": repr(transmitted),
              "spacing": ",".join(format(value, ".17g") for value in run.spacing),
              "dt": format(run.dt, ".17g")}
    values.update(sampling)
    command = [str(binary)]
    for key, value in values.items():
        command.extend(["--" + key, value])
    return command + (["--continuum-only"] if continuum_only else [])


def _binary_identity(binary: Path) -> dict:
    result = subprocess.run([str(binary), "--identity"], capture_output=True, text=True, check=True)
    identity = json.loads(result.stdout)
    if identity.get("program") != "interface-reference" or identity.get("reference_schema_version") != 1:
        raise RunDataError("Unsupported reference program or output schema")
    if not isinstance(identity.get("source_revision"), str):
        raise RunDataError("Reference executable did not report its source identity")
    return identity


def _parse_value(name: str, value: str):
    if name.endswith("_valid"):
        if value not in ("0", "1"):
            raise RunDataError(f"Invalid reference validity value: {value!r}")
        return value == "1"
    if name.endswith(("_reason", "_root")):
        return value
    if value == "":
        return None
    try:
        number = float(value)
    except ValueError as error:
        raise RunDataError(f"Invalid reference column {name}") from error
    if not math.isfinite(number):
        raise RunDataError(f"Non-finite reference column {name}")
    return number


def _reference_row(row: dict) -> dict:
    parsed = {name: _parse_value(name, value) for name, value in row.items()}
    for kind in ("continuum", "lattice"):
        valid, reason = parsed[kind + "_valid"], parsed[kind + "_reason"]
        if valid == bool(reason):
            raise RunDataError("Inconsistent reference applicability and reason")
        fields = [f"{kind}_{quantity}_{part}" for quantity in ("r", "t", "k1", "k2")
                  for part in ("real", "imag")]
        if valid and any(parsed[name] is None for name in fields):
            raise RunDataError("Applicable reference has missing coefficients or roots")
    return parsed


def _read_table(path: Path, expected: tuple[float, ...]) -> tuple[dict, ...]:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        names = reader.fieldnames or []
        rows = list(reader)
    required = {"omega_code", "continuum_R", "continuum_T", "continuum_angle_degrees", "continuum_transmission_weight"}
    for kind in ("continuum", "lattice"):
        required |= {kind + suffix for suffix in ("_valid", "_reason", "_incident_root", "_transmitted_root")}
        required |= {f"{kind}_{quantity}_{part}" for quantity in ("r", "t", "k1", "k2") for part in ("real", "imag")}
    if len(names) != len(set(names)) or not required.issubset(names) or len(rows) != len(expected):
        raise RunDataError("Reference table columns or row count disagree with the request")
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise RunDataError("Malformed reference table")
    result = tuple(_reference_row(row) for row in rows)
    if any(not math.isclose(row["omega_code"], omega, rel_tol=2e-14, abs_tol=0)
           for row, omega in zip(result, expected)):
        raise RunDataError("Reference frequency order differs from the measurement")
    return result


def _phase_coefficients(run: PlaneRun, sample, reference: dict) -> tuple[complex, complex]:
    """Apply Schneider §5.8.2 phase transfer with roots supplied by the C++ utility."""
    first, second = (coefficient(reference, "lattice", key) for key in ("k1", "k2"))
    if first is None or second is None or first.imag != 0 or second.imag != 0:
        raise RunDataError("Extraction audit needs two propagating lattice roots")
    if sample.incident_m is None or sample.incident_m == 0:
        raise RunDataError("A valid measured bin has no nonzero raw incident transform")
    plane = run.plane
    side = -1 if plane["incidence"] == "vacuum" else 1
    ki, kt = side * first.real, side * second.real
    dx = plane["face_position"] - plane["incident_probe_x"]
    other = plane["face_position"] - plane["other_probe_x"]
    reflection = (sample.total_m / sample.incident_m - 1) * cmath.exp(-2j * ki * dx)
    transmission = sample.transmitted_m / sample.incident_m * cmath.exp(1j * (kt * other - ki * dx))
    return reflection, transmission


def extraction_residuals(run: PlaneRun, rows: tuple[dict, ...]) -> list[dict]:
    """Return raw-to-published discrepancies; this is not a physical accuracy test."""
    if run.schema_m == 0:
        return []
    if run.schema_m == 2:
        return _tilted_extraction_residuals(run)
    result = []
    for sample, reference in zip(run.samples_m, rows):
        if not sample.valid_m:
            continue
        reflection, transmission = _phase_coefficients(run, sample, reference)
        result.append({"omega_code": sample.omega_m,
                       "reflection_residual": abs(reflection - sample.reflection_m),
                       "transmission_residual": abs(transmission - sample.transmission_m)})
    return result


def _fit_probe_values(probes: list, values: list, normal_k: float, kappa: float) -> tuple:
    """Repeat the declared one-column fit at the stored lifted phase coordinates."""
    basis = [cmath.exp(1j * (normal_k * probe.normal_position_m + kappa * probe.tangent_position_m))
             for probe in probes]
    amplitude = sum(value * phase.conjugate() for value, phase in zip(values, basis)) / len(values)
    residual = math.sqrt(sum(abs(value - amplitude * phase)**2 for value, phase in zip(values, basis)) / len(values))
    return amplitude, residual


def _tilted_bin_audit(run: PlaneRun, index: int, first: list, second: list) -> dict:
    sample = run.samples_m[index]
    ki, kr, kt = sample.normal_wavenumbers_m
    incoming = _fit_probe_values(first, [probe.incident_m for probe in first], ki, run.plane["kappa"])
    reflected = _fit_probe_values(first, [probe.total_m - probe.incident_m for probe in first], kr, run.plane["kappa"])
    transmitted = _fit_probe_values(second, [probe.total_m for probe in second], kt, run.plane["kappa"])
    if incoming[0] == 0:
        raise RunDataError("A valid tilted bin has zero reconstructed incident amplitude")
    result = {"omega_code": sample.omega_m, "reflection_residual": abs(reflected[0] / incoming[0] - sample.reflection_m),
              "transmission_residual": abs(transmitted[0] / incoming[0] - sample.transmission_m)}
    for name, fitted, stored, residual in zip(("incident", "reflected", "transmitted"),
                                             (incoming, reflected, transmitted), sample.face_amplitudes_m,
                                             sample.fit_residuals_m):
        result[name + "_face_residual"] = abs(fitted[0] - stored)
        result[name + "_fit_residual_difference"] = abs(fitted[1] / abs(incoming[0]) - residual)
    return result


def _tilted_extraction_residuals(run: PlaneRun) -> list[dict]:
    """Audit supplied directional roots and raw point fits; this is no staircase oracle."""
    groups = {}
    for probe in run.probes_m:
        groups.setdefault((probe.bin_m, probe.side_m), []).append(probe)
    return [_tilted_bin_audit(run, index, groups[index, "incident"], groups[index, "transmitted"])
            for index, sample in enumerate(run.samples_m) if sample.valid_m]


def _measurement_hashes(run: PlaneRun) -> dict:
    names = ("run.json", "spectra.csv") + (("probe-spectra.csv",) if run.schema_m == 2 else ())
    return {name: sha256_file(run.directory_m / name) for name in names}


def _record(run: PlaneRun, binary: Path, command: list[str]) -> dict:
    return {"reference_record_version": 1, "status": "started", "command": command,
            "started_utc": datetime.now(timezone.utc).isoformat(), "binary": str(binary),
            "binary_sha256": sha256_file(binary), "analysis_source": source_identity(),
            "measurement_directory": str(run.directory_m),
            "measurement_build": run.record_m.get("build", {}),
            "measurement_sha256": _measurement_hashes(run),
            "reference_document_sha256": sha256_file(Path(__file__).with_name("references.org")),
            "frequency_units": "normalized", "polarization": run.plane["polarization"],
            "component": run.plane["component"], "geometry": run.plane["geometry"],
            "phase_origin": "requested_face_anchor" if run.schema_m == 2 else "requested_face",
            "extraction_audit_basis": "recorded_directional_roots" if run.schema_m == 2 else "aligned_lattice_roots",
            "physics_validation": "not_performed"}


def _calculate(expected: tuple, binary: Path, directory: Path, command: list[str], record: dict) -> tuple:
    record["reference_build"] = _binary_identity(binary)
    with (directory / "reference.csv").open("x") as output, (directory / "stderr.txt").open("x") as errors:
        result = subprocess.run(command, stdout=output, stderr=errors, check=False)
    record["returncode"] = result.returncode
    if result.returncode:
        raise RunDataError(f"Reference calculator failed with exit code {result.returncode}")
    if sha256_file(binary) != record["binary_sha256"]:
        raise RunDataError("Reference executable changed during evaluation; use an immutable copy")
    return _read_table(directory / "reference.csv", expected)


def _evaluate(run: PlaneRun, binary: Path, directory: Path, command: list[str], record: dict,
              expected: tuple, audit: bool) -> ReferenceTable:
    try:
        rows = _calculate(expected, binary, directory, command, record)
        _check_angles(rows, record.get("angles_degrees"))
        if "angles_degrees" not in record and run.schema_m == 2:
            _check_tilted_angle(run, rows)
        record["extraction_residuals"] = extraction_residuals(run, rows) if audit else []
        record["status"] = "calculated"
        return ReferenceTable(directory, record, rows)
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        record.update(status="failed", error=str(error))
        raise
    finally:
        (directory / "reference.json").write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")


def _check_tilted_angle(run: PlaneRun, rows: tuple) -> None:
    central = rows[run.central_bin_m]
    if not central["continuum_valid"]:
        return
    signed = math.copysign(central["continuum_angle_degrees"], run.plane["kappa"])
    if not math.isclose(signed, run.plane["actual_angle_degrees"], rel_tol=2e-12, abs_tol=2e-12):
        raise RunDataError("Tilted recorded angle disagrees with the C++ continuum reference")


def _check_angles(rows: tuple, angles: tuple | None) -> None:
    if angles is None:
        return
    for row, angle in zip(rows, angles):
        if row.get("kz") is None or row["lattice_valid"] or row["lattice_reason"] != "not_requested":
            raise RunDataError("Angle curve requires a recorded wavenumber and continuum-only references")
        if row["continuum_valid"] and not math.isclose(row["continuum_angle_degrees"], angle,
                                                       rel_tol=1e-12, abs_tol=1e-12):
            raise RunDataError("Reference angle order differs from the request")


def evaluate(run: PlaneRun, binary: Path | str, output_directory: Path | str,
             continuum_only: bool = False) -> ReferenceTable:
    """Create a new reference directory; existing directories are never overwritten."""
    _lattice_applicability(run, continuum_only)
    binary, directory = Path(binary).resolve(strict=True), Path(output_directory).resolve()
    directory.mkdir(parents=True, exist_ok=False)
    frequencies = directory / "frequencies.txt"
    frequencies.write_text("".join(f"{value:.17g}\n" for value in run.frequencies))
    tangential_key = "kappa" if run.schema_m == 2 else "kz"
    sampling = {"kz": format(run.plane[tangential_key], ".17g"), "frequencies": str(frequencies)}
    command = _arguments(run, binary, sampling, continuum_only)
    record = _record(run, binary, command)
    record["tangential_wavenumber_source"] = "effective/plane/" + tangential_key
    return _evaluate(run, binary, directory, command, record, run.frequencies, not continuum_only or run.schema_m == 2)


def angle_curve(run: PlaneRun, binary: Path | str, angles_degrees,
                output_directory: Path | str) -> ReferenceTable:
    """Continuum curve at the recorded central frequency; C++ resolves each tangential wavevector."""
    angles = tuple(float(value) for value in angles_degrees)
    if not angles or any(not math.isfinite(value) or not 0 <= value < 90 for value in angles):
        raise ValueError("A nonempty angle vector in [0, 90) degrees is required")
    binary, directory = Path(binary).resolve(strict=True), Path(output_directory).resolve()
    directory.mkdir(parents=True, exist_ok=False)
    path = directory / "angles_degrees.txt"
    path.write_text("".join(f"{value:.17g}\n" for value in angles))
    omega = run.plane["omega0"]
    command = _arguments(run, binary, {"angles": str(path), "omega": format(omega, ".17g")}, True)
    record = _record(run, binary, command)
    record.update(sampling="continuum_angle_curve", angles_degrees=angles, omega_code=omega,
                  angle_definition="continuum phase angle relative to interface normal")
    return _evaluate(run, binary, directory, command, record, (omega,) * len(angles), False)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--continuum-only", action="store_true")
    args = parser.parse_args()
    try:
        result = evaluate(load_run(args.run), args.binary, args.output, args.continuum_only)
        print(result.directory_m / "reference.csv")
        return 0
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        parser.exit(2, f"references: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
