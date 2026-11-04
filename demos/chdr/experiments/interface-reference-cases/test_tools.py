#!/usr/bin/env python3
"""Focused checks of reference signs/branches and recorded-data interpretation.

The independent expected continuum values follow the boundary equations in
references.org (Ellingson §§5.6–5.7, 5.10). Synthetic records test transport of
units, masks and raw data; they are not solver validation results.
"""
import argparse
import cmath
import csv
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest

from read_run import RunDataError, load_run, sha256_file
from references import angle_curve, evaluate, extraction_residuals
from scan import ScanSettings, _run_attempt, output_input

REFERENCE = None


def calculate(polarization, epsilon1=1, epsilon2=4, kz=0, omega=.2, continuum_only=False):
    command = [str(REFERENCE), "--polarization", polarization, "--epsilon-incident", str(epsilon1),
               "--epsilon-transmitted", str(epsilon2), "--spacing", "1,1,1", "--dt", ".4",
               "--kz", str(kz), "--frequencies", "-"]
    if continuum_only:
        command.append("--continuum-only")
    result = subprocess.run(command, input=f"{omega:.17g}\n", text=True, capture_output=True, check=True)
    return next(csv.DictReader(io.StringIO(result.stdout)))


def synthetic_record():
    plane = dict(schema_version=1, central_bin=0, geometry="aligned_plane", incidence_plane="xz",
                 polarization="TE", component="E_y", incidence="vacuum", frequency_units="normalized",
                 coordinate_units="normalized", epsilon=4, omega0=.5, kx0=.5, kz=0, sigma_k=.05,
                 tail_sigmas=8, requested_cells_per_wavelength=20, cells_per_wavelength=20,
                 separation_factor=1, face_position=0, incident_probe_x=1, other_probe_x=-1,
                 probe_z=.5, requested_angle_degrees=0, actual_angle_degrees=0, cutoff_fraction=0,
                 gate_incident=10, gate_other=10, timing_samples=513)
    spectra = dict(amplitudes_units="dimensionless", angular_frequency_units="rad/s",
                   raw_units="normalized", transform="dt sum w f exp(+i omega t)",
                   raw_phase_origin="probe_coordinates", coefficient_phase_origin="requested_face")
    return dict(record_version=1, status="completed", application="planar-interface", input={"case": "plane"},
                output_factors={"system": "si", "frequency": 10}, build={"source_revision": "synthetic"},
                effective=dict(value_system="normalized", completed_solves=10, plane=plane,
                               box={"cells": [32, 2, 2], "spacing": [1, 1, 1]},
                               solver={"time_step": .1}, output={"spectra": spectra}))


def synthetic_sample():
    return dict(angular_frequency_rad_per_s=5, omega_code=.5, r_real=.2, r_imag=0,
                t_real=1.2, t_imag=0, valid=1, reason="", incident_real=2, incident_imag=0,
                total_real=2.2161209223472559, total_imag=.3365883939231586,
                transmitted_real=.16976928400248697, transmitted_imag=2.3939879678497307)


def write_synthetic(path, record, sample):
    (path / "run.json").write_text(json.dumps(record))
    with (path / "spectra.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, list(sample))
        writer.writeheader()
        writer.writerow(sample)


def tilted_record():
    record = synthetic_record()
    plane = record["effective"]["plane"]
    c, s = math.cos(math.radians(17)), math.sin(math.radians(17))
    plane.update(schema_version=2, geometry="periodic_tilted_plane", normal=[c, 0, s], tangent=[-s, 0, c],
                 anchor=[0, 0, 0], requested_cell_shape=[1, 1, 1], interface_tilt_degrees=17, kappa=.1,
                 requested_face_position=0, normal_period=32*c, tangent_period=32/s, primitive_tangent_period=1/s,
                 normal_dielectric_width=16*c, normal_incident_width=16*c, primitive_cells_x=1, primitive_cells_z=1,
                 period_cells_multiplier=32, tangent_mode=1, retained_mode_count=9, mode_cutoff_sigmas=9,
                 floquet_margin=.1, retained_floquet_margin=.1, fit_residual_limit=.01, band_sigmas=3,
                 incident_floor=1e-8, probe_cells=2, probe_half_width_cells=1,
                 probe_point_count_incident=3, probe_point_count_transmitted=3,
                 projection="unweighted_complex_plane_wave_fit", projection_basis="exp(i*(normal_k*u+kappa*v))",
                 projection_normalization="mean_over_points", actual_angle_degrees=math.degrees(math.asin(.2)),
                 angular_frequencies=[.5])
    for key in ("kx0", "kz", "face_position", "incident_probe_x", "other_probe_x", "probe_z"):
        del plane[key]
    record["effective"]["box"].update(cells=[32, 2, 32], spacing=[1, 1, c/s])
    record["effective"]["output"]["spectra"].update(face_amplitudes_units="normalized",
        raw_phase_origin="lifted_probe_coordinates", coefficient_phase_origin="requested_face_anchor")
    return record


def tilted_values():
    incoming = complex(2, .3)
    return incoming, incoming * complex(.2, .05), incoming * complex(1.2, -.1)


def tilted_sample():
    sample = dict(angular_frequency_rad_per_s=5, omega_code=.5, r_real=.2, r_imag=.05,
                  t_real=1.2, t_imag=-.1, valid=1, reason="", fit_residual_incident=0, fit_residual_r=0, fit_residual_t=0)
    for name, value, root, speed in zip(("incident", "reflected", "transmitted"), tilted_values(),
                                      (-.4, .55, -.7), (-.8, .9, -.6)):
        sample.update({name + "_face_real": value.real, name + "_face_imag": value.imag,
                       name + "_normal_k": root, name + "_normal_group_velocity": speed})
    return sample


def tilted_probe_rows(record):
    rows = []
    c, _, s = record["effective"]["plane"]["normal"]
    incoming, reflected, transmitted = tilted_values()
    for side, sign in (("incident", 1), ("transmitted", -1)):
        for point in range(3):
            u, v = sign * (1 + point / 2), 8 + point / 4
            incident = incoming * cmath.exp(1j * (-.4*u + .1*v)) if side == "incident" else 0j
            scattered = (reflected * cmath.exp(1j * (.55*u + .1*v)) if side == "incident"
                         else transmitted * cmath.exp(1j * (-.7*u + .1*v)))
            total = incident + scattered
            rows.append(dict(bin=0, omega_code=.5, side=side, point=point, x=c*u-s*v+32, y=.5, z=s*u+c*v,
                             lift_x=-1, lift_z=0, u=u, v=v, incident_real=incident.real, incident_imag=incident.imag,
                             total_real=total.real, total_imag=total.imag))
    return rows


def write_probe_rows(path, rows):
    with (path / "probe-spectra.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def shifted_tilted_probes(record, rows, offset):
    plane = record["effective"]["plane"]
    anchor = [offset*n + 1.5*plane["tangent_period"]*t for n, t in zip(plane["normal"], plane["tangent"])]
    plane.update(anchor=anchor, requested_face_position=offset)
    box = record["effective"]["box"]
    for row in rows:
        for axis, name in ((0, "x"), (2, "z")):
            position = anchor[axis] + row["u"]*plane["normal"][axis] + row["v"]*plane["tangent"][axis]
            length = box["cells"][axis]*box["spacing"][axis]
            winding = math.floor(position/length)
            row[name], row["lift_" + name] = position - winding*length, winding
    return rows


def negative_incidence_rows(record):
    plane = record["effective"]["plane"]
    plane.update(kappa=-.1, tangent_mode=-1, actual_angle_degrees=-plane["actual_angle_degrees"])
    rows = tilted_probe_rows(record)
    for row in rows:
        factor = cmath.exp(-.2j * row["v"])
        for name in ("incident", "total"):
            value = complex(row[name + "_real"], row[name + "_imag"]) * factor
            row[name + "_real"], row[name + "_imag"] = value.real, value.imag
    return rows


class ReferenceTests(unittest.TestCase):
    def test_angle_curve_uses_cpp_setup_and_preserves_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            record = synthetic_record()
            record["effective"]["plane"].update(polarization="TM", component="H_y")
            write_synthetic(path, record, synthetic_sample())
            angles = (0, math.degrees(math.atan(2)), 80)
            result = angle_curve(load_run(path), REFERENCE, angles, path / "curve")
            self.assertEqual(result.record_m["angles_degrees"], angles)
            self.assertTrue((result.directory_m / "angles_degrees.txt").is_file())
            self.assertEqual(len(result.rows_m), 3)
            self.assertAlmostEqual(result.rows_m[0]["continuum_r_real"], 1/3, places=14)
            self.assertAlmostEqual(result.rows_m[1]["continuum_r_real"], 0, places=14)
            for row, angle in zip(result.rows_m, angles):
                self.assertEqual(row["lattice_reason"], "not_requested")
                self.assertAlmostEqual(row["continuum_angle_degrees"], angle, places=12)
                self.assertAlmostEqual(row["continuum_R"] + row["continuum_T"], 1, places=14)

    def test_normal_component_sign_and_power(self):
        for polarization, reflection, transmission in [("TE", -1/3, 2/3), ("TM", 1/3, 4/3)]:
            row = calculate(polarization)
            self.assertEqual(row["continuum_valid"], "1")
            self.assertAlmostEqual(float(row["continuum_r_real"]), reflection, places=14)
            self.assertAlmostEqual(float(row["continuum_t_real"]), transmission, places=14)
            self.assertAlmostEqual(float(row["continuum_R"]) + float(row["continuum_T"]), 1, places=14)

    def test_brewster_and_complex_transmitted_branch(self):
        row = calculate("TM", kz=.2 * 2 / math.sqrt(5))
        self.assertAlmostEqual(float(row["continuum_r_real"]), 0, places=14)
        self.assertAlmostEqual(float(row["continuum_t_real"]), 1, places=14)
        row = calculate("TE", epsilon1=4, epsilon2=1, kz=.2 * math.sqrt(2))
        self.assertEqual(row["continuum_transmitted_root"], "evanescent")
        self.assertAlmostEqual(float(row["continuum_r_real"]), 1/3, places=14)
        self.assertAlmostEqual(float(row["continuum_r_imag"]), -2 * math.sqrt(2)/3, places=14)
        self.assertEqual(float(row["continuum_T"]), 0)
        self.assertEqual(row["lattice_valid"], "1")

    def test_applicability_is_separate(self):
        row = calculate("TE", epsilon2=2.13, kz=.5235987755982988, omega=.5234373912173248)
        self.assertEqual(row["continuum_valid"], "0")
        self.assertEqual(row["continuum_reason"], "incident_not_propagating")
        self.assertEqual(row["continuum_R"], "")
        self.assertEqual(row["lattice_valid"], "1")
        row = calculate("TE", continuum_only=True)
        self.assertEqual(row["lattice_reason"], "not_requested")
        self.assertEqual(row["lattice_r_real"], "")


class ReaderTests(unittest.TestCase):
    def test_tilted_signed_incidence_retains_negative_kappa(self):
        with tempfile.TemporaryDirectory() as directory:
            path, record = Path(directory), tilted_record()
            plane = record["effective"]["plane"]
            rows = negative_incidence_rows(record)
            write_synthetic(path, record, tilted_sample())
            write_probe_rows(path, rows)
            result = evaluate(load_run(path), REFERENCE, path / "reference", continuum_only=True)
            command = result.record_m["command"]
            self.assertEqual(float(command[command.index("--kz") + 1]), -.1)
            self.assertAlmostEqual(result.rows_m[0]["continuum_angle_degrees"], abs(plane["actual_angle_degrees"]))
            self.assertLess(result.record_m["extraction_residuals"][0]["reflection_residual"], 2e-14)
            plane.update(kappa=.1, tangent_mode=1)
            write_synthetic(path, record, tilted_sample())
            with self.assertRaisesRegex(RunDataError, "recorded angle"):
                evaluate(load_run(path), REFERENCE, path / "wrong-sign", continuum_only=True)

    def test_signed_tilt_offset_and_both_box_windings(self):
        with tempfile.TemporaryDirectory() as directory:
            path, record = Path(directory), tilted_record()
            rows = shifted_tilted_probes(record, tilted_probe_rows(record), -.3)
            write_synthetic(path, record, tilted_sample())
            write_probe_rows(path, rows)
            run = load_run(path)
            self.assertTrue(all(probe.winding_m[0] and probe.winding_m[1] for probe in run.probes_m))
            residual = extraction_residuals(run, ())[0]
            self.assertLess(residual["reflection_residual"], 2e-14)
            record["effective"]["plane"]["requested_face_position"] = .3
            write_synthetic(path, record, tilted_sample())
            with self.assertRaises(RunDataError):
                load_run(path)

    def test_tilted_asymmetric_projection_and_continuum_only(self):
        with tempfile.TemporaryDirectory() as directory:
            path, record, sample = Path(directory), tilted_record(), tilted_sample()
            write_synthetic(path, record, sample)
            write_probe_rows(path, tilted_probe_rows(record))
            run = load_run(path)
            with self.assertRaises(RunDataError):
                evaluate(run, REFERENCE, path / "forbidden")
            result = evaluate(run, REFERENCE, path / "reference", continuum_only=True)
            residuals = result.record_m["extraction_residuals"][0]
            self.assertLess(max(value for key, value in residuals.items() if key != "omega_code"), 2e-14)
            self.assertEqual(result.rows_m[0]["lattice_reason"], "not_requested")
            self.assertIn("probe-spectra.csv", result.record_m["measurement_sha256"])
            self.assertAlmostEqual(result.rows_m[0]["continuum_angle_degrees"], run.plane["actual_angle_degrees"])
            sample["reflected_normal_k"] = .4
            write_synthetic(path, record, sample)
            self.assertGreater(extraction_residuals(load_run(path), ())[0]["reflection_residual"], .01)

    def test_tilted_coordinate_and_qualification_rejections(self):
        with tempfile.TemporaryDirectory() as directory:
            path, record = Path(directory), tilted_record()
            write_synthetic(path, record, tilted_sample())
            rows = tilted_probe_rows(record)
            rows[0]["lift_x"] = 0
            write_probe_rows(path, rows)
            with self.assertRaises(RunDataError):
                load_run(path)
            write_probe_rows(path, tilted_probe_rows(record))
            record["effective"]["plane"]["retained_floquet_margin"] = -.1
            write_synthetic(path, record, tilted_sample())
            with self.assertRaises(RunDataError):
                load_run(path)

    def test_si_raw_units_and_phase_audit(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            record, sample = synthetic_record(), synthetic_sample()
            write_synthetic(path, record, sample)
            run = load_run(path)
            self.assertEqual(run.frequencies, (.5,))
            self.assertEqual(run.samples_m[0].incident_m, 2)
            roots = ({"lattice_k1_real": .5, "lattice_k1_imag": 0,
                      "lattice_k2_real": 1, "lattice_k2_imag": 0},)
            audit = extraction_residuals(run, roots)[0]
            self.assertLess(audit["reflection_residual"], 1e-15)
            self.assertLess(audit["transmission_residual"], 1e-15)
            sample["angular_frequency_rad_per_s"] = 6
            write_synthetic(path, record, sample)
            with self.assertRaises(RunDataError):
                load_run(path)

    def test_invalid_bin_retained_and_unknown_schema_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            record, sample = synthetic_record(), synthetic_sample()
            sample.update(valid=0, reason="nonfinite_spectrum", r_real="nan", r_imag="nan",
                          t_real="nan", t_imag="nan", total_real="nan")
            write_synthetic(path, record, sample)
            run = load_run(path)
            self.assertFalse(run.samples_m[0].valid_m)
            self.assertIsNone(run.samples_m[0].reflection_m)
            self.assertTrue(math.isnan(run.samples_m[0].total_m.real))
            record["effective"]["plane"]["schema_version"] = 999
            write_synthetic(path, record, sample)
            with self.assertRaises(RunDataError):
                load_run(path)

    @unittest.skipUnless(importlib.util.find_spec("yaml"), "PyYAML is needed only by the runner")
    def test_output_replacement_preserves_physical_text(self):
        import yaml
        forms = ["schema_version: 1\ncase: plane\nparameters: {courant: 9e-1}\n",
                 "{schema_version: 1, case: plane, parameters: {courant: 9e-1}}\n",
                 "{case: plane, parameters: {courant: 9e-1}, }\n",
                 "{case: plane, parameters: {courant: 9e-1} # trailing comment\n, # comma\n}\n",
                 "case: plane\nparameters: {courant: 9e-1}\noutput:\n  directory: old\n"]
        for original in forms:
            resolved = output_input(original, Path('/tmp/example output'))
            self.assertIn("parameters: {courant: 9e-1}", resolved)
            self.assertEqual(yaml.safe_load(resolved)["output"]["directory"], '/tmp/example output')
        with self.assertRaises(ValueError):
            output_input("case: plane\ncase: plane\n", Path('/tmp/example'))


def live_process(pid):
    """A killed orphan can remain a zombie until its system reaper runs."""
    try:
        return Path(f"/proc/{pid}/stat").read_text().split(")", 1)[1].split()[0] != "Z"
    except FileNotFoundError:
        return False


def cancellation_binary(path):
    fixture = path / "fixture"
    fixture.mkdir()
    write_synthetic(fixture, synthetic_record(), synthetic_sample())
    binary = path / "simulation"
    child = "import os,signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); print(os.getpid(),flush=True); time.sleep(60)"
    binary.write_text(f"#!{sys.executable}\nimport pathlib,shutil,subprocess,sys,yaml\n"
                      "data=yaml.safe_load(pathlib.Path(sys.argv[2]).read_text())\n"
                      "mode=data['parameters']['mode']\n"
                      f"if mode=='completed': shutil.copytree({str(fixture)!r},data['output']['directory'])\n"
                      "elif mode=='failed': sys.exit(7)\n"
                      "else:\n print('launcher',file=sys.stderr,flush=True)\n"
                      f" subprocess.Popen([sys.executable,'-c',{child!r}]).wait()\n")
    binary.chmod(0o755)
    return binary


def cancellation_inputs(path):
    inputs = []
    for name in ("completed", "failed", "running", "queued"):
        item = path / (name + ".yaml")
        item.write_text(f"case: plane\nparameters: {{mode: {name}}}\n")
        inputs.append(str(item))
    return inputs


def wait_for_descendant(root, process):
    deadline = time.monotonic() + 10
    while time.monotonic() < deadline and process.poll() is None:
        for output in root.glob("0002-running/*/stdout.txt"):
            text = output.read_text().strip()
            if text:
                return int(text)
        time.sleep(.02)
    raise AssertionError("Cancellation reproduction did not start its descendant")


def interrupt_scan(path):
    root, binary = path / "scan", cancellation_binary(path)
    command = [sys.executable, str(Path(__file__).with_name("scan.py")), "--binary", str(binary),
               "--output-root", str(root), "--jobs", "1", "--max-cores", "1", *cancellation_inputs(path)]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, start_new_session=True)
    descendant = None
    try:
        descendant = wait_for_descendant(root, process)
        process.send_signal(signal.SIGINT)
        stdout, stderr = process.communicate(timeout=8)
        manifest = json.loads(next(root.glob("scan-*.json")).read_text())
        return process.returncode, descendant, manifest, stdout, stderr
    finally:
        if descendant is not None and live_process(descendant):
            os.killpg(os.getpgid(descendant), signal.SIGKILL)
        if process.poll() is None:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait(timeout=2)


@unittest.skipUnless(sys.platform.startswith("linux"), "Process-state check uses Linux procfs")
class TimeoutTests(unittest.TestCase):
    @unittest.skipUnless(importlib.util.find_spec("yaml"), "PyYAML is needed by the runner")
    def test_cli_sigint_preserves_rows_and_stops_all_attempts(self):
        with tempfile.TemporaryDirectory() as directory:
            code, descendant, record, stdout, stderr = interrupt_scan(Path(directory))
            self.assertEqual(code, 130, stderr.decode())
            self.assertEqual(record["status"], "cancelled")
            self.assertEqual([item["status"] for item in record["attempts"]],
                             ["completed", "failed", "cancelled", "cancelled"])
            self.assertEqual(record["attempts"][1]["returncode"], 7)
            self.assertFalse(live_process(descendant), "SIGINT left a running descendant")
            self.assertTrue(record["cancellation_cleanup"])
            self.assertTrue(stdout.strip())
            for item in record["attempts"]:
                saved = json.loads((Path(item["directory"]) / "attempt.json").read_text())
                self.assertEqual(saved["status"], item["status"])
            queued = Path(record["attempts"][3]["directory"])
            self.assertNotIn("launcher", (queued / "stderr.txt").read_text())

    def test_timeout_stops_descendant_after_launcher_exits(self):
        child = "import os,signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); print(os.getpid(),flush=True); time.sleep(60)"
        parent = f"import subprocess,sys; print('launcher',file=sys.stderr,flush=True); subprocess.Popen([sys.executable,'-c',{child!r}]).wait()"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            settings = ScanSettings(Path(sys.executable), path, timeout_m=.4)
            record = {"command": [sys.executable, "-c", parent], "working_directory": directory}
            result = _run_attempt(path, record, settings, sha256_file(settings.binary_m))
            pid = int((path / "stdout.txt").read_text().strip())
            try:
                self.assertEqual(result["status"], "failed")
                self.assertEqual(result["error_type"], "TimeoutExpired")
                self.assertEqual(result["timeout_cleanup"], "group_signalled_launcher_reaped")
                self.assertIn("launcher", (path / "stderr.txt").read_text())
                self.assertFalse(live_process(pid), "Timed-out launch left a running descendant")
                self.assertLess(result["elapsed_seconds"], 8)
            finally:
                if live_process(pid):
                    os.kill(pid, signal.SIGKILL)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--reference", required=True, type=Path)
    options, remaining = parser.parse_known_args()
    REFERENCE = options.reference.resolve()
    unittest.main(argv=[__file__, *remaining])
