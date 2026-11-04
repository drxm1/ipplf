#!/usr/bin/env python3
"""Focused checks of reference signs/branches and recorded-data interpretation.

The independent expected continuum values follow the boundary equations in
references.org (Ellingson §§5.6–5.7, 5.10). Synthetic records test transport of
units, masks and raw data; they are not solver validation results.
"""
import argparse
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
from references import angle_curve, extraction_residuals
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
