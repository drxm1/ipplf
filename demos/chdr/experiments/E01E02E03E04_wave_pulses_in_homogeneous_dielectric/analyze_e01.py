#!/usr/bin/env python3
"""Analyze the completed E01 refinement study without invoking the solver.

Run this copied script in its result bundle, or pass --directory PATH.
Python 3 and Matplotlib are required. Inputs and measured quantities follow E01.typ.
Observed orders follow Salari and Knupp, SAND2000-1444, section 4.2, pp. 29--30
(library key salari2000CodeVerificationMms). They describe the recorded combined
longitudinal space/time refinement; no acceptance threshold is inferred here.
Each reference was evaluated at its numerical sample's actual position and time.
Different meshes retain their recorded nearest times; this script does not interpolate.
"""

import argparse
import csv
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import platform
import sys

MEASUREMENTS = ("A_L2", "E_L2", "B_L2", "E_Linf")
SAMPLE_TIME_KEYS = {"A_L2": "A", "E_L2": "E", "E_Linf": "E", "B_L2": "B"}
ORIGINAL_CHECKPOINTS = ("quarter_1", "quarter_2", "quarter_3", "final")
REFINEMENT_LADDER = (32, 64, 128, 256, 512, 1024, 2048, 4096)
COMPONENTS = {"Ax": "A_L2", "Ex": "E_L2", "By": "B_L2"}
COMPONENT_LABELS = {"Ax": r"$A_x$", "Ex": r"$E_x$", "By": r"$B_y$ (R2)"}
LABELS = {"A_L2": r"Relative $A$ $L^2$ error", "E_L2": r"Relative $E$ $L^2$ error",
          "B_L2": r"Relative R2 $B$ $L^2$ error", "E_Linf": r"Relative $E$ $L^\infty$ error"}


@dataclass
class LevelResults:
    """One completed level and its numerical data, keyed by quantity and checkpoint."""
    cells_per_pulse: int
    record: dict
    checkpoints: tuple
    errors: dict
    profiles: dict


def require(condition, message):
    """Reject incompatible or incomplete results before reporting convergence."""
    if not condition:
        raise ValueError(message)


def finite(value):
    """Parse a finite recorded number without accepting NaN or infinite errors."""
    result = float(value)
    require(math.isfinite(result), f"Nonfinite recorded value: {value}")
    return result


def read_json(path):
    """Read a manifest or provenance record from the supplied bundle."""
    return json.loads(path.read_text(encoding="utf-8"))


def planned_levels(root):
    """Require a contiguous prefix of the supported refinement ladder."""
    study = read_json(root / "study.json")
    require(study["schema_version"] == 1 and study["case"] == "E01-study", "Not an E01 study")
    require(study["status"] == "completed", "The study has not completed")
    entries = study["levels"]
    counts = [entry["cells_per_pulse"] for entry in entries]
    require(4 <= len(counts) <= len(REFINEMENT_LADDER)
            and counts == list(REFINEMENT_LADDER[:len(counts)]), "Unexpected refinement ladder")
    paths = [Path(entry["directory"]) for entry in entries]
    require(len(set(paths)) == len(paths), "Repeated level directory")
    require(all(not path.is_absolute() and ".." not in path.parts for path in paths),
            "Level directories must be relative paths inside the study")
    require(all(entry["status"] == "completed" for entry in entries), "Incomplete study level")
    return entries


def checkpoint_names(record):
    """Read named observations; recognize the original four-target records explicitly."""
    observations = record["effective"]["error_observations"]
    names = observations.get("protocol_checkpoint_names", list(ORIGINAL_CHECKPOINTS))
    require(isinstance(names, list) and names
            and all(isinstance(name, str) and name.isidentifier() for name in names),
            "Invalid protocol checkpoint names")
    require(len(set(names)) == len(names) and "initial" not in names, "Repeated checkpoint name")
    for key in ("protocol_target_times", "protocol_potential_steps", "protocol_field_steps"):
        require(len(observations[key]) == len(names), f"Checkpoint count differs from {key}")
    targets = [finite(value) for value in observations["protocol_target_times"]]
    require(targets[0] > 0 and all(b > a for a, b in zip(targets, targets[1:])),
            "Protocol target times must be positive and increasing")
    distances = observations.get("protocol_travel_distances_in_w")
    if distances is not None:
        require(len(distances) == len(names), "Checkpoint/travel-distance count mismatch")
        require(all(finite(value) > 0 for value in distances), "Invalid travel distance")
    return ("initial", *names)


def checkpoint_label(level, checkpoint):
    """Display the recorded travel distance when available, preserving its association."""
    if checkpoint == "initial":
        return "initial"
    observations = level.record["effective"]["error_observations"]
    distances = observations.get("protocol_travel_distances_in_w")
    if distances is None:
        return checkpoint.replace("_", " ")
    distance = distances[level.checkpoints.index(checkpoint) - 1]
    return f"after {distance:g}w"


def validate_record(record, count):
    """Check completion and the selected E01 case before reading its observations."""
    require(record["status"] == "completed" and record["input"]["case"] == "E01", "Incomplete E01 run")
    require(record["input"]["parameters"]["cells_per_pulse"] == count, "Level/input mismatch")
    effective = record["effective"]
    observations = effective["error_observations"]
    require(not observations["partial_final"], "Partial runs cannot enter the refinement study")
    final_step = max(observations["protocol_potential_steps"] + observations["protocol_field_steps"])
    require(effective["completed_solves"] == effective["requested_steps"] == final_step,
            "Run did not reach its prescribed endpoint")
    solver = effective["solver"]
    require(effective["value_system"] == "normalized", "E01 analysis expects normalized values")
    require((solver["formulation"], solver["reconstruction"], solver["spatial"], solver["stencil"])
            == ("O1VS", "R2", "staggered", "Standard"), "Unexpected E01 solver settings")
    require(all(effective["box"]["periodic"]), "E01 requires periodic boundaries")


def validate_resolution(record, count):
    """Match declared cell counts to the recorded mesh and physical pulse width."""
    effective = record["effective"]
    box, reference = effective["box"], effective["analytic_reference"]
    parameters = record["input"]["parameters"]
    if "transverse_cells_per_axis" in parameters:
        require(box["cells"][:2] == [parameters["transverse_cells_per_axis"]] * 2,
                "Transverse mesh differs from the selected input")
    hz = finite(box["spacing"][2])
    require(hz > 0 and math.isclose(count * hz, 2 * reference["support_half_width"],
                                    rel_tol=1e-12), "Pulse-cell count does not match its mesh")
    require(math.isclose(box["cells"][2] * hz, reference["period"], rel_tol=1e-12),
            "Longitudinal mesh length differs from the analytic period")


def error_row(row, record):
    """Retain the actual comparison time and cross-check the recorded mesh metadata."""
    result = {key: finite(row[key]) for key in
              ("target_time", "actual_time", "relative_error", "hx", "hy", "hz", "dt")}
    result.update({key: int(row[key]) for key in ("step", "nx", "ny", "nz")})
    effective = record["effective"]
    require([result[key] for key in ("nx", "ny", "nz")] == effective["box"]["cells"], "CSV mesh mismatch")
    require([result[key] for key in ("hx", "hy", "hz")] == effective["box"]["spacing"], "CSV spacing mismatch")
    require(result["dt"] == effective["solver"]["time_step"] > 0, "CSV timestep mismatch")
    require(result["relative_error"] >= 0, "Negative relative error")
    return result


def recorded_sample_time(observations, measurement, checkpoint, checkpoints):
    """Read the solver's resolved native label, including its initial time."""
    sample_times = observations.get("sample_times")
    require(isinstance(sample_times, dict),
            "run.json lacks LevelContract-derived error_observations.sample_times")
    quantity = SAMPLE_TIME_KEYS[measurement]
    times = sample_times.get(quantity)
    require(isinstance(times, list) and len(times) == len(checkpoints),
            f"run.json has incomplete sample_times for {quantity}")
    return finite(times[checkpoints.index(checkpoint)])


def validate_observation(row, measurement, checkpoint, record, checkpoints):
    """Match observations to the run record without reconstructing LevelContract offsets."""
    observations = record["effective"]["error_observations"]
    require(row["actual_time"] == recorded_sample_time(observations, measurement, checkpoint,
                                                       checkpoints),
            "Observation sample time differs from its recorded LevelContract label")
    if checkpoint == "initial":
        require(row["step"] == 0 and row["target_time"] == 0, "Invalid initial observation")
        return
    index = checkpoints.index(checkpoint) - 1
    steps = observations["protocol_potential_steps" if measurement == "A_L2" else "protocol_field_steps"]
    require(row["step"] == steps[index], "Observation step differs from run record")
    require(row["target_time"] == observations["protocol_target_times"][index], "Target time differs")


def load_errors(directory, record, checkpoints):
    """Require exactly the four measurements at every recorded observation."""
    errors = {}
    with (directory / "errors.csv").open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            key = row["measurement"], row["checkpoint"]
            require(key[0] in MEASUREMENTS and key[1] in checkpoints and key not in errors,
                    f"Unexpected or duplicate error row: {key}")
            errors[key] = error_row(row, record)
            validate_observation(errors[key], *key, record, checkpoints)
    require(len(errors) == len(MEASUREMENTS) * len(checkpoints), "Missing error observations")
    return errors


def append_profile(profiles, row, errors, checkpoints):
    """Validate one native profile sample against its matching error observation."""
    key = row["component"], row["checkpoint"]
    require(key[0] in COMPONENTS and key[1] in checkpoints, f"Unexpected profile: {key}")
    error = errors[COMPONENTS[key[0]], key[1]]
    require(int(row["step"]) == error["step"] and finite(row["actual_time"]) == error["actual_time"]
            and finite(row["target_time"]) == error["target_time"], "Profile/error observation mismatch")
    xy = finite(row["x"]), finite(row["y"])
    profile = profiles.setdefault(key, {"z": [], "numerical": [], "analytic": [], "xy": xy})
    require(xy == profile["xy"] and int(row["k"]) == len(profile["z"]), "Profile line/index mismatch")
    for column in ("z", "numerical", "analytic"):
        profile[column].append(finite(row[column]))


def load_profiles(directory, record, errors, checkpoints):
    """Preserve native z coordinates; never interpolate profiles between meshes or times."""
    profiles = {}
    with (directory / "profiles.csv").open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            append_profile(profiles, row, errors, checkpoints)
    require(len(profiles) == len(COMPONENTS) * len(checkpoints), "Missing profiles")
    for profile in profiles.values():
        z = profile["z"]
        require(len(z) == record["effective"]["box"]["cells"][2], "Incomplete profile")
        require(all(right > left for left, right in zip(z, z[1:])), "Unordered profile positions")
    return profiles


def physical_signature(record):
    """Compare fixed physical settings while allowing the prescribed longitudinal refinement."""
    effective = record["effective"]
    box = effective["box"]
    solver = {key: value for key, value in effective["solver"].items()
              if key not in ("time_step", "safety_factor")}
    lengths = [count * spacing for count, spacing in zip(box["cells"], box["spacing"])]
    return (effective["analytic_reference"], effective["material"], solver, box["periodic"],
            box["cells"][:2], box["spacing"][:2], lengths, effective["profiles"],
            effective["error_observations"]["protocol_target_times"],
            effective["error_observations"].get("protocol_travel_distances_in_w"))


def load_level(directory, count):
    """Resolve one level's observation schedule before reading its numerical data."""
    record = read_json(directory / "run.json")
    validate_record(record, count)
    validate_resolution(record, count)
    checkpoints = checkpoint_names(record)
    errors = load_errors(directory, record, checkpoints)
    profiles = load_profiles(directory, record, errors, checkpoints)
    return LevelResults(count, record, checkpoints, errors, profiles)


def load_study(root):
    """Load validated completed runs and ensure that they form one physical comparison."""
    levels = [load_level(root / entry["directory"], entry["cells_per_pulse"])
              for entry in planned_levels(root)]
    for coarse, fine in zip(levels, levels[1:]):
        require(coarse.checkpoints == fine.checkpoints, "Checkpoint schedule differs between levels")
        require(physical_signature(coarse.record) == physical_signature(fine.record),
                "Physical settings differ between levels")
        for key in ("hz", "dt"):
            require(coarse.errors["A_L2", "initial"][key] == 2 * fine.errors["A_L2", "initial"][key],
                    "Unexpected refinement ratio")
        require(all(coarse.profiles[key]["xy"] == fine.profiles[key]["xy"] for key in coarse.profiles),
                "Transverse profile line changed")
    return levels


def observed_order(coarse, fine):
    r"""Estimate \f$p=\log(e_c/e_f)/\log(h_c/h_f)\f$; Salari/Knupp section 4.2.

    Both errors must be positive. Zero-error pairs, including exactly sampled initial A,
    have no defined logarithmic order. E01 refines longitudinal spacing and time together.
    """
    if coarse["relative_error"] == 0 or fine["relative_error"] == 0:
        return None
    return ((math.log(coarse["relative_error"]) - math.log(fine["relative_error"]))
            / math.log(coarse["hz"] / fine["hz"]))


def convergence_rows(levels):
    """Keep both actual times beside every adjacent-level order estimate."""
    for index, level in enumerate(levels):
        for measurement in MEASUREMENTS:
            for checkpoint in level.checkpoints:
                row = level.errors[measurement, checkpoint]
                coarse = levels[index - 1].errors[measurement, checkpoint] if index else None
                yield {"cells_per_pulse": level.cells_per_pulse, "measurement": measurement,
                       "checkpoint": checkpoint, **{key: row[key] for key in
                       ("hz", "dt", "target_time", "actual_time", "relative_error")},
                       "coarser_cells_per_pulse": levels[index - 1].cells_per_pulse if index else "",
                       "coarser_actual_time": coarse["actual_time"] if coarse else "",
                       "observed_order": observed_order(coarse, row) if coarse else None}


def save_figure(figure, directory, name, plt):
    """Write shareable vector and raster copies without requiring a display server."""
    figure.savefig(directory / f"{name}.pdf", bbox_inches="tight")
    figure.savefig(directory / f"{name}.png", dpi=180, bbox_inches="tight")
    plt.close(figure)


def plot_convergence(levels, directory, plt):
    """Show each checkpoint's positive errors against recorded longitudinal spacing."""
    figure, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
    for axis, measurement in zip(axes.flat, MEASUREMENTS):
        for index, checkpoint in enumerate(levels[0].checkpoints):
            rows = [level.errors[measurement, checkpoint] for level in levels]
            positive = [row for row in rows if row["relative_error"] > 0]
            if positive:
                axis.loglog([row["hz"] for row in positive], [row["relative_error"] for row in positive],
                            "o-", color=f"C{index}", label=checkpoint_label(levels[0], checkpoint))
        spacings = [level.errors[measurement, "initial"]["hz"] for level in levels]
        axis.set_xticks(spacings, [f"{spacing:.3g}" for spacing in spacings],
                       rotation=30, ha="right", fontsize="small")
        axis.minorticks_off()
        axis.set(xlabel=r"$h_z$ (solver units)", ylabel=LABELS[measurement])
        axis.grid(True, which="both", alpha=0.25)
        axis.legend(fontsize="small")
    figure.suptitle("E01 longitudinal space/time refinement\n"
                   "Zero errors are omitted on log axes; values and orders are in convergence.csv")
    save_figure(figure, directory, "convergence", plt)


def plot_error_histories(levels, directory, plt):
    """Plot sparse observations at their actual times, including zero initial A error."""
    figure, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
    for axis, measurement in zip(axes.flat, MEASUREMENTS):
        for level in levels:
            rows = [level.errors[measurement, checkpoint] for checkpoint in level.checkpoints]
            axis.plot([row["actual_time"] for row in rows], [row["relative_error"] for row in rows],
                      "o-", label=f"N_pulse={level.cells_per_pulse}")
        axis.set(xlabel="Actual time (solver units)", ylabel=LABELS[measurement])
        axis.grid(True, alpha=0.25)
        axis.legend(fontsize="small", ncol=2 if len(levels) > 4 else 1)
    figure.suptitle("E01 sparse error observations")
    save_figure(figure, directory, "error_histories", plt)


def plot_profiles(levels, directory, plt):
    """Pair each numerical profile with its own reference at the same native samples."""
    for checkpoint in levels[0].checkpoints:
        figure, axes = plt.subplots(3, 1, figsize=(11, 9), layout="constrained")
        for axis, (component, measurement) in zip(axes, COMPONENTS.items()):
            for index, level in enumerate(levels):
                profile = level.profiles[component, checkpoint]
                time = level.errors[measurement, checkpoint]["actual_time"]
                label = f"N_pulse={level.cells_per_pulse}, t={time:.8g}"
                axis.plot(profile["z"], profile["numerical"], color=f"C{index}", label=label)
                axis.plot(profile["z"], profile["analytic"], "--", color=f"C{index}", linewidth=1)
            axis.set(xlabel="z (solver units)", ylabel=f"{COMPONENT_LABELS[component]} (solver units)")
            axis.grid(True, alpha=0.25)
            axis.legend(fontsize="small", ncol=2)
        figure.suptitle(f"E01 {checkpoint_label(levels[0], checkpoint)}: numerical (solid), continuum reference (dashed)\n"
                       "Each legend reports the profile's actual time; samples are not interpolated")
        save_figure(figure, directory, f"profiles_{checkpoint}", plt)


def analyze(root):
    """Validate first, then produce tables and figures from existing completed results."""
    (root / "analysis.json").unlink(missing_ok=True)
    levels = load_study(root)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    directory = root / "figures"
    directory.mkdir(exist_ok=True)
    rows = list(convergence_rows(levels))
    with (root / "convergence.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    plot_convergence(levels, directory, plt)
    plot_error_histories(levels, directory, plt)
    plot_profiles(levels, directory, plt)
    write_analysis_record(root, levels, matplotlib.__version__)


def write_analysis_record(root, levels, matplotlib_version):
    """Record completed analysis, dependency versions and the interpretation of empty orders."""
    record = {"status": "completed", "python_version": platform.python_version(),
              "matplotlib_version": matplotlib_version,
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "cells_per_pulse_levels": [level.cells_per_pulse for level in levels],
              "checkpoint_names": list(levels[0].checkpoints),
              "comparison": "recorded nearest native times; no time interpolation",
              "empty_orders": "first level or a zero error in an adjacent pair",
              "interpretation": "observed orders only; no automatic convergence verdict",
              "figures": [f"figures/{name}.{suffix}" for name in
                          ("convergence", "error_histories", *(f"profiles_{key}" for key in levels[0].checkpoints))
                          for suffix in ("png", "pdf")]}
    (root / "analysis.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")


def main():
    """Default to the copied script's directory, so the result bundle can be moved."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=Path(__file__).resolve().parent,
                        help="study directory containing study.json (default: this script's directory)")
    args = parser.parse_args()
    try:
        analyze(args.directory)
    except (OSError, ValueError, KeyError, ImportError) as error:
        parser.exit(1, f"E01 analysis failed: {error}\n")
    print(f"E01 analysis written to {args.directory / 'convergence.csv'} and figures/")


if __name__ == "__main__":
    main()
