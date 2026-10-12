"""Compare owned fields and positive reductions in one-rank and two-rank records.

Owned-cell differences use the field's maximum as their scale. Independent
positive sums and maxima use a relative comparison. Residuals arrive already
normalized by their absolute evaluation and therefore use an absolute comparison.
Continuum approximation errors are separate JSONL records and are not read here.

TODO: REVIEW.
"""
import json
import math
from pathlib import Path
import sys


# The producer writes Cartesian x/y/z vector components.
SpatialDim = 3


def finite_number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def expected_names():
    """List the registered producer cases, so matching omissions cannot pass."""
    states = [f"T1-{stencil}-{r}" for stencil in ("Std", "NS") for r in range(3)]
    states += [f"T2-{stencil}" for stencil in ("Std", "NS")]
    states += [f"{family}-{stencil}-{r}-0-{n}" for family in ("T3", "T4")
               for stencil in ("Std", "NS") for r in range(3) for n in (16, 32)]
    states += [f"{family}-Std-1-1-{n}" for family in ("T3", "T4") for n in (16, 32)]
    seams = [f"monitors-{prefix}{stencil}-{reading}" for prefix in ("", "nyquist-")
             for stencil, reading in (("Std", 0), ("NS", 0), ("Std", 1))]
    names = {f"{state}-{field}.json" for state in states + seams
             for field in ("current", "previous", "E", "B")}
    names.update(f"{state}-residuals.json" for state in seams)
    names.update(f"T5-Std-32-{direction}-{safety:.6f}.json"
                 for direction in range(5) for safety in (0, .98))
    names.update(f"T5-NS-axis-24-{safety:.6f}.json" for safety in (1, .99, .98))
    names.update(f"T5-NS-near-24-{length}-{safety:.6f}.json"
                 for length in (2, 4, 8) for safety in (1, .98))
    return names


def matching_names(left, right):
    """Require the same nonempty set of records from the two fresh fixture runs."""
    names = {path.name for path in left.glob("*.json")}
    other = {path.name for path in right.glob("*.json")}
    if not names or names != other:
        raise AssertionError(f"Missing or mismatched rank records: {names ^ other}")
    expected = expected_names()
    if names != expected:
        raise AssertionError(f"Producer coverage mismatch: {names ^ expected}")
    return sorted(names)


def compare_reduction(name, key, left, right):
    """A zero reference scale requires exact agreement; no epsilon floor is added."""
    if not finite_number(left) or not finite_number(right):
        raise AssertionError(f"Nonfinite or null reduction {name}: {key}")
    scale = max(abs(left), abs(right))
    if abs(left - right) > 1e-12 * scale:
        raise AssertionError(f"Reduction mismatch {name}: {key}")
    return abs(left - right) / scale if scale else 0


def compare_field(name, left, right):
    """Compare each globally ordered component, then independently reduced norms."""
    if set(left) != {"values", "max", "sum_abs", "sum_squares"}:
        raise AssertionError(f"Unexpected field schema: {name}")
    parts = name.removesuffix(".json").split("-")
    n = (int(parts[-2]) if parts[0] in ("T3", "T4") else
         32 if parts[0] == "T2" else 16)
    cells = n * n * (3 * n // 2 if "NS" in parts else n)
    expected = cells * (SpatialDim if parts[-1] in ("E", "B") else SpatialDim + 1)
    if len(left["values"]) != expected or len(left["values"]) != len(right["values"]):
        raise AssertionError(f"Shape mismatch: {name}")
    if any(not finite_number(record[key]) or record[key] < 0
           for record in (left, right) for key in ("max", "sum_abs", "sum_squares")):
        raise AssertionError(f"Invalid field statistic: {name}")
    if not all(finite_number(x) for x in left["values"] + right["values"]):
        raise AssertionError(f"Nonfinite or null field value: {name}")
    maxima = [max(abs(x) for x in record["values"]) for record in (left, right)]
    if any(record["max"] != maximum for record, maximum in zip((left, right), maxima)):
        raise AssertionError(f"Recorded maximum disagrees with owned values: {name}")
    scale = max(maxima)
    differences = (abs(x-y) for x, y in zip(left["values"], right["values"]))
    difference = max(differences, default=0)
    if difference > 1e-13 * scale:
        raise AssertionError(f"Per-cell mismatch {name}: {difference} / {scale}")
    reduction = max(compare_reduction(name, key, left[key], right[key])
                    for key in ("max", "sum_abs", "sum_squares"))
    return difference / scale if scale else 0, reduction


def compare_residuals(name, left, right):
    """Do not divide by a residual that is expected to be near zero."""
    a, b = left["normalized_residuals"], right["normalized_residuals"]
    if len(a) != 8 or len(a) != len(b):
        raise AssertionError(f"Residual count mismatch: {name}")
    for x, y in zip(a, b):
        if not finite_number(x) or not finite_number(y) or abs(x-y) > 1e-12:
            raise AssertionError(f"Normalized residual mismatch: {name}")
    a, b = left["residual_norms"], right["residual_norms"]
    if len(a) != 32 or len(a) != len(b):
        raise AssertionError(f"Residual norm count mismatch: {name}")
    return max(compare_reduction(name, str(i), x, y) for i, (x, y) in enumerate(zip(a, b)))


def compare_metrics(name, left, right):
    """The remaining numeric record contains a measured cumulative phase frequency."""
    for key in left:
        if isinstance(left[key], (float, int)):
            compare_reduction(name, key, left[key], right[key])
        elif left[key] != right[key]:
            raise AssertionError(f"Metadata mismatch {name}: {key}")


def compare(left, right):
    names = matching_names(left, right)
    max_cell = max_reduction = max_norm = 0
    bitwise_fields = fields = 0
    for name in names:
        a, b = (json.loads((folder / name).read_text()) for folder in (left, right))
        if set(a) != set(b):
            raise AssertionError(f"Schema mismatch: {name}")
        if "values" in a:
            cell, reduction = compare_field(name, a, b)
            max_cell, max_reduction = max(max_cell, cell), max(max_reduction, reduction)
            fields += 1
            bitwise_fields += all(float(x).hex() == float(y).hex()
                                  for x, y in zip(a["values"], b["values"]))
        elif "normalized_residuals" in a:
            max_norm = max(max_norm, compare_residuals(name, a, b))
        else:
            if set(a) != {"omega"}:
                raise AssertionError(f"Unexpected phase schema: {name}")
            compare_metrics(name, a, b)
    print(f"PASS: {len(names)} records; bitwise fields={bitwise_fields}/{fields}; "
          f"max per-cell relative={max_cell:.17g}; reductions={max_reduction:.17g}; "
          f"monitor norms={max_norm:.17g}")


if __name__ == "__main__":
    compare(Path(sys.argv[1]), Path(sys.argv[2]))
