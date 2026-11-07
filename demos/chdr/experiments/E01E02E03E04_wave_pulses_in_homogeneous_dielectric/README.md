# Homogeneous dielectric pulse experiments

The current executable implements E01's longitudinal refinement levels and its four
global error measurements. The reference and experiment protocol are in
[E01.typ](E01.typ). Longitudinal profiles and full-domain field snapshots remain
to be implemented.

From the project root, build and run:

```sh
bash scripts/build-ippl.sh --worktree main build --target experiments-wave-pulses-in-homogeneous-dielectric
OMP_NUM_THREADS=16 OMP_PLACES=cores OMP_PROC_BIND=spread \
  build/ippl-cpu/bin/experiments-wave-pulses-in-homogeneous-dielectric --input ippl/demos/chdr/experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01.yaml
```

The supplied input selects the full E01 target on the coarsest `16 x 16 x 704`
periodic mesh, advancing 1869 steps with timestep `1/32`.
`output.directory` must be new or empty; paths are relative to the process's
working directory. The run writes `errors.csv` and, after successful completion,
the common `run.json` provenance record. Omitting `output` writes the CSV to
standard output instead.

## Refinement selection

Set `parameters.cells_per_pulse` to the number of z-cells across the full pulse
width, as defined in E01.typ. This parameter is required. Each invocation runs one level:

```yaml
parameters: {cells_per_pulse: 64}
output: {directory: build/results/e01-n64}
```

| Cells per pulse | Z-cells | Z-spacing | Timestep | Full-run steps |
|---|---:|---|---|---:|
| 32 | 704 | 1/16 | 1/32 | 1869 |
| 64 | 1408 | 1/32 | 1/64 | 3737 |
| 128 | 2816 | 1/64 | 1/128 | 7473 |
| 256 | 5632 | 1/128 | 1/256 | 14945 |
| 512 (optional) | 11264 | 1/256 | 1/512 | 29890 |

The transverse mesh remains `16 x 16` with spacing `1/16`, and the physical box
remains `1 x 1 x 44`. Material, pulse parameters and target times are fixed.
Choose a new output directory for each run. The CSV and run record contain the
actual mesh, timestep and selected observation indices.

## Short runs

To inspect initialization or make a short run, copy the input, set `run.steps`
and change `output.directory`, for example:

```yaml
run: {steps: 8}
output: {directory: build/results/e01-short}
```

`run.steps: 0` records initialization only. Limits from zero through the selected
level's full-run step count are supported. `run.duration`, SI units, other
application parameters and observation overrides are rejected.

## Error observations

The output contains one row per selected measurement:

| Measurement | Definition in E01.typ |
|---|---|
| `E_L2` | E01.meas01: global relative electric-field L2 error |
| `B_L2` | E01.meas02: global relative R2 magnetic-field L2 error |
| `E_Linf` | E01.meas03: global relative componentwise electric-field maximum error |
| `A_L2` | E01.meas04: global relative vector-potential L2 error |

Each row includes the mesh identifier, cell counts, spacings, timestep,
checkpoint name, target time, completed-step index and actual comparison time.
References are evaluated at the numerical component's position and actual time.
All owned cells and all three vector components contribute; ghost cells do not.

The full run saves initialization and the four targets `T/4`, `T/2`, `3T/4`, `T`
specified in E01.typ. On the coarsest level, potential observations use steps
467, 934, 1401 and 1868; electric/magnetic observations use steps 468, 935, 1402
and 1869. Each refined level selects its nearest available native times using
its actual timestep.
The initial potential time is zero and the initial R2 field time is `-dt/2`.

A short run keeps those same protocol targets and skips any it has not reached.
It additionally measures its last stored state with checkpoint `partial_final`,
unless that quantity group was already measured at a protocol checkpoint there.
For this extra endpoint, `target_time` is the integer potential time; the E/B
comparison time is half a step earlier. Initialization-only runs have no extra
endpoint rows. A complete run has twenty measurement rows.

The first results are single-grid comparisons, not a convergence study.
