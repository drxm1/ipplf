# Homogeneous dielectric pulse experiments

The current executable implements E01's longitudinal refinement levels,
four global error measurements, matching longitudinal profiles and full-domain
snapshots of E, R2-averaged B and raw B.
The reference and experiment protocol are in [E01.typ](E01.typ).

From the project root, build and run:

```sh
bash scripts/build-ippl.sh --worktree main build --target experiments-wave-pulses-in-homogeneous-dielectric
OMP_NUM_THREADS=16 OMP_PLACES=cores OMP_PROC_BIND=spread \
  build/ippl-cpu/bin/experiments-wave-pulses-in-homogeneous-dielectric --input ippl/demos/chdr/experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01.yaml
```

The supplied input selects the full E01 target on the coarsest `16 x 16 x 704`
periodic mesh, advancing 1869 steps with timestep `1/32`.
`output.directory` must be new or empty; paths are relative to the process's
working directory. The run writes `errors.csv`, `profiles.csv`, field frames under
`fields/` and,
after successful completion, the common `run.json` provenance record.
Omitting `output` writes only the error CSV to standard output.

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
application parameters and observation names other than `fields` are rejected.

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

## Longitudinal profiles

`profiles.csv` contains numerical and analytic `Ax`, `Ex` and R2 `By`
at their respective error observations, including initialization and any
`partial_final` endpoint. Each row records the component, checkpoint,
target time, completed-step index, actual time, global z-cell index `k`,
native coordinates `x,y,z` and the two values.

All refinements sample the fixed transverse cell indices `(i,j) = (7,7)`.
For these components, this gives the same physical line `x = 1/2, y = 15/32`.
Use the recorded `z` coordinates for plotting: `By` is half a longitudinal
cell beyond `Ax` and `Ex`, following the staggered placement in E01.typ.
Its last owned sample is at `z = 44` and is retained; ghost cells are excluded.
Only the selected line is transferred to rank zero for writing.

A complete run contains five observations per component, giving `15 * Nz`
profile rows (10,560 on the coarsest mesh), in addition to the CSV header.

## Full-domain snapshots

Open `fields/fields.pvti.series` in ParaView to load the full field series.
Each selected solver state contains three vector arrays: `E`, `B_R2` and
`B_raw`. The first two are the fields used in the error measurements;
raw B is reconstructed into separate storage only when saving a frame.
All ranks write their owned cells through the shared `core/output/FieldOutput.h`.

The default interval is the nearest integer step to the first quarter target,
as specified in E01.typ: 467, 934, 1868 and 3736 steps on the four main meshes.
Snapshots include initialization and the final completed state,
without duplicating a frame when the final state is already on the cadence.
Thus the full coarsest run saves steps 0, 467, 934, 1401, 1868 and 1869.
A short run retains the same interval and includes its own final state.

Use the common observation options to override this schedule, for example:

```yaml
observations:
  fields:
    every: 100
    include_final: true
```

`first`, `last` and `enabled` are supported as described in
[the common input documentation](../../core/input/README.md#observation-schedules).
`enabled: false` disables snapshots while keeping errors and profiles.
By default, snapshots are enabled when `output.directory` is supplied;
explicitly enabling them requires that directory.
These options do not change the error or longitudinal-profile observations.
The resolved field schedule is recorded in `run.json`.

Each array has its own `<name>_time` and `<name>_component_offsets` metadata,
derived from the solver's `LevelContract`. The offsets form three tuples,
one per Cartesian component, in units of mesh spacing from the cell centre.
`<name>_scale_from_solver` records the value conversion factor (one for E01).
The series uses the potential time; E and R2 B carry their half-step time,
while raw B carries the integer time. Preserve these individual time labels
when using the snapshots for comparisons; see the reconstruction convention in E01.typ.
