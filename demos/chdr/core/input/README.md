# Common experiment input and run records

The field applications retain their executable names. Each accepts `--input FILE`
for the same strict YAML schema. Historical command-line arguments are translated
into the same `RunInput` structure before setup and execution. YAML input cannot
be mixed with case or output overrides on the command line.

`Application.h` handles this selection, IPPL initialization/finalization and fatal
error reporting. `ChDRSimulationManager` owns and advances the simulation. The
application supplies its material, initial data, measurements and interpretation
of their results. Measurements use the manager's existing callback; a shared
`ObservationSchedule` defines when each measurement is called.

## Input schema

For example, the existing coarse vacuum-wave case can run for four steps:

```yaml
schema_version: 1
case: coarse
units:
  system: normalized
run:
  steps: 4
output:
  directory: vacuum-output
observations:
  fields:
    every: 2
    include_final: true
```

The required keys are `schema_version`, `case` and `units`. Version 1 is supported.
The optional keys are:

| Key | Meaning |
|---|---|
| `parameters` | A mapping interpreted and validated by the selected application. |
| `run.steps` | Nonnegative integer number of advances. Omission retains the case default. |
| `run.duration` | Duration in the selected time units; mutually exclusive with `steps`. |
| `output.directory` | Nonempty output path relative to the process working directory, or an absolute path. |
| `observations` | Mapping from supported measurement names to their schedule overrides. |

Unknown keys are errors. Each application also rejects unsupported case names,
parameter names and measurement names, including disabled measurements. A
parameter must have its expected type; a fractional or overflowing step count
is rejected. Boolean options accept `true`/`false` and `0`/`1`; Conduit also
represents quoted `"true"` and `"false"` as the same text values.
The schema supplies a common route for existing cases, not arbitrary
changes to a case's physics or numerical contract. See each application's README
for its supported parameters and restrictions.

Omitting `output` retains stdout-only operation where the application supports it.
A case that requires output can reject that omission. Output requires a new or empty
directory, checked before advancing the simulation. A successful recorded run writes `run.json` and
will not replace an existing record.

## Observation schedules

Each observation accepts `enabled`, `first`, `every`, `last` and `include_final`.
Omitted fields retain that application's defaults independently. The standardized
field-output observation name is `fields`; other names reflect the actual measured
quantities and can differ between applications.

Indices count completed solver advances. State zero is the initialized state.
`first` and `last` are inclusive; `last: -1` leaves the upper end open. `every`
must be positive, and its cadence starts at `first`. `include_final: true`
additionally selects the final run state if it lies inside the interval. That
state is not sampled twice when it already lies on the cadence.

Some existing measurements require every state. For example, a Fourier
accumulator requires consecutive samples inside its gate. Such a case can reject
incompatible schedule overrides. Scheduling does not change a measurement's
definition or its field time labels, and does not add a new monitor algorithm.

## Normalized values and SI

`units: {system: normalized}` preserves the existing code-unit values without
assigning physical units to a preset. SI input/output requires both positive
normalization scales:

```yaml
units:
  system: si
  length_m: 0.001
  charge_density_C_m3: 8.8541878188e-9
```

Here `length_m` is the SI length corresponding to one normalized length and
`charge_density_C_m3` is the SI charge density corresponding to one normalized
charge density. The preset's normalized geometry and numerical time step remain
unchanged; these scales specify their physical interpretation. The current presets
do not expose dimensional geometry or initial-field parameters in YAML.
`run.duration` is in seconds for SI input and in code time for normalized input.
The duration must equal an integer number of the case's actual time steps within
floating-point roundoff. The reader neither truncates the duration nor adjusts
the time step to fit it.

With scales \(L_0\) and \(\rho_0\), the conversion factors are

\[
t_0=\frac{L_0}{c},\qquad
E_0=\frac{\rho_0L_0}{\varepsilon_0},\qquad
B_0=\frac{E_0}{c},\qquad
\phi_0=E_0L_0,\qquad A_0=B_0L_0.
\]

These follow by substituting scaled variables into SI Maxwell's equations and
the potential definitions [fitzpatrick2014ClassicalElectromagnetism,
Eqs. (1.1), (1.3), (1.11)–(1.16); deduction]. The velocity and frequency factors
are \(c\) and \(1/t_0\). Physical field energy in the linear nondispersive,
nonmagnetic model has the factor
\(U_0=\varepsilon_0E_0^2L_0^3\), obtained from the electric and magnetic
energy density [same source, Eq. (1.109), p. 20; dielectric electric-energy
weighting: Eq. (4.70); deduction]. The implementation
uses the published CODATA-2022 values of \(c\) and \(\varepsilon_0\)
[mohr2025Codata2022, Table XXXIII, p. 025002-45].

In SI mode the output adapters convert physical coordinates, time and fields.
Each application's quantity labels specify which diagnostic columns are physical
SI quantities. A potential invariant is not relabelled as physical energy merely
because it is quadratic. Existing dimensionless checks remain in code units;
conversion takes place at input/output boundaries.

## Records and extension points

`core/output/RunRecord.h` writes the common record after a successful run and
output finalization. Its main entries are:

| Entry | Content |
|---|---|
| `input` | Supplied common options, explicit unit selection and parameter overrides. |
| `effective` | Actual box, solver time step, reconstruction, formulation, completed advances and resolved observation schedules; numerical values are labelled normalized. |
| `output_factors` | Standard conversion factors; output descriptions identify which are used. |
| `build` | Source revision/dirty state at build time, compiler, build type and Kokkos backend. |
| `execution` | MPI rank count, actual owned index ranges, rank-zero host and selected thread settings. |

Case metadata and output-column descriptions supplement these shared entries.
Boolean entries in the Conduit JSON record use numeric `0` and `1`.
Floating values use JSON decimals that round-trip to the original double values.
Source identity is refreshed at each build; a failed Git query is labelled
`unrecorded`. It is not a hash of every source or output file. `status: completed` means the
application returned successfully; it is not a physics-accuracy claim.

To connect an existing case, translate its arguments to `RunInput`, validate its
supported names, resolve its duration after selecting the actual time step, and
attach its existing measurements to the manager callback. Call each measurement
only when its resolved schedule selects the state, finalize its output, then
record the effective settings. Core code does not depend on experiment-specific
or optional diagnostic implementations.

The fast `chdr.input.unit.serial` test covers strict parsing, partial schedule
overrides, SI factors and integer-duration conversion. It does not establish the
accuracy of an electromagnetic experiment.
