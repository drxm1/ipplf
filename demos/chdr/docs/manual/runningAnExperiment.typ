#import "./preamble.typ": *

#let chapter = [
  #revision[
    #set raw(theme: none)

    = Running an experiment <running-an-experiment>
    The executables `chdr-vacuum-wave`, `planar-interface` and `interface-movies` require `--input FILE`.
    They use the same YAML schema.
    The required keys are `schema_version: 1`, `case` and `units`.
    The optional keys are listed below.
    The following input runs the coarse vacuum case with `chdr-vacuum-wave --input experiment.yaml`.

    ```yaml
    schema_version: 1
    case: coarse
    units: {system: normalized}
    run: {steps: 4}
    output: {directory: vacuum-output}
    observations:
      fields: {every: 2, include_final: true}
    ```

    #table(
      columns: (auto, 1fr),
      inset: 5pt,
      stroke: 0.4pt,
      table.header([*Optional key*], [*Meaning*]),
      [`parameters`], [Case-specific parameters; omitted values retain the case defaults.],
      [`run.steps`], [Nonnegative integer number of advances.],
      [`run.duration`], [Duration in the selected time units; mutually exclusive with `steps`. It must resolve to an integer number of the case's actual time steps. Omitting both retains the case duration.],
      [`output.directory`], [New or empty directory, absolute or relative to the working directory. Required by `interface-movies`; the other programs also support runs without saved output.],
      [`observations`], [Supported measurement names with individual schedule overrides; omitted fields retain their defaults.],
    )

    IPPL handles its own options before the YAML selector.
    Each application's README lists its supported cases and parameters.

    A case supplies the domain, material, solver settings, initial data, measurements and interpretation of their results.
    `ChDRSimulationManager` owns the field storage, material and solver,
    advances the state, and calls the supplied observation function at the initial state and after each advance.
    `runApplication` handles input selection, IPPL initialization and finalization, and fatal errors.
    This interface is implemented for the existing presets.
    It does not expose arbitrary geometry or initial fields through YAML;
    in particular, the `plane` case retains its prescribed duration and Fourier gates.
    The SI scales give the physical interpretation described in @units-si-input-energy.

    == Observations and output <experiment-observations-output>
    Observation indices count completed advances. State zero is the initialized state, and a run with zero advances can still observe it. Each observation accepts `enabled`, `first`, `every`, `last` and `include_final`. The interval from `first` to `last` is inclusive; `last: -1` leaves its upper end open. The cadence starts at `first`, with positive spacing `every`. `include_final: true` additionally selects the final run state if it lies inside this interval. A final state already on the cadence is sampled once. For example, `first: 1`, `every: 3`, `last: 8` and `include_final: true` select states 1, 4, 7 and 8 in a run of eight advances.

    `fields` controls field output independently of the other measurements. Some measurements require every state: vacuum diagnostics must remain enabled from state zero through the final state, and the `plane` Fourier observation must retain `first: 0`, `every: 1` and `last: -1`. Its own gate selects the accumulation interval. Incompatible overrides are rejected. Changing an observation schedule changes when a measurement is called; it does not change the measurement's definition.

    Field frames retain separate time labels for the potentials and the raw electric and magnetic arrays. The series time and `TimeValue` use the potential time; `E_time` and `B_time` identify the array times supplied by the solver. They must be read separately. In SI mode, `FieldOutput` converts coordinates to metres, times and `dt` to seconds, electric fields to volts per metre and magnetic fields to teslas. The solver arrays, ownership indices and dimensionless component offsets are unchanged. The metadata below records the convention. Energy columns follow @units-si-input-energy.

    #table(
      columns: (auto, 1fr),
      inset: 5pt,
      stroke: 0.4pt,
      table.header([*Field metadata*], [*Meaning*]),
      [`step`], [Completed-advance index, not output-frame number.],
      [`TimeValue`, `potential_time`], [Potential time and the time used for the dataset and series.],
      [`E_time`, `B_time`, `dt`], [Separate raw-array times and the actual numerical time step, all in the selected output units.],
      [`normalization`], [`1` for normalized output; `0` for SI output.],
      [`length_scale_from_solver`, `time_scale_from_solver`], [Multipliers applied to coordinates and times.],
      [`E_scale_from_solver`, `B_scale_from_solver`], [Multipliers applied independently to the raw field arrays.],
      [`component_offsets`], [When present, component positions as dimensionless fractions of the mesh spacing.],
    )

    == Run records <experiment-run-records>
    A run with an output directory writes `run.json` after the simulation and output have completed successfully. An existing record is not replaced. Runs without an output directory have no record. The entries describe the supplied input, the settings actually used and the execution environment:

    #table(
      columns: (auto, 1fr),
      inset: 5pt,
      stroke: 0.4pt,
      table.header([*Entry*], [*Content*]),
      [`record_version`, `status`, `application`], [Record format version `1`, status `completed`, and the executable name.],
      [`input`], [Parsed supplied values, including units, parameters and explicit schedule overrides; not the original YAML text.],
      [`effective`], [Normalized box and solver settings, actual time step, stencil, formulation and reconstruction, requested and completed advances, resolved schedules, case metadata and output-unit descriptions.],
      [`output_factors`], [Selected unit system and conversion factors for length, time, fields, potentials, physical field energy, velocity and frequency.],
      [`build`], [Source revision and dirty state captured at build time, compiler, build type, Kokkos execution space and host concurrency.],
      [`execution`], [MPI rank count, each rank's owned index bounds, their ordering, rank-zero host, and the values of `OMP_NUM_THREADS` and `OMP_PROC_BIND` when set.],
    )

    The schedule flags `enabled` and `include_final` use numeric `0` and `1`. Floating values are written with decimal representations that recover the original doubles on reading. Source identity is refreshed at build time; an unavailable Git value is recorded as `unrecorded`. The record identifies the source revision and whether the checkout was dirty, but contains no copy of the uncommitted changes and no hash of every source or output file. A completed record establishes that the simulation and output reached successful completion. It is not a physics-accuracy claim, a convergence result or a conservation proof.

    The input tests verify strict parsing, partial schedule overrides, SI factors and integer-duration conversion. The normalized/SI run comparisons check field data, separate time labels, energy conversion, plane spectra and agreement between the recorded and written time step. These affected checks passed on the CPU build during integration. They verify the interface and conversions for the tested cases; they do not establish new electromagnetic accuracy results or GPU behavior.
  ]
]
