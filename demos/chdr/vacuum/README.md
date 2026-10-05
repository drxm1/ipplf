# Periodic Standard vacuum wave

This source-free laboratory-frame benchmark advances a periodic sinusoid for one
full period with IPPL's existing Standard solver. It reports diagnostic CSV on
rank zero and fails with nonzero status when a fixed acceptance check fails.
Optional VTK snapshots preserve the raw E/B arrays and their distinct timestamps;
companion checks read them with actual VTK/ParaView readers.

## Build and run

Use your existing IPPL compiler/MPI/Kokkos configuration. From the IPPL source
root, replace `build` with your build directory:

```sh
cmake -S . -B build -DIPPL_ENABLE_CHDR=ON -DIPPL_ENABLE_UNIT_TESTS=ON \
  -DIPPL_ENABLE_FFT=ON -DIPPL_USE_STANDARD_FOLDERS=ON
cmake --build build --target chdr-vacuum-wave test-chdr-vacuum test-chdr-vacuum-output --parallel 2
OMP_NUM_THREADS=1 OMP_PROC_BIND=false build/bin/chdr-vacuum-wave > coarse.csv
OMP_NUM_THREADS=1 OMP_PROC_BIND=false build/bin/chdr-vacuum-wave --fine > fine.csv
ctest --test-dir build -R '^chdr\.vacuum\.' --no-tests=error --output-on-failure
```

See the [IPPL build resources](../../../README.md#resources) for dependencies.
The existing global unit-test setup requires FFT support and GoogleTest; this
test main uses the flag API available in GoogleTest 1.12 and newer.
CHDR also configures the mesh example's Catalyst/Conduit dependency.
With `IPPL_USE_STANDARD_FOLDERS=OFF`, the executable is under
`build/demos/chdr/vacuum/` instead of `build/bin/`.

| Option | Meaning |
| --- | --- |
| Default | 8 × 8 × 64 cells; 128 steps to newest-potential time 1 |
| `--fine` | 16 × 16 × 128 cells; 256 steps to newest-potential time 1 |
| `--steps N` | Partial run with 0 ≤ N ≤ the chosen grid's full-period step count |
| `--output DIRECTORY` | Save owned E/B fields at quarter periods and the final state |

`--steps 0` measures initialization; `--steps 1` checks one update.
The timestep comes from the solver. There is no adjustable CFL option.

Omitting `--output` keeps the run CSV-only. The output directory must be new or
empty, shared by all MPI ranks, and supplied identically on each rank. For example:

```sh
OMP_NUM_THREADS=1 OMP_PROC_BIND=false build/bin/chdr-vacuum-wave --output wave-run > wave.csv
```

Open `wave-run/wave.pvti.series` in ParaView; [ViewingGuide.org](ViewingGuide.org)
explains array selection, timestamps, camera orientation and the optional headless
rendering script. A full run saves five frames. `--steps 0` saves initialization;
`--steps 1` saves initialization and the final one-step state.

IPPL processes its options first. Its current initialization can exit with status
zero before the benchmark runs, including for help/version and some invalid option
values; require the expected CSV rows as well as a successful exit. Place `--debug`
last because IPPL's option scan skips the following token. See
[`ippl::initialize`](../../../src/Ippl.cpp) for these inherited behaviors.

## Reference, storage and time labels

The box is `[0,1/8] × [0,1/8] × [0,1]`, periodic on all faces and decomposed
along z. Light speed and wavelength are normalized to one. All components are
cell-centered. Our reference travels along +z with electric polarization x and
magnetic polarization y: the four-potential is
`(0, sin(2*pi*(z-time))/(2*pi), 0, 0)`, and continuum Ex/By are both
`cos(2*pi*(z-time))`. Other components and the source are zero.

This is our substitution into Fallahi,
[*MITHRA 2.0*, arXiv:2009.13645v1](https://arxiv.org/abs/2009.13645v1),
“Wave Equation,” (3.6)–(3.9). The Standard stencil/dispersion starting relations
are “FDTD for Wave Equation” and “Numerical Dispersion in FDTD,” (3.11)–(3.20).

We sample the exact wave at 0 and -dt to supply the two histories required by
the three-level update. The [Standard solver](../../../src/MaxwellSolvers/StandardFDTDSolver.hpp)
chooses `dt = min(spacing)/2`: 1/128 coarse and 1/256 fine.

| After n completed solves | Physical reference time |
| --- | --- |
| Newest potential, raw B | `n*dt` |
| Previous potential | `(n-1)*dt` |
| Raw E | `(n-1/2)*dt` |

The initial raw E time is therefore negative. These labels follow
[`FDTDSolverBase::solve/timeShift/evaluate_EB`](../../../src/MaxwellSolvers/FDTDSolverBase.hpp):
for this zero-scalar-potential case E differences two histories, while B curls
the newest history. They use MITHRA's continuum identities (3.8)–(3.9);
IPPL does not perform the paper's magnetic-field time averaging.

## Diagnostics and fixed acceptance

Measurements cover owned cells only. A local Kokkos reduction combines error
squares, maxima, Fourier sums and counts; MPI combines scalar summaries.
No diagnostic changes the fields or copies a full field to the host.
Continuum relative L2 error is the square root of the global error-square sum
divided by the global reference-square sum; a nonpositive/nonfinite denominator
is a failure. Reference-zero components use absolute errors.

The discrete reference is independent of the initializer and solver reconstruction.
Our axial reduction of MITHRA (3.11)–(3.17) gives
`u[n+1] = 2*cos(Omega)*u[n] - u[n-1]`, with `u[0]=1` and
`u[-1]=exp(i*k*dt)`, where `k=2*pi` and
`Omega=2*asin((dt/hz)*sin(k*hz/2))`. The closed-form oracle keeps both temporal
branches implied by those histories. Comparison to the continuum remains a
separate accuracy check.

For direction and speed we measure the normalized coefficient
`C=(2*k/M)*sum(Ax*exp(-i*k*z))` over M owned cells. For the continuum reference
`C=-i*exp(-i*k*t)`, so +z propagation decreases phase. Successive principal
phase increments are accumulated; mean speed is `-sum(increments)/(k*t)`.
This is our Fourier substitution in MITHRA (3.6)–(3.9), with its discrete
dispersion relation (3.17) supplying the comparison.

| Check | Fixed requirement |
| --- | --- |
| Continuum Ex/By relative L2 at quarter periods and final state | ≤0.01 coarse; ≤0.003 fine |
| Ratio of coarse/fine maximum checkpoint errors | Between 3 and 5, for each field |
| Potential/E/B discrete-reference maximum error, every step | ≤1e-10; potentials scaled by amplitude 1/k |
| Inactive components and zero source, every step | ≤1e-12; all owned values must be finite |
| Fourier amplitude / successive phase increment | >0.5 / strictly between -pi and 0 |
| Cumulative phase speed, after step zero | Within 0.002 of one |
| Serial/MPI2 global diagnostics | Agree within 1e-11, with exact owned-cell counts |

These are benchmark requirements, not tunable solver parameters or general
accuracy guarantees. The driver checks each run; the Python benchmark also
checks refinement and cross-run agreement. The actual-reader check additionally
compares every saved E/B component with its reference and compares serial/MPI2
fields by global cell index within the fixed 1e-11 tolerance.

## Field files and readback

Each checkpoint writes one `.vti` piece per rank and one rank-zero `.pvti` wrapper,
including for a serial run. `wave.pvti.series` lists the completed wrappers and
their explicit potential times. `E_raw` and `B_raw` are three-component Cell Data;
halos are excluded and tuples are serialized with x varying fastest. Geometry
uses the fields' actual global origin, spacing and owned extents.

Field Data records `potential_time`, `E_time`, `B_time`, `dt`, `step` and
`normalization`. Code 1 denotes the dimensionless convention above, with no SI
scale assigned. ASCII doubles use round-trip precision and the classic locale.
These choices follow [VTK XML formats](https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html),
[VTK time metadata](https://docs.vtk.org/en/latest/design_documents/IOXMLTimeInFieldData.html)
and [ParaView temporal series](https://docs.paraview.org/en/latest/UsersGuide/dataIngestion.html#handling-temporal-file-series).

Writes are synchronous. Ordinary rank-local I/O failures are communicated before
the next publication stage: all pieces finish before a wrapper is published,
and only completed frames enter the series. On failure, temporary or unlisted
files can remain; the previous complete series is retained. Use a fresh directory
for a new run. This does not provide recovery from process or MPI failure.

## CSV and tests

There is one row for every measured state, including step zero:

| Columns | Meaning |
| --- | --- |
| `step,dt,timeA,timeE,timeB` | Completed solves, actual step size and distinct array times |
| `relativeE,relativeB` | Continuum Ex/By relative L2 errors |
| `discreteError,zeroError,sourceError` | Normalized maximum discrete, inactive/zero and source errors |
| `phaseReal,phaseImag,phaseAmplitude` | Normalized potential Fourier coefficient and magnitude |
| `phaseIncrement,phaseSpeed` | Signed increment in radians and cumulative speed |
| `count` | Global owned-cell count: 4096 coarse, 32768 fine |

Initial `phaseSpeed` is empty because no time has elapsed. Rank zero attempts to
print each sample before its own validation; another rank's MPI abort can prevent
the failing row from being written or delivered. A CSV alone does not imply success.

The 29 GoogleTests cover first-step/setup behavior, synthetic diagnostic faults
(including non-root NaNs and poisoned halos), independent scalar recurrence,
phase tracking and acceptance gates. They run serially, on MPI2 and shuffled MPI2.
Unseeded shuffling uses seed1; explicit nonzero seeds are preserved. All ranks
must receive matching test-selection/order options. Halo guards remain debug-only;
the test target explicitly enables them.

Two driver smoke tests use `--steps 1`. The propagation test runs coarse serial,
fine serial and coarse MPI2 through a full period, checking every CSV row with
an independent scalar recurrence. The CLI test covers bounded/invalid arguments.
Python3 and an MPI launcher are needed for these integration checks. Actual-reader
tests additionally need `pvpython`; CMake discovers it or accepts
`-DCHDR_VACUUM_PVPYTHON=/path/to/pvpython`. The driver itself does not link VTK or
ParaView. Unavailable checks are explicitly skipped and do not count as acceptance.

Six output tests exercise a distinguishable 4×3×6 fixture and four coordinated
failure cases. The actual reader checks every component, coordinates, halos and
metadata; the full-wave reader checks fifteen frames and serial/MPI field agreement.
Short output-enabled runs also exercise initialization-only and one-step readback.

CTest logs, CSV and JSON manifests default to the test build directory's `runs/`;
set `CHDR_VACUUM_TEST_OUTPUT_DIR` to relocate them. The wrapper scripts use
`--results` for this record directory; it is not a driver field-output option.

## Source reading order

Start with [the driver](vacuum-wave.cpp) and [shared setup](VacuumSetup.cpp).
[VacuumDiagnostics.h](VacuumDiagnostics.h) defines measurement and phase state;
[VacuumChecks.h](VacuumChecks.h) applies the fixed policy.
[VacuumOutput.h](VacuumOutput.h) handles synchronous collective file output. The
[field container](../../fel/FELFieldContainer.hpp) owns mesh/layout and E/B/J
and must outlive the solver. Tests and source comments document the independent
references and numerical assumptions.
