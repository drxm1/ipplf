# Tests of the ChDR code

This directory holds the tests of the ChDR code, one directory per unit of the
code they test: `core/<unit>`, `opt-in/diagnostics`, `apps/<program>`, `1d`,
and `src/<component>` for the IPPL library changes made for ChDR. The tests of
the mesh example live in `../mesh/`. Shared test support is in `../testing/`.

This file lists every CTest entry: what it checks, against which reference,
from which source file, and how long it ran. The tables describe the code at
the time of writing; when a test changes, its row changes with it.

## Running the tests

Configure IPPL with `-DIPPL_ENABLE_CHDR=ON -DIPPL_ENABLE_FEL=ON
-DIPPL_ENABLE_FFT=ON -DIPPL_ENABLE_SOLVERS=ON -DBUILD_TESTING=ON
-DIPPL_ENABLE_UNIT_TESTS=ON` (see [`../README.md`](../README.md)). Then:

```sh
ctest --test-dir build -R '^chdr\.' -LE 'long|record' --output-on-failure   # fast set
ctest --test-dir build -R '^chdr\.' -L long --output-on-failure             # ladders and long runs
ctest --test-dir build -R '^chdr\.' -L record -V                            # reference functions and records
```

The entries are named `chdr.<unit>.<suite>[.<variant>]`. The variant `serial`
runs one process; `mpi2` runs two MPI ranks with one OpenMP thread each (the
mesh entries `chdr.mesh.mpi2` and `chdr.mesh.brick.mpi2` use two threads per
rank and declare four processors).
CTest counts the declared processors of each entry against `--parallel`.

The current registration has 211 entries, including two disabled former
holdouts; 164 enabled entries belong to the fast selection above. The timing
columns retain the baseline run of 2026-11-03, before the registration cleanup:
222 of 226 entries ran and passed in 682 s with `--parallel 16` on a 16-core
AMD Ryzen 9 9950X3D (RelWithDebInfo, GCC 15.3, Open MPI 5.0.10, Kokkos 5.2.2).
The former fast set had 176 enabled entries and ran in 57 s. These are
historical measurements, not timings of the current selection. A pair `a / b`
gives the serial and `mpi2` variant. Entry times include contention from the
other tests running at the same time and will grow on slower machines.

The `long` label describes work: roughly \f$10^8\f$ cell-steps (cells times
solves, summed over all runs of an entry), or long trajectories of 10000 steps
or ten wave periods. Existing classifications are retained rather than changed
to match timings on this machine. The `record` label selects record generation
and checks of reference functions that do not evolve the solver; these entries
can still check finiteness or the reference functions. The `holdout` and
`reserved` labels identify the two disabled former acceptance cases.
`unavailable` identifies skipped placeholders when a required interpreter,
MPI launcher or reader is absent. The `mpi` label identifies entries that
launch multiple MPI ranks, including launches inside Python scripts; a
one-process comparator does not carry it. A grouped row lists the union of
its variants' labels. The timeout is 180 s for `medium.Convergence.*` and
`mur.reflection.*`.

## How to read the tables

"Reference" names what the test compares the code with:

| Reference | Meaning |
|---|---|
| exact lattice value | the exact solution of the discrete equations, computed independently of the solver; a correct code meets it to round-off (typically `1e-12` or `1e-13`) |
| continuum analytic | a closed-form solution of the continuous equations; the tolerance is a discretization error or an observed order of convergence |
| self-consistency | an identity, an invariant, a contract, or the equality of two code paths; "refusal" marks checks that a wrong input throws |
| one against two ranks | the same run on one rank and on two ranks |
| regression record | exact numerical or bitwise equality with a frozen code path, or a record that the test only prints |
| smoke | the program runs to the end and writes the expected files |

Most evolving solver tests also require, after their last solve, that the
normalized residuals of the divergence of the output B, the divergence of the
raw B, the three Faraday components and the Gauss formula are each at most
`1e-12` (`requireFinalIdentities` in `core/solver/MonitorTestSupport.h`). A
normalized residual is the residual of a cell divided by its uncancelled
absolute scale, maximized over the cells. The tables call these checks the final
identities. "Default step" means the inherited
step: \f$\Delta t=h/2\f$ on cubic Standard cells and \f$\Delta t=h_z\f$ for
NonStandard. R1, R2 and R2Prime are the E and B reconstruction conventions of
`ChdrSolver`; "collocated" and "staggered" are its two spatial readings.

Two entries are disabled: `chdr.solver.mur.holdout` and
`chdr.interface.reserved`, the former held-out final acceptance checks of
earlier work. Four further disabled GoogleTests, also former holdouts, have
no CTest entry:
`DISABLED_Holdout.Collocated`, `DISABLED_Holdout.StaggeredRecord`,
`MediumHoldout.DISABLED_Fields` and `MediumHoldout.DISABLED_Phase`.

## Overview

| Directory | Entries | Test code (files, lines) | Time sum (s), all / fast set |
|---|---|---|---|
| `src/FieldLayout` | 2 | 1, 188 | 3 / 3 |
| `src/MaxwellSolvers` | 12 | 14, 1726 | 30 / 30 |
| `core/material` | 6 | 5, 421 | 10 / 10 |
| `core/simulation` | 2 | 1, 124 | 3 / 3 |
| `core/input` | 2 | 2, 352 | 14 / 14 |
| `core/solver` | 146 | 55, 8394 | 5251 / 458 |
| `core/solver/formulation` | 2 | 1, 167 | 4 / 4 |
| `opt-in/diagnostics` | 11 | 7, 1373 | 23 / 23 |
| `apps/vacuum-wave` | 12 | 13, 1917 | 74 / 74 |
| `apps/planar-interface` | 5 | 2, 372 | 364 / 8 |
| `apps/interface-movies` | 3 | 1, 119 | 17 / 17 |
| `../mesh` | 7 | 2, 714 | 7 / 7 |
| `1d` | 1 | 1 (in `../1d/`), 199 | 1 / 1 |

"Test code" counts the C++ and Python files of the directory itself (CMake
files excluded). These counts, and the line counts in the "Source (lines)"
columns below, refer to commit `ade845ffe`, before the comment pass that
followed this index. The time sums add the baseline wall times of the retained
entries, grouped by the current labels; they are not a new run or the time of
a parallel run.

## `src/FieldLayout`: per-axis periodic layouts

`PeriodicAxes.cpp` tests IPPL field layouts that carry one periodicity flag per
axis. CMake enters `tests/src` only when IPPL is the top-level project. Run it
with `ctest --test-dir build -R '^chdr\.src\.PeriodicAxes\.' --output-on-failure`.

| Entry | Source (lines) · suite | What it checks | Reference | Labels | Time (s) |
|---|---|---|---|---|---|
| `chdr.src.PeriodicAxes.{serial,mpi2}` | `PeriodicAxes.cpp` (188) · `PeriodicAxes` (5 tests) | On a 6 by 5 by 4 box, for all 8 periodicity masks and 3 split axes, `fillHalo()` gives every allocated cell the label that the oracle `haloSourceIndex` predicts (ghosts beyond an open face keep the sentinel −1), and `accumulateHalo()` gives every owned cell the predicted deposit sum; the multigrid preconditioner refuses a mixed layout and accepts uniform ones. | exact lattice value; refusal | `unit`, `mpi` | 1.6 / 1.6 |

## `src/MaxwellSolvers`: Mur faces of the FDTD potential solvers

These tests check the IPPL changes to `StandardFDTDSolver` and
`NonStandardFDTDSolver` for open faces: the face description `FDTDBoundary`
and its validation, the cells that the interior stencil writes, the Mur
weights at a general speed and their dispatch to faces, edges and corners, the
halo exchange after the Mur update, and the potential ghosts beyond absorbing
faces. `BoundaryOracles.h` and `OuterLayerOracles.h` hold expected values
written independently of the solvers; `ShippedOrderSolver.h` restores the
original finishing order (halo exchange, then Mur) as a negative control.
[`HaloOrder.md`](src/MaxwellSolvers/HaloOrder.md) and
[`OpenFaceClosure.md`](src/MaxwellSolvers/OpenFaceClosure.md) describe two
suites in more detail. Run the tests with
`ctest --test-dir build -R '^chdr\.src\.' --output-on-failure`.

| Entry | Source (lines) · suite | What it checks | Reference | Labels | Time (s) |
|---|---|---|---|---|---|
| `chdr.src.HaloOrder.mpi2` | `HaloOrder.cpp` (70), `HaloOrderMixed.cpp` (97), `HaloOrderPulse.cpp` (97) · `HaloOrder.*` (6 tests) | The affine drift checks on one and two ranks: the scaled error is at most `1e-12` over 64 steps, and the old-order control misses in the mixed box. A Gaussian pulse over 256 steps crosses the rank seam and matches a one-rank solver to `1e-13` of its maximum; the old-order control differs by at least `1e-3`. | exact lattice value; one against two ranks | `unit`, `mpi` | 5.0 |
| `chdr.src.OpenBoundaryRegression.{serial,mpi2}` | `OpenBoundaryRegression.cpp` (143) · `OriginalOrder` | Over 256 steps from fixed pseudo-random histories, the potentials, E, raw B and averaged B are bitwise equal to those of the old finishing order with the same kernels, in periodic and all-absorbing Standard and NonStandard boxes; the absorbing boxes exclude the outermost owned E and B layer and run on `MPI_COMM_SELF` in both entries, so only the periodic cases use two ranks. | regression record (bitwise) | `unit`, `mpi` | 5.7 / 5.3 |
| `chdr.src.OpenFaceClosure.{serial,mpi2}` | `OpenFaceClosure.cpp` (130) · `OpenFaceClosure` (4 tests) | After `prepareHistories()`, E and raw B on the outermost layer match a 40-coefficient polynomial oracle to `1e-13` of the analytic maximum, zero ghosts miss by at least 0.1, and one step of either stencil closes the ghosts of the new level; `mpi2` also compares two ranks with a one-rank replay. | exact lattice value; one against two ranks | `unit`, `mpi` | 1.6 / 1.6 |
| `chdr.src.FDTDBoundary.{serial,mpi2}` | `FDTDBoundary.cpp` (120) · `FDTDBoundary` (4 tests) | `fromLayout` marks exactly the faces of open axes as Mur faces of speed 1 for all 8 axis patterns, and `validateBoundary` refuses invalid speeds, mismatched face kinds, wrong ghost widths, too few cells and (two ranks) a rank with one cell. | refusal; exact lattice value (cell sets) | `unit`, `mpi` | 1.6 / 1.6 |
| `chdr.src.MurSpeed.serial` | `MurSpeed.cpp` (194) · `MurSpeed` (2 tests) | The face, edge and corner Mur weights at speeds 0.6, 1 and 1.2 match the closed forms to \f$4\epsilon\f$ times the largest weight, and at speed 1 their bits equal the original unit-speed arithmetic. | exact lattice value; regression record (bitwise) | `unit` | 1.6 |
| `chdr.src.InteriorBounds.{serial,mpi2}` | `InteriorBounds.cpp` (163) · `InteriorBounds` (5 tests, 2 compile-time checks) | One step of either stencil without the boundary update writes \f$\Delta t^2\f$ into exactly the interior cells of the oracle mask for all 8 axis patterns; the solver constructors refuse invalid descriptions before any communication. | exact lattice value; refusal | `unit`, `mpi` | 1.6 / 1.6 |
| `chdr.src.MurDispatch.{serial,mpi2}` | `MurDispatch.cpp` (139) · `EveryOpenFaceEdgeAndCorner` | After `applyBCs()`, every face, edge and corner cell equals the output of the Mur struct that its mask selects and interior cells keep their values; an all-open dispatch must fail on the mixed patterns. | self-consistency (routing against an independent mask) | `unit`, `mpi` | 1.6 / 1.6 |

## `core/material`: material owner

`test-chdr-material` tests `MaterialOwner`: the rasterization of a permittivity
map onto cells, the exact global statistics, the ghost cells at rank seams and
at periodic or open faces, and the cached coefficient fields of the O1-VS
formulation. `OwnerOracles.h` gives the expected value of every allocated cell
independently of the owner. The rasterization, halo and complete-map cases
repeat with each axis as the split axis; the two tests of homogeneous values
use the split along x only. Run the tests with
`ctest --test-dir build -R '^chdr\.material\.' --output-on-failure`.

| Entry | Source (lines) · suite | What it checks | Reference | Labels | Time (s) |
|---|---|---|---|---|---|
| `chdr.material.Owner.{serial,mpi2}` | `test_owner.cpp` (104) · `MaterialOwner` (2 tests) | On a 12 by 14 by 16 grid, every allocated permittivity cell is bitwise equal to the map at its clamped source cell, and count, checksums, minimum and maximum equal a serial host scan; a NaN map makes `homogeneousValue()` throw. | exact lattice value; one against two ranks | `unit`, `mpi` | 1.6 / 1.6 |
| `chdr.material.Halo.{serial,mpi2}` | `test_halo.cpp` (67) · `MaterialHalo` (3 tests) | With a distinct integer per cell, no allocated cell (faces, edges, corners, rank seams) differs from its source cell, for periodic, open and the 6 mixed axis patterns. | exact lattice value | `unit`, `mpi` | 1.6 / 1.6 |
| `chdr.material.Derived.{serial,mpi2}` | `test_derived.cpp` (73) · `MaterialDerived` (2 tests) | Every allocated cell of \f$Q=\varepsilon^{-2}\f$ and of the three harmonic face fields is bitwise equal to its closed form; the cache returns the same storage and refuses use before `derive`. | exact lattice value; self-consistency (cache, refusal) | `unit`, `mpi` | 1.6 / 1.6 |

## `core/simulation`: simulation manager

`test-chdr-simulation` tests `ChDRSimulationManager`, the class that builds and
advances every ChDR run. The expected call lists follow from the contract of
`ippl::BaseManager::run`; no field value is compared. Run the tests with
`ctest --test-dir build -R '^chdr\.simulation\.' --output-on-failure`.

| Entry | Source (lines) · suite | What it checks | Reference | Labels | Time (s) |
|---|---|---|---|---|---|
| `chdr.simulation.unit.{serial,mpi2}` | `test_simulation.cpp` (124) · `Simulation` (5 tests) | The manager records the expected observer states and solve counts for `advance` and `run`, keeps the open-periodic-open layout and the safety factor, selects the vacuum update without material, uses an explicit step, keeps independent observation schedules, and refuses invalid steps and (two ranks only) a rank with one cell. | self-consistency (call contract, refusal) | `unit`, `mpi` | 1.7 / 1.6 |

## `core/input`: common input and SI units

`test-chdr-input` tests the common YAML input of the three programs: strict
keys and types, partial overrides, the SI scale factors and the conversion of
an SI duration to steps. `check_runs.py` runs `chdr-vacuum-wave`,
`interface-movies` and `planar-interface` once in code units and once in SI
units and compares every published number; without Python 3, CMake keeps
the entry as an `unavailable` placeholder. Run both with
`ctest --test-dir build -R '^chdr\.input\.' --output-on-failure`.

| Entry | Source (lines) · suite | What it checks | Reference | Labels | Time (s) |
|---|---|---|---|---|---|
| `chdr.input.unit.serial` | `test_input.cpp` (152) · `RunInput` (2 tests) | The reader refuses unknown keys, malformed step counts and schedules, invalid units, mixed selectors and broken YAML; the SI factors for 1 mm equal fixed CODATA 2022 values, and an SI duration converts to an exact number of steps. | refusal; closed-form SI factors | `unit` | 1.6 |
| `chdr.input.runs.si` | `check_runs.py` (200) | Each program runs one case in code units and again for the equal SI duration (`chdr-vacuum-wave` and `interface-movies` 4 steps; `planar-interface` its complete pulse of 1180 steps); every SI number (record factors, VTK origin, spacing, times, fields, CSV columns, spectra) equals the code-unit number times its factor to relative `2e-12`. | self-consistency (SI run against code-unit run) | `integration` | 12.1 |

## `core/solver`: the potential solver

`test-chdr-solver` is one executable for all suites of this directory (146
entries). The sections below group the entries by subject. Run a subject with
the pattern given in its section, or all fast entries with
`ctest --test-dir build -R '^chdr\.solver\.' -LE 'long|record' --output-on-failure`.

### Construction, operators and periodic halos

These entries test the construction contract of `ChdrSolver`, the discrete
operators in `core/solver/operators/` and the periodic box that
`makeBoxFieldContainer` builds. The name `compact` is the earlier name of the
staggered difference. Run them with
`ctest --test-dir build -R '^chdr\.solver\.(contract|compact|staggered|setup|operators)\.' --output-on-failure`.

| Entry | Source (lines) · suite | What it checks | Reference | Labels | Time (s) |
|---|---|---|---|---|---|
| `chdr.solver.contract.{serial,mpi2}` | `test_solver_contract.cpp` (173) · `ChdrSolverContract` (7 tests) | With R1, `ChdrSolver` gives E and B exactly equal to IPPL's Standard solver (maximum difference zero) after each of 128 steps of a vacuum plane wave; it declares the label times, honours the selected reconstruction and refuses a halo width of 2, rebinding, distinct layout owners and (two ranks) a rank with one cell. | regression record (exact numerical equality); refusal | `unit`, `mpi` | 1.7 / 1.6 |
| `chdr.solver.compact.serial` | `test_compact.cpp` (72) · `CompactOperators` | The forward and backward `staggeredDifference` of one sampled sine mode equal the signed half-cell symbol on each axis to `1e-12`. | exact lattice value | `unit` | 1.6 |
| `chdr.solver.staggered.{serial,mpi2}` | `test_staggered.cpp` (345) · `ChdrStaggered` (7 tests) | In the staggered reading with R2, the two-pass Laplacian equals the 7-point Laplacian on rough data, the discrete identities and the continuity monitor hold to `1e-12`, discrete gauge waves give E and B at most `1e-12` of the gradient scale, and R1, R2Prime and NonStandard are refused. | exact lattice value; self-consistency (identities, refusal) | `unit`, `mpi` | 2.2 / 2.0 |
| `chdr.solver.setup.{serial,mpi2}` | `test_setup.cpp` (95) · `BoxSetup` (2 tests) | After `preparePotentialHalo`, every allocated cell of the four-potential (faces, edges, corners) equals the periodic copy of its integer identifier, for each split axis; invalid geometries throw. | exact lattice value; refusal | `unit`, `mpi` | 1.6 / 1.6 |
| `chdr.solver.operators.{serial,mpi2}` | `test_operators.cpp` (111) · `Operators` | The Standard and NonStandard `laplacian` and the `centeredDifference` of one mode on unequal spacings equal an independently coded Fourier symbol to `1e-12`. | exact lattice value | `unit`, `mpi` | 1.6 / 1.6 |

### Time step

These entries test the time-step setters of `ChdrSolver` against the stability
bound, and the amplitude of random potential histories below and above it. Run
them with
`ctest --test-dir build -R '^chdr\.solver\.(long\.)?timestep\.' --output-on-failure`.

| Entry | Source (lines) · suite | What it checks | Reference | Labels | Time (s) |
|---|---|---|---|---|---|
| `chdr.solver.timestep.{serial,mpi2}` | `test_timestep.cpp` (294) · `TimeStep` (3 tests) | Safety factors 0.2, 0.9 and 0.98 give \f$\Delta t=S\,\Delta t_{\max}\f$ to `1e-15`, invalid factors and step changes after a solve throw, and at \f$1.02\,\Delta t_{\max}\f$ random histories grow at least 1000-fold within 60 solves. | exact lattice value (bound); refusal | `unit`, `mpi` | 1.7 / 1.7 |
| `chdr.solver.long.timestep.{StandardBounded,AnisotropicBounded,NonStandardBounded}.{serial,mpi2}` | `test_timestep.cpp` · `TimeStepLong` (one test each) | At \f$0.98\,\Delta t_{\max}\f$, over 10000 solves, the maximum of random histories in the last 1000 post-solve samples is at most twice the maximum in the first 1000. | self-consistency (envelope bound) | `long`, `mpi` | 3.9 / 3.6, 6.1 / 4.6, 14.4 / 9.1 |

### Manufactured solution and rank comparison

These entries test the evolution under a manufactured discrete source and the
agreement of one-rank and two-rank runs. The producers (`manufactured`,
`rank_cases`, the 12 ladders below, `dispersion` and `convergence.staggered`)
write owned fields below `CHDR_SOLVER_RECORD_ROOT`; the CTest fixtures run
`chdr.solver.ranks.prepare` first and `chdr.solver.ranks` last. Selecting
`chdr.solver.ranks` alone therefore runs 34 entries. Without Python 3 or an MPI
launcher, CMake registers the two script entries with the label `unavailable`.

| Entry | Source (lines) · suite | What it checks | Reference | Labels | Time (s) |
|---|---|---|---|---|---|
| `chdr.solver.gather.{serial,mpi2}` | `test_gather.cpp` (66) · `RankGather` (2 tests) | `gatherOwned` returns on rank 0 exactly the owned values in global order without halo cells, for each split axis; a repeated index throws. | exact lattice value; refusal | `unit`, `mpi` | 1.6 / 1.6 |
| `chdr.solver.manufactured.{serial,mpi2}` | `test_manufactured.cpp` (307) · `Manufactured` (2 tests) | For R1, R2 and R2Prime, 256 solves driven by a sampled discrete manufactured source keep all four potentials within `1e-12` of the exact samples; a charge sampled half a step late misses by at least `1e-3`. | exact lattice value | `unit`, `mpi` | 11.6 / 7.1 |
| `chdr.solver.rank_cases.{serial,mpi2}` | `test_ranks.cpp` (150) · `RankCases` (3 tests) | Writes the fields of two ten-period gauge-wave runs (Standard, NonStandard) and of six monitor runs on rough data (three configurations, with and without Nyquist content) for the rank comparison; in process it asserts only the final identities. | one against two ranks (producer) | `unit`, `mpi` | 3.9 / 3.0 |
| `chdr.solver.ranks.prepare` | `prepare_rank_records.py` (18) | Moves earlier records aside and creates empty record directories; it checks nothing. | smoke (fixture setup) | `unit`, `ranks` | 0.0 |
| `chdr.solver.ranks` | `compare_ranks.py` (156) | The one-rank and two-rank directories hold the same 193 records, and they agree per cell to `1e-13` of the field maximum and in maxima, sums, residuals and phases to `1e-12`. | one against two ranks | `unit`, `ranks` | 10.3 |

### Refinement ladders and wave norms

Each of the 12 `convergence` entries runs one refinement ladder (\f$N=16,32,64\f$,
one period, split x) for one stencil, one wave family and one reconstruction,
and asserts the observed orders \f$\log_2(e_{32}/e_{64})\f$ against the
continuum wave. The pair (64, 128) and the first-order R1 control run in the
`long.Convergence` entries. Run them with
`ctest --test-dir build -R '^chdr\.solver\.(convergence\.(Standard|NonStandard)|wave_norms)\.' --output-on-failure`.

| Entry | Source (lines) · suite | What it checks | Reference | Labels | Time (s) |
|---|---|---|---|---|---|
| `chdr.solver.convergence.{Standard,NonStandard}.Mixed.R1.{serial,mpi2}` | `test_convergence.cpp` (152), `Convergence.h` (387) · `Convergence.*MixedR1` | For the source-free mixed wave with R1, the B order is at least 1.8 on (32, 64); final identities. | continuum analytic | `unit`, `mpi` | 3.3 / 2.6, 5.6 / 3.7 |
| `chdr.solver.convergence.{Standard,NonStandard}.Mixed.{R2,R2Prime}.{serial,mpi2}` | as above · `Convergence.*Mixed{R2,R2Prime}` | For the mixed wave, the E, B and Gauss-residual orders are at least 1.8 on (32, 64); final identities. | continuum analytic | `unit`, `mpi` | 3.4 to 5.7 / 2.6 to 3.9 |
| `chdr.solver.convergence.{Standard,NonStandard}.Sourced.R1.{serial,mpi2}` | as above · `Convergence.*SourcedR1` | For the sourced Lorenz-gauge wave with R1, the B order is at least 1.8 on (32, 64), and the gauge-residual order is at least 1.8 when the fine gauge error exceeds `1e-12`; final identities. | continuum analytic | `unit`, `mpi` | 3.8 / 2.9, 7.4 / 4.7 |
| `chdr.solver.convergence.{Standard,NonStandard}.Sourced.{R2,R2Prime}.{serial,mpi2}` | as above · `Convergence.*Sourced{R2,R2Prime}` | For the sourced wave, the E, B and Gauss orders are at least 1.8 on (32, 64), and the gauge order is at least 1.8 when the fine gauge error exceeds `1e-12`; final identities. | continuum analytic | `unit`, `mpi` | 4.1 to 7.6 / 2.9 to 4.8 |
| `chdr.solver.wave_norms.{serial,mpi2}` | `test_convergence.cpp`, `WaveNorms.h` (219) · `WaveNorms` | The gauge ratio that `waveErrors` measures equals an independent complex-amplitude calculation to `1e-12` for R1, R2 and R2Prime. | exact lattice value | `unit`, `mpi` | 1.6 / 1.6 |

### Long ladders and gauge-residual histories

`test_convergence_long.cpp` and `test_staggered_convergence.cpp` extend the
ladders to \f$N=128\f$. `test_gauge_growth.cpp` evolves a wave for ten periods
and records the global \f$L^2\f$ norm of the discrete gauge residual after every
solve. Run them with
`ctest --test-dir build -R '^chdr\.solver\.(long\.(Convergence|staggered|GaugeGrowth)|convergence\.staggered)\.' --output-on-failure`.

| Entry | Source (lines) · suite | What it checks | Reference | Labels | Time (s) |
|---|---|---|---|---|---|
| `chdr.solver.long.Convergence.Standard.{serial,mpi2}` | `test_convergence_long.cpp` (30) · `LongConvergence.Standard` | The six collocated ladders (mixed and sourced wave, R1, R2, R2Prime) on \f$N=16\f$ to 128: B order, and for R2 and R2Prime the E and Gauss orders, at least 1.8 on (32, 64) and (64, 128); R1 E and Gauss order between 0.8 and 1.2 on (64, 128); for the sourced wave the gauge order at least 1.8 when the fine gauge error exceeds `1e-12`; final identities at every level. | continuum analytic | `long`, `mpi` | 311.1 / 168.0 |
| `chdr.solver.long.Convergence.NonStandard.{serial,mpi2}` | as above · `LongConvergence.NonStandard` | The same six ladders and rules on NonStandard meshes \f$N\times N\times 3N/2\f$. | continuum analytic | `long`, `mpi` | 548.7 / 339.1 |
| `chdr.solver.convergence.staggered.{serial,mpi2}` | `test_staggered_convergence.cpp` (27) · `StaggeredConvergence.Quick` | The mixed and the sourced wave with R2 in the staggered reading on \f$N=16,32,64\f$: E, B and Gauss orders at least 1.8 on (32, 64), and for the sourced wave the gauge order at least 1.8 when the fine gauge error exceeds `1e-12`; it also writes rank records. | continuum analytic | `unit`, `mpi` | 6.0 / 3.9 |
| `chdr.solver.long.staggered.{serial,mpi2}` | as above · `LongStaggeredConvergence.Standard` | The staggered ladders on \f$N=16\f$ to 128 with the same rules on both pairs. | continuum analytic | `long`, `mpi` | 95.9 / 54.0 |
| `chdr.solver.long.GaugeGrowth.{Standard,NonStandard}.{serial,mpi2}` | `test_gauge_growth.cpp` (186) · `LongGaugeGrowth` | Over ten periods with R2, collocated, the largest gauge-residual norm in the last period is at most 1.05 times the largest in the first; for the exact discrete gauge wave the first norm equals its closed form and the ratio equals 1 to `1e-12`. | self-consistency; exact lattice value | `long`, `mpi` | 3.1 / 2.5, 5.7 / 3.8 |
| `chdr.solver.long.GaugeGrowth.StaggeredRecord.{serial,mpi2}` | as above · `LongGaugeGrowth.StaggeredRecord` | Writes the gauge-residual history of the discrete gauge wave in the staggered reading; it asserts only finite norms and the final identities. | regression record (prints only) | `record`, `mpi` | 2.4 / 2.1 |

### Gauge-wave cancellation and level equivalence

`test_levels.cpp` loads the exact discrete gauge wave of `Waves.h`; for this
wave R2 and R2Prime give zero E and B, and R1 gives a closed-form E. The test
divides the largest field error by the largest sampled
\f$\lvert D_c\phi\rvert\f$. `test_equivalence.cpp` compares R2 and R2Prime
bit by bit. Run them with
`ctest --test-dir build -R '^chdr\.solver\.(levels|equivalence|gauge)\.' --output-on-failure`.

| Entry | Source (lines) · suite | What it checks | Reference | Labels | Time (s) |
|---|---|---|---|---|---|
| `chdr.solver.levels.{serial,mpi2}` | `test_levels.cpp` (258) · `Levels` (2 tests) | Over 128 steps of an axial and an oblique gauge wave, the normalized E and B errors of R1, R2 and R2Prime stay at most `1e-12`. | exact lattice value | `unit`, `mpi` | 1.9 / 1.9 |
| `chdr.solver.equivalence.{serial,mpi2}` | `test_equivalence.cpp` (209) · `LevelEquivalence` (2 tests) | R2 and R2Prime on separate containers give bitwise equal histories, A, E and B after every one of 60 solves of the sourced wave. | self-consistency (bitwise) | `unit`, `mpi` | 1.7 / 1.7 |
| `chdr.solver.gauge.{ObliqueStandard,ObliqueNonStandard}.{serial,mpi2}` | `test_levels.cpp` · `GaugeWave.Oblique*` | The same rule over ten periods of an oblique gauge wave on 32 by 32 by 32 (Standard) or 32 by 32 by 48 (NonStandard) cells. | exact lattice value | `unit`, `mpi` | 3.4 / 2.9, 7.5 / 5.3 |
| `chdr.solver.gauge.AxialNonStandardStable.{serial,mpi2}` | `test_levels.cpp` · `GaugeWave.AxialNonStandardStable` | The same rule for an axial NonStandard gauge wave at safety factor 0.98 over 640 steps. | exact lattice value | `unit`, `mpi` | 3.5 / 3.0 |
| `chdr.solver.gauge.AxialNonStandardBoundaryRecord.{serial,mpi2}` | `test_levels.cpp` · `GaugeWave.AxialNonStandardBoundaryRecord` | The axial case at the stability bound (\f$S=1\f$): prints the normalized E and B and the z-parity projection; it asserts only finite values and the final identities. | regression record (prints only) | `unit`, `mpi` | 3.5 / 3.1 |

### Dispersion

`test_dispersion.cpp` starts one travelling mode of \f$A_x\f$, projects the
evolved array on that mode after every step for one period, and compares the
mean phase advance per step with the discrete dispersion relation
\f$\Omega\Delta t=2\arcsin\big(\tfrac{\Delta t}{2}\sqrt{S_z(c_x+c_y)+c_z}\big)\f$,
\f$c_d=(2\sin(k_dh_d/2)/h_d)^2\f$, with \f$S_z=1\f$ for Standard. All cases use
R1 and the collocated reading. Run them with
`ctest --test-dir build -R '^chdr\.solver\.(record\.)?dispersion\.' --output-on-failure`.

| Entry | Source (lines) · suite | What it checks | Reference | Labels | Time (s) |
|---|---|---|---|---|---|
| `chdr.solver.dispersion.{serial,mpi2}` | `test_dispersion.cpp` (422) · `Dispersion` (2 tests) | For 19 modes (Standard \f$N=32\f$, NonStandard \f$N_z=24\f$, two or three safety factors), the measured phase per step equals the discrete relation to `1e-12`, and the NonStandard axial phase velocity at \f$S=1\f$ equals 1; it writes 19 rank records. | exact lattice value | `unit`, `mpi` | 2.4 / 2.1 |
| `chdr.solver.record.dispersion.serial` | as above · `RecordDispersion` (2 tests) | The same rule on 37 finer cases (Standard \f$N=64\f$, NonStandard \f$N_z=24,48,96\f$); writes one phase-velocity record per case. | exact lattice value; regression record | `long` | 34.7 |
| `chdr.solver.record.dispersion.mpi2` | as above | The same 37 cases on two ranks. | exact lattice value; regression record | `long`, `mpi` | 18.3 |
| `chdr.solver.record.endpoint.{serial,mpi2}` | `test_endpoint_records.cpp` (174) · `TimeStepRecords` (4 tests) | At the stability bound \f$S=1\f$, prints the amplitude of the highest periodic scalar mode for 256 steps beside the host recurrence; it asserts only finite values and the final identities. | regression record (prints only) | `unit`, `mpi` | 2.4 / 2.1 |

### Rank-layout instruments and co-moving reconstruction

`SerialReference.h` repeats a distributed run on `MPI_COMM_SELF` at rank 0, and
`OwnedStatistics.h` computes the maximum, the absolute sum and the square sum
of the owned cells with pairwise sums; `test_serial_reference.cpp` and
`test_statistics.cpp` test these tools. `test_comoving.cpp` loads the
potentials of a charge distribution that moves rigidly along z at
\f$\gamma=118.417071\f$ (`ComovingPotentials.h`) and calls `reconstructEB()`
once, without time stepping. Run them with
`ctest --test-dir build -R '^chdr\.solver\.(serial_reference|statistics|comoving)\.' --output-on-failure`.

| Entry | Source (lines) · suite | What it checks | Reference | Labels | Time (s) |
|---|---|---|---|---|---|
| `chdr.solver.serial_reference.{serial,mpi2}` | `test_serial_reference.cpp` (41), `SerialReference.h` (193) · `SerialReference` | Eight solves of the sourced wave agree with a one-rank replay per cell to `1e-13` of the field maximum and in maximum and sums to `1e-12`; only the `mpi2` entry compares two decompositions. | one against two ranks | `unit`, `mpi` | 1.6 / 1.6 |
| `chdr.solver.statistics.{serial,mpi2}` | `test_statistics.cpp` (101), `OwnedStatistics.h` (184) · `StableStatistics` (4 tests) | Pairwise sums are exact for leaves \f$2^{53},1,1,1,1\f$ and for integer fields with NaN ghosts; an owned NaN makes the maximum infinite and the sums non-finite. | exact lattice value; self-consistency (NaN guard) | `unit`, `mpi` | 1.6 / 1.6 |
| `chdr.solver.comoving.{serial,mpi2}` | `test_comoving.cpp` (298), `ComovingPotentials.h` (425) · `ComovingT7` (7 tests) | For one mode and a Gaussian bunch train, the measured \f$L^2\f$ and maximum errors of the reconstructed \f$E_z\f$ and \f$E_x\f$ equal the predictions of the exact difference symbols to relative `1e-9`, and for the bunch train also their leading-order series to 0.10, for R1, R2, R2Prime and the staggered reading. | exact lattice value; continuum analytic | `unit`, `mpi` | 2.9 / 3.5 |

### Homogeneous dielectric

These entries test `ChdrSolver` with `MaterialFormulation::Homogeneous`, the
update for one uniform relative permittivity \f$\varepsilon\ge1\f$. Most medium
suites live in the same source files as the vacuum suites. The update runs the
vacuum kernel with \f$\Delta t/\sqrt\varepsilon\f$, so several medium entries
measure the same numbers as a vacuum entry. Run the short entries with
`ctest --test-dir build -R '^chdr\.solver\.medium\.' -LE 'long|record' --output-on-failure`.
Selecting `chdr.solver.medium.ranks` alone also runs its preparation and its two
producers. Without Python 3 or an MPI launcher, the two medium rank-script
entries remain registered as `unavailable` placeholders.

| Entry | Source (lines) · suite | What it checks | Reference | Labels | Time (s) |
|---|---|---|---|---|---|
| `chdr.solver.medium.Validity.{serial,mpi2}` | `test_medium_guard.cpp` (93) · `MediumValidity` | Five invalid material rows (one bad cell, NaN, \f$\varepsilon<1\f$, \f$\alpha\ne1\f$, NonStandard) throw from `validateMaterialUpdate`, and the stable step for \f$\varepsilon=1,2.13,4\f$ equals \f$\sqrt{\varepsilon/\sum_d h_d^{-2}}\f$ to `1e-15`. | refusal; exact lattice value (bound) | `unit`, `mpi` | 1.6 / 1.7 |
| `chdr.solver.medium.Guard.{serial,mpi2}` | `test_medium_guard.cpp` · `MediumGuard` | The `ChdrSolver` constructor refuses the invalid rows on every rank and accepts the valid ones, and the vacuum monitors refuse a medium. | refusal | `unit`, `mpi` | 1.6 / 1.7 |
| `chdr.solver.medium.Equivalence.{serial,mpi2}` | `test_equivalence.cpp` (209) · `MediumEquivalence` | At \f$\varepsilon\equiv1\f$ the Homogeneous update gives bitwise the same potentials, E and B as the vacuum `ChdrSolver` at each of 256 steps, for R1, R2, R2Prime and the staggered reading. | regression record (bitwise) | `unit`, `mpi` | 2.4 / 2.4 |
| `chdr.solver.medium.Manufactured.{serial,mpi2}` | `test_manufactured.cpp` (307) · `MediumManufactured` | With physical sources \f$\rho=\varepsilon^2s_0\f$, \f$J_d=\varepsilon s_d\f$, the potentials follow a discrete three-mode history to `1e-12` over 256 steps for \f$\varepsilon=2.13,4\f$; four wrong source conventions miss by at least `1e-3`. | exact lattice value | `unit`, `mpi` | 6.1 / 4.1 |
| `chdr.solver.medium.TimeStep.{serial,mpi2}` | `test_timestep.cpp` (294) · `MediumTimeStep` | The step setter gives \f$S\sqrt{\varepsilon/\sum_d h_d^{-2}}\f$ to `1e-15` and refuses invalid factors; above the bound random histories grow at least 1000-fold. | exact lattice value (bound); refusal | `unit`, `mpi` | 1.7 / 1.7 |
| `chdr.solver.medium.Dispersion.{serial,mpi2}` | `test_dispersion.cpp` (422) · `MediumDispersion.RankCases` | The phase per step of 20 single modes (\f$\varepsilon=2.13,4\f$, \f$N=32\f$) equals \f$2\arcsin((\Delta t/2)\sqrt{K^2/\varepsilon})\f$ to `1e-12`; writes 20 rank records. | exact lattice value | `unit`, `mpi` | 2.5 / 2.1 |
| `chdr.solver.medium.Source.{serial,mpi2}` | `test_convergence.cpp` (152) · `MediumSource` (2 tests) | A driven Lorenz-gauge wave at \f$\varepsilon=2.13\f$ reaches E and B orders of at least 1.8 on (32, 64) in both readings. | continuum analytic; one against two ranks | `unit`, `mpi` | 7.9 / 5.9 |
| `chdr.solver.medium.Convergence.{Collocated,Staggered}.{serial,mpi2}` | `test_convergence.cpp` · `MediumConvergence` | Six source-free plane waves at speed \f$1/n\f$ reach E and B orders of at least 1.8 on (32, 64); the \f$N=32\f$ run is replayed on one rank. | continuum analytic; one against two ranks | `unit`, `mpi` | 31.2 / 22.5, 33.7 / 25.5 |
| `chdr.solver.medium.long.Dispersion.{serial,mpi2}` | `test_dispersion.cpp` · `MediumDispersionLong` | The dispersion rule of `MediumDispersion` at \f$N=64\f$. | exact lattice value | `long`, `mpi` | 26.1 / 18.7 |
| `chdr.solver.medium.record.Dispersion.{serial,mpi2}` | `test_dispersion.cpp` · `MediumDispersionRecord` | Prints phase-error-per-wavelength records for 48 cases (\f$N=16\f$ to 128, three directions, two step rules); it asserts only a nonzero amplitude. | regression record (prints only) | `record`, `mpi` | 506.6 / 342.8 |
| `chdr.solver.medium.long.TimeStep.{serial,mpi2}` | `test_timestep.cpp` · `MediumTimeStepLong` | At safety 0.98 over 10000 steps, the late maximum of random histories is at most twice the early maximum. | self-consistency (envelope bound) | `long`, `mpi` | 17.9 / 15.9 |
| `chdr.solver.medium.long.Convergence.{Collocated,Staggered}.{serial,mpi2}` | `test_convergence.cpp` · `MediumConvergenceLong` | The six plane-wave ladders on \f$N=16\f$ to 128 with orders at least 1.8 on (32, 64) and (64, 128). | continuum analytic; one against two ranks | `long`, `mpi` | 482.4 / 322.9, 500.5 / 278.1 |
| `chdr.solver.medium.long.Source.{Collocated,Staggered}.{serial,mpi2}` | `test_convergence.cpp` · `MediumSourceLong` | The driven-wave ladder on \f$N=16\f$ to 128 with orders at least 1.8 on both pairs. | continuum analytic; one against two ranks | `long`, `mpi` | 61.0 / 40.0, 68.0 / 34.3 |
| `chdr.solver.medium.ranks.prepare` | `prepare_rank_records.py` (18) | Creates empty record directories under `runs/medium-ranks`; it checks nothing. | smoke (fixture setup) | `unit`, `ranks` | 0.0 |
| `chdr.solver.medium.ranks` | `compare_ranks.py` (156) | The 36 medium records of one and two ranks agree per cell to `1e-13` of the field maximum and in sums and phases to `1e-12`. | one against two ranks | `unit`, `ranks` | 0.3 |

### Mur faces and open boxes

These entries test the Mur faces of `ChdrSolver` on boxes with one or more
nonperiodic axes, in vacuum and in a homogeneous dielectric with matched Mur
speed: the face description and the Mur speed, the reflection coefficient
against the discrete theory of Mur's condition, the fields on the outermost
layer and at edges and corners, long-run growth, and one against two ranks.
Run the short entries with
`ctest --test-dir build -R '^chdr\.solver\.(mur|open)\.' -LE 'long|record' --output-on-failure`.

| Entry | Source (lines) · suite | What it checks | Reference | Labels | Time (s) |
|---|---|---|---|---|---|
| `chdr.solver.mur.faces.{serial,mpi2}` | `test_mur_faces.cpp` (67) · `MurFaces` (2 tests) | Lower-face kinds follow the layout, all face kinds and speeds survive `initialize()`, NonStandard is refused on layouts with a Mur axis, and the matched and physical Mur speeds at \f$\varepsilon=2.13\f$ are exact. | refusal; exact lattice value (speed) | `unit`, `mpi` | 1.6 / 1.7 |
| `chdr.solver.mur.medium.{serial,mpi2}` | `test_mur_medium.cpp` (120) · `MurMedium` (2 tests) | In an all-Mur box, Homogeneous at \f$\Delta t\f$ equals the vacuum kernel at \f$\Delta t/\sqrt\varepsilon\f$ bit for bit over 256 steps; a test-only kernel speed leaves step, face kinds and speeds unchanged. | regression record (bitwise); self-consistency | `unit`, `mpi` | 2.4 / 2.1 |
| `chdr.solver.mur.reflection.{serial,mpi2}` | `test_mur_reflection.cpp` (151) · `MurReflection.ExactTwoWaveRows` | After 30 steps from an exact two-wave state, the fitted reflection coefficient and the near-face cells match the discrete Mur prediction to `1e-11` for 13 rows (vacuum, \f$\varepsilon=2.13\f$, anisotropic cells); a perturbed speed and the old halo order miss. | exact lattice value; one against two ranks | `unit`, `mpi` | 28.9 / 17.3 |
| `chdr.solver.mur.field-reflection.{serial,mpi2}` | `test_mur_field_reflection.cpp` (102) · `MurFieldReflection` | TE \f$E_y\f$ and TM averaged \f$B_y\f$ near the face match the incident-plus-reflected oracle to `1e-10`, with the correct signs. | exact lattice value; one against two ranks | `unit`, `mpi` | 2.0 / 1.9 |
| `chdr.solver.mur.refinement.serial` | `test_mur_refinement.cpp` (124) · `MurReflectionRefinement.CentrePlaneOrder` | Without a solver, the discrete Mur reflection converges to the continuum value with order in [1.9, 2.1] for 11 rows. | continuum analytic (oracle only) | `record` | 1.7 |
| `chdr.solver.mur.records.serial` | `test_mur_refinement.cpp` · `MurReflectionRefinement.CodeConfirmedRows` | Without a solver, prints 117 reflection records (discrete and continuum coefficients); asserts nothing. | regression record (prints only) | `record` | 1.6 |
| `chdr.solver.mur.outer-layer.{serial,mpi2}` | `test_mur_outer_layer.cpp` (195) · `MurOuterLayer` | On the outermost layers of open axes, E, averaged B and raw B of a polynomial fixture match their analytic values to `1e-13` in both readings with R2; zero ghosts miss by at least 0.1. | exact lattice value; one against two ranks | `unit`, `mpi` | 1.6 / 1.7 |
| `chdr.solver.mur.families.{serial,mpi2}` | `test_mur_families.cpp` (181) · `MurEdgeCornerFamilies` (2 tests) | After 3 steps, the potentials in the causal region of a Mur edge or corner equal an exact sum of lattice plane waves to `1e-13` for 72 fixtures; the old finishing order must miss by at least `1e-3` on the mixed edge and on two ranks, and passes for the serial corner-zero control. | exact lattice value; one against two ranks | `unit`, `mpi` | 1.9 / 2.0 |
| `chdr.solver.mur.stability.serial` | `test_mur_stability.cpp` (190) · `MurStability.InnerDataEnvelope` | For 16 open-box cases over 10000 steps, the late maximum of random potentials is at most twice the early maximum. | self-consistency (envelope bound) | `long` | 44.1 |
| `chdr.solver.mur.stability.mpi2` | as above | Four of the cases for 1000 steps on two ranks, with the same rule and a one-rank replay to `1e-13`. | one against two ranks; self-consistency | `long`, `mpi` | 3.3 |
| `chdr.solver.mur.stability-control.serial` | as above · `MurStability.MismatchedSpeedControl` | With a mismatched Mur speed in a dielectric, the late maximum is at least \f$10^6\f$ times the early maximum (negative control). | self-consistency (negative control) | `long` | 5.6 |
| `chdr.solver.open.beam-ranks.mpi2` | `test_open_beam.cpp` (188) · `OpenBeam.ReducedRankWitness` | A 300-step beam run in an all-Mur 40 by 80 by 200 box equals a one-rank replay to `1e-13`, and its face Fourier sums to `1e-12`. | one against two ranks | `long`, `mpi` | 70.1 |
| `chdr.solver.open.beam-full.mpi2` | `test_open_beam.cpp` · `OpenBeam.FullBoxAgainstPeriodicWitness` | The reflection coefficient of the central plane-wave part of a Gaussian beam in a 40 by 160 by 400 all-Mur box differs from a transversely periodic witness by at most `1e-4`. | self-consistency (beam against periodic witness) | `long`, `mpi` | 312.6 |
| `chdr.solver.open.record.gauge` | `test_open_gauge_records.cpp` (142) · `OpenGaugeRecords` | Prints the gauge and Gauss residuals at a Mur face next to the plane-wave prediction; asserts nothing. | regression record (prints only) | `record` | 1.7 |
| `chdr.solver.open.record.zero-ghost` | `test_mur_outer_layer.cpp` · `MurOuterLayerRecords` | Prints the raw-B error at an open face with zero ghosts at four resolutions; asserts nothing. | regression record (prints only) | `record` | 1.6 |
| `chdr.solver.open.record.all-cell-drift` | `test_mur_stability.cpp` · `MurStabilityRecords` | Prints the envelope ratio of random data on all cells over 10000 steps; asserts nothing. | regression record (prints only) | `record` | 4.9 |
| `chdr.solver.mur.holdout` | `test_mur_reflection.cpp` · `MurReflection.DISABLED_Holdout` | Disabled; the former held-out acceptance check: two reflection rows at \f$\varepsilon_r=3.1\f$ on anisotropic cells, two ranks. | — | `holdout`, `mpi` | disabled |

## `core/solver/formulation`: O1-VS

The executable `test-chdr-formulation` tests the O1-VS material formulation of
`ChdrSolver`, the update for planar dielectric interfaces in the staggered
reading. Each entry runs all five tests of the suite `O1VS`. Run it with
`ctest --test-dir build -R '^chdr\.solver\.formulation\.' --output-on-failure`.

| Entry | Source (lines) · suite | What it checks | Reference | Labels | Time (s) |
|---|---|---|---|---|---|
| `chdr.solver.formulation.o1vs.{serial,mpi2}` | `test_o1vs_regression.cpp` (167) · `O1VS` (5 tests) | At \f$\varepsilon\equiv1\f$ O1-VS gives the vacuum staggered R2 step bit for bit over 256 steps in a periodic and in a Mur box, at \f$\varepsilon=2.13\f$ it matches Homogeneous to `1e-13`, invalid configurations are refused, and the harmonic face mean of large values stays finite. | regression record (bitwise); refusal | `formulation`, `unit`, `mpi` | 1.9 / 1.8 |

## `opt-in/diagnostics`: monitors, probe and energy

These tests check the residual monitors, the Gauss formula check, the Fourier
probe and the potential invariant of `opt-in/diagnostics/`. One executable,
`test-chdr-diagnostics`, holds all suites. `chdr.diagnostics.probe.ranks` runs
the long probe test on one and on two ranks through `compare_probe_ranks.py`
and compares the two output files. Run them with
`ctest --test-dir build -R '^chdr\.diagnostics\.' --output-on-failure`.

| Entry | Source (lines) · suite | What it checks | Reference | Labels | Time (s) |
|---|---|---|---|---|---|
| `chdr.diagnostics.energy.mpi2` | `test_interface_energy.cpp` (118) · `InterfaceEnergy` (2 tests) | An O1-VS slab (\f$\varepsilon_r=2.13\f$) evolves random histories for 50 solves, split across the slab face and along it: the vector and scalar potential invariants drift by at most `1e-12`, histories are bitwise equal to a one-rank replay, and energy records agree within `1e-14`. | one against two ranks; self-consistency | `unit`, `mpi` | 1.7 |
| `chdr.diagnostics.monitors.{serial,mpi2}` | `test_monitors.cpp` (224) · `ChdrMonitors` (5 tests) | On rough data with Nyquist content, the normalized divergence and Faraday residuals stay at most `1e-12`, a planted E error is detected, the continuity residual matches a planted value \f$p\f$ (maximum \f$\lvert p\rvert\f$, L2 norm \f$\lvert p\rvert\sqrt{960}\f$, signed cell values) to `1e-12` of its scale, and wrong call order is refused. | exact lattice value (identities); refusal | `unit`, `mpi` | 1.7 / 1.7 |
| `chdr.diagnostics.monitor-records.{serial,mpi2}` | `test_monitor_records.cpp` (95) · `ChdrMonitorRecords` (2 tests) | The normalized mismatch between the Gauss defect and its discrete formula is at most `1e-12`, and the report keeps this mismatch and the L2 norm of the defect after the live source changes. | exact lattice value; self-consistency | `unit`, `mpi` | 1.7 / 1.7 |
| `chdr.diagnostics.probe.{serial,mpi2}` | `test_probe.cpp` (154) · `FourierProbe` (5 tests) | `FourierProbe` refuses invalid input, its phase equals an integer reference to `2e-15` at steps near \f$2^{52}\f$, Hann weights match \f$\sin^2(\pi n/N)\f$, and (two ranks) an unowned box fails on all ranks. | exact lattice value; refusal | `unit`, `mpi` | 1.7 / 1.7 |
| `chdr.diagnostics.mur-layer-refusals.serial` | `test_mur_layer.cpp` (36) · `MurLayer` | E, B and raw-B probe boxes on Mur faces and the monitors on a Mur layout throw; a box one cell inside and a \f$\Phi\f$ box on the face do not. | refusal | `unit` | 1.6 |
| `chdr.diagnostics.probe_wiring.{serial,mpi2}` | `test_probe_full.cpp` (628) · `FourierProbeSpectra` (4 tests) | Channel labels and coordinates equal fixed tables, a staggered plane equals the finite-window closed form to `1e-12` of its scale, and probes on an evolved mode equal the discrete samples. | exact lattice value | `unit`, `mpi` | 1.9 / 1.8 |
| `chdr.diagnostics.probe.ranks` | `compare_probe_ranks.py` (118) | Runs `FrozenGatesAndNaNGhosts` on one and two ranks; all 17292 complex outputs agree to `1e-14` of the scale. | one against two ranks | `unit`, `probe`, `mpi` | 5.6 |

## `apps/vacuum-wave`: the vacuum-wave program and the field writer

These tests check `chdr-vacuum-wave` (a plane wave in a periodic vacuum box,
Standard update, one period) and the field writer `core/output/FieldOutput.h`,
which all three programs use. `test-chdr-vacuum` holds the GoogleTests; Python
scripts run the program, check its CSV output and, with ParaView's `pvpython`,
read its field files. Without `pvpython` the reader entries report a skip.
The oracle follows `fallahi2020mithra20fullwavesimulation`, Eqs. (3.6)–(3.20).
Run them with `ctest --test-dir build -R '^chdr\.vacuum\.' --output-on-failure`.

| Entry | Source (lines) · suite | What it checks | Reference | Labels | Time (s) |
|---|---|---|---|---|---|
| `chdr.vacuum.unit.{serial,mpi2}` | `test_first_step.cpp` (317), `test_setup.cpp` (181), `test_diagnostics.cpp` (333), `test_policy.cpp` (115) · 22 tests | Initial histories and one solve match an independent host oracle to `1e-10`, the relative \f$L^2\f$ error against the continuum wave is at most 0.01, a 128-step scalar recurrence matches to `1e-12`, and domain errors and injected faults are refused on all ranks. The two-rank entry shuffles the test order with the same fixed seed on both ranks. | exact lattice value; continuum analytic; refusal | `unit`, `mpi` | 1.9 / 1.7 |
| `chdr.vacuum.output.{serial,mpi2}` | `output_fixture.cpp` (187), `read_output_fixture.py` (163) | A three-frame synthetic fixture written under a comma-decimal locale reads back through VTK and ParaView bit-exact, with exactly the owned cells and the series times. | exact lattice value (fixture formula) | `output`, `reader`, `mpi` | 6.0 / 5.4 |
| `chdr.vacuum.output.failure.{piece,wrapper,series,nonfinite}` | as above | A publication failure or a NaN on one of two ranks reaches every rank and leaves a readable series of the earlier frame. | refusal | `output`, `reader`, `mpi` | 5.5, 5.2, 5.2, 5.8 |
| `chdr.vacuum.propagation` | `run_wave_checks.py` (122), `wave_checks.py` (84) | Three runs (coarse, fine, coarse on two ranks): CSV errors at most `1e-10`, phase speed within 0.002 of 1, relative \f$L^2\f$ errors at most 0.01 and 0.003 with coarse-to-fine ratio in [3, 5], two ranks equal to one rank within `1e-11`. | continuum analytic (second order); exact lattice value; one against two ranks | `integration`, `mpi` | 6.1 |
| `chdr.vacuum.cli` | `run_cli_checks.py` (90) | Six valid command lines exit 0; seven invalid ones (one on two ranks) exit with the expected reason. | refusal; smoke | `integration`, `mpi` | 21.4 |
| `chdr.vacuum.readback` | `verify_wave_output.py` (156) | The 15 field frames of `chdr.vacuum.propagation` read back with the expected extents and times, match the discrete fields to `1e-10` and the CSV errors to `1e-12`, and two ranks equal one rank within `1e-11`. | exact lattice value; continuum analytic; one against two ranks | `reader` | 6.6 |
| `chdr.vacuum.readback.partial` | `verify_partial_output.py` (59) | The zero-step and one-step outputs of `chdr.vacuum.cli` list exactly their frames and pass the frame checks. | exact lattice value; continuum analytic | `reader` | 3.7 |

## `apps/planar-interface`: the planar-interface program

These tests measure the reflection and transmission of a plane pulse at a
planar dielectric face with the functions of `apps/planar-interface/PlanePulse.h`,
and run the program `planar-interface` once. The expected values come from
`InterfaceOracles.h`: the exact lattice coefficients of the staggered Standard
scheme and the continuum Fresnel coefficients, both written independently of
the solver. Run the fast entries with
`ctest --test-dir build -R '^chdr\.interface\.(normal|smoke)' --output-on-failure`
and the Brewster entry with `ctest --test-dir build -R '^chdr\.interface\.long\.brewster$'`.

| Entry | Source (lines) · suite | What it checks | Reference | Labels | Time (s) |
|---|---|---|---|---|---|
| `chdr.interface.normal` | `test_plane_interface.cpp` (131), `InterfaceOracles.h` (241) · `PlaneInterface.Normal` | A TE pulse at normal incidence onto \f$\varepsilon_r=2.13\f$: \f$r\f$ and \f$t\f$ equal the exact lattice values to `1e-11` and Fresnel to `5e-3`, the invariant drifts by at most `1e-12` and equals Yee's energy at the start. | exact lattice value; continuum analytic (Fresnel) | — | 2.5 |
| `chdr.interface.normal.mpi2` | as above · `PlaneInterface.NormalSeam` | The same pulse and checks on two ranks with the rank seam on the material face. | exact lattice value; continuum analytic | `mpi` | 2.3 |
| `chdr.interface.long.brewster` | as above · `PlaneInterface.Brewster` | A TM pulse at the Brewster angle onto \f$\varepsilon_r=2.13\f$: \f$r\f$ and \f$t\f$ equal the exact lattice values to `1e-11` and \f$\lvert r\rvert\le10^{-2}\f$ against the Fresnel zero. | exact lattice value; continuum analytic | `long` | 356.1 |
| `chdr.interface.reserved` | as above · `PlaneInterface.DISABLED_Reserved` | Disabled; the former held-out acceptance check: a TM pulse from \f$\varepsilon_r=3.4\f$ at 18 degrees on two ranks. | — | `long`, `reserved`, `mpi` | disabled |
| `chdr.interface.smoke` | `smoke.cmake` (44) | A YAML half-space case gives the expected frames, energy rows and solve count, and the legacy command gives the same settings and byte-identical field files. | smoke; regression record (bitwise) | `smoke` | 3.4 |

## `apps/interface-movies`: the interface-movies program

These tests check the starts of the eight movie cases of `interface-movies`,
the per-face Mur speeds of `MovieFaceSolver` and the two energy definitions at
load. Run them with
`ctest --test-dir build -R '^chdr\.interface\.movies\.' --output-on-failure`.

| Entry | Source (lines) · suite | What it checks | Reference | Labels | Time (s) |
|---|---|---|---|---|---|
| `chdr.interface.movies.unit.{serial,mpi2}` | `test_movie_start.cpp` (119) · `Movies` (4 tests) | All eight starts pass their placement check and misplaced starts throw; with unit face speeds the potential \f$A_n\f$ after 48 steps equals that of the library manager exactly (maximum difference zero); the TE start loads the expected components and the TM histories have weighted divergence below `1e-13`. | self-consistency; regression record (exact numerical equality) | `unit`, `mpi` | 2.1 / 1.8 |
| `chdr.interface.movies.smoke` | `movies-smoke.cmake` (25) | Each of the eight cases runs four steps and writes 2 field frames and an energy row for state 4. | smoke | `smoke` | 13.6 |

## `../mesh`: the mesh example

These entries belong to the mesh example `mesh-1`, which the project received
with the code. `chdr.mesh.geometry` checks the geometry predicates and the YAML
input without an IPPL mesh; `verify_mesh_output.py` checks each `mesh-1`
output against an independent Cartesian description of the radiator and
compares the two rank counts. Run them with
`ctest --test-dir build -R '^chdr\.mesh\.' --output-on-failure`.

| Entry | Source (lines) · suite | What it checks | Reference | Labels | Time (s) |
|---|---|---|---|---|---|
| `chdr.mesh.geometry` | `test_geometry.cpp` (322) | Prism and brick volumes to `1e-13`, point classification on the host and in a Kokkos reduction, refusal of degenerate shapes, and the parsing of `config.yaml`. | continuum analytic (volumes); refusal | — | 0.5 |
| `chdr.mesh.{serial,mpi2}` | `mesh-1` with `config.yaml` | Builds the prism radiator on one and on two ranks and writes the output that `chdr.mesh.compare` reads. | smoke (fixture setup) | `mpi` | 1.6 / 1.7 |
| `chdr.mesh.compare` | `verify_mesh_output.py` (392) | Cell values and the dielectric count match an independent description of the prism, and the two rank counts agree exactly in counts and checksum. | one against two ranks; exact lattice value | — | 0.2 |
| `chdr.mesh.brick.{serial,mpi2}` | `mesh-1` with `config-brick.yaml` | Builds the brick radiator on one and on two ranks. | smoke (fixture setup) | `mpi` | 1.6 / 1.7 |
| `chdr.mesh.brick.compare` | `verify_mesh_output.py` | The same checks for the brick. | one against two ranks; exact lattice value | — | 0.2 |

## `1d`: the one-dimensional prototype

`chdr.1d.unittest` runs the 13 `unittest` methods of `../1d/test_chdr_1d.py`
against the one-dimensional spectral prototype `../1d/chdr_1d.py`, which the
project received with the code. The entry needs a Python with NumPy, SciPy and
Matplotlib (`CHDR_1D_PYTHON`); without them it reports a skip. Earlier results
depend on the exact bytes of both files, so change them only on purpose. Run
it with `ctest --test-dir build -R '^chdr\.1d\.' --output-on-failure`.

| Entry | Source (lines) · suite | What it checks | Reference | Labels | Time (s) |
|---|---|---|---|---|---|
| `chdr.1d.unittest` | `test_chdr_1d.py` (199) · `HalfSpaceTests` (13 methods) | Interface conditions and Maxwell's equations per mode to `2e-11`, radiated flux against the work on the electron to `2e-9`, vacuum fields against the Bessel solution to `2e-7`, zero radiation below threshold, and parameter refusals. | continuum analytic; refusal | `unit` | 1.2 |

## Citation keys

The keys below are cited in the test files and their support headers. Each
entry gives the reference and the files that cite it; the files name the
section or equation they use.

- `bamberger1987SecondOrderCorner`: A. Bamberger, P. Joly and J. E. Roberts
  (1987), *Second Order Absorbing Boundary Conditions for the Wave Equation: A
  Solution for the Corner Problem*, <https://inria.hal.science/inria-00075909>.
  Cited in `src/MaxwellSolvers/BoundaryOracles.h`, `core/solver/ReflectionOracles.h`,
  `core/solver/OpenRunSupport.h`, `core/solver/test_mur_families.cpp`.
- `chew2014generalizedgauge`: W. C. Chew (2014), *Vector Potential
  Electromagnetic Theory with Generalized Gauge for Inhomogeneous Anisotropic
  Media*, [arXiv:1406.4780](https://arxiv.org/abs/1406.4780). Cited in
  `core/material/OwnerOracles.h`, `apps/planar-interface/InterfaceOracles.h`,
  `opt-in/diagnostics/test_interface_energy.cpp`,
  `core/solver/formulation/test_o1vs_regression.cpp` and ten files of `core/solver`.
- `christlieb2024GaugeConserving`: A. J. Christlieb, W. A. Sands and S. R.
  White (2025), *A Particle-in-Cell Method for Plasmas with a Generalized
  Momentum Formulation, Part III: A Family of Gauge Conserving Methods*,
  Journal of Scientific Computing 104(2), 38,
  [doi:10.1007/s10915-025-02953-7](https://doi.org/10.1007/s10915-025-02953-7);
  cited version [arXiv:2410.18414v1](https://arxiv.org/abs/2410.18414). Cited in
  `opt-in/diagnostics/test_monitors.cpp`, `test_monitor_records.cpp` and eight
  files of `core/solver`.
- `clemens2001FiniteIntegration`: M. Clemens and T. Weiland (2001), *Discrete
  Electromagnetism with the Finite Integration Technique*, Progress In
  Electromagnetics Research 32, 65–87,
  [doi:10.2528/PIER00080103](https://doi.org/10.2528/PIER00080103). Cited in
  `core/solver/test_staggered.cpp`, `core/solver/MonitorTestSupport.h`,
  `opt-in/diagnostics/test_monitors.cpp`.
- `ellingson2020ElectromagneticsVol2`: S. W. Ellingson (2020),
  *Electromagnetics, Vol. 2*, Virginia Tech Publishing,
  [doi:10.21061/electromagnetics-vol-2](https://doi.org/10.21061/electromagnetics-vol-2).
  Cited in `apps/planar-interface/InterfaceOracles.h`.
- `engquist1977AbsorbingBoundaryWaves`: B. Engquist and A. Majda (1977),
  *Absorbing Boundary Conditions for the Numerical Simulation of Waves*,
  Mathematics of Computation 31(139), 629–651,
  [doi:10.1090/S0025-5718-1977-0436612-4](https://doi.org/10.1090/S0025-5718-1977-0436612-4).
  Cited in `src/MaxwellSolvers/BoundaryOracles.h`, `core/solver/ReflectionOracles.h`.
- `eymard2000FiniteVolume`: R. Eymard, T. Gallouët and R. Herbin (2000),
  *Finite Volume Methods*, in Handbook of Numerical Analysis, vol. 7, Elsevier,
  713–1020, doi:10.1016/S1570-8659(00)07005-8. Cited in
  `core/material/OwnerOracles.h`, `apps/planar-interface/InterfaceOracles.h`.
- `fallahi2020mithra20fullwavesimulation`: A. Fallahi (2020), *MITHRA 2.0: A
  Full-Wave Simulation Tool for Free Electron Lasers*,
  [arXiv:2009.13645](https://arxiv.org/abs/2009.13645). Cited in
  `src/MaxwellSolvers/OuterLayerOracles.h`, `apps/planar-interface/InterfaceOracles.h`,
  `opt-in/diagnostics/test_probe_full.cpp`, seven files of `apps/vacuum-wave` and
  twenty files of `core/solver`.
- `fitzpatrick2014ClassicalElectromagnetism`: R. Fitzpatrick (2014), *Classical
  Electromagnetism*, <https://farside.ph.utexas.edu/teaching/jk1/jk1.html>.
  Cited in `core/input/check_runs.py`, `core/input/test_input.cpp`,
  `core/solver/test_comoving.cpp`, `core/solver/ComovingPotentials.h`.
- `fornberg1988FiniteDifferenceFormulas`: B. Fornberg (1988), *Generation of
  Finite Difference Formulas on Arbitrarily Spaced Grids*,
  [doi:10.1090/S0025-5718-1988-0935077-0](https://doi.org/10.1090/S0025-5718-1988-0935077-0).
  Cited in `src/MaxwellSolvers/OuterLayerOracles.h`, `core/solver/test_mur_outer_layer.cpp`.
- `gedney2011IntroductionFdtd`: S. D. Gedney (2011), *Introduction to the
  Finite-Difference Time-Domain (FDTD) Method for Electromagnetics*,
  [doi:10.1007/978-3-031-01712-4](https://doi.org/10.1007/978-3-031-01712-4).
  Cited in `core/solver/MediumOracles.h`.
- `gustafsson1972StabilityTheoryII`: B. Gustafsson, H.-O. Kreiss and A.
  Sundström (1972), *Stability Theory of Difference Approximations for Mixed
  Initial Boundary Value Problems. II*, Mathematics of Computation 26(119),
  649–686, [doi:10.1090/S0025-5718-1972-0341888-3](https://doi.org/10.1090/S0025-5718-1972-0341888-3).
  Cited in `apps/planar-interface/InterfaceOracles.h`.
- `heinzel2002SpectrumSpectralDensity`: G. Heinzel, A. Rüdiger and R.
  Schilling (2002), *Spectrum and spectral density estimation by the Discrete
  Fourier transform (DFT), including a comprehensive list of window functions
  and some new flat-top windows*, Max-Planck-Institut für Gravitationsphysik,
  [hdl:11858/00-001M-0000-0013-557A-5](https://hdl.handle.net/11858/00-001M-0000-0013-557A-5).
  Cited in `opt-in/diagnostics/test_probe.cpp`, `test_probe_full.cpp`,
  `core/solver/test_comoving.cpp`, `test_convergence.cpp`, `test_gauge_growth.cpp`.
- `higdon1986AbsorbingDifferenceMultidim`: R. L. Higdon (1986), *Absorbing
  Boundary Conditions for Difference Approximations to the Multi-Dimensional
  Wave Equation*, Mathematics of Computation 47, 437–459,
  [doi:10.1090/S0025-5718-1986-0856696-4](https://doi.org/10.1090/S0025-5718-1986-0856696-4).
  Cited in `core/solver/ReflectionOracles.h`.
- `jeannerod2013ImprovedErrorBounds`: C.-P. Jeannerod and S. M. Rump (2013),
  *Improved Error Bounds for Inner Products in Floating-Point Arithmetic*,
  [doi:10.1137/120894488](https://doi.org/10.1137/120894488). Cited in
  `core/solver/OwnedStatistics.h`.
- `jones1988MultiPointImpedance`: M. G. Jones and T. L. Parrott (1988),
  *Evaluation of a Multi-Point Method for Determining Acoustic Impedance*, NASA
  TM-100637, [NTRS 19880017529](https://ntrs.nasa.gov/citations/19880017529).
  Cited in `core/solver/test_open_beam.cpp`, `core/solver/ReflectionRun.h`,
  `core/solver/TwoWaveState.h`.
- `kreiss2002SecondOrderWave`: H.-O. Kreiss, N. A. Petersson and J. Yström
  (2002), *Difference Approximations for the Second Order Wave Equation*, SIAM
  J. Numer. Anal. 40, 1940–1967,
  [doi:10.1137/S0036142901397435](https://doi.org/10.1137/S0036142901397435).
  Cited in `core/solver/MediumOracles.h`, `core/solver/RandomFields.h`,
  `core/solver/test_mur_stability.cpp`, `opt-in/diagnostics/test_interface_energy.cpp`.
- `mohr2025Codata2022`: P. J. Mohr, D. B. Newell, B. N. Taylor and E. Tiesinga
  (2025), *CODATA Recommended Values of the Fundamental Physical Constants:
  2022*, Reviews of Modern Physics 97(2), 025002,
  [doi:10.1103/RevModPhys.97.025002](https://doi.org/10.1103/RevModPhys.97.025002).
  Cited in `core/input/check_runs.py`, `core/input/test_input.cpp`.
- `mur1981AbsorbingBoundaryFdtd`: G. Mur (1981), *Absorbing Boundary
  Conditions for the Finite-Difference Approximation of the Time-Domain
  Electromagnetic-Field Equations*,
  [doi:10.1109/TEMC.1981.303970](https://doi.org/10.1109/TEMC.1981.303970).
  Cited in `src/MaxwellSolvers/BoundaryOracles.h`, `src/MaxwellSolvers/HaloOrder.cpp`,
  `core/material/OwnerOracles.h`, `core/solver/ReflectionOracles.h`,
  `core/solver/OpenRunSupport.h`, `core/solver/test_mur_faces.cpp`,
  `core/solver/test_mur_refinement.cpp`.
- `oskooi2010meep`: A. F. Oskooi, D. Roundy, M. Ibanescu, P. Bermel, J. D.
  Joannopoulos and S. G. Johnson (2010), *Meep: A flexible free-software package
  for electromagnetic simulations by the FDTD method*, Computer Physics
  Communications 181, 687–702,
  [doi:10.1016/j.cpc.2009.11.008](https://doi.org/10.1016/j.cpc.2009.11.008).
  Cited in `core/solver/test_open_beam.cpp`, `opt-in/diagnostics/test_probe_full.cpp`,
  `src/MaxwellSolvers/HaloOrder.cpp`.
- `ryu2016potentialfdtd`: C. J. Ryu, A. Y. Liu, W. E. I. Sha and W. C. Chew
  (2016), *Finite-Difference Time-Domain Simulation of the Maxwell–Schrödinger
  System*, IEEE Journal on Multiscale and Multiphysics Computational Techniques
  1, 40–47, [doi:10.1109/JMMCT.2016.2605378](https://doi.org/10.1109/JMMCT.2016.2605378).
  Cited in `src/MaxwellSolvers/OuterLayerOracles.h`, `core/material/OwnerOracles.h`,
  `apps/planar-interface/InterfaceOracles.h`, `opt-in/diagnostics/test_probe_full.cpp`
  and six files of `core/solver`.
- `sakkaplangkul2021YeeDispersiveMedia`: P. Sakkaplangkul and V. A. Bokil
  (2021), *Convergence Analysis of Yee-FDTD Schemes for 3D Maxwell's Equations
  in Linear Dispersive Media*, International Journal of Numerical Analysis and
  Modeling 18(4), 524–568,
  [doi:10.4208/ijnam.v18.n4.p524](https://doi.org/10.4208/ijnam.v18.n4.p524).
  Cited in `apps/planar-interface/test_plane_interface.cpp`,
  `apps/interface-movies/test_movie_start.cpp`.
- `salari2000CodeVerificationMms`: K. Salari and P. Knupp (2000), *Code
  Verification by the Method of Manufactured Solutions*, Sandia National
  Laboratories, SAND2000-1444, [doi:10.2172/759450](https://doi.org/10.2172/759450).
  Cited in `core/solver/test_manufactured.cpp`, `core/solver/Convergence.h`.
- `trefethen1982GroupVelocity`: L. N. Trefethen (1982), *Group Velocity in
  Finite Difference Schemes*, SIAM Review 24(2), 113–136,
  [doi:10.1137/1024038](https://doi.org/10.1137/1024038). Cited in
  `core/solver/ReflectionOracles.h`, `apps/planar-interface/InterfaceOracles.h`.
- `trefethen1996FiniteDifferenceSpectral`: L. N. Trefethen (1996), *Finite
  Difference and Spectral Methods for Ordinary and Partial Differential
  Equations*, unpublished text,
  <https://people.maths.ox.ac.uk/trefethen/pdetext.html>. Cited in
  `core/solver/test_operators.cpp`, `core/solver/test_mur_reflection.cpp`.
- `yee1966MaxwellIsotropicMedia`: K. S. Yee (1966), *Numerical Solution of
  Initial Boundary Value Problems Involving Maxwell's Equations in Isotropic
  Media*, IEEE Transactions on Antennas and Propagation 14(3), 302–307,
  [doi:10.1109/TAP.1966.1138693](https://doi.org/10.1109/TAP.1966.1138693).
  Cited in `apps/planar-interface/InterfaceOracles.h`.
