# Periodic vacuum wave: first step

This small laboratory-frame demo initializes a source-free sinusoidal wave,
reconstructs its electric and magnetic fields, advances **one timestep** with
IPPL's `StandardFDTDSolver`, and exits. A normal successful run is silent and
writes no field files. Numerical checks live in the companion tests.

## Build and run

Use an existing IPPL build configured with your compiler, MPI and Kokkos
dependencies; see the [IPPL build documentation](../../../README.md#resources).
From the **IPPL source root**, replacing `build` with your build directory:

```sh
cmake -S . -B build -DIPPL_ENABLE_CHDR=ON -DIPPL_ENABLE_SOLVERS=ON \
  -DIPPL_ENABLE_UNIT_TESTS=ON -DIPPL_ENABLE_FFT=ON \
  -DIPPL_USE_STANDARD_FOLDERS=ON
cmake --build build --target chdr-vacuum-wave test-chdr-vacuum --parallel 2
OMP_NUM_THREADS=1 OMP_PROC_BIND=false build/bin/chdr-vacuum-wave
ctest --test-dir build -R '^chdr\.vacuum\.' \
  --no-tests=error --output-on-failure
```

The unit-test option requires GoogleTest and enables testing; IPPL's existing
unit-test configuration also requires FFT support. Enabling CHDR additionally
configures the mesh demo's Catalyst/Conduit dependency; the FEL executable need
not be enabled. With `IPPL_USE_STANDARD_FOLDERS=OFF`, the wave executable is
`build/demos/chdr/vacuum/chdr-vacuum-wave` instead of `build/bin/chdr-vacuum-wave`.

## Fixed setup and field times

The grid has **8 × 8 × 64 cell-centered cells** in
`[0, 1/8] × [0, 1/8] × [0, 1]`, periodic on every face and decomposed along z.
Light speed and wavelength are normalized to one. The continuum reference
travels along +z, with electric polarization along x and magnetic polarization
along y; see [the analytic potential](VacuumWave.h). Its starting relations are
Fallahi, [*MITHRA 2.0*, arXiv:2009.13645v1](https://arxiv.org/abs/2009.13645v1),
“Wave Equation” and “FDTD for Wave Equation,” equations (3.6)–(3.15).

The [Standard solver](../../../src/MaxwellSolvers/StandardFDTDSolver.hpp) selects
`dt = min(spacing)/2`, here `1/128`. We sample the exact wave at `0` and `-dt`
to supply the two histories required by MITHRA's three-level update (3.14)–(3.15).
For this zero-scalar-potential case, the raw fields have different times:

| State | Newest potential | Raw E | Raw B |
| --- | --- | --- | --- |
| After initialization/reconstruction | `0` | `-dt/2` | `0` |
| After one `solve()` | `dt` | `dt/2` | `dt` |

These timestamps follow IPPL's actual
[`evaluate_EB()` reconstruction](../../../src/MaxwellSolvers/FDTDSolverBase.hpp),
using the potential identities in MITHRA (3.8)–(3.9). E differences two histories;
B curls the newest one. The initial negative E time is intentional.

## Tests and source reading

The **14 GoogleTests** run serially and, when an MPI launcher is available, on
two ranks. They check initial histories/E/B, one update and exact history copies,
periodic face halos, shifted coordinates and unequal spacing, domain restrictions,
and debug guards for invalid halo configurations. They compare against independent
references, check zero sources/inactive components and reject nonfinite values.
[The first-step tests](tests/test_first_step.cpp) define the regression thresholds:
`1e-10` normalized discrete-reference error, `1e-12` inactive-component error,
1% relative continuum E/B error, and exact source-zero/history-copy checks.
These are chosen requirements for this fixed grid. Applying the implemented
reconstruction to the reference wave gives initial relative errors about `0.0016`
for B and `0.00010` for E, leaving margin below 1% (our calculation from MITHRA
(3.8)–(3.9) and `evaluate_EB()` cited above). The discrete tolerance tests agreement
with the update independently of this continuum truncation error.
These tests cover initialization and one step; longer propagation and field output
are later work.

CTest also runs the actual driver as serial/MPI2 smoke checks and repeats the
MPI2 unit suite with shuffled order. All ranks must receive the same test-selection,
repetition and shuffle options. Unseeded shuffling uses seed `1` on every rank;
explicit nonzero seeds are preserved. Halo validation remains debug-only, and the
test target explicitly enables those checks.

Start with [the driver](vacuum-wave.cpp), then [shared setup](VacuumSetup.cpp).
The reused [field container](../../fel/FELFieldContainer.hpp) owns mesh/layout
and E/B/J storage and must outlive the solver. Mathematical assumptions and
derivations are documented beside the setup and test helpers.
