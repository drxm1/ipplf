# ChDR: Cherenkov diffraction radiation in IPPL

This directory holds the ChDR code of my master's project:
- a potential-based FDTD solver (wrapper) for dielectric media on IPPL's fields in `core/solver/`,
- some optional measurement utilities in `opt-in/diagnostics/`,
- some utilities to write vtk outputs (`core/output/`),
- currently three little programs under `apps`:
  - `vacuum-wave`: test evolving a vacuum wave for one period,
  - `planar-interface`: a few simple interface cases,
  - `interface-movies`: some videos of ring/annulus-shaped initialized fields hitting a few different dielectric geometries,
- tests (to be organized and documented better),
- the mesh geometry for comparison against the meep simulation under `mesh/` that was provided to me,
- a one-dimensional Python code for later 1d comparisons (`1d/`) and a literature review (`ideas-lit-research/`) which were also provided to me.

`demos/CMakeLists.txt` adds this directory with `IPPL_ENABLE_CHDR=ON`; the tests need `BUILD_TESTING` and `IPPL_ENABLE_UNIT_TESTS`.

## Layout

| Directory | Content | Documentation |
|---|---|---|
| `core/` | The code that no experiment, test or optional measurement owns, one unit per directory: the base shared by everything, the field container, the sampling contract and the abort helper (`base/`); the potential solver `ChdrSolver` with the box constructor and its own parts, the discrete operators (`solver/operators/`) and the material formulations (`solver/formulation/`); the material description (`material/`); the simulation manager (`simulation/`); the field writer (`output/`) | [`core/README.md`](core/README.md), [`core/solver/README.md`](core/solver/README.md) |
| `opt-in/diagnostics/` | The residual monitors, the Gauss formula check, the Fourier probe and the two field energies; nothing in `core/` depends on them | [`opt-in/diagnostics/README.md`](opt-in/diagnostics/README.md) |
| `opt-in/output/` | `FieldRunOutput.h`: frames, material and `energy.csv` for a run of the simulation manager | The header |
| `apps/vacuum-wave/` | `chdr-vacuum-wave`: the vacuum plane wave in a periodic box, the benchmark of the inherited vacuum update | [`apps/vacuum-wave/README.md`](apps/vacuum-wave/README.md), [`ViewingGuide.org`](apps/vacuum-wave/ViewingGuide.org) |
| `apps/planar-interface/` | `planar-interface`: one planar dielectric in a periodic box; the plane-pulse reflection measurement and the packet cases | [`apps/planar-interface/README.md`](apps/planar-interface/README.md) |
| `apps/interface-movies/` | `interface-movies`: an open box with Mur faces, ring starts and eight movie cases; `render_movie.py`, which renders the output of both interface programs | [`apps/interface-movies/README.md`](apps/interface-movies/README.md) |
| `tests/` | One directory per unit of the code it tests: `core/<unit>`, `opt-in/diagnostics`, `apps/<program>`, `1d`, and `src/<component>` for the IPPL library changes made for ChDR (per-axis periodic layouts in `src/FieldLayout`, Mur faces in `src/MaxwellSolvers`) | The section Tests below; `tests/src/MaxwellSolvers/HaloOrder.md` and `OpenFaceClosure.md` |
| `testing/` | Shared test support: the GoogleTest main `TestMain.cpp` and `Util.h`, included as `testing/Util.h` | |
| `mesh/` | `mesh-1`: the supervisor's distributed mesh and dielectric radiator example with Catalyst output, and its geometry test | [`mesh/README.md`](mesh/README.md), [`mesh/StudentGuide.md`](mesh/StudentGuide.md) |
| `1d/` | The one-dimensional spectral prototype in Python: the field and radiated spectrum of a point electron near an infinite planar interface, by analytic matching per transverse mode | [`1d/README.md`](1d/README.md) |
| `ideas-lit-reseach/` | The literature review `chdr_literature_review.tex` and reading notes | |

## Includes and targets

Every header is included by its path from this directory, for example
`#include "core/solver/ChdrSolver.h"`, `#include "opt-in/diagnostics/Monitors.h"`,
`#include "apps/vacuum-wave/VacuumSetup.h"` or
`#include "tests/core/solver/CoreTestSupport.h"`. A bare include such as
`#include "ChdrSolver.h"` does not compile. The CMake targets export this
directory as their include root:

| Target | Kind | Content |
|---|---|---|
| `chdr-core` | static library | `core/`; links `ippl`, `chdr-material` and `chdr-fel-headers` |
| `chdr-material` | interface | `core/material/` (header only) |
| `chdr-fel-headers` | interface | The field container headers of `demos/fel` |
| `chdr-diagnostics` | interface | `opt-in/diagnostics/`; links `chdr-core` |
| `chdr-run-output` | interface | `opt-in/output/` |
| `chdr-testing` | interface | `testing/`; compiles the GoogleTest main into each test executable |
| `chdr-vacuum-wave`, `planar-interface`, `interface-movies`, `mesh-1` | executables | The programs |

The namespaces follow the directories. `core/base/` declares its names in
`chdr` (`Util.h` in `chdr::util`); `core/solver/` with its subdirectories in
`chdr::solver`; `core/material/` in `chdr::material`; `core/simulation/` in
`chdr::simulation`; `core/output/` and `opt-in/output/` in `chdr::output`;
`opt-in/diagnostics/` in `chdr::diagnostics`; the programs in
`chdr::apps::vacuumWave`, `chdr::apps::planarInterface` and
`chdr::apps::interfaceMovies`. The test support headers of `tests/core/solver`
and `tests/core/material` use `chdr::solver::test` and `chdr::material::test`
(`CoreTestSupport.h` uses `chdr::solver::coreTest`), the tests of a program
use the program's namespace, `tests/opt-in/diagnostics/test_mur_layer.cpp`
uses `chdr::diagnostics::test`, and `testing/Util.h` uses
`chdr::testing::util`.

No file under `core/` includes a file from `apps/`, `opt-in/`, `tests/` or
`testing/`.

## Build and tests

From the IPPL checkout, with the usual toolchain and dependency options:

```sh
cmake -S . -B build -DIPPL_ENABLE_CHDR=ON -DIPPL_ENABLE_FEL=ON \
  -DIPPL_ENABLE_FFT=ON -DIPPL_ENABLE_SOLVERS=ON \
  -DBUILD_TESTING=ON -DIPPL_ENABLE_UNIT_TESTS=ON
cmake --build build --parallel 2
ctest --test-dir build -R '^chdr\.' -LE 'long|record' --output-on-failure
```

The CTest entries are named `chdr.<unit>.<suite>[.<variant>]`, where the unit is
the directory under `tests/` and the variant is `serial` or `mpi2` (two MPI
ranks, one OpenMP thread each):

| Prefix | Tests of | Directory |
|---|---|---|
| `chdr.solver.` | the potential solver, the operators and the sampling contract | `tests/core/solver` |
| `chdr.solver.formulation.` | the O1-VS formulation | `tests/core/solver/formulation` |
| `chdr.material.`, `chdr.simulation.` | the material owner, the simulation manager | `tests/core/<unit>` |
| `chdr.diagnostics.` | the opt-in diagnostics | `tests/opt-in/diagnostics` |
| `chdr.vacuum.`, `chdr.interface.` | the three programs; `chdr.interface.movies.*` belong to `interface-movies` | `tests/apps/<program>` |
| `chdr.mesh.` | the mesh example | `mesh/` |
| `chdr.1d.` | the one-dimensional prototype (`unittest`; needs a Python with NumPy, SciPy and Matplotlib, `CHDR_1D_PYTHON`) | `tests/1d` |
| `chdr.src.` | the IPPL library changes made for ChDR: per-axis periodic layouts (`PeriodicAxes`) and the Mur faces (`FDTDBoundary`, `MurSpeed`, `MurDispatch`, `InteriorBounds`, `HaloOrder`, `OpenFaceClosure`, `OpenBoundaryRegression`) | `tests/src/<component>` |

Labels: `long` marks ladders and trajectories above one minute, `record`
marks entries that only print records, `holdout` marks the disabled holdout
entries, and `unavailable` marks entries registered as skipped because Python,
an MPI launcher or ParaView is missing. Every entry carries the label `chdr`,
so `ctest -L chdr` selects the same set as `ctest -R '^chdr\.'`.

## Documentation

`doxygen Doxyfile` in this directory writes the HTML reference of `core/`,
`opt-in/`, `apps/` and `testing/`, with this file and the READMEs as pages, to
`docs/html/index.html`. It needs Doxygen, Graphviz, and LaTeX with Ghostscript
for the formulas. The tests are not included; the comment at the top of
`Doxyfile` says how to add them. `mesh/` has its own `Doxyfile`.
