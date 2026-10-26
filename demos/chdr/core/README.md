# ChDR core

The core is the part of the ChDR code that no experiment, test or optional
measurement owns: the shared field container and its sampling contract, the
material description, the discrete operators, the material formulations, the
potential solver, the simulation manager and the field writer. No file under
`core/` includes a file from `apps/`, `opt-in/`, `tests/` or `testing/`.

Include a core header by its path from `demos/chdr`, for example
`#include "core/solver/ChdrSolver.h"`. The CMake target `chdr-core` is a static
library of `base/FieldContainer.cpp` and `solver/BoxSetup.cpp`; it exports the
`demos/chdr` root as include directory and links `ippl`, `chdr-material` and the
FEL field-container headers. Bare includes such as `#include "ChdrSolver.h"` do
not compile.

| Unit | Headers | Content | Tests |
|---|---|---|---|
| `base/` | `FieldContainer`, `SpatialDimension`, `HalfStepOffset`, `LevelContract`, `Util` | The field container shared with the FEL demo, its ownership check and storage allocation; the spatial dimension; the space and time sampling contract of every solver quantity; the abort helper `chdr::core::util::abortFailure` | `Levels`, `LevelEquivalence` in `tests/core/solver` |
| `material/` | `MaterialOwner`, `O1VSCoefficients`, `PlanarSlab`, `RectangleXZ` | The distributed cell permittivity with its derived face coefficients, and the two material regions | `tests/core/material` |
| `operators/` | `Operators`, `StaggeredOperators`, `Reconstruction`, `StaggeredReconstruction` | Collocated and staggered differences, and the reconstruction of E and B from the potential histories in each reading | `Operators`, `CompactOperators`, `ChdrStaggered` in `tests/core/solver` |
| `formulation/` | `MaterialFormulation`, `MaterialStepContext`, `HomogeneousUpdate`, `O1VSOperator`, `O1VSUpdate` | The material update that the solver selects: the homogeneous rescaling of IPPL's vacuum kernel and the O1-VS variable-scalar gauge | `tests/core/formulation` |
| `solver/` | `ChdrSolver`, `BoxSetup` | The potential solver over IPPL's Standard and NonStandard FDTD solvers with periodic or Mur faces, and the constructor of a box field container with selectable periodic axes | `tests/core/solver` |
| `simulation/` | `ChDRSimulationManager`, `SimulationConfig` | One lifecycle for a run: box, material map, solver settings, stepping with observers | `tests/core/simulation` |
| `output/` | `FieldOutput` | Synchronous owned-cell VTI/PVTI output of fields and material | The output fixture and reader checks in `tests/apps/vacuum-wave` |

Both spatial readings, collocated and staggered, are kept on purpose; they are
two discretizations, not two copies of one. `ChDRSimulationManager` builds its
field container with `solver::makeBoxFieldContainer`, so the box constructor is
the one place that defines mesh and MPI ownership.

[`solver/README.md`](solver/README.md) documents the solver in detail: the
readings and time labels, the homogeneous material, the periodic and Mur faces,
the difference operators, the time step, units, and the solver tests. The
residual monitors, the Fourier probe and the field energies are optional
measurements under [`../opt-in/diagnostics`](../opt-in/diagnostics/README.md).
The programs that use the core are under `../apps`.
