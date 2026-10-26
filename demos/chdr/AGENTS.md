# AGENTS.md

## Repository purpose
- IPPL is a C++ Particle and fields Framework. Preserve physical correctness over stylistic cleanup.
- In this directory we model Cherenkov Diffraction Radiation. The solver extends IPPL's Standard and NonStandard FDTD potential solvers and reuses the field container of ../fel.
- Add meaningful unit tests for classes, and for functions where it is reasonable to do so. One directory under `tests/` per unit of the code (README.md maps them).
- Use Doxygen extensively, also to document the mathematics and physics: inline `\f$ ... \f$`, displayed `\f[ ... \f]`, standard LaTeX only, and cite the library key of every non-trivial equation or argument, as the headers do. Never just guess where your equation or argument might have come from or infer from second-hand notes. Instead double-check the sources latex and pdf files directly and ensure that the equations match exactly. If the match is not exact, state it clearly as a deduction and state how to deduce from the citation. If that deduction requires too many steps, check if there is a better source to cite, e.g. one that is closer to the statement you need to support.

## Build
- Configure with CMake using the project’s normal toolchain and dependency prefixes. `IPPL_ENABLE_CHDR=ON` adds this directory. The commands are in README.md.
- Build outside of source, in a directory under a canonical `build/` dir.
- Run the relevant test subset after changes (`ctest -R '^chdr\.'`, see also README.md for a list). Do not run test sets that are irrelevant to the changes you made, I mean run tests lazily - only when needed.
- External dependencies are in the _deps directory fetched with cmake (the exact setup might differ slightly in my master's projects nix environment, doesn't really matter).

## External dependencies
- heFFte and Kokkos are external dependencies.
- Prefer changing IPPL adapter/wrapper code over patching upstream dependencies.
- When a change touches parallel kernels, consider thoroughly and explain when necessary the execution space, memory space, and data movement implications.

## Physics / numerics rules
- For algorithmic changes, consider thoroughly and state expected impact on conservation, stability, and reproducibility.
- Flag any change that may alter floating-point behavior or MPI/GPU execution order.
- Add or update (usually at least one) regression/sanity test for physics-facing changes, but in a reasonable way. Do not add tests that provide no value.
- Keep an eye on parallel efficiency at least on the openmp level

## Layout and standards of this directory
- README.md maps the directory: `core/` (the code that no experiment, test or optional measurement owns), `opt-in/` (measurements and output that nothing in `core/` depends on), `apps/` (the programs), `tests/` (one directory per unit of the code it tests), `testing/` (shared test support).
- Include every header by its path from this directory, for example `#include "core/solver/ChdrSolver.h"`. Bare includes do not compile.
- No file under `core/` includes a file from `apps/`, `opt-in/`, `tests/` or `testing/`.
- Name exactly what is meant: the executable, the file, the class, the quantity. Avoid general words such as "the driver", "demo" or "the shell" for a named thing.
- Internal identifiers that you are using belong to the project records, not to code, comments, CMake files, READMEs or anywhere else in ippl's source control.
- Prefer functions of 20 lines or fewer; extract coherent helpers from longer routines. Do not sacrifice good code-comments or physical correctness to achieve under 20 LOC functions. Functions that call kokkos loops/reductions should be under 20 LOC if the kokkos reduction call's lambda function's lines are not counted, and any other function must be 20 LOC or less in almost any case.
- Write READMEs and comments in precise scientific language and keep sentences well readable (in the spirit of ASD-STE100).
- The pre-commit hook runs clang-format and cmake-format; keep every file clean under both.

# Naming Conventions
- Variables should use camel casing
- Compile time constants should use capital casing
- Member variables must be suffixed with `_m`

# Math

Mathematical constants should be obtained from `Kokkos::numbers` (moved out of `experimental` in Kokkos 4). Any instance of a mathematical function in host-only code should use symbols from the standard library, e.g. `std::sqrt`. Any instances occurring in device-code or code that might be run on a device (such as those marked with `KOKKOS_INLINE_FUNCTION`) should use the symbols from Kokkos, e.g. `Kokkos::sqrt`, to ensure performance and portability on GPUs.
