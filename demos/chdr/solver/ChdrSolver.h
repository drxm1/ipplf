/** @file ChdrSolver.h
 * @brief Periodic demo-level shell; wave updates remain in IPPL's base solvers.
 */
#pragma once

// Solver internals use the communication declarations provided by the IPPL umbrella.
#include "Ippl.h"

#include <type_traits>

#include "FELFieldContainer.hpp"
#include "MaxwellSolvers/NonStandardFDTDSolver.h"
#include "VacuumSetup.h"

namespace chdr::solver {

    /// Dimension of a cauchy slice.
    /// Must match the amount of components that \f$\vec{E}, \vec{B}, \vec{A}, \vec{J}\f$ have.
    constexpr unsigned SpatialDim = 3;

    /// Caller-owned mesh, layout, source and output storage container.
    using ChDRFieldContainer = FELFieldContainer<double, SpatialDim>;

    /// Three-component electric or magnetic field storage typeon an IPPL mesh.
    using EMFieldsType = VField_t<double, SpatialDim>;

    /// Four-component scalar/vector potential or charge/current storage type.
    /// Three spacial components (A or J) + one temporal component (\f$\phi\f$ or \f$\rho\f$).
    using AJFourFieldsType = SourceField_t<double, SpatialDim>;

    /// Periodic Standard FDTD Solver, like in the vacuum example.
    /// Later on we will enable different boundary conditions,
    /// maybe Mur or PML_ADI.
    using StandardSolverPeriodic =
        ippl::StandardFDTDSolver<EMFieldsType, AJFourFieldsType, ippl::periodic>;

    /// Periodic longitudinally(===z)-smoothed NonStandard vacuum update.
    /// Later on we will enable different boundary conditions,
    /// maybe Mur or PML_ADI.
    using NonStandardSolverPeriodic =
        ippl::NonStandardFDTDSolver<EMFieldsType, AJFourFieldsType, ippl::periodic>;

}  // namespace chdr::solver
