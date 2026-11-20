/** @file SolverConfiguration.h
 * @brief The solver configurations that the pulse experiments can select.
 */
#pragma once

#include <array>
#include <stdexcept>
#include <string>

#include "core/base/LevelContract.h"
#include "core/input/RunInput.h"
#include "core/solver/formulation/MaterialFormulation.h"

namespace chdr::experiments {
    /** @brief Material formulation, spatial placement and E/B reconstruction of one run.
     * @details E01 to E04 use O1VS. Three further configurations permit a comparison with it:
     * the Homogeneous formulation on the same staggered layout, and the Homogeneous
     * formulation on the collocated layout with the R1 or the R2Prime time levels.
     * In the collocated ones ChdrSolver runs the inherited step(), timeShift() and
     * evaluate_EB() of IPPL's StandardFDTDSolver, with the time step
     * \f$h_t/\sqrt{\varepsilon_r}\f$ in the update. R1 loads \f$\phi\f$ at the time of
     * \f$\mathbf A\f$, as the authors of the solver define the time levels
     * [winkler2024ChargeConservingFDTD, Eqs. (3.1)-(3.4), pp. 16-17]. R2Prime loads
     * \f$\phi\f$ half a step earlier, so that evaluate_EB() gives \f$\mathbf E\f$ at the
     * single time \f$t_{n-1/2}\f$; ChdrSolver then averages \f$\mathbf B\f$ over both slots
     * (LevelContract).
     */
    struct SolverConfiguration {
        /// Name in the inputs and the run records: "O1VS", "Homogeneous", "Standard" or
        /// "StandardR2Prime".
        const char* name_m;
        solver::MaterialFormulation formulation_m;
        SpatialDiscretization spatial_m;
        EBReconstructionConvention reconstruction_m;

        /// Whether E and B use the staggered positions of the R2 contract.
        bool staggered() const { return spatial_m == SpatialDiscretization::Staggered; }
    };

    /// The configuration of E01 to E04: O1VS update, staggered layout, R2.
    inline constexpr SolverConfiguration O1VSConfiguration{
        "O1VS", solver::MaterialFormulation::O1VS, SpatialDiscretization::Staggered,
        EBReconstructionConvention::R2};
    /// The Homogeneous update (scaled step) on the same staggered layout with R2.
    inline constexpr SolverConfiguration HomogeneousConfiguration{
        "Homogeneous", solver::MaterialFormulation::Homogeneous, SpatialDiscretization::Staggered,
        EBReconstructionConvention::R2};
    /// IPPL's StandardFDTDSolver path: Homogeneous update, collocated layout, R1.
    inline constexpr SolverConfiguration StandardConfiguration{
        "Standard", solver::MaterialFormulation::Homogeneous, SpatialDiscretization::Collocated,
        EBReconstructionConvention::R1};
    /// The same path with \f$\phi\f$ half a step behind \f$\mathbf A\f$ (R2Prime).
    inline constexpr SolverConfiguration StandardR2PrimeConfiguration{
        "StandardR2Prime", solver::MaterialFormulation::Homogeneous,
        SpatialDiscretization::Collocated, EBReconstructionConvention::R2Prime};

    /// Read the required `parameters.solver_configuration` of an input that selects it.
    inline SolverConfiguration solverConfiguration(const input::ChdrYamlRunInput& request) {
        if (request.parameters_m.has_path("solver_configuration")
            && request.parameters_m["solver_configuration"].dtype().is_string()) {
            const std::string name = request.parameters_m["solver_configuration"].as_string();
            for (const auto& configuration :
                 std::array{O1VSConfiguration, HomogeneousConfiguration, StandardConfiguration,
                            StandardR2PrimeConfiguration})
                if (name == configuration.name_m)
                    return configuration;
        }
        throw std::invalid_argument(
            "parameters.solver_configuration must be O1VS, Homogeneous, "
            "Standard or StandardR2Prime");
    }
}  // namespace chdr::experiments
