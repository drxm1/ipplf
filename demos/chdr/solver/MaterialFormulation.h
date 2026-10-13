/** @file MaterialFormulation.h
 * @brief Material selection, validity and the single host factory registration point.
 *
 * Responsibility split:
 * - chdr::material stores and samples permittivity, manages ghosts and provides derived coefficient
 * fields.
 * - chdr::solver (here) selects the equations and the update algorithm that use the material data.
 */
#pragma once

#include <memory>

#include "LevelContract.h"
#include "MaterialOwner.h"
#include "MaterialStepContext.h"
#include "Operators.h"

namespace chdr::solver {
    /// Available material equations; each formulation owns its validity rule.
    enum class MaterialFormulation {
        Homogeneous
    };

    /// Stable host label written alongside formulation parameters in run records.
    inline const char* formulationName(MaterialFormulation formulation) {
        switch (formulation) {
            case MaterialFormulation::Homogeneous:
                return "Homogeneous";
        }
        return "Unknown";
    }

    /** @brief Select the material equations and their generalized-gauge coefficient.
     * @details alpha_m is the dimensionless gauge parameter of
     * [chew2014generalizedgauge, Eq. (14) and text after Eq. (18)]. The implemented
     * homogeneous Lorenz update accepts unity only; it is not a time-step factor.
     */
    struct FormulationParameters {
        MaterialFormulation formulation_m = MaterialFormulation::Homogeneous;
        double alpha_m                    = 1.0;
    };

    /** @brief Record the solver stencil, component placement and E/B reconstruction choices.
     * @details The stencil advances potentials, spatialDiscretization_m places
     * components on the grid, and ebReconstructionConvention_m selects the field
     * reconstruction's time labels. These are separate choices, checked together
     * before base construction or formulation workspace allocation.
     */
    struct SolverDiscretizationSelection {
        Stencil stencil_m;
        SpatialDiscretization spatialDiscretization_m;
        EBReconstructionConvention ebReconstructionConvention_m;
    };

    /** @brief Refuse unsupported material equations before any solver collective.
     * @details Reads only cached, rank-identical statistics. Every participating
     * rank must supply the same parameters and selection.
     * The statistics summarize owned relative-permittivity samples over all ranks.
     * No field scan or MPI call occurs here, so a rejected common selection throws
     * on all ranks before a caller can enter base-constructor communication.
     */
    inline void validateMaterialUpdate(const FormulationParameters& parameters,
                                       const material::MapStatistics& statistics,
                                       const SolverDiscretizationSelection& selection);

    /** @brief Validate and allocate formulation workspace on the prototype's mesh/layout.
     * @details The returned object owns its scratch storage. The prototype supplies
     * the allocation geometry; it is not a destination for this factory. Its mesh
     * and layout objects must outlive the returned update. Validation is repeated
     * here so direct factory callers cannot bypass the formulation's restrictions.
     * @pre The owner and prototype use the same mesh and layout objects.
     */
    inline std::unique_ptr<MaterialUpdate> makeMaterialUpdate(
        const FormulationParameters& parameters, const material::MaterialOwner& owner,
        const SolverDiscretizationSelection& selection, AJFourFieldsType& prototype);
}  // namespace chdr::solver
