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

#include "HomogeneousUpdate.h"
#include "LevelContract.h"
#include "MaterialOwner.h"
#include "MaterialStepContext.h"
#include "O1VSUpdate.h"
#include "Operators.h"

namespace chdr::solver {
    /// Available material equations; each formulation owns its validity rule.
    enum class MaterialFormulation {
        Homogeneous,
        O1VS
    };

    /// Stable host label written alongside formulation parameters in run records.
    inline const char* formulationName(MaterialFormulation formulation) {
        switch (formulation) {
            case MaterialFormulation::Homogeneous:
                return "Homogeneous";
            case MaterialFormulation::O1VS:
                return "O1VS";
        }
        return "Unknown";
    }

    /** @brief Select the material equations and their generalized-gauge coefficient.
     * @details alpha_m is the dimensionless gauge parameter of
     * [chew2014generalizedgauge, Eqs. (15) and (18), p. 3]. The implemented
     * homogeneous Lorenz update accepts unity only; O1-VS accepts finite alpha
     * at least one with representable chi. It is not a time-step factor.
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
     * HomogeneousUpdate implements the constant-coefficient, unit-gauge equations
     * documented in that class. Its supported range is finite relative permittivity
     * at least one and the Standard stencil; a heterogeneous map cannot use this
     * formulation. homogeneous() checks finiteness and bitwise-equal extrema.
     * Both spatial placements and all reconstruction conventions are accepted here;
     * they do not change the formulation's potential recurrence.
     */
    inline void validateMaterialUpdate(const FormulationParameters& parameters,
                                       const material::MapStatistics& statistics,
                                       const SolverDiscretizationSelection& selection) {
        if (parameters.formulation_m != MaterialFormulation::Homogeneous)
            throw IpplException("validateMaterialUpdate", "Unknown material formulation");
        if (selection.stencil_m != Stencil::Standard || parameters.alpha_m != 1.0
            || !statistics.homogeneous() || statistics.minimum_m < 1.0)
            throw IpplException(
                "validateMaterialUpdate",
                "Homogeneous requires Standard, finite constant permittivity >= 1 and alpha = 1");
    }

    /** @brief Validate layout-dependent equations without communicating.
     * @details O1-VS supports periodic staggered R2 with Standard differences,
     * finite epsilon and alpha at least one, and representable \f$\chi=\alpha\varepsilon^2\f$
     * [chew2014generalizedgauge, Eq. (18)]. Cached extrema make refusals
     * identical on all ranks. The three-argument entry remains homogeneous-only.
     */
    inline void validateMaterialUpdate(const FormulationParameters& parameters,
                                       const material::MapStatistics& statistics,
                                       const SolverDiscretizationSelection& selection,
                                       const material::Layout& layout) {
        if (parameters.formulation_m != MaterialFormulation::O1VS) {
            validateMaterialUpdate(parameters, statistics, selection);
            return;
        }
        if (selection.stencil_m != Stencil::Standard
            || selection.spatialDiscretization_m != SpatialDiscretization::Staggered
            || selection.ebReconstructionConvention_m != EBReconstructionConvention::R2
            || !statistics.allFinite_m || statistics.minimum_m < 1
            || !std::isfinite(parameters.alpha_m) || parameters.alpha_m < 1
            || !std::isfinite(parameters.alpha_m * statistics.maximum_m * statistics.maximum_m)
            || !layout.isPeriodic(0) || !layout.isPeriodic(1) || !layout.isPeriodic(2))
            throw IpplException("validateMaterialUpdate",
                                "O1VS requires periodic Standard staggered R2, finite epsilon >= "
                                "1, alpha >= 1 and chi");
    }

    /** @brief Validate and allocate formulation workspace on the prototype's mesh/layout.
     * @details The returned object owns its scratch storage. The prototype supplies
     * the allocation geometry; it is not a destination for this factory. Its mesh
     * and layout objects must outlive the returned update. Validation is repeated
     * here so direct factory callers cannot bypass the formulation's restrictions.
     * @pre The owner and prototype use the same mesh and layout objects.
     */
    inline std::unique_ptr<MaterialUpdate> makeMaterialUpdate(
        const FormulationParameters& parameters, const material::MaterialOwner& owner,
        const SolverDiscretizationSelection& selection, AJFourFieldsType& prototype) {
        validateMaterialUpdate(parameters, owner.statistics(), selection, owner.layout());
        if (parameters.formulation_m == MaterialFormulation::O1VS)
            return std::make_unique<O1VSUpdate>(owner, parameters.alpha_m, prototype);
        return std::make_unique<HomogeneousUpdate>(owner.homogeneousValue(), prototype);
    }
}  // namespace chdr::solver
