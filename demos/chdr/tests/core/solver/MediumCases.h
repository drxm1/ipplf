/** @file MediumCases.h
 * @brief One shared table of accepted and refused homogeneous material selections.
 */
#pragma once
#include <limits>

#include "core/formulation/MaterialFormulation.h"

namespace chdr::solver::test {
    /** @brief Expected acceptance for one material map and solver selection.
     * @details epsilon_m is the background relative permittivity. defect_m is zero
     * for a constant map, one for one finite mismatch, or two for one NaN sample.
     * The rejected rows isolate these map defects, a subunit value, a nonunit gauge
     * parameter and the NonStandard stencil. The final rows are valid controls.
     */
    struct MediumCase {
        double epsilon_m, alpha_m;
        bool nonstandard_m;
        int defect_m;
        bool accepted_m;
    };
    inline constexpr MediumCase MediumCases[] = {
        {2.13, 1, false, 1, false}, {2.13, 1, false, 2, false}, {0.5, 1, false, 0, false},
        {2.13, 2, false, 0, false}, {2.13, 1, true, 0, false},  {1, 1, false, 0, true},
        {2.13, 1, false, 0, true},  {4, 1, false, 0, true}};

    /** @brief Place any defect at the global origin, independent of rank decomposition.
     * @details Exactly one rank owns the defective sample. All ranks must nevertheless
     * receive the same refusal from the owner's cached global statistics; testing a
     * local-only validity check would miss the defect on the other ranks.
     */
    struct GuardMap {
        MediumCase row_m;
        KOKKOS_INLINE_FUNCTION double operator()(material::CellIndex global, material::Position,
                                                 bool& marked) const {
            marked = false;
            if (global[0] == 0 && global[1] == 0 && global[2] == 0) {
                if (row_m.defect_m == 1)
                    return row_m.epsilon_m + 1;
                if (row_m.defect_m == 2)
                    return std::numeric_limits<double>::quiet_NaN();
            }
            return row_m.epsilon_m;
        }
    };
}  // namespace chdr::solver::test
