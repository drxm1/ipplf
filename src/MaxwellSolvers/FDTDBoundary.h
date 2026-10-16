/**
 * @file FDTDBoundary.h
 * @brief Describe periodic and Mur faces and partition writes between the FDTD updates.
 *
 * The periodic and absorbing solver templates select one boundary kind for the whole box.
 * In the absorbing case, the interior stencil skips the outermost owned layer, and
 * second_order_mur_boundary_conditions updates those owned cells through separate face,
 * edge and corner passes. Those boundary cells are not exterior ghosts.
 *
 * A mixed box needs a per-face description shared by the interior and boundary updates.
 * For example, with periodic y and Mur x faces, a cell on an x face at the y wrap belongs
 * to one Mur face, not a two-Mur-face edge. The description identifies the active Mur
 * faces and their speeds; InteriorBounds excludes their owned layers from the interior
 * stencil while retaining cells at periodic wraps. FieldLayout remains the source of
 * periodic exchange topology, and FDTDBoundary::fromLayout derives paired face kinds.
 *
 * This file provides configuration and write-region bounds, not the Mur update formulas,
 * halo exchange or exterior-ghost completion. Write ownership here means which update
 * computes a cell, not which MPI rank stores it.
 */
#ifndef IPPL_FDTD_BOUNDARY_H
#define IPPL_FDTD_BOUNDARY_H

#include <array>
#include <cmath>

#include "Types/Vector.h"

#include "Utility/IpplException.h"

#include "FieldLayout/FieldLayout.h"

namespace ippl {

    /** @brief Boundary operator at one face; periodic faces occur in opposing pairs.
     *
     * Answers: Is the considered face periodic or does it use Mur BCs?
     */
    enum class FaceKind {
        periodic,
        mur
    };

    /**
     * @brief Host description of the six faces of a three-dimensional potential solver.
     *
     * Face 2*d is the lower face of axis d; face 2*d+1 is its upper face.
     * Speeds use the interior update's units, in which its wave speed is one.
     * Periodic faces ignore their speed. Kernels capture the mask and scalar values,
     * rather than a host object or a virtual boundary operation.
     */
    struct FDTDBoundary {
        std::array<FaceKind, 6> kind{};
        std::array<double, 6> speed{1, 1, 1, 1, 1, 1};

        /** @brief Derive paired faces from the layout, using one speed on every open face. */
        static FDTDBoundary fromLayout(const FieldLayout<3>& layout, double waveSpeed = 1) {
            FDTDBoundary result;
            for (unsigned face = 0; face < 6; ++face) {
                result.kind[face] =
                    layout.isPeriodic(face / 2) ? FaceKind::periodic : FaceKind::mur;
                result.speed[face] = waveSpeed;
            }
            return result;
        }

        /** @brief Return the face bits used to distinguish open faces from periodic wraps.
         *
         * Which faces the boundary kernels should treat as absorbing.
         *
         * The 6 bits answer: Does this face use Mur BCs?
         * { (bit, face) } = { (0, lower x), (1, upper x), (2, lower y), (3, upper y), (4, lower z),
         * (5, upper z) }
         *
         * One can then ask whether a particular face uses Mur:
         * bool usesMur = (mask & (1u << face)) != 0
         */
        unsigned murMask() const {
            unsigned mask = 0;
            for (unsigned face = 0; face < 6; ++face)
                if (kind[face] == FaceKind::mur)
                    mask |= 1u << face;
            return mask;
        }
    };

    /**
     * @brief Strict bounds of the interior update in ghost-offset global indices.
     *
     * An owned sample has index 1 through n. Mur owns the first and last samples
     * of an open axis; the interior update owns all samples of a periodic axis.
     * This value type can be captured by a Kokkos device kernel.
     */
    struct InteriorBounds {
        Vector<unsigned, 3> lower_m, upper_m;

        InteriorBounds(const NDIndex<3>& domain, const FDTDBoundary& boundary) {
            for (unsigned d = 0; d < 3; ++d) {
                lower_m[d] = boundary.kind[2 * d] == FaceKind::mur ? 1 : 0;
                upper_m[d] =
                    domain[d].length() + (boundary.kind[2 * d + 1] == FaceKind::mur ? 0 : 1);
            }
        }

        /** @brief Test exclusive ownership by the interior stencil, excluding every Mur layer. */
        KOKKOS_INLINE_FUNCTION bool contains(unsigned i, unsigned j, unsigned k) const {
            return lower_m[0] < i && i < upper_m[0] && lower_m[1] < j && j < upper_m[1]
                   && lower_m[2] < k && k < upper_m[2];
        }
    };

}  // namespace ippl
#endif
