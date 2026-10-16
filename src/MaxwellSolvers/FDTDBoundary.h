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

    /**
     * @brief Reject incompatible face descriptions before any MPI communication.
     *
     * Inputs must be identical on all ranks. Every rank checks the stored domains
     * of all ranks, so an undersized partition cannot leave peers in a halo exchange.
     * Exactly one ghost layer is required. At least two owned cells per local axis keep
     * the inward A_np1 samples read by Mur on the same rank. Three global cells on an
     * open axis support the quadratic ghost extension.
     * The speed bound follows Trefethen and Halpern, "Well-posedness of one-way wave
     * equations and absorbing boundary conditions", Math. Comp. 47 (1986), 421-435,
     * Theorem 2, p. 424; it is not a discrete stability guarantee.
     * Intersecting open faces require the matched speed because the compatibility
     * equations eliminate normal derivatives with the interior wave equation
     * [bamberger1987SecondOrderCorner, section II, equations (2.10)-(2.11)].
     *
     * Ghost widths count layers of cells on each side of the local owned domain, not
     * field components or the total number of ghost cells. This function receives only
     * the widths, not the field objects. The three width checks are identical and may
     * also be used to validate the three stored potential levels.
     *
     * @param boundary Kind and Mur speed of each of the six physical faces. Face kinds
     * must match the layout's periodic axes; speed entries on periodic faces are ignored.
     * @param layout Global domain, periodic exchange topology and stored rank-local domains
     * to check. Every rank inspects all partitions without communicating.
     * @param sourceGhosts Source field's ghost-layer width, obtained from getNghost();
     * must be exactly one.
     * @param electricGhosts Electric field E's ghost-layer width, obtained from getNghost();
     * must be exactly one, independently of the source field's width.
     * @param magneticGhosts Magnetic field B's ghost-layer width, obtained from getNghost();
     * must be exactly one, independently of the source and electric fields' widths.
     * @throws IpplException If any storage, domain, face-kind or Mur-speed check fails.
     */
    inline void validateBoundary(const FDTDBoundary& boundary, const FieldLayout<3>& layout,
                                 int sourceGhosts, int electricGhosts, int magneticGhosts) {
        auto require = [](bool valid, const char* message) {
            if (!valid)
                throw IpplException("validateBoundary", message);
        };
        require(sourceGhosts == 1 && electricGhosts == 1 && magneticGhosts == 1,
                "FDTD fields require exactly one ghost cell");
        unsigned openAxes = 0;
        for (unsigned d = 0; d < 3; ++d)
            openAxes += boundary.kind[2 * d] == FaceKind::mur;
        for (unsigned d = 0; d < 3; ++d) {
            require(layout.getDomain()[d].first() == 0, "FDTD global indices must start at zero");
            const auto expected = layout.isPeriodic(d) ? FaceKind::periodic : FaceKind::mur;
            require(boundary.kind[2 * d] == expected && boundary.kind[2 * d + 1] == expected,
                    "FDTD face kinds must agree with layout periodicity");
            if (expected == FaceKind::mur) {
                require(layout.getDomain()[d].length() >= 3, "An open axis needs three cells");
                for (unsigned face = 2 * d; face < 2 * d + 2; ++face) {
                    const double speed = boundary.speed[face];
                    require(std::isfinite(speed) && speed > 0 && speed < std::sqrt(2.0),
                            "Mur speed must be finite, positive and below sqrt(2)");
                    require(openAxes < 2 || speed == 1,
                            "Intersecting Mur faces need matched speed");
                }
            }
            const auto domains = layout.getHostLocalDomains();
            for (unsigned rank = 0; rank < domains.extent(0); ++rank)
                require(domains(rank)[d].length() >= 2, "Every local axis needs two owned cells");
        }
    }
}  // namespace ippl
#endif
