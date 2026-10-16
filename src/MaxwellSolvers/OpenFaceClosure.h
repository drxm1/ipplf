/** @file OpenFaceClosure.h
 * @brief Complete potential ghosts beyond the physical faces of an absorbing box.
 */
#pragma once

#include <Kokkos_Core.hpp>

#include "Utility/IpplException.h"

namespace ippl {
    /** @brief Check the storage needed by the three-point physical ghost closure.
     * @details All ranks inspect all local domains, so an undersized remote partition is also
     * rejected before a caller exchanges halos. Two owned samples and one exchanged sample
     * suffice on a boundary rank; the global axis must contain at least three samples.
     */
    template <typename Field>
    void validateOpenFaceClosure(const Field& field) {
        static_assert(Field::dim == 3, "The FDTD potential closure requires three dimensions.");
        if (field.getNghost() != 1) {
            throw IpplException("validateOpenFaceClosure", "Exactly one ghost layer is required.");
        }
        const auto& layout = field.getLayout();
        for (unsigned d = 0; d < Field::dim; ++d) {
            if (layout.getDomain()[d].length() < 3) {
                throw IpplException("validateOpenFaceClosure",
                                    "At least three global cells are required.");
            }
            for (int rank = 0; rank < layout.comm.size(); ++rank) {
                if (layout.getLocalNDIndex(rank)[d].length() < 2) {
                    throw IpplException("validateOpenFaceClosure",
                                        "At least two local cells are required.");
                }
            }
        }
    }

    namespace detail {
        /** @brief Extrapolate one physical ghost plane, including its transverse ghost cells.
         * @details The caller selects only a locally owned physical face. The inward third
         * sample can be a halo sample when the rank owns only two cells along the normal.
         */
        template <typename Field>
        void extrapolateOpenFace(Field& field, unsigned axis, bool upper) {
            const auto view  = field.getView();
            const unsigned a = (axis + 1) % 3, b = (axis + 2) % 3;
            const int ghost = upper ? view.extent(axis) - 1 : 0, inward = upper ? -1 : 1;
            using Policy = Kokkos::MDRangePolicy<typename Field::execution_space, Kokkos::Rank<2>>;
            Kokkos::parallel_for(
                "Extrapolate open face", Policy({0, 0}, {view.extent(a), view.extent(b)}),
                KOKKOS_LAMBDA(int i, int j) {
                    int p[3];
                    p[a]          = i;
                    p[b]          = j;
                    p[axis]       = ghost + inward;
                    const auto u0 = view(p[0], p[1], p[2]);
                    p[axis] += inward;
                    const auto u1 = view(p[0], p[1], p[2]);
                    p[axis] += inward;
                    const auto u2          = view(p[0], p[1], p[2]);
                    p[axis]                = ghost;
                    view(p[0], p[1], p[2]) = 3 * u0 - 3 * u1 + u2;
                });
            Kokkos::fence();
        }
    }  // namespace detail

    /** @brief Extend every component through the physical ghost layer of an absorbing box.
     * @details Call after halo exchange. Only ranks owning a physical face write it. Axes are
     * completed in x, y, z order over the full transverse extent, including ghosts. Thus later
     * axes extend earlier extensions at edges and corners (a tensor-product deduction).
     *
     * The quadratic extension \f$g = 3u_0 - 3u_1 + u_2\f$ converts the unchanged centred first
     * difference into \f$(-3u_0+4u_1-u_2)/(2h)\f$ at a lower face. This is the second-order
     * one-sided formula in Fornberg (1988), "Generation of Finite Difference Formulas on
     * Arbitrarily Spaced Grids", Table 3, p. 704, doi:10.1090/S0025-5718-1988-0935077-0.
     * No owned value is changed. The closure does not impose an additional wave condition.
     * @pre validateOpenFaceClosure(field) succeeds. Owned values and exchanged interior
     * samples are current. An exchange before Mur can leave its seam copies stale;
     * this closure does not refresh them. Complete corner ghosts require the halo
     * exchange to follow Mur's updates as well.
     */
    template <typename Field>
    void fillOpenFaceGhosts(Field& field) {
        const auto& layout = field.getLayout();
        const auto& local  = layout.getLocalNDIndex();
        const auto& global = layout.getDomain();
        for (unsigned d = 0; d < Field::dim; ++d) {
            if (local[d].first() == global[d].first()) {
                detail::extrapolateOpenFace(field, d, false);
            }
            if (local[d].last() == global[d].last()) {
                detail::extrapolateOpenFace(field, d, true);
            }
        }
    }
}  // namespace ippl
