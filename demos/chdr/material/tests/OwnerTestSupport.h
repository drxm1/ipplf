/** @file OwnerTestSupport.h
 * @brief Non-cubic grids and independent host observations of complete material allocations.
 */
#pragma once
#include <gtest/gtest.h>

#include "FELFieldContainer.hpp"
#include "MaterialOwner.h"
#include "OwnerOracles.h"

namespace chdr::material::test {
    inline const CellIndex Cells(12, 14, 16);

    /// Select one decomposition axis while exercising nonzero origins and unequal spacings.
    template <class Periodicity>
    inline FELFieldContainer<double, 3> makeGrid(unsigned axis, Periodicity periodic) {
        Position h(0.125, 0.25, 0.5), lower(-0.5, 1.0, 2.0), upper;
        ippl::NDIndex<3> domain;
        std::array<bool, 3> decomposition{false, false, false};
        decomposition[axis] = true;
        for (unsigned d = 0; d < 3; ++d) {
            domain[d] = ippl::Index(Cells[d]);
            upper[d]  = lower[d] + h[d] * Cells[d];
        }
        return FELFieldContainer<double, 3>(h, lower, upper, decomposition, domain, lower,
                                            periodic);
    }

    /** @brief Deterministic index map for halo and coefficient observations.
     * @details mode 0 produces unique integer IDs, mode 1 a reproducible two-value
     * map, mode 2 values in [1,16). No production indexing or hashing helper is used.
     */
    struct IndexMap {
        int mode_m = 0;
        KOKKOS_INLINE_FUNCTION double operator()(CellIndex g) const {
            const unsigned id     = (g[0] * 14 + g[1]) * 16 + g[2];
            const unsigned random = (id * 1664525U + 1013904223U) ^ (id * 22695477U);
            if (mode_m == 0)
                return 1 + id;
            if (mode_m == 1)
                return (random >> 12) % 2 ? 2.13 : 4.0;
            return 1 + 15.0 * (random % 65536) / 65536;
        }
        KOKKOS_INLINE_FUNCTION double operator()(CellIndex g, Position, bool& marked) const {
            const double value = (*this)(g);
            marked             = value < 3;
            return value;
        }
    };

    /** @brief Check every allocated cell on the host against an independently supplied value.
     * @details Expected takes the unwrapped global index; the oracle chooses its
     * source cell. Return one rank-global count so all ranks observe test failures.
     */
    template <class Expected>
    std::uint64_t mismatches(ScalarFieldConstView view, Layout& layout, Expected expected) {
        const auto host      = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), view);
        const auto local     = layout.getLocalNDIndex();
        std::uint64_t errors = 0;
        for (int i = 0; i < int(view.extent(0)); ++i)
            for (int j = 0; j < int(view.extent(1)); ++j)
                for (int k = 0; k < int(view.extent(2)); ++k) {
                    const CellIndex global(i + local[0].first() - 1, j + local[1].first() - 1,
                                           k + local[2].first() - 1);
                    errors += Kokkos::bit_cast<std::uint64_t>(host(i, j, k))
                              != Kokkos::bit_cast<std::uint64_t>(expected(global));
                }
        MPI_Allreduce(MPI_IN_PLACE, &errors, 1, MPI_UINT64_T, MPI_SUM, layout.comm);
        return errors;
    }
}  // namespace chdr::material::test
