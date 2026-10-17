/** @file RandomFields.h
 * @brief Decomposition-independent history and source samples for interface checks.
 */
#pragma once
#include "CoreTestSupport.h"

namespace chdr::interfaceTest {
    using namespace solver;
    /** @brief Hash global indices, then remove each component's global mean.
     * @details Negative seeds select signed-zero data without mean subtraction.
     * Independent seeds excite all three axes, four components and both levels.
     * Removing the constant mode avoids the double-root secular solution of
     * leapfrog [kreiss2002SecondOrderWave, Eq. (2.8)] in the short medium comparison.
     */
    inline void randomize(AJFourFieldsType& field, int seed, bool removeMean = true) {
        const auto v     = field.getView();
        const auto first = field.getLayout().getLocalNDIndex().first();
        Kokkos::parallel_for(
            field.getFieldRangePolicy(), KOKKOS_LAMBDA(int i, int j, int k) {
                for (unsigned c = 0; c < 4; ++c) {
                    const auto bits = material::detail::hash(
                        (i + first[0]) * 104729ULL + (j + first[1]) * 257ULL
                        + (k + first[2]) * 7919ULL + c * 65537ULL + seed * 99991ULL);
                    v(i, j, k)[c] = seed < 0 ? -0.0 : double(bits % 10007) / 10007 - 0.5;
                }
            });
        if (seed < 0 || !removeMean)
            return;
        const auto domain  = field.getLayout().getDomain();
        const double count = domain[0].length() * domain[1].length() * domain[2].length();
        for (unsigned c = 0; c < 4; ++c) {
            double mean = 0;
            Kokkos::parallel_reduce(
                field.getFieldRangePolicy(),
                KOKKOS_LAMBDA(int i, int j, int k, double& sum) { sum += v(i, j, k)[c]; }, mean);
            MPI_Allreduce(MPI_IN_PLACE, &mean, 1, MPI_DOUBLE, MPI_SUM, field.getLayout().comm);
            mean /= count;
            Kokkos::parallel_for(
                field.getFieldRangePolicy(),
                KOKKOS_LAMBDA(int i, int j, int k) { v(i, j, k)[c] -= mean; });
        }
    }
}  // namespace chdr::interfaceTest
