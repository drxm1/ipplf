/** @file test_compact.cpp
 * @brief Natural-position symbols of the forward and backward staggeredDifference.
 * @details Signed half-cell symbols follow by Fourier substitution into
 * [ryu2016potentialfdtd, Appendix A, Eqs. (19)-(20)] (deduction).
 *
 * The test samples \f$\sin(\mathbf k\cdot\mathbf x)\f$ at cell centres, with
 * \f$\mathbf k=2\pi(1,2,-1)\f$, on 12 by 10 by 8 cells of the unit box (split along x). On
 * each axis d, the forward and the backward difference must equal
 * \f$\frac{2}{h_d}\sin(k_dh_d/2)\cos(\theta\pm k_dh_d/2)\f$, where \f$\theta\f$ is the phase
 * at the cell centre, to \f$10^{-12}\f$ relative to the symbol. A nonfinite value fails. The
 * reference is an exact lattice value. Each rank samples the mode at its own local indices
 * and reads no field or halo value; on two ranks each rank samples a shorter local x range,
 * and MPI_Allreduce combines the maximum errors. The suite name keeps "compact", the earlier
 * name of the staggered difference. CTest runs the suite CompactOperators as
 * chdr.solver.compact.serial.
 */
#include <Kokkos_MathematicalConstants.hpp>
#include <gtest/gtest.h>
#include <limits>

#include "core/solver/BoxSetup.h"
#include "core/solver/operators/StaggeredOperators.h"

namespace {
    using namespace chdr;
    using namespace chdr::solver;
    struct SampledMode {
        ippl::Vector<double, 3> waveNumber_m, spacing_m;
        KOKKOS_INLINE_FUNCTION double operator()(int i, int j, int k) const {
            return Kokkos::sin(waveNumber_m[0] * (i + 0.5) * spacing_m[0]
                               + waveNumber_m[1] * (j + 0.5) * spacing_m[1]
                               + waveNumber_m[2] * (k + 0.5) * spacing_m[2]);
        }
    };

    /// Compare both signed staggered derivative symbols at their respective half-cell positions.
    struct CompactSymbolError {
        SampledMode sample_m;
        ippl::Vector<double, 3> symbol_m;
        static constexpr double INFINITY_VALUE = std::numeric_limits<double>::infinity();
        KOKKOS_INLINE_FUNCTION void operator()(int i, int j, int k, double& maximum) const {
            const auto waveNumber = sample_m.waveNumber_m, spacing = sample_m.spacing_m;
            const double phase = waveNumber[0] * (i + 0.5) * spacing[0]
                                 + waveNumber[1] * (j + 0.5) * spacing[1]
                                 + waveNumber[2] * (k + 0.5) * spacing[2];
            for (unsigned d = 0; d < 3; ++d) {
                const double forward = staggeredDifference(sample_m, i, j, k, d, spacing[d], true);
                const double backward =
                    staggeredDifference(sample_m, i, j, k, d, spacing[d], false);
                const double expectedForward =
                    symbol_m[d] * Kokkos::cos(phase + waveNumber[d] * spacing[d] / 2.0);
                const double expectedBackward =
                    symbol_m[d] * Kokkos::cos(phase - waveNumber[d] * spacing[d] / 2.0);
                const double forwardError =
                    Kokkos::abs(forward - expectedForward) / Kokkos::abs(symbol_m[d]);
                const double backwardError =
                    Kokkos::abs(backward - expectedBackward) / Kokkos::abs(symbol_m[d]);
                maximum = Kokkos::max(
                    maximum, Kokkos::isfinite(forwardError) ? forwardError : INFINITY_VALUE);
                maximum = Kokkos::max(
                    maximum, Kokkos::isfinite(backwardError) ? backwardError : INFINITY_VALUE);
            }
        }
    };

    /// The forward and backward staggeredDifference of one sampled mode equal the signed
    /// half-cell symbol on all three axes to \f$10^{-12}\f$.
    TEST(CompactOperators, NaturalPositionSingleModeSymbols) {
        auto fields = makeBoxFieldContainer({12, 10, 8}, {1.0, 1.0, 1.0}, 0);
        initializeFieldStorage(fields);
        const auto spacing = fields.getMesh().getMeshSpacing();
        const auto waveNumber =
            2.0 * Kokkos::numbers::pi_v<double> * ippl::Vector<double, 3>(1, 2, -1);
        const SampledMode sample{waveNumber, spacing};
        ippl::Vector<double, 3> symbol;
        for (unsigned d = 0; d < 3; ++d)
            symbol[d] = 2.0 * std::sin(waveNumber[d] * spacing[d] / 2.0) / spacing[d];
        double error = 0.0;
        Kokkos::parallel_reduce("Compact single-mode symbols", fields.getJ().getFieldRangePolicy(),
                                CompactSymbolError{sample, symbol}, Kokkos::Max<double>(error));
        MPI_Allreduce(MPI_IN_PLACE, &error, 1, MPI_DOUBLE, MPI_MAX, fields.getFL().comm);
        EXPECT_LE(error, 1e-12);
    }

}  // namespace
