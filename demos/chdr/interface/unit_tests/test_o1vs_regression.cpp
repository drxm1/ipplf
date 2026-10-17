/** @file test_o1vs_regression.cpp
 * @brief Vacuum bit patterns, uniform-medium source scaling and admissibility.
 */
#include <gtest/gtest.h>

#include "Monitors.h"
#include "RandomFields.h"

namespace {
    using namespace chdr;
    using namespace chdr::interfaceTest;
    using Solver             = ChdrSolver<StandardSolverPeriodic>;
    constexpr auto R2        = EBReconstructionConvention::R2;
    constexpr auto Staggered = SpatialDiscretization::Staggered;

    /// Compare all owned component bits, or global max error divided by global reference max.
    template <class Field>
    double difference(const Field& left, const Field& right, bool bits) {
        const auto a = left.getView(), b = right.getView();
        double error = coreTest::localMaximum(
            left, KOKKOS_LAMBDA(int i, int j, int k) {
                double value = 0;
                for (unsigned c = 0; c < Field::value_type::dim; ++c)
                    value = Kokkos::max(
                        value, bits ? double(Kokkos::bit_cast<std::uint64_t>(a(i, j, k)[c])
                                             != Kokkos::bit_cast<std::uint64_t>(b(i, j, k)[c]))
                                    : Kokkos::abs(a(i, j, k)[c] - b(i, j, k)[c]));
                return value;
            });
        error = coreTest::globalMaximum(left, error);
        if (bits)
            return error;
        const double scale = coreTest::globalMaximum(
            right, coreTest::localMaximum(
                       right, KOKKOS_LAMBDA(int i, int j, int k) {
                           double value = 0;
                           for (unsigned c = 0; c < Field::value_type::dim; ++c)
                               value = Kokkos::max(value, Kokkos::abs(b(i, j, k)[c]));
                           return value;
                       }));
        return scale == 0 ? error : error / scale;
    }
    /// Fresh sources share their physical values; both methods receive the same physical dt.
    void equivalent(double epsilon, bool signedZero) {
        auto fields = makeBoxFieldContainer({8, 6, 10}, {8, 7.8, 8});
        initializeFieldStorage(fields);
        material::MaterialOwner owner(
            fields.getMesh(), fields.getFL(),
            KOKKOS_LAMBDA(material::CellIndex, material::Position, bool& marked) {
                marked = false;
                return epsilon;
            });
        EMFieldsType e(fields.getMesh(), fields.getFL()), b(fields.getMesh(), fields.getFL());
        Solver candidate(fields.getJ(), fields.getE(), fields.getB(), R2, Staggered, owner,
                         {MaterialFormulation::O1VS, 1});
        auto reference = coreTest::makeTestSolver<StandardSolverPeriodic>(
            fields.getJ(), e, b, R2, Staggered, epsilon == 1 ? 0 : epsilon);
        candidate.setSafetyFactor(0.9);
        reference->setTimeStepUnchecked(candidate.getDt());
        randomize(candidate.A_n, signedZero ? -1 : 23);
        randomize(candidate.A_nm1, signedZero ? -1 : 71);
        // Opposite zero signs make every first-step vacuum summand negative zero;
        // changing the effective source's zero sign then changes the result bits.
        if (signedZero)
            candidate.A_nm1 = 0.0;
        Kokkos::deep_copy(reference->A_n.getView(), candidate.A_n.getView());
        Kokkos::deep_copy(reference->A_nm1.getView(), candidate.A_nm1.getView());
        candidate.prepareHistories();
        reference->prepareHistories();
        double error = difference(candidate.A_nm1, reference->A_nm1, epsilon == 1);
        for (int step = 0; step < (epsilon == 1 ? 256 : 20); ++step) {
            randomize(fields.getJ(), signedZero ? -1 : 100 + step);
            candidate.solve();
            reference->solve();
            error = std::max(error, difference(candidate.A_n, reference->A_n, epsilon == 1));
            if (epsilon == 1)
                error = std::max({error, difference(fields.getE(), e, true),
                                  difference(fields.getB(), b, true)});
        }
        std::cout << "epsilon=" << epsilon << " signed_zero=" << signedZero << " error=" << error
                  << '\n';
        EXPECT_LE(error, epsilon == 1 ? 0 : 1e-13);
    }
    TEST(O1VS, VacuumBits) {
        equivalent(1, false);
        equivalent(1, true);
    }
}  // namespace
