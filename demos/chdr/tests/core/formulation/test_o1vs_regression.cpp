/** @file test_o1vs_regression.cpp
 * @brief Vacuum bit patterns, uniform-medium source scaling and admissibility.
 */
#include <array>
#include <gtest/gtest.h>

#include "opt-in/diagnostics/Monitors.h"
#include "tests/core/solver/RandomFields.h"

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
    void equivalent(double epsilon, bool signedZero,
                    std::array<bool, 3> periodic = {true, true, true}) {
        auto fields = makeBoxFieldContainer({8, 6, 10}, {8, 7.8, 8}, 0, periodic);
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
    TEST(O1VS, OpenVacuumBits) {
        equivalent(1, false, {false, true, false});
        equivalent(1, true, {false, true, false});
    }
    TEST(O1VS, UniformMedium) {
        equivalent(2.13, false);
    }

    /// Exercise each refusal independently using cached metadata before any solver construction.
    TEST(O1VS, Refusals) {
        auto fields = makeBoxFieldContainer({8, 6, 10}, {8, 7.8, 8});
        material::MapStatistics stats;
        stats.minimum_m = 1;
        stats.maximum_m = 2.13;
        const SolverDiscretizationSelection selected{Stencil::Standard, Staggered, R2};
        auto open = makeBoxFieldContainer({8, 6, 10}, {8, 7.8, 8}, 0, {false, true, false});
        EXPECT_THROW(
            validateMaterialUpdate({MaterialFormulation::O1VS, 2}, stats, selected, open.getFL()),
            IpplException);
        for (double alpha : {0.5, std::numeric_limits<double>::quiet_NaN()})
            EXPECT_THROW(validateMaterialUpdate({MaterialFormulation::O1VS, alpha}, stats, selected,
                                                fields.getFL()),
                         IpplException);
        for (auto invalid : {SolverDiscretizationSelection{Stencil::NonStandard, Staggered, R2},
                             SolverDiscretizationSelection{Stencil::Standard,
                                                           SpatialDiscretization::Collocated, R2},
                             SolverDiscretizationSelection{Stencil::Standard, Staggered,
                                                           EBReconstructionConvention::R1}})
            EXPECT_THROW(validateMaterialUpdate({MaterialFormulation::O1VS, 1}, stats, invalid,
                                                fields.getFL()),
                         IpplException);
        for (double bad : {0.5, 1e155}) {
            stats.minimum_m = stats.maximum_m = bad;
            EXPECT_THROW(validateMaterialUpdate({MaterialFormulation::O1VS, 1}, stats, selected,
                                                fields.getFL()),
                         IpplException);
        }
        stats.allFinite_m = false;
        EXPECT_THROW(
            validateMaterialUpdate({MaterialFormulation::O1VS, 1}, stats, selected, fields.getFL()),
            IpplException);
        initializeFieldStorage(fields);
        for (double epsilon : {1.0, 2.13}) {
            material::MaterialOwner owner(
                fields.getMesh(), fields.getFL(),
                KOKKOS_LAMBDA(material::CellIndex, material::Position, bool& marked) {
                    marked = false;
                    return epsilon;
                });
            Solver solver(fields.getJ(), fields.getE(), fields.getB(), R2, Staggered, owner,
                          {MaterialFormulation::O1VS, 2});
            Monitors monitor(solver.A_n);
            EXPECT_THROW(monitor.beforeSolve(solver), std::logic_error);
        }
    }
    TEST(O1VS, FiniteHarmonicFace) {
        auto fields = makeBoxFieldContainer({8, 6, 10}, {8, 6, 10});
        initializeFieldStorage(fields);
        material::MaterialOwner owner(
            fields.getMesh(), fields.getFL(),
            KOKKOS_LAMBDA(material::CellIndex cell, material::Position, bool& marked) {
                marked = cell[0] < 4;
                return marked ? 9e153 : 1e154;
            });
        const auto face = owner.derive(material::HarmonicFace<0>{});
        EXPECT_EQ(coreTest::globalMaximum(
                      fields.getJ(), coreTest::localMaximum(
                                         fields.getJ(),
                                         KOKKOS_LAMBDA(int i, int j, int k) {
                                             return Kokkos::isfinite(face(i, j, k)) ? 0.0 : 1.0;
                                         })),
                  0);
    }
}  // namespace
