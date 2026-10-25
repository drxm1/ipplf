/** @file test_equivalence.cpp
 * @brief Check bitwise equivalence of shifted histories and unit-permittivity evolution.
 * @details The R2/R2Prime comparison uses consistently shifted histories.
 * The extra scalar slot is produced by the same base update, not by an
 * analytic approximation. A one-slot scalar shift evolves independently of
 * the common vector histories (deduction from the componentwise wave update
 * of [fallahi2020mithra20fullwavesimulation, Eqs. (3.11)-(3.15)]).
 * Bitwise comparison assumes the two kernels use matching floating-point
 * contraction settings; the run record retains their actual compile commands.
 * The material comparison requires the unit-permittivity update to reproduce
 * the owner-free vacuum solver's potential histories and E/B fields exactly.
 *
 * TODO: REVIEW.
 */
#include <bit>
#include <cstdint>
#include <gtest/gtest.h>

#include "CoreTestSupport.h"
#include "MonitorTestSupport.h"
#include "Waves.h"
#include "core/base/SpatialDimension.h"
#include "testing/Util.h"

static_assert(
    chdr::solver::SpatialDim == 3,
    "test_equivalence.cpp requires three spatial dimensions for its Cartesian implementation");

namespace {
    using namespace chdr::solver;
    using namespace chdr::solver::coreTest;
    using chdr::testing::util::preparePotentialHalo;

    /** @brief Produce R2's extra scalar history while retaining the common vector slots.
     * @details prime.step() fills scratch A_np1 without shifting its histories or
     * increasing the shell counter. Only scalar components are copied into R2;
     * both solvers therefore keep identical A histories and consistently shifted
     * phi histories. Copied owned values receive true-copy halos before any solve.
     */
    template <class SolverType>
    void seedScalarShift(SolverType& r2, SolverType& prime, LorenzWave wave) {
        writeSource(prime, wave.sources(), 0);
        prime.step();
        const auto old = r2.A_nm1.getView(), current = r2.A_n.getView();
        const auto primeNow = prime.A_n.getView(), primeNext = prime.A_np1.getView();
        Kokkos::parallel_for(
            r2.A_n.getFieldRangePolicy(), KOKKOS_LAMBDA(int i, int j, int k) {
                old(i, j, k)[0]     = primeNow(i, j, k)[0];
                current(i, j, k)[0] = primeNext(i, j, k)[0];
            });
        preparePotentialHalo(r2.A_n);
        preparePotentialHalo(r2.A_nm1);
    }

    /// @brief Independent host data needed for one equal-label comparison; no mesh owners move.
    struct EquivalenceSnapshot {
        AJFourFieldsType::host_mirror_type current_m, previous_m, prime_m;
        EMFieldsType::host_mirror_type electric_m, primeElectric_m, magnetic_m, primeMagnetic_m;
    };

    /// @brief Retain current bytes in the original snapshot order after both solves finish.
    template <class SolverType>
    EquivalenceSnapshot snapshot(SolverType& r2, SolverType& prime) {
        return {copyToHost(r2.A_n),
                copyToHost(r2.A_nm1),
                copyToHost(prime.A_n),
                copyToHost(r2.electricField()),
                copyToHost(prime.electricField()),
                copyToHost(r2.magneticField()),
                copyToHost(prime.magneticField())};
    }

    /// @brief Compare representations, including the sign bit of zero, without a tolerance.
    bool differentBits(double first, double second) {
        return std::bit_cast<std::uint64_t>(first) != std::bit_cast<std::uint64_t>(second);
    }

    /** @brief Count every differing owned scalar/vector/field component at equal labels.
     * @details Previous R2 phi matches current R2Prime phi; current vector
     * potentials and reconstructed E/B match directly. Halos are not counted.
     */
    unsigned long long differingComponents(const EquivalenceSnapshot& s,
                                           const ippl::NDIndex<SpatialDim>& local) {
        unsigned long long count = 0;
        for (int i = 1; i <= static_cast<int>(local[0].length()); ++i)
            for (int j = 1; j <= static_cast<int>(local[1].length()); ++j)
                for (int k = 1; k <= static_cast<int>(local[2].length()); ++k) {
                    count += differentBits(s.previous_m(i, j, k)[0], s.prime_m(i, j, k)[0]);
                    for (unsigned c = 0; c < SpatialDim; ++c) {
                        count +=
                            differentBits(s.current_m(i, j, k)[c + 1], s.prime_m(i, j, k)[c + 1]);
                        count +=
                            differentBits(s.electric_m(i, j, k)[c], s.primeElectric_m(i, j, k)[c]);
                        count +=
                            differentBits(s.magnetic_m(i, j, k)[c], s.primeMagnetic_m(i, j, k)[c]);
                    }
                }
        return count;
    }

    /// @brief Keep the source timing distinct while comparing all sixty completed solves.
    template <class SolverType>
    void compareEvolution(SolverType& r2, SolverType& prime, LorenzWave wave) {
        unsigned long long differences = 0;
        const auto local               = r2.A_n.getLayout().getLocalNDIndex();
        Monitors r2Monitors(r2.A_n), primeMonitors(prime.A_n);
        for (int step = 0; step < 60; ++step) {
            writeSource(r2, wave.sources(), step);
            writeSource(prime, wave.sources(), step);
            chdr::solver::test::advanceWithFinalMonitor(r2, r2Monitors, step == 59);
            chdr::solver::test::advanceWithFinalMonitor(prime, primeMonitors, step == 59);
            differences += differingComponents(snapshot(r2, prime), local);
        }
        MPI_Allreduce(MPI_IN_PLACE, &differences, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                      r2.A_n.getLayout().comm);
        EXPECT_EQ(differences, 0ULL);
    }

    /// @brief Construct both owners in their final locations before allocating dependent fields.
    template <class PotentialSolverBase>
    void equivalent(std::array<int, SpatialDim> cells) {
        auto fields      = makeBoxFieldContainer(cells, {1, 1, 1});
        auto primeFields = makeBoxFieldContainer(cells, {1, 1, 1});
        initializeFieldStorage(fields);
        initializeFieldStorage(primeFields);
        ChdrSolver<PotentialSolverBase> r2(fields.getJ(), fields.getE(), fields.getB(),
                                           EBReconstructionConvention::R2);
        ChdrSolver<PotentialSolverBase> prime(primeFields.getJ(), primeFields.getE(),
                                              primeFields.getB(),
                                              EBReconstructionConvention::R2Prime);
        const auto wave = makeSourcedWave();
        loadHistory(r2, wave);
        loadHistory(prime, wave);
        seedScalarShift(r2, prime, wave);
        compareEvolution(r2, prime, wave);
    }

    TEST(LevelEquivalence, Standard) {
        equivalent<StandardSolverPeriodic>({8, 8, 16});
    }
    TEST(LevelEquivalence, NonStandard) {
        equivalent<NonStandardSolverPeriodic>({8, 8, 12});
    }

    /** @brief Count unequal component representations on this rank's owned cells.
     * @details Matching layouts are a fixture precondition. Host copies preserve the
     * values for inspection; differentBits distinguishes signed zero without a
     * numerical tolerance. The caller combines counts across fields, steps and ranks.
     */
    template <class Field>
    unsigned long long fieldBitDifferences(const Field& first, const Field& second) {
        const auto a = copyToHost(first), b = copyToHost(second);
        const auto local         = first.getLayout().getLocalNDIndex();
        unsigned long long count = 0;
        for (int i = 1; i <= int(local[0].length()); ++i)
            for (int j = 1; j <= int(local[1].length()); ++j)
                for (int k = 1; k <= int(local[2].length()); ++k)
                    for (unsigned c = 0; c < Field::value_type::dim; ++c)
                        count += differentBits(a(i, j, k)[c], b(i, j, k)[c]);
        return count;
    }

    /** @brief Require exact recovery of the owner-free solver at unit permittivity.
     * @details Identical labelled source-free histories evolve for 256 steps.
     * Both potential histories and reconstructed E/B are compared every step;
     * the full global difference count must be zero.
     * The two solvers share mesh/layout and the unchanged zero source, but own
     * separate potential histories and write separate E/B fields. Both retain their
     * constructor-default physical step; this checks the material route without
     * hiding a default-step difference behind an explicit safety-factor setter.
     * This source-free test checks propagation and reconstruction equivalence;
     * nonzero charge/current scaling requires the manufactured-source tests.
     */
    void vacuumMaterialEquivalent(EBReconstructionConvention reconstruction,
                                  SpatialDiscretization reading) {
        auto fields = makeBoxFieldContainer({16, 16, 16}, {1, 1, 1}, 0);
        initializeFieldStorage(fields);
        auto vacuum = makeTestSolver<StandardSolverPeriodic>(
            fields.getJ(), fields.getE(), fields.getB(), reconstruction, reading);
        EMFieldsType electric, magnetic;
        electric.initialize(fields.getMesh(), fields.getFL());
        magnetic.initialize(fields.getMesh(), fields.getFL());
        auto medium     = makeTestSolver<StandardSolverPeriodic>(fields.getJ(), electric, magnetic,
                                                                 reconstruction, reading, 1);
        const auto wave = makeMixedWave();
        loadHistory(*vacuum, wave);
        loadHistory(*medium, wave);
        unsigned long long differences = 0;
        for (int step = 0; step < 256; ++step) {
            vacuum->solve();
            medium->solve();
            differences += fieldBitDifferences(vacuum->A_n, medium->A_n);
            differences += fieldBitDifferences(vacuum->A_nm1, medium->A_nm1);
            differences += fieldBitDifferences(vacuum->electricField(), medium->electricField());
            differences += fieldBitDifferences(vacuum->magneticField(), medium->magneticField());
        }
        MPI_Allreduce(MPI_IN_PLACE, &differences, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                      MPI_COMM_WORLD);
        EXPECT_EQ(differences, 0ULL);
    }
    /// Exercise all collocated time reconstructions and the staggered R2 placement.
    TEST(MediumEquivalence, EveryReconstruction) {
        for (auto reconstruction : {EBReconstructionConvention::R1, EBReconstructionConvention::R2,
                                    EBReconstructionConvention::R2Prime})
            vacuumMaterialEquivalent(reconstruction, SpatialDiscretization::Collocated);
        vacuumMaterialEquivalent(EBReconstructionConvention::R2, SpatialDiscretization::Staggered);
    }
}  // namespace
