/** @file test_comoving.cpp
 * @brief Compare reconstructed co-moving fields with independent analytic predictions.
 * @details Output rows tagged T7a describe one Fourier mode; T7b describes the
 * neutralized Gaussian train. P1 columns contain exact-symbol predictions and
 * P2 columns contain their leading-order expansions, defined by comovingSymbols().
 *
 * TODO: REVIEW.
 */
#include <algorithm>
#include <gtest/gtest.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <type_traits>

#include "ComovingPotentials.h"
#include "MonitorTestSupport.h"
#include "RankRecords.h"
#include "core/base/SpatialDimension.h"

static_assert(
    chdr::solver::SpatialDim == 3,
    "test_comoving.cpp requires three spatial dimensions for its Cartesian implementation");

namespace {
    using namespace chdr::solver;
    using namespace chdr::solver::test;

    struct ErrorNorms {
        double l2_m, maximum_m;
    };
    struct ComovingComparison {
        ErrorNorms measured_m, exact_m, leading_m;
    };

    /** @brief Accumulate owned-cell errors without allowing a NaN to hide in a maximum.
     * @details The reference is the analytic E component at the reconstruction's
     * comparison time and position. Long-double local positive sums are converted
     * to plain doubles for MPI; halos never contribute to any metric.
     */
    struct MeasuredNorms {
        long double errorSum_m = 0.0L, exactSum_m = 0.0L;
        double maxima_m[2]{};

        void add(double actual, double reference) {
            const double error = actual - reference;
            if (!std::isfinite(actual)) {
                errorSum_m  = std::numeric_limits<long double>::infinity();
                maxima_m[0] = std::numeric_limits<double>::infinity();
            } else
                errorSum_m += static_cast<long double>(error) * error;
            exactSum_m += static_cast<long double>(reference) * reference;
            maxima_m[0] = std::max(maxima_m[0], std::abs(error));
            maxima_m[1] = std::max(maxima_m[1], std::abs(reference));
        }
        ErrorNorms global(MPI_Comm comm) {
            double sums[2]{double(errorSum_m), double(exactSum_m)};
            MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_DOUBLE, MPI_SUM, comm);
            MPI_Allreduce(MPI_IN_PLACE, maxima_m, 2, MPI_DOUBLE, MPI_MAX, comm);
            return {std::sqrt(sums[0] / sums[1]), maxima_m[0] / maxima_m[1]};
        }
    };

    /// Measure the actual solver kernel on owned cells; the reference arrays are independent.
    ErrorNorms measuredError(EMFieldsType& field, unsigned component,
                             const ComovingPotentials& wave, const std::vector<double>& exact) {
        const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, field.getView());
        const auto local = field.getLayout().getLocalNDIndex();
        const int ghost  = field.getNghost();
        MeasuredNorms norms;
        for (int i = local[0].first(); i <= local[0].last(); ++i)
            for (int j = local[1].first(); j <= local[1].last(); ++j)
                for (int k = local[2].first(); k <= local[2].last(); ++k)
                    norms.add(host(i - local[0].first() + ghost, j - local[1].first() + ghost,
                                   k - local[2].first() + ghost)[component],
                              exact[wave.index(i, j, k)]);
        return norms.global(field.getLayout().comm);
    }

    /** @brief Use Fourier-coefficient norms for relative L2 and sampled series for maxima.
     * @details Orthogonality of distinct integer box modes cancels the common
     * transform normalization in the ratio of squared coefficient norms. This is
     * the discrete Parseval identity, deduced from the DFT conventions in
     * [heinzel2002SpectrumSpectralDensity, Sec. 5, Eqs. (3)-(5)]. The maximum ratio
     * instead uses the independently evaluated error and reference series at the
     * actual sample positions. comovingSymbols() defines both prediction arrays;
     * no measured field enters this calculation.
     */
    ErrorNorms predictedError(const std::vector<Complex>& error,
                              const std::vector<Complex>& reference,
                              const std::vector<double>& errorGrid,
                              const std::vector<double>& referenceGrid) {
        long double errorNorm = 0.0L, referenceNorm = 0.0L;
        double errorMaximum = 0.0, referenceMaximum = 0.0;
        for (std::size_t i = 0; i < error.size(); ++i) {
            errorNorm += std::norm(error[i]);
            referenceNorm += std::norm(reference[i]);
            errorMaximum     = std::max(errorMaximum, std::abs(errorGrid[i]));
            referenceMaximum = std::max(referenceMaximum, std::abs(referenceGrid[i]));
        }
        return {std::sqrt(double(errorNorm / referenceNorm)), errorMaximum / referenceMaximum};
    }

    /// Preserve global-component traversal within each owned block for the bitwise R2 check.
    std::vector<double> ownedElectric(EMFieldsType& field) {
        const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, field.getView());
        const auto local = field.getLayout().getLocalNDIndex();
        const int ghost  = field.getNghost();
        std::vector<double> values;
        for (int i = 0; i < int(local[0].length()); ++i)
            for (int j = 0; j < int(local[1].length()); ++j)
                for (int k = 0; k < int(local[2].length()); ++k)
                    for (unsigned c = 0; c < SpatialDim; ++c)
                        values.push_back(host(i + ghost, j + ghost, k + ghost)[c]);
        return values;
    }

    /** @brief Check equality under the same floating-point contraction setting.
     * @details The two reconstruction kernels must be compiled with compatible
     * contraction flags; the evidence records the actual compile commands.
     */
    void checkEquivalentFields(EBReconstructionConvention ebReconstructionConvention,
                               EMFieldsType& field, std::vector<double>& r2Values) {
        if (ebReconstructionConvention == EBReconstructionConvention::R2)
            r2Values = ownedElectric(field);
        if (ebReconstructionConvention == EBReconstructionConvention::R2Prime) {
            EXPECT_EQ(ownedElectric(field), r2Values) << "R2 and R2Prime differ";
        }
    }

    /** @brief Evaluate one component at its own spatial offset and chosen output time.
     * @details R1 has two recorded comparison times; R2/R2Prime use the common
     * half step. For the staggered reading the E component is evaluated on its
     * positive face. ComovingPotentials documents the analytic field derivation
     * from the boosted potentials, including the cancellation in the longitudinal
     * component [fitzpatrick2014ClassicalElectromagnetism, Sec. 12.17, Eq. (12.190)].
     */
    template <class SolverType>
    ComovingComparison compareComponent(SolverType& solver, ComovingPotentials& wave,
                                        unsigned component, bool r1Half, double label) {
        const auto componentOffset =
            solver.levelContract().positionOffset(Quantity::Electric, component);
        const std::array<double, SpatialDim> offset{componentOffset[0], componentOffset[1],
                                                    componentOffset[2]};
        const auto symbols =
            comovingSymbols(wave, solver.getDt(), solver.ebReconstructionConvention(),
                            solver.spatialDiscretization(), component, r1Half);
        const auto exact  = wave.evaluate(symbols.analytic_m, label, offset);
        const auto p1Grid = wave.evaluate(symbols.exactError_m, label, offset);
        const auto p2Grid = wave.evaluate(symbols.leadingError_m, label, offset);
        return {measuredError(solver.electricField(), component, wave, exact),
                predictedError(symbols.exactError_m, symbols.analytic_m, p1Grid, exact),
                predictedError(symbols.leadingError_m, symbols.analytic_m, p2Grid, exact)};
    }

    /** @brief Check the exact prediction for every case and the leading expansion for the train.
     * @details The single-mode cases resolve the exact trigonometric symbols. The
     * smooth Gaussian spectrum additionally tests their leading-order approximation;
     * its finite-mesh truncation error is assessed separately from round-off agreement.
     */
    void checkComparison(const ComovingComparison& result, bool singleMode) {
        const auto measured = result.measured_m, p1 = result.exact_m, p2 = result.leading_m;
        EXPECT_TRUE(std::isfinite(measured.l2_m) && std::isfinite(measured.maximum_m));
        EXPECT_LE(std::abs(measured.l2_m / p1.l2_m - 1.0), 1e-9);
        EXPECT_LE(std::abs(measured.maximum_m / p1.maximum_m - 1.0), 1e-9);
        if (!singleMode) {
            EXPECT_LE(std::abs(measured.l2_m / p2.l2_m - 1.0), 0.10);
            EXPECT_LE(std::abs(measured.maximum_m / p2.maximum_m - 1.0), 0.10);
        }
    }

    template <class SolverType>
    void recordComparison(SolverType& solver, const ComovingComparison& result, bool singleMode,
                          unsigned component, double label) {
        if (ippl::Comm->rank() != 0)
            return;
        const std::string name = std::string(singleMode ? "T7a-" : "T7b-")
                                 + std::to_string(int(SolverType::stencil())) + "-"
                                 + std::to_string(int(solver.ebReconstructionConvention())) + "-"
                                 + std::to_string(int(solver.spatialDiscretization())) + "-"
                                 + std::to_string(component) + (label == 0 ? "-whole" : "-half");
        std::filesystem::create_directories(recordDirectory());
        std::ofstream out(recordDirectory() / (name + ".jsonl"));
        out << nlohmann::json{{"label", label},
                              {"dt", solver.getDt()},
                              {"measured_l2", result.measured_m.l2_m},
                              {"measured_max", result.measured_m.maximum_m},
                              {"exact_symbol_l2", result.exact_m.l2_m},
                              {"exact_symbol_max", result.exact_m.maximum_m},
                              {"leading_series_l2", result.leading_m.l2_m},
                              {"leading_series_max", result.leading_m.maximum_m}}
                   .dump()
            << '\n';
    }

    template <class SolverType>
    void printComparison(SolverType& solver, const ComovingComparison& result, bool singleMode,
                         unsigned component, double label) {
        if (ippl::Comm->rank() != 0)
            return;
        const auto measured = result.measured_m, p1 = result.exact_m, p2 = result.leading_m;
        std::cout << std::setprecision(17) << (singleMode ? "T7a" : "T7b") << " stencil="
                  << (solver.stencil() == Stencil::NonStandard ? "NonStandard" : "Standard")
                  << " reconstruction=" << int(solver.ebReconstructionConvention())
                  << " reading=" << int(solver.spatialDiscretization())
                  << " component=" << component << " label/dt=" << label / solver.getDt()
                  << " L2=" << measured.l2_m << " P1_L2=" << p1.l2_m << " P2_L2=" << p2.l2_m
                  << " max=" << measured.maximum_m << " P1_max=" << p1.maximum_m
                  << " P2_max=" << p2.maximum_m << '\n';
    }

    template <class SolverType>
    void checkOutput(SolverType& solver, ComovingPotentials& wave, bool singleMode, bool r1Half) {
        const double label = solver.ebReconstructionConvention() == EBReconstructionConvention::R1
                                 ? (r1Half ? -0.5 * solver.getDt() : 0.0)
                                 : solver.labelTime(Quantity::Electric, solver.completedSolves());
        for (unsigned component : {2u, 0u}) {
            SCOPED_TRACE(::testing::Message()
                         << (singleMode ? "T7a" : "T7b") << " stencil=" << int(solver.stencil())
                         << " reconstruction=" << int(solver.ebReconstructionConvention())
                         << " reading=" << int(solver.spatialDiscretization())
                         << " component=" << component << " label/dt=" << label / solver.getDt());
            const auto comparison = compareComponent(solver, wave, component, r1Half, label);
            checkComparison(comparison, singleMode);
            printComparison(solver, comparison, singleMode, component, label);
            recordComparison(solver, comparison, singleMode, component, label);
        }
    }

    /// Keep the field owner stable: helpers receive it by reference and the solver stays local.
    template <class PotentialSolverBase>
    void runReconstruction(ChDRFieldContainer& fields, ComovingPotentials& wave, bool singleMode,
                           SpatialDiscretization spatialDiscretization,
                           EBReconstructionConvention ebReconstructionConvention,
                           std::vector<double>& r2Values) {
        ChdrSolver<PotentialSolverBase> solver(fields.getJ(), fields.getE(), fields.getB(),
                                               ebReconstructionConvention, spatialDiscretization);
        constexpr bool NONSTANDARD = std::is_same_v<PotentialSolverBase, NonStandardSolverPeriodic>;
        EXPECT_DOUBLE_EQ(solver.getDt(), NONSTANDARD ? 1.0 / 64.0 : 0.5 / 64.0);
        wave.load(solver);
        solver.reconstructEB();
        requireReconstructedDivB(solver);
        checkEquivalentFields(ebReconstructionConvention, fields.getE(), r2Values);
        for (bool r1Half : {false, true})
            if (!r1Half || ebReconstructionConvention == EBReconstructionConvention::R1)
                checkOutput(solver, wave, singleMode, r1Half);
    }

    template <class PotentialSolverBase>
    void comovingCase(bool singleMode, SpatialDiscretization spatialDiscretization =
                                           SpatialDiscretization::Collocated) {
        constexpr bool NONSTANDARD = std::is_same_v<PotentialSolverBase, NonStandardSolverPeriodic>;
        const double hz = 1.0 / 64.0, hp = (NONSTANDARD ? 1.5 : 1.0) * hz;
        const std::array<int, SpatialDim> cells = singleMode
                                                      ? std::array<int, SpatialDim>{8, 8, 64}
                                                      : std::array<int, SpatialDim>{32, 32, 64};
        ComovingPotentials wave(cells, {hp, hp, hz}, singleMode);
        auto fields = makeBoxFieldContainer(cells, wave.lengths(), 2);
        initializeFieldStorage(fields);
        std::vector<double> r2Values;
        for (auto ebReconstructionConvention :
             {EBReconstructionConvention::R1, EBReconstructionConvention::R2,
              EBReconstructionConvention::R2Prime})
            if (spatialDiscretization == SpatialDiscretization::Collocated
                || ebReconstructionConvention == EBReconstructionConvention::R2)
                runReconstruction<PotentialSolverBase>(fields, wave, singleMode,
                                                       spatialDiscretization,
                                                       ebReconstructionConvention, r2Values);
    }

    TEST(ComovingT7, RejectsNonperiodicSingleMode) {
        constexpr double PI = ComovingPotentials::PI;
        EXPECT_THROW(
            ComovingPotentials::requirePeriodicWave({16 * PI, 0, 6 * PI}, {0.1875, 0.1875, 1.0}),
            std::invalid_argument);
        EXPECT_NO_THROW(ComovingPotentials::requirePeriodicWave({32 * PI / 3, 0, 6 * PI},
                                                                {0.1875, 0.1875, 1.0}));
    }
    TEST(ComovingT7, StandardSingleMode) {
        comovingCase<StandardSolverPeriodic>(true);
    }
    TEST(ComovingT7, NonStandardSingleMode) {
        comovingCase<NonStandardSolverPeriodic>(true);
    }
    TEST(ComovingT7, StandardBunchTrain) {
        comovingCase<StandardSolverPeriodic>(false);
    }
    TEST(ComovingT7, NonStandardBunchTrain) {
        comovingCase<NonStandardSolverPeriodic>(false);
    }
    TEST(ComovingT7, StaggeredSingleMode) {
        comovingCase<StandardSolverPeriodic>(true, SpatialDiscretization::Staggered);
    }
    TEST(ComovingT7, StaggeredBunchTrain) {
        comovingCase<StandardSolverPeriodic>(false, SpatialDiscretization::Staggered);
    }
}  // namespace
