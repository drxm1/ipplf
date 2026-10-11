/** @file test_manufactured.cpp
 * @brief Discrete manufactured solutions check the source sign, component order and time labels.
 * @details The source is manufactured from sampled discrete differences, not
 * inferred from a continuum equation. The stencils follow
 * [fallahi2020mithra20fullwavesimulation, Eqs. (3.11)-(3.15), (3.23)-(3.27)],
 * but the forcing uses IPPL's positive source sign: the printed MITHRA update
 * has the opposite sign. The sign is fixed by the wave equations (3.6)-(3.7).
 * Every term is exercised [salari2000CodeVerificationMms, section 3.1].
 * T1 identifies these manufactured-source records; the negative control
 * shifts charge sampling by half a step.
 *
 * TODO: REVIEW.
 */
#include <Kokkos_MathematicalConstants.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <utility>

#include "ChdrSolver.h"
#include "MonitorTestSupport.h"
#include "Operators.h"
#include "PeriodicSetup.h"
#include "RankRecords.h"
#include "SpatialDimension.h"

static_assert(
    chdr::solver::SpatialDim == 3,
    "test_manufactured.cpp requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver::test {
    namespace mmsDetail {
        /// Four-component storage shared by sampled potentials and their manufactured source.
        using FourComponentView = decltype(std::declval<AJFourFieldsType&>().getView());

        /// Frozen three-mode potential; each component has distinct time dependence.
        struct Manufactured {
            double shift_m = 0;
            KOKKOS_FUNCTION double mode(int m, double x, double y, double z, double t,
                                        unsigned c) const {
                const double pi = Kokkos::numbers::pi_v<double>;
                if (m == 0)
                    return Kokkos::cos(2 * pi * (x + 2 * y - z - (1 + 0.31 * c) * t) + c);
                if (m == 1)
                    return 0.5
                           * Kokkos::sin(2 * pi * (3 * x - y + 2 * z + (0.5 + 0.17 * c) * t)
                                         + 0.3 * c);
                return 0.25 * Kokkos::cos(2 * pi * (y + z - 0.23 * (1 + c) * t) + 0.7 * c);
            }
            KOKKOS_FUNCTION double operator()(double x, double y, double z, double t,
                                              unsigned c) const {
                if (c == 0)
                    t += shift_m;
                double value = 0;
                for (int m = 0; m < 3; ++m)
                    value += mode(m, x, y, z, t, c);
                return value;
            }
        };

        /** @brief Apply the literal discrete MMS source to exact sampled slots.
         * @details \f$s=(\psi_{n+1}-2\psi_n+\psi_{n-1})/\Delta t^2-L_h\psi_n\f$.
         * This is a deduction from the cited wave update. Only the middle spatial
         * sample needs a halo; time neighbours are read at the same owned cell.
         * old_m, current_m and next_m contain consecutive exact potential samples:
         * scalar potential at component 0, vector potential at components 1 through
         * SpatialDim. source_m receives the corresponding four-component forcing.
         * @tparam Scheme Solver's spatial wave stencil, fixed before kernel launch.
         */
        template <Stencil Scheme>
        struct Forcing {
            FourComponentView old_m, current_m, next_m, source_m;
            ippl::Vector<double, SpatialDim> spacing_m;
            double dt_m;
            KOKKOS_INLINE_FUNCTION void operator()(int i, int j, int k) const {
                for (unsigned c = 0; c < 4; ++c)
                    source_m(i, j, k)[c] =
                        (next_m(i, j, k)[c] - 2 * current_m(i, j, k)[c] + old_m(i, j, k)[c])
                            / (dt_m * dt_m)
                        - laplacian<Scheme>(component(current_m, c), i, j, k, spacing_m);
            }
        };

        /// Temporary samples borrow the existing mesh/layout; this object is never moved.
        class Sources {
        public:
            explicit Sources(AJFourFieldsType& prototype) {
                for (auto& field : samples_m)
                    field.initialize(prototype.get_mesh(), prototype.getLayout());
            }
            template <class SolverType>
            void write(SolverType& solver, Manufactured sample, int step) {
                for (int slot = 0; slot < 3; ++slot) {
                    writeSource(solver, sample, step + slot - 1);
                    Kokkos::deep_copy(samples_m[slot].getView(), solver.sourceField().getView());
                }
                samples_m[1].fillHalo();
                const Forcing<SolverType::stencil()> forcing{samples_m[0].getView(),
                                                             samples_m[1].getView(),
                                                             samples_m[2].getView(),
                                                             solver.sourceField().getView(),
                                                             solver.A_n.get_mesh().getMeshSpacing(),
                                                             solver.getDt()};
                Kokkos::parallel_for("discrete MMS forcing", solver.A_n.getFieldRangePolicy(),
                                     forcing);
            }

        private:
            AJFourFieldsType samples_m[3];
        };

        /** @brief Compare one evolved slot with its exact component-labelled sample.
         * @details Normalization uses the maximum actual reference sample, not
         * an analytic amplitude bound. A nonfinite difference becomes infinity,
         * so a maximum reduction cannot accidentally conceal a NaN.
         */
        struct Error {
            FourComponentView current_m;
            Manufactured exact_m;
            ippl::Vector<int, SpatialDim> first_m;
            ippl::Vector<double, SpatialDim> spacing_m;
            double dt_m;
            int step_m, phiLabel_m;
            bool normalization_m;
            KOKKOS_INLINE_FUNCTION void operator()(int i, int j, int k, double& result) const {
                const double x = (i + first_m[0] - 0.5) * spacing_m[0],
                             y = (j + first_m[1] - 0.5) * spacing_m[1];
                const double z = (k + first_m[2] - 0.5) * spacing_m[2];
                for (unsigned c = 0; c < 4; ++c) {
                    const double expected =
                        exact_m(x, y, z, (step_m + 1 + (c == 0 ? 0.5 * phiLabel_m : 0)) * dt_m, c);
                    const double difference = current_m(i, j, k)[c] - expected;
                    const double error      = Kokkos::isfinite(difference)
                                                  ? Kokkos::abs(difference)
                                                  : std::numeric_limits<double>::infinity();
                    result = Kokkos::fmax(result, normalization_m ? Kokkos::abs(expected) : error);
                }
            }
        };

        template <class SolverType>
        double maximumError(SolverType& solver, const Manufactured& exact, int step) {
            Error sample{solver.A_n.getView(),
                         exact,
                         solver.A_n.getLayout().getLocalNDIndex().first(),
                         solver.A_n.get_mesh().getMeshSpacing(),
                         solver.getDt(),
                         step,
                         solver.levelContract().halfStepLabel(Quantity::Phi).count(),
                         false};
            double values[2] = {};
            Kokkos::parallel_reduce("T1 error", solver.A_n.getFieldRangePolicy(), sample,
                                    Kokkos::Max<double>(values[0]));
            sample.normalization_m = true;
            Kokkos::parallel_reduce("T1 normalization", solver.A_n.getFieldRangePolicy(), sample,
                                    Kokkos::Max<double>(values[1]));
            MPI_Allreduce(MPI_IN_PLACE, values, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            return values[0] / values[1];
        }

        template <class SolverType>
        void report(SolverType& solver, bool wrongLabel, double worst) {
            const std::string stencil =
                SolverType::stencil() == Stencil::NonStandard ? "NS" : "Std";
            const int ebReconstructionConvention =
                static_cast<int>(solver.ebReconstructionConvention());
            std::cout << "T1 stencil=" << stencil
                      << " reconstruction=" << ebReconstructionConvention
                      << " wrong_label=" << wrongLabel << " max_error=" << worst << '\n';
            ASSERT_TRUE(std::isfinite(worst));
            if (wrongLabel) {
                EXPECT_GE(worst, 1e-3);
            } else {
                EXPECT_LE(worst, 1e-12);
                recordState("T1-" + stencil + "-" + std::to_string(ebReconstructionConvention),
                            solver);
            }
        }
    }  // namespace mmsDetail

    /// Run freely for 256 solves; the reference never overwrites an evolved history.
    template <class PotentialSolverBase>
    void runManufactured(EBReconstructionConvention ebReconstructionConvention, bool wrongLabel) {
        constexpr bool NS = std::is_same_v<PotentialSolverBase, NonStandardSolverPeriodic>;
        auto fields       = makePeriodicFieldContainer({16, 16, NS ? 24 : 16}, {1, 1, 1}, 0);
        initializeFieldStorage(fields);
        ChdrSolver<PotentialSolverBase> solver(fields.getJ(), fields.getE(), fields.getB(),
                                               ebReconstructionConvention);
        const mmsDetail::Manufactured exact, source{wrongLabel ? solver.getDt() / 2 : 0};
        loadHistory(solver, exact);
        mmsDetail::Sources samples(solver.A_n);
        Monitors monitors(solver.A_n);
        double worst = 0;
        for (int n = 0; n < (wrongLabel ? 1 : 256); ++n) {
            samples.write(solver, source, n);
            advanceWithFinalMonitor(solver, monitors, n == (wrongLabel ? 0 : 255));
            worst = std::max(worst, mmsDetail::maximumError(solver, exact, n));
        }
        mmsDetail::report(solver, wrongLabel, worst);
    }

    template <class PotentialSolverBase>
    void manufacturedLabels() {
        for (auto ebReconstructionConvention :
             {EBReconstructionConvention::R1, EBReconstructionConvention::R2,
              EBReconstructionConvention::R2Prime}) {
            runManufactured<PotentialSolverBase>(ebReconstructionConvention, false);
            runManufactured<PotentialSolverBase>(ebReconstructionConvention, true);
        }
    }
    TEST(Manufactured, StandardLabelsAndNegativeControls) {
        manufacturedLabels<StandardSolverPeriodic>();
    }
    TEST(Manufactured, NonStandardLabelsAndNegativeControls) {
        manufacturedLabels<NonStandardSolverPeriodic>();
    }
}  // namespace chdr::solver::test
