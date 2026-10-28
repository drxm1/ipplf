/** @file test_open_gauge_records.cpp
 * @brief Record gauge and Gauss residuals after componentwise Mur reflection.
 */
#include <gtest/gtest.h>

#include "OpenRunSupport.h"
#include "TwoWaveState.h"
#include "core/solver/BoxSetup.h"

namespace chdr::solver::test {
    namespace {
        using reflection::Case;
        using reflection::Complex;

        /** @brief Staggered TMT or pure-gauge incidence, with each component's own Mur plane.
         * @details TMT has \f$\phi=0,\ A_x=K_z,\ A_z=K_x\f$; pure gauge has
         * \f$\phi=iW,\ A_x=-iK_x,\ A_z=iK_z\f$ for unit gauge function amplitude.
         * These satisfy the incident vacuum gauge on the Standard lattice
         * [Chew (2014), section 2, Eq. (11); Ryu et al. (2016), Appendix A].
         * TwoWaveState supplies the unchanged oracle roots and componentwise reflection.
         */
        struct GaugeState {
            reflection::TwoWaveState real_m, imaginary_m;
            Kokkos::Array<double, 4> realAmplitude_m{}, imaginaryAmplitude_m{};
            double hx_m;
            int part_m;
            GaugeState(const Case& row, bool pureGauge, int part)
                : real_m(row, SpatialDiscretization::Staggered, false, 0)
                , imaginary_m(row, SpatialDiscretization::Staggered, false, 1)
                , hx_m(row.h_m[0])
                , part_m(part) {
                const double theta = row.oracle().reflectedRoot.real();
                const double kx    = 2 * std::sin(theta / 2) / hx_m;
                const double kz    = 2 * std::sin(row.tangent()[2] * row.h_m[2] / 2) / row.h_m[2];
                const double w     = 2 * std::sin(row.omega_m * row.dt_m / 2) / row.dt_m;
                const std::array<Complex, 4> amplitude =
                    pureGauge
                        ? std::array<Complex, 4>{Complex(0, w), Complex(0, -kx), 0., Complex(0, kz)}
                        : std::array<Complex, 4>{0., kz, 0., kx};
                for (unsigned q = 0; q < 4; ++q) {
                    const auto a       = amplitude[q] * std::polar(1., -theta * (q == 1 ? 1. : .5));
                    realAmplitude_m[q] = a.real();
                    imaginaryAmplitude_m[q] = a.imag();
                }
            }
            KOKKOS_FUNCTION double operator()(double x, double y, double z, double t,
                                              unsigned q) const {
                const double shifted = x - (q == 1 ? .5 * hx_m : 0.);
                const double r = real_m(shifted, y, z, t, 2), i = imaginary_m(shifted, y, z, t, 2);
                return part_m == 0 ? realAmplitude_m[q] * r - imaginaryAmplitude_m[q] * i
                                   : realAmplitude_m[q] * i + imaginaryAmplitude_m[q] * r;
            }
        };

        /** @brief Gauge at integer time and Gauss at half time from the completed array ghosts.
         * @details Backward face-to-centre differences give \f$G=\delta_t\phi+D_-\cdot A\f$
         * and \f$D_-\cdot E\f$ [Ryu et al. (2016), Appendix A; Chew (2014), Eqs. (9),(11)].
         * Evaluate the row documented by staggeredDetail::electricAt at owned and backward
         * neighbor indices: \f$E_d=-(A_d[n]-A_d[n-1])/\Delta t-D_{+,d}\phi[n-1]\f$.
         * It has time label \f$t_{n-1/2}\f$ and needs only the same closed potential ghosts.
         * No electric-field ghost extrapolation or additional boundary operator is introduced.
         */
        std::array<double, 2> residuals(OpenStandardSolver& solver, int layer) {
            const auto a = coreTest::copyToHost(solver.A_n), b = coreTest::copyToHost(solver.A_nm1);
            const auto h        = solver.A_n.get_mesh().getMeshSpacing();
            const auto electric = [&](std::array<int, 3> p, unsigned d) {
                auto next = p;
                ++next[d];
                return -(a(p[0], p[1], p[2])[d + 1] - b(p[0], p[1], p[2])[d + 1]) / solver.getDt()
                       - (b(next[0], next[1], next[2])[0] - b(p[0], p[1], p[2])[0]) / h[d];
            };
            const int i = layer + 1;
            std::array<double, 2> result{(a(i, 1, 1)[0] - b(i, 1, 1)[0]) / solver.getDt(), 0.};
            for (unsigned d = 0; d < 3; ++d) {
                int p[3]{i, 1, 1};
                --p[d];
                result[0] += (a(i, 1, 1)[d + 1] - a(p[0], p[1], p[2])[d + 1]) / h[d];
                result[1] += (electric({i, 1, 1}, d) - electric({p[0], p[1], p[2]}, d)) / h[d];
            }
            return result;
        }

        /// The physical face includes closure truncation; layer three retains the plane-wave
        /// oracle.
        void gaugeRecord(const Case& row, bool pureGauge) {
            auto fields = makeBoxFieldContainer(row.cells(), row.lengths(), 2, {false, true, true});
            initializeFieldStorage(fields);
            std::array<std::array<Complex, 2>, 2> measured{};
            for (int part = 0; part < 2; ++part) {
                auto solver = makeOpenTestSolver(fields.getJ(), fields.getE(), fields.getB(), 0,
                                                 EBReconstructionConvention::R2,
                                                 SpatialDiscretization::Staggered);
                solver->setTimeStepUnchecked(row.dt_m);
                loadHistory<StandardSolverPeriodic>(*solver, GaugeState(row, pureGauge, part));
                for (int step = 0; step < 30; ++step)
                    solver->solve();
                for (unsigned p = 0; p < 2; ++p) {
                    const auto values = residuals(*solver, p == 0 ? 0 : 3);
                    for (unsigned q = 0; q < 2; ++q)
                        measured[p][q] +=
                            Complex(part == 0 ? values[q] : 0., part == 1 ? values[q] : 0.);
                }
            }
            const auto r       = row.oracle();
            const double theta = r.reflectedRoot.real(), kz = row.tangent()[2];
            const double w   = 2 * std::sin(row.omega_m * row.dt_m / 2) / row.dt_m;
            const double sx  = 2 * std::sin(theta / 2) / row.h_m[0],
                         sz  = 2 * std::sin(kz * row.h_m[2] / 2) / row.h_m[2];
            const Complex ax = pureGauge ? Complex(0, -sx) : Complex(sz, 0);
            const auto gauge = reflectedGaugeResidual(r, row.h_m[0], 0., ax, true);
            // Positive incident scale from the absolute temporal and divergence terms.
            const double scale = pureGauge ? w * w + sx * sx + sz * sz : 2 * sx * sz;
            for (unsigned p = 0; p < 2; ++p) {
                const int j      = p == 0 ? 0 : 3;
                const auto phase = std::polar(
                    1., theta * (j + .5) + kz * row.h_m[2] / 2 - row.omega_m * 30 * row.dt_m);
                const auto expected = gauge * phase;
                // Interior Gauss = -delta_t G, hence iW G at half time (deduction from
                // the staggered symbols; Christlieb et al. (2025), Theorem 3.5).
                const auto gauss =
                    Complex(0, w) * expected * std::polar(1., row.omega_m * row.dt_m / 2);
                std::cout << "gauge-record " << row.name_m
                          << " representation=" << (pureGauge ? "pure-gauge" : "TMT")
                          << " layer=" << j << " closed-gauge=" << measured[p][0]
                          << " plane-gauge=" << expected << " closed-Gauss=" << measured[p][1]
                          << " plane-Gauss=" << gauss << " gauge-scale=" << scale
                          << " Gauss-scale=" << w * scale << " gauge-difference/scale="
                          << std::abs(measured[p][0] - expected) / scale
                          << " Gauss-difference/scale="
                          << std::abs(measured[p][1] - gauss) / (w * scale) << '\n';
            }
        }
    }  // namespace

    TEST(OpenGaugeRecords, VacuumRepresentatives) {
        if (ippl::Comm->size() != 1)
            GTEST_SKIP() << "This record uses one rank";
        for (const auto& row : {reflection::cases[0], reflection::cases[1]})
            for (bool pureGauge : {false, true})
                gaugeRecord(row, pureGauge);
    }
}  // namespace chdr::solver::test
