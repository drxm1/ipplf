/** @file E04Control.h
 * @brief Matched periodic replay for E04's startup and pre-contact comparison.
 */
#pragma once

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03E04Setup.h"

namespace chdr::experiments::e04 {
    /** @brief Retain the two native periodic field states needed by E04.
     * @details The replay shares the open run's mesh, rank decomposition, physical timestep
     * and analytic histories. It stops at the E03 final field state or the requested short
     * limit. Snapshots stay in the field memory space for direct differences on the same rank.
     */
    class E04PeriodicControl {
        using Solver = e03::SphericalSimulation::Solver;
        using View   = EMFieldsType::view_type;

        struct Frame {
            View e_m, b_m;
            std::array<int, 6> bounds_m{};
            long long step_m = -1;
        };

    public:
        E04PeriodicControl(simulation::ChdrCartesianGridSettings grid, const e03::E03Reference& ref,
                           int stopStep, int controlStep) {
            grid.periodic_m = {true, true, true};
            auto periodic   = e03::makeSphericalSimulation(grid, ref);
            periodic.load(ref);
            periodic.advance(stopStep, [&](long long step, Solver& state, const auto&) {
                if (step == 0 || step == controlStep)
                    frames_m[step == 0 ? 0 : 1] = capture(state, step);
            });
        }

        /// Return the absolute native L2 differences in E and R2 B, excluding ghosts.
        std::array<double, 2> difference(Solver& state, long long step) const {
            const auto& frame = frames_m[step == 0 ? 0 : 1];
            if (frame.step_m != step || frame.bounds_m != localBounds(state))
                throw std::runtime_error("E04 periodic control does not match the current state");
            auto squares = differenceSquares(state, frame);
            MPI_Allreduce(MPI_IN_PLACE, squares.data(), 2, MPI_DOUBLE, MPI_SUM,
                          state.A_n.getLayout().comm);
            const auto h = state.A_n.get_mesh().getMeshSpacing();
            for (double& value : squares) {
                if (!std::isfinite(value))
                    throw std::runtime_error("Nonfinite E04 periodic difference");
                value = std::sqrt(value * h[0] * h[1] * h[2]);
            }
            return squares;
        }

    private:
        static std::array<int, 6> localBounds(Solver& state) {
            const auto domain = state.A_n.getLayout().getLocalNDIndex();
            return {domain[0].first(), domain[0].last(),  domain[1].first(),
                    domain[1].last(),  domain[2].first(), domain[2].last()};
        }

        static View copyField(const char* name, View field) {
            View copy(Kokkos::view_alloc(Kokkos::WithoutInitializing, std::string(name)),
                      field.extent(0), field.extent(1), field.extent(2));
            Kokkos::deep_copy(copy, field);
            return copy;
        }

        static Frame capture(Solver& state, long long step) {
            return {copyField("E04 periodic E", state.electricField().getView()),
                    copyField("E04 periodic R2 B", state.magneticField().getView()),
                    localBounds(state), step};
        }

        /// Sum native squared differences [salari2000CodeVerificationMms, §4.1, Eq. (15)].
        static std::array<double, 2> differenceSquares(Solver& state, const Frame& frame) {
            const auto e = state.electricField().getView(), b = state.magneticField().getView();
            const auto periodicE = frame.e_m, periodicB = frame.b_m;
            std::array<double, 2> squares{};
            Kokkos::parallel_reduce(
                "E04 periodic field differences", state.A_n.getFieldRangePolicy(),
                KOKKOS_LAMBDA(int i, int j, int k, double& eSum, double& bSum) {
                    for (unsigned c = 0; c < 3; ++c) {
                        const double de = e(i, j, k)[c] - periodicE(i, j, k)[c];
                        const double db = b(i, j, k)[c] - periodicB(i, j, k)[c];
                        eSum += de * de;
                        bSum += db * db;
                    }
                },
                squares[0], squares[1]);
            return squares;
        }

        std::array<Frame, 2> frames_m;
    };
}  // namespace chdr::experiments::e04
