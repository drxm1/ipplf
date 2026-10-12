/** @file WaveNorms.h
 * @brief Continuum error measurements, kept separate from exact-identity monitors.
 * @details ChDRFieldContainer and analytic references are compared at the same component
 * position and declared time. The references follow
 * [fallahi2020mithra20fullwavesimulation, Eqs. (3.6)-(3.9)]; half-step labels
 * follow [christlieb2024GaugeConserving, Eqs. (25)-(29)]. Each returned norm states its own
 * denominator.
 *
 * TODO: REVIEW.
 */
#pragma once
#include <array>
#include <cmath>
#include <utility>

#include "ChdrSolver.h"
#include "SpatialDimension.h"
#include "Waves.h"

static_assert(chdr::solver::SpatialDim == 3,
              "WaveNorms.h requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver::test {
    struct WaveErrors {
        double electric_m, magnetic_m, gauss_m, gauge_m;
    };

    /** @brief Difference used by the continuum norm, independently of Monitors.
     * @details Collocated divergence uses the centred difference; the staggered
     * face-to-centre divergence uses the backward compact difference. The latter
     * samples each face component at the centre where charge is defined.
     * The centred difference is the direct symmetric Taylor quotient
     * [fallahi2020mithra20fullwavesimulation, Eqs. (3.13)-(3.14)]. The compact
     * quotient is a deduction from the face/centre placement and discrete
     * divergence in [ryu2016potentialfdtd, Appendix A, Eqs. (19)-(24)].
     * @tparam FieldView Three-dimensional view with a component-indexable value at each index.
     * @param fieldView E or four-potential samples with valid neighbours for the chosen difference.
     * @param i,j,k Local array indices including the halo offset.
     * @param component Stored component: E uses 0 through SpatialDim-1; A uses 1 through
     * SpatialDim.
     * @param direction Spatial derivative direction.
     * @param spacing Positive mesh spacing in that direction.
     * @param backward Select face-to-centre differences; otherwise use centred differences.
     */
    template <class FieldView>
    KOKKOS_FUNCTION double normDifference(FieldView fieldView, int i, int j, int k, int component,
                                          int direction, double spacing, bool backward) {
        int high[SpatialDim] = {i, j, k}, low[SpatialDim] = {i, j, k};
        if (backward)
            --low[direction];
        else {
            ++high[direction];
            --low[direction];
        }
        return (fieldView(high[0], high[1], high[2])[component]
                - fieldView(low[0], low[1], low[2])[component])
               / (spacing * (backward ? 1 : 2));
    }

    namespace waveNormDetail {
        using EMView        = decltype(std::declval<EMFieldsType&>().getView());
        using PotentialView = decltype(std::declval<AJFourFieldsType&>().getView());

        /** @brief Device-copyable measurement inputs; no solver pointer is captured.
         * @details A view copy retains storage ownership. The object is assembled
         * on the host after E halos are current, then copied into each reduction.
         * R1 has no single E time: its prescribed control uses the vector-potential
         * midpoint. Other reconstructions use the same midpoint as their actual label.
         * current_m and previous_m view A_n and A_nm1. Each element stores scalar
         * potential at component 0 and vector potential at components 1 through
         * SpatialDim. LevelContract supplies their time labels and component positions.
         */
        struct Samples {
            EMView electric_m, magnetic_m;
            PotentialView current_m, previous_m;
            ippl::Vector<int, SpatialDim> first_m;
            ippl::Vector<double, SpatialDim> spacing_m;
            LevelContract contract_m;
            LorenzWave wave_m;
            double dt_m, electricTime_m, magneticTime_m;

            template <class SolverType>
            Samples(SolverType& solver, const LorenzWave& wave)
                : electric_m(solver.electricField().getView())
                , magnetic_m(solver.magneticField().getView())
                , current_m(solver.A_n.getView())
                , previous_m(solver.A_nm1.getView())
                , first_m(solver.A_n.getLayout().getLocalNDIndex().first())
                , spacing_m(solver.A_n.get_mesh().getMeshSpacing())
                , contract_m(solver.levelContract())
                , wave_m(wave)
                , dt_m(solver.getDt())
                , electricTime_m(
                      solver.ebReconstructionConvention() == EBReconstructionConvention::R1
                          ? (solver.completedSolves() - 0.5) * dt_m
                          : solver.labelTime(Quantity::Electric, solver.completedSolves()))
                , magneticTime_m(solver.labelTime(Quantity::Magnetic, solver.completedSolves())) {
                const auto origin = solver.A_n.get_mesh().getOrigin();
                if (solver.A_n.getNghost() != 1 || origin[0] != 0 || origin[1] != 0
                    || origin[2] != 0)
                    throw std::invalid_argument(
                        "Wave norms require origin zero and one ghost layer");
            }

            KOKKOS_INLINE_FUNCTION ippl::Vector<double, SpatialDim> centre(int i, int j,
                                                                           int k) const {
                return {(i + first_m[0] - 0.5) * spacing_m[0],
                        (j + first_m[1] - 0.5) * spacing_m[1],
                        (k + first_m[2] - 0.5) * spacing_m[2]};
            }

            /// Compare each vector component at its own spatial offset and time.
            KOKKOS_INLINE_FUNCTION double fieldTerm(int metric, int c, int i, int j, int k) const {
                const bool electric = metric < 2;
                const auto x        = centre(i, j, k);
                const auto offset   = contract_m.positionOffset(
                    electric ? Quantity::Electric : Quantity::Magnetic, c);
                const double reference =
                    electric ? wave_m.electric(x[0] + offset[0] * spacing_m[0],
                                               x[1] + offset[1] * spacing_m[1],
                                               x[2] + offset[2] * spacing_m[2], electricTime_m, c)
                             : wave_m.magnetic(x[0] + offset[0] * spacing_m[0],
                                               x[1] + offset[1] * spacing_m[1],
                                               x[2] + offset[2] * spacing_m[2], magneticTime_m, c);
                const double value = electric ? electric_m(i, j, k)[c] : magnetic_m(i, j, k)[c];
                return metric % 2 == 0 ? value - reference : reference;
            }

            /// Keep the original x/y/z accumulation order for each divergence.
            KOKKOS_INLINE_FUNCTION Kokkos::Array<double, 3> divergences(int i, int j, int k) const {
                double electric = 0, current = 0, previous = 0;
                const bool backward =
                    contract_m.spatialDiscretization_m == SpatialDiscretization::Staggered;
                for (int d = 0; d < static_cast<int>(SpatialDim); ++d) {
                    electric += normDifference(electric_m, i, j, k, d, d, spacing_m[d], backward);
                    current += normDifference(current_m, i, j, k, d + 1, d, spacing_m[d], backward);
                    previous +=
                        normDifference(previous_m, i, j, k, d + 1, d, spacing_m[d], backward);
                }
                return {electric, current, previous};
            }

            /** @brief Gauss and gauge samples at their respective contract labels.
             * @details Charge is sampled at E's comparison time, including the R1
             * control. Gauge uses the contract's vector interpolation weight.
             * This is a continuum-error diagnostic, not the type-A Gauss identity.
             * Reference: christlieb2024GaugeConserving, Eqs. (25)-(29).
             */
            KOKKOS_INLINE_FUNCTION double constraintTerm(int metric, int i, int j, int k) const {
                const auto x        = centre(i, j, k);
                const auto div      = divergences(i, j, k);
                const double weight = contract_m.gaugeWeight();
                const double gauss  = div[0] - wave_m.source(x[0], x[1], x[2], electricTime_m, 0);
                const double gauge  = (current_m(i, j, k)[0] - previous_m(i, j, k)[0]) / dt_m
                                     + weight * div[1] + (1 - weight) * div[2];
                return metric == 4 ? gauss : metric == 5 ? div[1] : metric == 6 ? gauge : div[2];
            }
        };

        /// One scalar reduction per metric, using the original component sum order.
        struct Reduction {
            Samples samples_m;
            int metric_m;
            KOKKOS_INLINE_FUNCTION void operator()(int i, int j, int k, double& result) const {
                if (metric_m < 4) {
                    for (int c = 0; c < static_cast<int>(SpatialDim); ++c) {
                        const double term = samples_m.fieldTerm(metric_m, c, i, j, k);
                        result += term * term;
                    }
                } else {
                    const double term = samples_m.constraintTerm(metric_m, i, j, k);
                    result += term * term;
                }
            }
        };

        template <class SolverType>
        std::array<double, 8> globalSums(SolverType& solver, const Samples& samples) {
            std::array<double, 8> sums{};
            for (int metric = 0; metric < 8; ++metric)
                Kokkos::parallel_reduce("continuum L2", solver.A_n.getFieldRangePolicy(),
                                        Reduction{samples, metric}, sums[metric]);
            MPI_Allreduce(MPI_IN_PLACE, sums.data(), 8, MPI_DOUBLE, MPI_SUM,
                          solver.A_n.getLayout().comm);
            return sums;
        }
    }  // namespace waveNormDetail

    /** @brief Relative E/B L2 errors and the frozen Gauss/gauge norm ratios.
     * @details The common cell volume cancels on each uniform grid. Gauss is
     * normalized by the continuum field norm times the wave-number magnitude;
     * gauge by the numerical vector-potential divergence norm.
     * Only owned cells enter any sum; ghost values are used solely by differences.
     */
    template <class SolverType>
    WaveErrors waveErrors(SolverType& solver, const LorenzWave& wave) {
        solver.electricField().fillHalo();
        const auto sums = waveNormDetail::globalSums(solver, waveNormDetail::Samples(solver, wave));
        double kn2      = 0;
        for (int d = 0; d < static_cast<int>(SpatialDim); ++d)
            kn2 += wave.waveNumber_m[d] * wave.waveNumber_m[d];
        return {std::sqrt(sums[0] / sums[1]), std::sqrt(sums[2] / sums[3]),
                std::sqrt(sums[4] / (kn2 * sums[1])), std::sqrt(sums[6] / sums[5])};
    }
}  // namespace chdr::solver::test
