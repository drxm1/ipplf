/** @file E15Measurements.h
 * @brief Electric and magnetic projections and the field left after the fits (E15.meas01-03).
 */
#pragma once
#include <complex>
#include <mpi.h>
#include <nlohmann/json.hpp>

#include "core/simulation/ChDRSimulationManager.h"
#include "experiments/E15_homogeneous_long_run_stability/E15Wave.h"

namespace chdr::experiments::e15 {
    using Simulation = simulation::ChDRSimulationManager<>;
    using Solver     = Simulation::Solver;

    /** @brief Complex spatial basis at the native component locations.
     * @details Equal wave-vector components permit summing the three coordinates before
     * evaluating the phase. Factor out the common temporal phase for both fits.
     * This avoids subtracting a large accumulated phase from a small spatial phase in
     * the residual. The temporal factor is restored in the reported coefficients.
     * See E15.typ, "Measured quantities", and Schneider (2010), section 7.3.
     */
    struct Basis {
        Wave wave_m;
        int indexShift_m;  ///< Sum of global lower indices less the three ghost offsets.
        Kokkos::Array<double, 3> offsets_m;
        bool electric_m;

        KOKKOS_INLINE_FUNCTION auto operator()(int i, int j, int k, unsigned d) const {
            const double angle =
                TwoPi / wave_m.cells_m * (i + j + k + indexShift_m + 1.5 + offsets_m[d]);
            const double amplitude = wave_m.amplitude(electric_m, d);
            return Kokkos::complex<double>(-amplitude * Kokkos::sin(angle),
                                           amplitude * Kokkos::cos(angle));
        }
    };

    template <class Field>
    Basis makeBasis(Field& field, LevelContract contract, const Wave& wave, Quantity quantity) {
        const auto first = field.getLayout().getLocalNDIndex().first();
        Basis basis{wave,
                    first[0] + first[1] + first[2] - 3 * field.getNghost(),
                    {},
                    quantity == Quantity::Electric};
        for (unsigned d = 0; d < 3; ++d) {
            const auto offset  = contract.positionOffset(quantity, d);
            basis.offsets_m[d] = offset[0] + offset[1] + offset[2];
        }
        return basis;
    }

    /// Project owned components onto the complex spatial mode and normalize after MPI summation.
    template <class Field>
    Kokkos::complex<double> fit(Field& field, Basis basis) {
        const auto view = field.getView();
        double sums[3]{};
        Kokkos::parallel_reduce(
            "E15 mode projection", field.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k, double& re, double& im, double& norm) {
                for (unsigned d = 0; d < 3; ++d) {
                    const auto mode    = basis(i, j, k, d);
                    const double value = view(i, j, k)[d];
                    re += value * mode.real();
                    im -= value * mode.imag();
                    norm += Kokkos::norm(mode);
                }
            },
            sums[0], sums[1], sums[2]);
        MPI_Allreduce(MPI_IN_PLACE, sums, 3, MPI_DOUBLE, MPI_SUM, field.getLayout().comm);
        return {2 * sums[0] / sums[2], 2 * sums[1] / sums[2]};
    }

    /// Sum the pointwise residual squared directly to retain errors near machine precision.
    template <class Field>
    double residualSquares(Field& field, Basis basis, Kokkos::complex<double> coefficient) {
        const auto view = field.getView();
        double sum      = 0;
        Kokkos::parallel_reduce(
            "E15 residual", field.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k, double& squares) {
                for (unsigned d = 0; d < 3; ++d) {
                    const double difference =
                        view(i, j, k)[d] - (coefficient * basis(i, j, k, d)).real();
                    squares += difference * difference;
                }
            },
            sum);
        return sum;
    }

    /** @brief Apply the two fits and restore the phase of the continuum reference.
     * @details In E15.meas03 the reference squared norm sums to
     * \f$N^3\epsilon_r\omega^2\f$: both reference fields have equal energy,
     * with mean squared sinusoid one half. See E15.typ, "Residual outside the wave".
     */
    inline nlohmann::json observe(Solver& solver, Wave wave, long long step) {
        auto& e       = solver.electricField();
        auto& b       = solver.magneticField();
        const auto eb = makeBasis(e, solver.levelContract(), wave, Quantity::Electric);
        const auto bb = makeBasis(b, solver.levelContract(), wave, Quantity::Magnetic);
        const auto ce = fit(e, eb), cb = fit(b, bb);
        double squares = Epsilon * residualSquares(e, eb, ce) + residualSquares(b, bb, cb);
        MPI_Allreduce(MPI_IN_PLACE, &squares, 1, MPI_DOUBLE, MPI_SUM, e.getLayout().comm);
        const int n = wave.cells_m;
        const double residual =
            std::sqrt(squares / (n * n * n * Epsilon * wave.omega_m * wave.omega_m));
        if (!std::isfinite(residual))
            throw std::runtime_error("Nonfinite E15 field or fit");
        const double time   = solver.labelTime(Quantity::Electric, step);
        const auto rotation = std::polar(1.0, wave.omega_m * time);
        const auto electric = std::complex<double>(ce.real(), ce.imag()) * rotation;
        const auto magnetic = std::complex<double>(cb.real(), cb.imag()) * rotation;
        return {{"n", step},
                {"t", time},
                {"c_E_re", electric.real()},
                {"c_E_im", electric.imag()},
                {"c_B_re", magnetic.real()},
                {"c_B_im", magnetic.imag()},
                {"r", residual}};
    }
}  // namespace chdr::experiments::e15
