/** @file ComovingPotentials.h
 * @brief Analytic periodic co-moving potentials and independent field-error symbols.
 */
#pragma once

#include <Kokkos_MathematicalConstants.hpp>
#include <array>
#include <cmath>
#include <complex>
#include <map>
#include <tuple>
#include <vector>

#include "core/base/SpatialDimension.h"
#include "core/solver/BoxSetup.h"
#include "testing/Util.h"

static_assert(
    chdr::solver::SpatialDim == 3,
    "ComovingPotentials.h requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver::test {
    using Complex = std::complex<double>;
    using chdr::testing::util::preparePotentialHalo;

    /** @brief Finite Fourier series of one mode or a transversely neutralized Gaussian train.
     * @pre The test meshes have even cell counts and resolve the selected single mode.
     * @details Use normalized units \f$c=\epsilon_0=\mu_0=1\f$ and modes
     * \f$e^{i(\mathbf k\cdot\mathbf x-\omega t)}\f$. Rigid motion along z gives
     * \f$\omega=\beta k_z\f$, \f$A_z=\beta\phi\f$ and
     * \f$\widehat\phi=\widehat\rho/(k_x^2+k_y^2+k_z^2/\gamma^2)\f$.
     * This is a deduction from the four-potential boost
     * [fitzpatrick2014ClassicalElectromagnetism, Sec. 12.18, Eqs. (12.196)-(12.201)]
     * and the wave equations [fallahi2020mithra20fullwavesimulation, Eqs. (3.6)-(3.7)].
     * Substitution in \f$\mathbf E=-\partial_t\mathbf A-\nabla\phi\f$ gives
     * \f$E_x=-\partial_x\phi\f$ and \f$E_z=-\gamma^{-2}\partial_z\phi\f$:
     * the longitudinal component is the remainder of two much larger terms
     * [fallahi2020mithra20fullwavesimulation, Eq. (3.9); deduction].
     * In the Gaussian branch, the density has coefficients
     * \f$(Q/V)e^{-\sum_d k_d^2\sigma_d^2/2}e^{-i\mathbf k\cdot\mathbf x_c}\f$,
     * where Q is the unit bunch charge before neutralization and the bunch is
     * centred in the periodic box. Removing every
     * mode with \f$k_x=k_y=0\f$ subtracts the transverse mean of each longitudinal
     * slice: the bunch is accompanied by a co-moving neutralizing sheet. This
     * defines the test density rather than approximating an isolated free-space bunch.
     * Nyquist planes are omitted so retained modes have distinct conjugate partners
     * at both centre and face sample positions.
     */
    class ComovingPotentials {
    public:
        static constexpr double GAMMA = 118.417071;
        static constexpr double PI    = Kokkos::numbers::pi_v<double>;

        ComovingPotentials(std::array<int, SpatialDim> cells,
                           std::array<double, SpatialDim> spacing, bool singleMode)
            : cells_m(cells)
            , spacing_m(spacing)
            , beta_m(std::sqrt(1.0 - 1.0 / (GAMMA * GAMMA))) {
            for (unsigned d = 0; d < SpatialDim; ++d)
                initializeAxis(d);
            initializeCoefficients(singleMode);
            if (singleMode)
                requirePeriodicWave({2.0 * PI / lengths_m[0], 0.0, 6.0 * PI / lengths_m[2]},
                                    lengths_m);
        }

        /// Require integer box modes so periodic halo wrapping reproduces the same analytic wave.
        static void requirePeriodicWave(std::array<double, SpatialDim> wave,
                                        std::array<double, SpatialDim> lengths) {
            for (unsigned d = 0; d < SpatialDim; ++d) {
                const double mode = wave[d] * lengths[d] / (2.0 * PI);
                if (!std::isfinite(mode) || std::abs(mode - std::round(mode)) > 1e-12)
                    throw std::invalid_argument("Co-moving mode does not fit the periodic box");
            }
        }

        std::size_t index(int i, int j, int k) const {
            return (std::size_t(i) * cells_m[1] + j) * cells_m[2] + k;
        }
        std::array<double, SpatialDim> waveVector(int i, int j, int k) const {
            return {waveNumbers_m[0][i], waveNumbers_m[1][j], waveNumbers_m[2][k]};
        }
        const std::array<int, SpatialDim>& cells() const { return cells_m; }
        const std::array<double, SpatialDim>& spacing() const { return spacing_m; }
        const std::array<double, SpatialDim>& lengths() const { return lengths_m; }
        const std::vector<Complex>& coefficients() const { return phi_m; }
        double beta() const { return beta_m; }

        /** @brief Visit Fourier-array indices in a fixed order shared by preparation and checks.
         * @tparam ModeVisitor Host callable accepting three integer Fourier-array indices.
         * @param visitMode Operation called once for every stored mode, including zero
         * coefficients.
         */
        template <class ModeVisitor>
        void forEachMode(ModeVisitor visitMode) const {
            for (int i = 0; i < cells_m[0]; ++i)
                for (int j = 0; j < cells_m[1]; ++j)
                    for (int k = 0; k < cells_m[2]; ++k)
                        visitMode(i, j, k);
        }

        /** @brief Evaluate the finite series by three separable direct host sums.
         * @details The precomputed basis uses cell centres; time and component offsets
         * supply additional Fourier phases. This generates analytic data for the test,
         * while the measured fields always come from the actual reconstruction kernel.
         * The z/y/x transform order and each inner summation order are fixed.
         */
        std::vector<double> evaluate(const std::vector<Complex>& coefficients, double time,
                                     std::array<double, SpatialDim> offset = {}) const {
            auto values = phasedCoefficients(coefficients, time, offset);
            for (int axis = 2; axis >= 0; --axis)
                transformAxis(values, axis);
            std::vector<double> real(values.size());
            for (std::size_t i = 0; i < real.size(); ++i)
                real[i] = values[i].real();
            return real;
        }

        /// Repeated slot labels reuse identical doubles, making R2/R2Prime equality a bitwise test.
        const std::vector<double>& potential(double time,
                                             std::array<double, SpatialDim> offset = {}) {
            const auto key = std::make_tuple(time, offset[0], offset[1], offset[2]);
            auto found     = samples_m.find(key);
            if (found == samples_m.end())
                found = samples_m.emplace(key, evaluate(phi_m, time, offset)).first;
            return found->second;
        }

        /** @brief Load both actual history slots at their declared times and positions.
         * @details A_z is sampled on its face for the staggered reading. Scalar phi
         * stays at the centre and takes the reconstruction's own half-step label.
         * Both histories receive periodic halos before reconstructEB() reads them.
         * @tparam SolverType Solver exposing A_n, A_nm1 and their LevelContract.
         * @param solver Destination with initialized potential storage and a valid timestep.
         */
        template <class SolverType>
        void load(SolverType& solver) {
            const auto levels = solver.levelContract();
            const double dt   = solver.getDt();
            const auto scalar = levels.positionOffset(Quantity::Phi, 0);
            const auto vector = levels.positionOffset(Quantity::VectorPotential, 2);
            const std::array<double, SpatialDim> scalarOffset{scalar[0], scalar[1], scalar[2]};
            const std::array<double, SpatialDim> vectorOffset{vector[0], vector[1], vector[2]};
            const double phiTime = levels.labelTime(Quantity::Phi, 0, dt);
            loadSlot(solver.A_n, potential(phiTime, scalarOffset), potential(0.0, vectorOffset));
            loadSlot(solver.A_nm1, potential(phiTime - dt, scalarOffset),
                     potential(-dt, vectorOffset));
            preparePotentialHalo(solver.A_n);
            preparePotentialHalo(solver.A_nm1);
            solver.markHistoryLoaded();
        }

    private:
        int signedMode(int index, unsigned axis) const {
            return index <= (cells_m[axis] - 1) / 2 ? index : index - cells_m[axis];
        }

        /// Store the signed wave numbers and centre-sampled basis for one Cartesian direction.
        void initializeAxis(unsigned d) {
            const int count = cells_m[d];
            lengths_m[d]    = count * spacing_m[d];
            waveNumbers_m[d].resize(count);
            basis_m[d].resize(count * count);
            for (int frequency = 0; frequency < count; ++frequency) {
                const int mode              = signedMode(frequency, d);
                waveNumbers_m[d][frequency] = 2.0 * PI * mode / lengths_m[d];
                for (int index = 0; index < count; ++index)
                    basis_m[d][index * count + frequency] =
                        std::polar(1.0, 2.0 * PI * mode * (index + 0.5) / count);
            }
        }

        /// The centred Gaussian's spatial phase is a real sign for integer box modes.
        double gaussianCoefficient(std::array<double, SpatialDim> wave,
                                   std::array<int, SpatialDim> mode) const {
            const std::array<double, SpatialDim> sigma{3.0 * spacing_m[0], 3.0 * spacing_m[1],
                                                       4.0 * spacing_m[2]};
            const double volume = lengths_m[0] * lengths_m[1] * lengths_m[2];
            double exponent     = 0.0;
            for (unsigned d = 0; d < SpatialDim; ++d)
                exponent += wave[d] * wave[d] * sigma[d] * sigma[d];
            const double denominator =
                wave[0] * wave[0] + wave[1] * wave[1] + wave[2] * wave[2] / (GAMMA * GAMMA);
            const double centreSign = (mode[0] + mode[1] + mode[2]) % 2 == 0 ? 1.0 : -1.0;
            return centreSign * std::exp(-0.5 * exponent) / (volume * denominator);
        }

        /// Apply the neutralizing-sheet and conjugate-mode conventions described by this class.
        double coefficient(int i, int j, int k, bool singleMode) const {
            const std::array<int, SpatialDim> mode{signedMode(i, 0), signedMode(j, 1),
                                                   signedMode(k, 2)};
            if (singleMode)
                return j == 0
                               && ((mode[0] == 1 && mode[2] == 3)
                                   || (mode[0] == -1 && mode[2] == -3))
                           ? 0.5
                           : 0.0;
            if ((i == 0 && j == 0) || i == cells_m[0] / 2 || j == cells_m[1] / 2
                || k == cells_m[2] / 2)
                return 0.0;
            return gaussianCoefficient(waveVector(i, j, k), mode);
        }
        void initializeCoefficients(bool singleMode) {
            phi_m.resize(std::size_t(cells_m[0]) * cells_m[1] * cells_m[2], 0.0);
            forEachMode([&](int i, int j, int k) {
                phi_m[index(i, j, k)] = coefficient(i, j, k, singleMode);
            });
        }

        /// Apply time and position phases to coefficients before any spatial summation.
        std::vector<Complex> phasedCoefficients(const std::vector<Complex>& coefficients,
                                                double time,
                                                std::array<double, SpatialDim> offset) const {
            auto values = coefficients;
            forEachMode([&](int i, int j, int k) {
                const auto wave = waveVector(i, j, k);
                double phase    = -wave[2] * beta_m * time;
                for (unsigned d = 0; d < SpatialDim; ++d)
                    phase += wave[d] * spacing_m[d] * offset[d];
                values[index(i, j, k)] *= std::polar(1.0, phase);
            });
            return values;
        }

        /// Preserve the direct sum's frequency order and complex multiply/add grouping.
        Complex transformedValue(const std::vector<Complex>& values, int axis,
                                 std::array<int, SpatialDim> position) const {
            const int target = position[axis];
            Complex sum(0.0, 0.0);
            for (int frequency = 0; frequency < cells_m[axis]; ++frequency) {
                position[axis] = frequency;
                sum += values[index(position[0], position[1], position[2])]
                       * basis_m[axis][target * cells_m[axis] + frequency];
            }
            return sum;
        }
        void transformAxis(std::vector<Complex>& values, int axis) const {
            std::vector<Complex> scratch(values.size());
            forEachMode([&](int i, int j, int k) {
                scratch[index(i, j, k)] = transformedValue(values, axis, {i, j, k});
            });
            values.swap(scratch);
        }

        /** @brief Fill scalar and longitudinal vector components in a host potential slot.
         * @param potentialA Field supplying the owned-cell layout and halo width.
         * @param hostPotential Writable host allocation for that potential slot.
         * @param phi Scalar-potential samples at the scalar component's time and position.
         * @param phiAtAz Scalar-series samples at the longitudinal vector component's
         * time and position; multiplication by beta follows the class's boost relation.
         */
        void fillOwnedSlot(const AJFourFieldsType& potentialA,
                           AJFourFieldsType::host_mirror_type hostPotential,
                           const std::vector<double>& phi,
                           const std::vector<double>& phiAtAz) const {
            const auto local = potentialA.getLayout().getLocalNDIndex();
            const int ghost  = potentialA.getNghost();
            for (int i = local[0].first(); i <= local[0].last(); ++i)
                for (int j = local[1].first(); j <= local[1].last(); ++j)
                    for (int k = local[2].first(); k <= local[2].last(); ++k) {
                        auto& value = hostPotential(i - local[0].first() + ghost,
                                                    j - local[1].first() + ghost,
                                                    k - local[2].first() + ghost);
                        value[0]    = phi[index(i, j, k)];
                        value[3]    = beta_m * phiAtAz[index(i, j, k)];
                    }
        }
        /// Initialize one potential slot; fillOwnedSlot defines the two sample arrays.
        void loadSlot(AJFourFieldsType& potentialA, const std::vector<double>& phi,
                      const std::vector<double>& phiAtAz) const {
            auto host = Kokkos::create_mirror_view(potentialA.getView());
            Kokkos::deep_copy(host, ippl::Vector<double, 4>(0.0));
            fillOwnedSlot(potentialA, host, phi, phiAtAz);
            Kokkos::deep_copy(potentialA.getView(), host);
        }

        std::array<int, SpatialDim> cells_m;
        std::array<double, SpatialDim> spacing_m, lengths_m;
        double beta_m;
        std::array<std::vector<double>, SpatialDim> waveNumbers_m;
        std::array<std::vector<Complex>, SpatialDim> basis_m;
        std::vector<Complex> phi_m;
        std::map<std::tuple<double, double, double, double>, std::vector<double>> samples_m;
    };

    struct ComovingSymbols {
        std::vector<Complex> analytic_m, exactError_m, leadingError_m;
    };
    struct ComovingFactors {
        Complex time_m, scalar_m, exSeries_m;
        double c1_m, c2_m;
    };

    /// R2's time term is centred at the output label; its spatial coefficient depends on reading.
    inline ComovingFactors centredFactors(std::array<double, SpatialDim> mode,
                                          std::array<double, SpatialDim> spacing, double dt,
                                          double beta,
                                          SpatialDiscretization spatialDiscretization) {
        const double gamma2 = ComovingPotentials::GAMMA * ComovingPotentials::GAMMA;
        const double spatialSeries =
            spatialDiscretization == SpatialDiscretization::Staggered ? 1.0 / 24.0 : 1.0 / 6.0;
        const double frequency = beta * mode[2];
        return {
            {2.0 * std::sin(frequency * dt / 2.0) / dt, 0.0},
            {1.0, 0.0},
            {-spatialSeries * std::pow(mode[0] * spacing[0], 2), 0.0},
            0.0,
            gamma2
                * (std::pow(beta, 4) * dt * dt / 24.0 - spatialSeries * spacing[2] * spacing[2])};
    }

    /** @brief Relabel the two R1 contributions independently.
     * @details At the whole step the backward A difference carries a half-step
     * phase; at the half step the scalar gradient carries the opposite phase.
     * With the negative time phase used by ComovingPotentials, these factors are
     * \f$e^{+i\omega\Delta t/2}\f$ for the time difference at the whole step and
     * \f$e^{-i\omega\Delta t/2}\f$ for the scalar contribution at the half step.
     * Expanding these factors together with the difference symbols produces the
     * leading-error coefficients below. For transverse R1 at the half step,
     * exSeries_m keeps the first-order phase term and the second-order spatial
     * term, but omits the second-order temporal phase term; it is a mixed-order
     * leading approximation, not a full second-order Taylor polynomial.
     * This is a deduction from the reconstructed
     * fields [fallahi2020mithra20fullwavesimulation, Eqs. (3.50)-(3.55)].
     */
    inline void applyR1Label(ComovingFactors& factors, double kz, double hz, double dt, double beta,
                             bool halfLabel) {
        const Complex imaginary(0.0, 1.0);
        const double gamma2    = ComovingPotentials::GAMMA * ComovingPotentials::GAMMA,
                     frequency = beta * kz;
        if (halfLabel) {
            factors.scalar_m = std::polar(1.0, -frequency * dt / 2.0);
            factors.c1_m     = -gamma2 * beta * dt / 2.0;
            factors.c2_m -= gamma2 * beta * beta * dt * dt / 8.0;
            factors.exSeries_m -= imaginary * frequency * dt / 2.0;
        } else {
            factors.time_m *= std::polar(1.0, frequency * dt / 2.0);
            factors.c1_m = -gamma2 * std::pow(beta, 3) * dt / 2.0;
            factors.c2_m = gamma2 * (std::pow(beta, 4) * dt * dt - hz * hz) / 6.0;
        }
    }

    struct ComovingModeSymbols {
        Complex analytic_m, exactError_m, leadingError_m;
    };

    /** @brief Evaluate exact-symbol and leading-order errors for one Fourier coefficient.
     * @pre Component is zero for the transverse electric field or two for the longitudinal field.
     * @details These predictions use analytic coefficients and difference symbols;
     * no measured field value enters the calculation. comovingSymbols documents
     * the symbols and the interpretation of both error arrays.
     */
    inline ComovingModeSymbols modeSymbols(const ComovingPotentials& wave, std::size_t index,
                                           std::array<double, SpatialDim> mode, double dt,
                                           EBReconstructionConvention ebReconstructionConvention,
                                           SpatialDiscretization spatialDiscretization,
                                           unsigned component, bool r1HalfLabel) {
        const Complex imaginary(0.0, 1.0);
        const auto spacing  = wave.spacing();
        const double beta   = wave.beta(),
                     gamma2 = ComovingPotentials::GAMMA * ComovingPotentials::GAMMA;
        auto factors        = centredFactors(mode, spacing, dt, beta, spatialDiscretization);
        if (ebReconstructionConvention == EBReconstructionConvention::R1)
            applyR1Label(factors, mode[2], spacing[2], dt, beta, r1HalfLabel);
        const double spaceSymbol =
            spatialDiscretization == SpatialDiscretization::Staggered
                ? 2.0 * std::sin(mode[component] * spacing[component] / 2.0) / spacing[component]
                : std::sin(mode[component] * spacing[component]) / spacing[component];
        const auto phi = wave.coefficients()[index];
        const Complex analytic =
            -imaginary * mode[component] * phi / (component == 2 ? gamma2 : 1.0);
        const Complex reconstructed =
            component == 2
                ? imaginary * (beta * factors.time_m - spaceSymbol * factors.scalar_m) * phi
                : -imaginary * spaceSymbol * factors.scalar_m * phi;
        return {analytic, reconstructed - analytic,
                analytic
                    * (component == 2
                           ? imaginary * factors.c1_m * mode[2] + factors.c2_m * mode[2] * mode[2]
                           : factors.exSeries_m)};
    }

    /** @brief Produce analytic fields and two independently predicted error series.
     * @details The collocated derivative factor is \f$q_d=\sin(k_dh_d)/h_d\f$;
     * the compact factor at the component's face is
     * \f$q_d=2\sin(k_dh_d/2)/h_d\f$. The one-step time factor at the half step is
     * \f$W=2\sin(\omega\Delta t/2)/\Delta t\f$
     * [fallahi2020mithra20fullwavesimulation, Eqs. (3.50)-(3.55);
     * ryu2016potentialfdtd, App. A, Eqs. (19)-(24)].
     * For R2, the reconstructed coefficients are \f$-iq_x\widehat\phi\f$ and
     * \f$i(\beta W-q_z)\widehat\phi\f$; applyR1Label supplies the phase factors
     * when the two contributions use different times. Subtracting the analytic
     * coefficients gives exactError_m, called P1 in the test records.
     * leadingError_m, called P2, keeps the leading Taylor terms of these symbols.
     * In particular, \f$q_d/k_d-1\f$ has leading term \f$-(k_dh_d)^2/6\f$ at
     * centres and \f$-(k_dh_d)^2/24\f$ at faces, while \f$W/\omega-1\f$ has
     * leading term \f$-(\omega\Delta t)^2/24\f$. The ratios use their continuous
     * limits at zero wave number or frequency. These are expansions of the
     * individual difference symbols; applyR1Label also expands the R1 phase factors
     * before the field-error coefficients are formed. They are deductions from the
     * displayed symbols; no reconstructed field or numerical difference is reused.
     */
    inline ComovingSymbols comovingSymbols(const ComovingPotentials& wave, double dt,
                                           EBReconstructionConvention ebReconstructionConvention,
                                           SpatialDiscretization spatialDiscretization,
                                           unsigned component, bool r1HalfLabel) {
        ComovingSymbols result;
        const auto size = wave.coefficients().size();
        result.analytic_m.resize(size);
        result.exactError_m.resize(size);
        result.leadingError_m.resize(size);
        wave.forEachMode([&](int i, int j, int k) {
            const auto index = wave.index(i, j, k);
            const auto mode =
                modeSymbols(wave, index, wave.waveVector(i, j, k), dt, ebReconstructionConvention,
                            spatialDiscretization, component, r1HalfLabel);
            result.analytic_m[index]     = mode.analytic_m;
            result.exactError_m[index]   = mode.exactError_m;
            result.leadingError_m[index] = mode.leadingError_m;
        });
        return result;
    }
}  // namespace chdr::solver::test
