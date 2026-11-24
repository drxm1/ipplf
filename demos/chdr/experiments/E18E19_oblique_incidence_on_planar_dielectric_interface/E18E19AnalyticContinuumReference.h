/** @file E18E19AnalyticContinuumReference.h
 * @brief Continuum TE and TM packets at the dielectric face of E18 and E19.
 */
#pragma once

#include <Kokkos_MathematicalConstants.hpp>
#include <array>
#include <catalyst_conduit.hpp>
#include <cmath>
#include <complex>
#include <string>
#include <vector>

namespace chdr::experiments::e18e19 {
    /** @brief One term of the compact spectrum prescribed in E18.typ.
     * @details The reflection coefficient multiplies the TE potential or the TM
     * magnetic amplitude. The transmission coefficient multiplies the potential.
     * The distinction follows E18-coefficients and E19-coefficients, derived from
     * [ellingson2020ElectromagneticsVol2, Eqs. (5.147), (5.184), (5.187)].
     */
    struct PacketMode {
        double q_m, omega_m, qSecond_m, amplitude_m, reflection_m, transmission_m;
    };

    /** @brief Evaluate the specified finite plane-wave sum on the host.
     * @details All waves have the same tangential wavenumber. Their complex
     * amplitudes therefore depend only on z and time, with the common factor
     * \f$\exp(i k_x x)\f$ applied after summation. This is an algebraic
     * factorization of E18-phases and E19-analytic-potentials.
     * The constructor normalizes the prescribed spectrum to unit total amplitude.
     * E and B follow by differentiation using [chew2014generalizedgauge,
     * Eqs. (5)--(6)]. No spatial or temporal discretization enters this reference.
     */
    class E18E19AnalyticContinuumReference {
    public:
        static constexpr double INTERFACE_Z         = 28;
        static constexpr double DIELECTRIC_EPSILON  = 2.13;
        static constexpr double BOX_LENGTH_Z        = 48;
        static constexpr double INITIAL_CENTER_Z    = 23;
        static constexpr double KX                  = 2 * Kokkos::numbers::pi;
        static constexpr double SPECTRUM_HALF_WIDTH = 0.85 * KX;
        static constexpr double SPECTRUM_SHARPNESS  = 8;
        static constexpr double TRAIN_PERIOD        = 4 * BOX_LENGTH_Z;
        static constexpr double MODE_SPACING        = 2 * Kokkos::numbers::pi / TRAIN_PERIOD;

        const bool tm_m;

        /// Construct E18 for false and E19 for true, with their common prescribed spectrum.
        explicit E18E19AnalyticContinuumReference(bool tm)
            : tm_m(tm) {
            const int halfCount = static_cast<int>(std::lround(SPECTRUM_HALF_WIDTH / MODE_SPACING));
            double normalization = 0;
            for (int index = -halfCount; index <= halfCount; ++index) {
                modes_m.push_back(makeMode(KX + index * MODE_SPACING));
                normalization += modes_m.back().amplitude_m;
            }
            for (auto& mode : modes_m)
                mode.amplitude_m /= normalization;
        }

        /// Return the normalized continuum modes used by initialization and modes.csv.
        const std::vector<PacketMode>& modes() const { return modes_m; }

        /// Return the eleven specified spectral probes followed by the Brewster probe.
        std::array<double, 12> measuredWavenumbers() const {
            std::array<double, 12> result;
            for (int index = -5; index <= 5; ++index)
                result[index + 5] = KX + index * SPECTRUM_HALF_WIDTH / 10;
            result.back() = KX / std::sqrt(DIELECTRIC_EPSILON);
            return result;
        }

        /** @brief Fresnel coefficient with the polarization convention of the specification.
         * @details Substitute \f$\eta_j=1/n_j\f$ and the normal wavenumbers into
         * [ellingson2020ElectromagneticsVol2, Eqs. (5.147), (5.184)]. E19 uses
         * \f$R=-\Gamma_{TM}\f$ because its reflected magnetic basis points along
         * positive y, as explained below E19-analytic-EB.
         */
        double reflection(double q) const {
            const double qSecond = std::sqrt(DIELECTRIC_EPSILON * (KX * KX + q * q) - KX * KX);
            const double scaledQ = tm_m ? DIELECTRIC_EPSILON * q : q;
            return (scaledQ - qSecond) / (scaledQ + qSecond);
        }

        /** @brief Energy-weighted reflected fraction of the entire continuum packet.
         * @details E18-energy-fractions and E19-energy-fractions weight each
         * coefficient squared by \f$a_m^2\omega_m^2\f$. The plane-wave power
         * fraction is given in [ellingson2020ElectromagneticsVol2, Example 5.7].
         */
        double reflectedEnergyFraction() const {
            double incident = 0, reflected = 0;
            for (const auto& mode : modes_m) {
                const double weight =
                    mode.amplitude_m * mode.amplitude_m * mode.omega_m * mode.omega_m;
                incident += weight;
                reflected += weight * mode.reflection_m * mode.reflection_m;
            }
            return reflected / incident;
        }

        /// Host potential sampler in storage order phi, Ax, Ay, Az.
        double referenceFourPotentialA(double x, double, double z, double time,
                                       unsigned component) const {
            if (component == 0)
                return 0;
            return sample(x, z, time, component - 1, Sample::Potential);
        }

        /// Host electric-field sampler at the caller's own component position and time.
        double referenceE(double x, double, double z, double time, unsigned component) const {
            return sample(x, z, time, component, Sample::Electric);
        }

        /// Host magnetic-field sampler at the caller's own component position and time.
        double referenceB(double x, double, double z, double time, unsigned component) const {
            return sample(x, z, time, component, Sample::Magnetic);
        }

        /** @brief Complex potential amplitude before applying the common x phase.
         * @details Initialization tabulates this value at the contract's z positions
         * and both stored times. Component is Cartesian, from zero through two.
         */
        std::complex<double> potentialAmplitude(double z, double time, unsigned component) const {
            return amplitude(z, time, component, Sample::Potential);
        }

        /// Preserve the prescribed geometry and continuum prediction beside the solver input.
        void record(conduit_cpp::Node node, const std::string& experiment) const {
            node["name"]                       = experiment + " oblique-incidence continuum packet";
            node["polarization"]               = tm_m ? "TM" : "TE";
            node["epsilon_side1"]              = 1.0;
            node["epsilon_side2"]              = DIELECTRIC_EPSILON;
            node["interface_z"]                = INTERFACE_Z;
            node["initial_center_z"]           = INITIAL_CENTER_Z;
            node["tangential_wavenumber"]      = KX;
            node["normal_wavenumber_center"]   = KX;
            node["spectrum_half_width"]        = SPECTRUM_HALF_WIDTH;
            node["spectrum_sharpness"]         = SPECTRUM_SHARPNESS;
            node["train_period"]               = TRAIN_PERIOD;
            node["mode_spacing"]               = MODE_SPACING;
            node["mode_count"]                 = static_cast<int>(modes_m.size());
            node["reflection_energy_fraction"] = reflectedEnergyFraction();
            node["transmission_energy_fraction"] = 1 - reflectedEnergyFraction();
        }

    private:
        enum class Sample {
            Potential,
            Electric,
            Magnetic
        };
        enum class Branch {
            Incident,
            Reflected,
            Transmitted
        };
        std::vector<PacketMode> modes_m;

        /// Prescribed compact spectral profile, evaluated before normalization.
        static double spectralWeight(double q) {
            const double u = (q - KX) / SPECTRUM_HALF_WIDTH;
            if (std::abs(u) >= 1)
                return 0;
            return std::exp(SPECTRUM_SHARPNESS * (1 - 1 / (1 - u * u)));
        }

        /// Dispersion and refraction of E18-modes, with the chosen Fresnel coefficients.
        PacketMode makeMode(double q) const {
            const double omega       = std::sqrt(KX * KX + q * q);
            const double qSecond     = std::sqrt(DIELECTRIC_EPSILON * omega * omega - KX * KX);
            const double coefficient = reflection(q);
            const double transmission =
                (1 + coefficient) / (tm_m ? std::sqrt(DIELECTRIC_EPSILON) : 1);
            return {q, omega, qSecond, spectralWeight(q), coefficient, transmission};
        }

        /// TE basis and the three TM bases of E19-polarization.
        std::array<double, 3> polarization(const PacketMode& mode, Branch branch) const {
            if (!tm_m)
                return {0, 1, 0};
            const double kz = branch == Branch::Transmitted ? mode.qSecond_m
                              : branch == Branch::Reflected ? -mode.q_m
                                                            : mode.q_m;
            const double waveNumber =
                mode.omega_m * (branch == Branch::Transmitted ? std::sqrt(DIELECTRIC_EPSILON) : 1);
            return {kz / waveNumber, 0, -KX / waveNumber};
        }

        /// Normal phase of E18-phases, with the common tangential phase removed.
        static double phase(const PacketMode& mode, Branch branch, double z, double time) {
            if (branch == Branch::Incident)
                return mode.q_m * (z - INITIAL_CENTER_Z) - mode.omega_m * time;
            if (branch == Branch::Reflected)
                return mode.q_m * (2 * INTERFACE_Z - z - INITIAL_CENTER_Z) - mode.omega_m * time;
            return mode.qSecond_m * (z - INTERFACE_Z) + mode.q_m * (INTERFACE_Z - INITIAL_CENTER_Z)
                   - mode.omega_m * time;
        }

        /** @brief Apply the continuum derivative to one real polarization vector.
         * @details With the phase of E18, the multipliers are \f$i\omega\f$
         * for E and \f$i\mathbf k\times\f$ for B, from
         * [chew2014generalizedgauge, Eqs. (5)--(6)].
         */
        std::complex<double> componentFactor(const PacketMode& mode, Branch branch,
                                             unsigned component, Sample quantity) const {
            const auto basis = polarization(mode, branch);
            if (quantity == Sample::Potential)
                return basis[component];
            if (quantity == Sample::Electric)
                return {0, mode.omega_m * basis[component]};
            const double kz = branch == Branch::Transmitted ? mode.qSecond_m
                              : branch == Branch::Reflected ? -mode.q_m
                                                            : mode.q_m;
            const std::array<double, 3> curl{-kz * basis[1], kz * basis[0] - KX * basis[2],
                                             KX * basis[1]};
            return {0, curl[component]};
        }

        /// One normalized incident, reflected or transmitted plane-wave contribution.
        std::complex<double> contribution(const PacketMode& mode, Branch branch, double z,
                                          double time, unsigned component, Sample quantity) const {
            const double coefficient = branch == Branch::Reflected     ? mode.reflection_m
                                       : branch == Branch::Transmitted ? mode.transmission_m
                                                                       : 1;
            return mode.amplitude_m * coefficient
                   * componentFactor(mode, branch, component, quantity)
                   * std::polar(1.0, phase(mode, branch, z, time));
        }

        /** @brief Sum the selected branches at one normal position and physical time.
         * @details E19 assigns A_z on the face to the vacuum branch for initialization.
         * The field formulas and all other components use the transmitted branch
         * on the face, as stated in E19-analytic-EB.
         */
        std::complex<double> amplitude(double z, double time, unsigned component,
                                       Sample quantity) const {
            // TE has A and E along y. TM has B along y.
            const bool polarizedAlongY = quantity == Sample::Magnetic ? tm_m : !tm_m;
            if ((component == 1) != polarizedAlongY)
                return 0;
            const bool vacuum =
                z < INTERFACE_Z
                || (z == INTERFACE_Z && component == 2 && quantity == Sample::Potential);
            std::complex<double> result = 0;
            for (const auto& mode : modes_m) {
                if (vacuum)
                    result += contribution(mode, Branch::Incident, z, time, component, quantity)
                              + contribution(mode, Branch::Reflected, z, time, component, quantity);
                else
                    result += contribution(mode, Branch::Transmitted, z, time, component, quantity);
            }
            return result;
        }

        /// Restore the common tangential phase and return the physical real-valued field.
        double sample(double x, double z, double time, unsigned component, Sample quantity) const {
            return std::real(std::polar(1.0, KX * x) * amplitude(z, time, component, quantity));
        }
    };
}  // namespace chdr::experiments::e18e19
