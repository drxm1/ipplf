/** @file E16E17AnalyticContinuumReference.h
 * @brief Continuum pulse at a planar dielectric interface for E16 and E17.
 */
#pragma once

#include <Kokkos_Macros.hpp>
#include <Kokkos_MathematicalConstants.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
#include <catalyst_conduit.hpp>
#include <cmath>
#include <string>

namespace chdr::experiments::e16e17 {
    /** @brief Compact potential envelope with the carrier specified in E16.typ.
     * @details With \f$q=1-(s/w)^2\f$, the profile inside its support is
     * \f$a(s)=A_0\exp[\beta(1-1/q)]\cos(k_c s)\f$.
     * E16 defines \f$A_0=w=1\f$, \f$\beta=8\f$ and
     * \f$k_c=2\pi/\lambda_c\f$. The profile and its derivative vanish outside
     * \f$|s|<w\f$. This is a prescribed initial condition, with no periodic extension.
     */
    struct CompactCarrierProfile {
        static constexpr double Amplitude        = 1;
        static constexpr double SupportHalfWidth = 1;
        static constexpr double Sharpness        = 8;

        double carrierWavelength_m = 0.25;

        /** @brief Evaluate the prescribed potential profile at displacement \f$s\f$.
         * @details Test the support before dividing by \f$q\f$, which vanishes at its edges.
         */
        KOKKOS_INLINE_FUNCTION double value(double s) const {
            if (Kokkos::abs(s) >= SupportHalfWidth)
                return 0;
            return envelope(s) * Kokkos::cos(waveNumber() * s);
        }

        /** @brief Differentiate the profile of E16.typ by the product and chain rules.
         * @details The envelope derivative is
         * \f$p_\beta'=-2\beta s p_\beta/(w^2q^2)\f$. Thus
         * \f$a'=p_\beta'\cos(k_c s)-k_cp_\beta\sin(k_c s)\f$.
         */
        KOKKOS_INLINE_FUNCTION double derivative(double s) const {
            if (Kokkos::abs(s) >= SupportHalfWidth)
                return 0;
            const double q             = 1 - s * s / (SupportHalfWidth * SupportHalfWidth);
            const double envelopeValue = envelope(s);
            const double slope =
                -2 * Sharpness * s * envelopeValue / (SupportHalfWidth * SupportHalfWidth * q * q);
            const double k = waveNumber();
            return slope * Kokkos::cos(k * s) - k * envelopeValue * Kokkos::sin(k * s);
        }

        /// Record the prescribed profile parameters for reproduction of the run.
        void record(conduit_cpp::Node& node) const {
            node["potential_amplitude"] = Amplitude;
            node["support_half_width"]  = SupportHalfWidth;
            node["profile_sharpness"]   = Sharpness;
            node["carrier_wavelength"]  = carrierWavelength_m;
        }

    private:
        /// Carrier wavenumber of E16.typ, in the incident medium.
        KOKKOS_INLINE_FUNCTION double waveNumber() const {
            return 2 * Kokkos::numbers::pi / carrierWavelength_m;
        }

        /// Evaluate the compact envelope only inside its support.
        KOKKOS_INLINE_FUNCTION double envelope(double s) const {
            const double q = 1 - s * s / (SupportHalfWidth * SupportHalfWidth);
            return Amplitude * Kokkos::exp(Sharpness * (1 - 1 / q));
        }
    };

    /** @brief Sample the continuum potentials and fields on both sides of one interface.
     * @details E16.typ gives the incident, reflected and transmitted phases. The
     * nonmagnetic, nondispersive media have refractive indices
     * \f$n_j=\sqrt{\epsilon_j}\f$ and the vacuum speed is one. The coefficients
     * \f$\Gamma=(n_1-n_2)/(n_1+n_2)\f$ and \f$\tau=1+\Gamma\f$ follow from
     * [ellingson2020ElectromagneticsVol2, Eqs. (5.16), (5.17), (5.27)]. Applying
     * \f$\mathbf E=-\partial_t\mathbf A-\nabla\phi\f$ and
     * \f$\mathbf B=\nabla\times\mathbf A\f$ to these potentials gives the fields
     * below [chew2014generalizedgauge, Eqs. (5)--(6)]. Only \f$A_x,E_x,B_y\f$
     * are nonzero. The reference is valid in the periodic box while the pulses stay
     * away from its second material interface, as specified in E16.typ and E17.typ.
     */
    struct E16E17AnalyticContinuumReferenceSampler {
        double epsilon1_m;
        double epsilon2_m;
        double z0_m;
        double zf_m;
        double v1_m;
        double indexRatio_m;
        double reflection_m;
        double transmission_m;
        CompactCarrierProfile profile_m;

        /** @brief Set the media, geometry and carrier wavelength of the reference.
         * @pre Both permittivities and the carrier wavelength are positive and finite.
         * @details Square roots are evaluated on the host. The resulting sampler is
         * copied by value into Kokkos kernels and has no owned memory or device transfers.
         */
        E16E17AnalyticContinuumReferenceSampler(double epsilon1, double epsilon2, double z0 = 2.5,
                                                double zf = 4, double carrierWavelength = 0.25)
            : epsilon1_m(epsilon1)
            , epsilon2_m(epsilon2)
            , z0_m(z0)
            , zf_m(zf)
            , v1_m(1 / std::sqrt(epsilon1))
            , indexRatio_m(std::sqrt(epsilon2) / std::sqrt(epsilon1))
            , reflection_m((1 - indexRatio_m) / (1 + indexRatio_m))
            , transmission_m(1 + reflection_m)
            , profile_m{carrierWavelength} {}

        /// Sample the four-potential in storage order phi, Ax, Ay, Az.
        KOKKOS_INLINE_FUNCTION double operator()(double, double, double z, double time,
                                                 unsigned component) const {
            return component == 1 ? Ax(z, time) : 0;
        }

        /// Use the same continuum sampler for initialization and potential observations.
        KOKKOS_INLINE_FUNCTION double referenceFourPotentialA(double x, double y, double z,
                                                              double time,
                                                              unsigned component) const {
            return (*this)(x, y, z, time, component);
        }

        /// Sample E at the caller's component position and physical time.
        KOKKOS_INLINE_FUNCTION double referenceE(double, double, double z, double time,
                                                 unsigned component) const {
            return component == 0 ? Ex(z, time) : 0;
        }

        /// Sample B at the caller's component position and physical time.
        KOKKOS_INLINE_FUNCTION double referenceB(double, double, double z, double time,
                                                 unsigned component) const {
            return component == 1 ? By(z, time) : 0;
        }

        /// Potential of E16.typ, with the transmitted branch at the interface itself.
        KOKKOS_INLINE_FUNCTION double Ax(double z, double time) const {
            if (z < zf_m)
                return profile_m.value(incidentPhase(z, time))
                       + reflection_m * profile_m.value(reflectedPhase(z, time));
            return transmission_m * profile_m.value(transmittedPhase(z, time));
        }

        /// Time differentiation of Ax, using Chew (6) with zero scalar potential.
        KOKKOS_INLINE_FUNCTION double Ex(double z, double time) const {
            if (z < zf_m)
                return v1_m
                       * (profile_m.derivative(incidentPhase(z, time))
                          + reflection_m * profile_m.derivative(reflectedPhase(z, time)));
            return v1_m * transmission_m * profile_m.derivative(transmittedPhase(z, time));
        }

        /// Spatial differentiation of Ax, using Chew (5) and the phases of E16.typ.
        KOKKOS_INLINE_FUNCTION double By(double z, double time) const {
            if (z < zf_m)
                return profile_m.derivative(incidentPhase(z, time))
                       - reflection_m * profile_m.derivative(reflectedPhase(z, time));
            return indexRatio_m * transmission_m * profile_m.derivative(transmittedPhase(z, time));
        }

        /** @brief Record the continuum coefficients and profile beside the solver input.
         * @param node Borrowed `analytic_reference` child returned by node indexing.
         * @param experiment Experiment name, E16, E17 or E16-S.
         */
        void record(conduit_cpp::Node node, const std::string& experiment) const {
            node["name"]                         = experiment + " normal-incidence continuum pulse";
            node["epsilon_side1"]                = epsilon1_m;
            node["epsilon_side2"]                = epsilon2_m;
            node["initial_center_z"]             = z0_m;
            node["interface_z"]                  = zf_m;
            node["incident_speed"]               = v1_m;
            node["refractive_index_ratio"]       = indexRatio_m;
            node["reflection_coefficient"]       = reflection_m;
            node["transmission_coefficient"]     = transmission_m;
            node["reflection_energy_fraction"]   = reflection_m * reflection_m;
            node["transmission_energy_fraction"] = 1 - reflection_m * reflection_m;
            profile_m.record(node);
        }

    private:
        /// Incident phase of E16.typ, increasing along z.
        KOKKOS_INLINE_FUNCTION double incidentPhase(double z, double time) const {
            return z - z0_m - v1_m * time;
        }

        /// Reflected phase of E16.typ, decreasing along z.
        KOKKOS_INLINE_FUNCTION double reflectedPhase(double z, double time) const {
            return 2 * zf_m - z - z0_m - v1_m * time;
        }

        /// Transmitted phase of E16.typ, continuous with the incident phase at the face.
        KOKKOS_INLINE_FUNCTION double transmittedPhase(double z, double time) const {
            return zf_m - z0_m - v1_m * time + indexRatio_m * (z - zf_m);
        }
    };
}  // namespace chdr::experiments::e16e17
