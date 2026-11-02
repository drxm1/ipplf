/** @file UnitSystem.h
 * @brief Convert normalized simulation values at the common input/output boundary.
 */
#pragma once

#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <string>

namespace chdr::input {
    /** @brief Multipliers from solver values to either normalized output or SI output.
     * @details Choose length and charge-density scales \f$L_0\f$ and \f$\rho_0\f$.
     * Substitution into SI Gauss's law fixes \f$E_0=\rho_0 L_0/\varepsilon_0\f$.
     * Setting \f$t_0=L_0/c\f$ in Faraday's law fixes \f$B_0=E_0/c\f$;
     * the potential definitions give \f$\phi_0=E_0L_0\f$ and
     * \f$A_0=B_0L_0\f$. These are deductions from
     * [fitzpatrick2014ClassicalElectromagnetism, Eqs. (1.1), (1.3), (1.11)-(1.16)].
     * For the linear, nondispersive nonmagnetic model, physical field energy uses
     * \f$U_0=\varepsilon_0 E_0^2L_0^3\f$, by substituting the field and volume
     * scales into [fitzpatrick2014ClassicalElectromagnetism, Eq. (1.109), p. 20;
     * dielectric electric-energy weighting: Eq. (4.70)].
     * This multiplier does not identify a potential invariant as physical energy.
     * The SI constants below are from [mohr2025Codata2022, Table XXXIII, p. 025002-45].
     * Normalized mode returns unit multipliers and assigns no physical scale to a preset.
     */
    struct UnitSystem {
        static constexpr double SPEED_OF_LIGHT_SI      = 299792458.0;
        static constexpr double VACUUM_PERMITTIVITY_SI = 8.8541878188e-12;

        std::string system_m = "normalized";
        double length_m      = 1.0;  ///< Metres per normalized length in SI mode.
        double chargeDensityScale_m =
            1.0;  ///< Coulombs per cubic metre per normalized charge density.

        /// Whether numerical input and output quantities are expressed in SI.
        bool isSI() const { return system_m == "si"; }
        /// Multiply code length by this factor to obtain metres in SI mode.
        double lengthFactor() const { return isSI() ? length_m : 1.0; }
        /// Multiply code time by this factor to obtain seconds in SI mode.
        double timeFactor() const { return isSI() ? length_m / SPEED_OF_LIGHT_SI : 1.0; }
        /// Multiply code velocity by this factor to obtain metres per second in SI mode.
        double velocityFactor() const { return isSI() ? SPEED_OF_LIGHT_SI : 1.0; }
        /// Multiply code frequency by this factor to obtain inverse seconds in SI mode.
        double frequencyFactor() const { return 1.0 / timeFactor(); }

        /// Multiply code electric field by this factor to obtain volts per metre in SI mode.
        double electricFactor() const {
            return isSI() ? static_cast<double>(static_cast<long double>(chargeDensityScale_m)
                                                * length_m / VACUUM_PERMITTIVITY_SI)
                          : 1.0;
        }
        /// Multiply code magnetic field by this factor to obtain teslas in SI mode.
        double magneticFactor() const {
            return isSI() ? electricFactor() / SPEED_OF_LIGHT_SI : 1.0;
        }
        /// Multiply code scalar potential by this factor to obtain volts in SI mode.
        double scalarPotentialFactor() const {
            return isSI()
                       ? static_cast<double>(static_cast<long double>(electricFactor()) * length_m)
                       : 1.0;
        }
        /// Multiply code vector potential by this factor to obtain tesla-metres in SI mode.
        double vectorPotentialFactor() const {
            return isSI()
                       ? static_cast<double>(static_cast<long double>(magneticFactor()) * length_m)
                       : 1.0;
        }
        /// Multiply a physical normalized field energy by this factor to obtain joules.
        double energyFactor() const {
            const long double electric = electricFactor(), length = length_m;
            return isSI() ? static_cast<double>(VACUUM_PERMITTIVITY_SI * electric * electric
                                                * length * length * length)
                          : 1.0;
        }

        /** @brief Reject invalid or unrepresentable conversion scales before allocating fields.
         * @details Long-double intermediates reduce unnecessary overflow; every resulting
         * double multiplier must still be finite and strictly positive for output conversion.
         */
        void validate() const {
            if (system_m != "normalized" && system_m != "si")
                throw std::invalid_argument("Unit system must be normalized or si");
            if (!isSI())
                return;
            for (double factor : {length_m, chargeDensityScale_m, timeFactor(), frequencyFactor(),
                                  electricFactor(), magneticFactor(), scalarPotentialFactor(),
                                  vectorPotentialFactor(), energyFactor()})
                if (!std::isfinite(factor) || factor <= 0)
                    throw std::invalid_argument(
                        "SI normalization scales must yield finite positive factors");
        }
    };
}  // namespace chdr::input
