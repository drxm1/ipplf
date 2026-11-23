/** @file test_analytic_reference.cpp
 * @brief Check field reconstruction, interface continuity and pulse energy fractions.
 */
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "experiments/E16E17_normal_incidence_on_planar_dielectric_interface/E16E17AnalyticContinuumReference.h"

namespace {
    using Reference = chdr::experiments::e16e17::E16E17AnalyticContinuumReferenceSampler;

    /// Fail in release builds too, and retain the name of the physical check.
    void expectClose(double actual, double expected, double tolerance, const std::string& name) {
        if (!std::isfinite(actual)
            || std::abs(actual - expected) > tolerance * std::max(1.0, std::abs(expected)))
            throw std::runtime_error(name + ": got " + std::to_string(actual) + ", expected "
                                     + std::to_string(expected));
    }

    /** @brief Independently differentiate the sampled potential at a point in either medium.
     * @details Chew [chew2014generalizedgauge, Eqs. (5)--(6)] requires
     * \f$E_x=-\partial_t A_x\f$ and \f$B_y=\partial_z A_x\f$ here.
     * Centered numerical differences test the analytic derivatives and the reflected sign.
     */
    void checkFieldDerivatives(const Reference& reference, double z, double time) {
        constexpr double DifferenceStep = 1e-6;
        const double electric =
            -(reference.Ax(z, time + DifferenceStep) - reference.Ax(z, time - DifferenceStep))
            / (2 * DifferenceStep);
        const double magnetic =
            (reference.Ax(z + DifferenceStep, time) - reference.Ax(z - DifferenceStep, time))
            / (2 * DifferenceStep);
        expectClose(reference.referenceE(0, 0, z, time, 0), electric, 2e-7, "electric derivative");
        expectClose(reference.referenceB(0, 0, z, time, 1), magnetic, 2e-7, "magnetic derivative");
    }

    /** @brief Check the tangential fields at the same physical interface from both sides.
     * @details Their continuity is the construction used in
     * [ellingson2020ElectromagneticsVol2, Eqs. (5.9), (5.13)]. With equal permeability,
     * continuity of H also gives continuity of B. The potential is continuous by
     * [chew2014generalizedgauge, Eq. (23)].
     */
    void checkInterface(const Reference& reference) {
        const double left =
            std::nextafter(reference.zf_m, -std::numeric_limits<double>::infinity());
        for (double phase : {-0.31, -0.07, 0.17, 0.43}) {
            const double time = (reference.zf_m - reference.z0_m - phase) / reference.v1_m;
            expectClose(reference.Ax(left, time), reference.Ax(reference.zf_m, time), 1e-12,
                        "continuous potential");
            expectClose(reference.Ex(left, time), reference.Ex(reference.zf_m, time), 1e-12,
                        "continuous electric field");
            expectClose(reference.By(left, time), reference.By(reference.zf_m, time), 1e-12,
                        "continuous magnetic field");
        }
    }

    /** @brief Integrate the continuum field energy with midpoint quadrature.
     * @details The density \f$(\epsilon E_x^2+B_y^2)/2\f$ follows from
     * [ellingson2020ElectromagneticsVol2, Eqs. (3.20), (3.22)] in E16's units.
     * This independent integral tests both the transmitted amplitude and its spatial scaling.
     */
    double pulseEnergy(const Reference& reference, double lower, double upper, double epsilon,
                       double time) {
        constexpr unsigned Cells = 20000;
        const double spacing     = (upper - lower) / Cells;
        double energy            = 0;
        for (unsigned cell = 0; cell < Cells; ++cell) {
            const double z        = lower + (cell + 0.5) * spacing;
            const double electric = reference.Ex(z, time), magnetic = reference.By(z, time);
            energy += 0.5 * spacing * (epsilon * electric * electric + magnetic * magnetic);
        }
        return energy;
    }

    /** @brief Check the integrated reflection and transmission fractions after separation.
     * @details The expected fractions are \f$\Gamma^2\f$ and \f$1-\Gamma^2\f$,
     * derived for pulses in E16.typ from
     * [ellingson2020ElectromagneticsVol2, Eqs. (5.22), (5.24)].
     */
    void checkEnergyFractions(const Reference& reference) {
        const double initial = pulseEnergy(reference, 0, reference.zf_m, reference.epsilon1_m, 0);
        const double time    = 3 / reference.v1_m;
        const double reflected =
            pulseEnergy(reference, 0, reference.zf_m, reference.epsilon1_m, time);
        const double transmitted =
            pulseEnergy(reference, reference.zf_m, 10, reference.epsilon2_m, time);
        const double fraction = reference.reflection_m * reference.reflection_m;
        expectClose(reflected / initial, fraction, 2e-12, "reflected energy fraction");
        expectClose(transmitted / initial, 1 - fraction, 2e-12, "transmitted energy fraction");
    }

    /// Test compact support exactly at both endpoints and derivative values in the interior.
    void checkProfile() {
        const chdr::experiments::e16e17::CompactCarrierProfile profile;
        constexpr double DifferenceStep = 1e-6;
        for (double s : {-1.1, -1.0, 1.0, 1.1}) {
            expectClose(profile.value(s), 0, 0, "compact potential support");
            expectClose(profile.derivative(s), 0, 0, "compact derivative support");
        }
        for (double s : {-0.7, -0.13, 0.0, 0.27, 0.8}) {
            const double derivative =
                (profile.value(s + DifferenceStep) - profile.value(s - DifferenceStep))
                / (2 * DifferenceStep);
            expectClose(profile.derivative(s), derivative, 2e-7, "profile derivative");
        }
    }

    /// Exercise each active pulse branch and both directions through the material contrast.
    void checkExperiment(double epsilon1, double epsilon2) {
        const Reference reference(epsilon1, epsilon2);
        checkFieldDerivatives(reference, reference.z0_m + 0.13, 0);
        const double time = 3 / reference.v1_m;
        checkFieldDerivatives(reference, 2 * reference.zf_m - reference.z0_m - 3 - 0.13, time);
        checkFieldDerivatives(reference, reference.zf_m + 1.67 / reference.indexRatio_m, time);
        checkInterface(reference);
        checkEnergyFractions(reference);
        expectClose(reference.referenceFourPotentialA(0, 0, 2.5, 0, 0), 0, 0, "zero phi");
        expectClose(reference.referenceE(0, 0, 2.5, 0, 1), 0, 0, "zero transverse E component");
        expectClose(reference.referenceB(0, 0, 2.5, 0, 0), 0, 0, "zero transverse B component");
    }

    /// Equal media recover a single traveling pulse throughout the domain.
    void checkEqualMedia() {
        const Reference reference(2.13, 2.13);
        const double time = 3 / reference.v1_m;
        for (double z : {2.3, 3.8, 4.0, 5.17, 5.63}) {
            const double phase = z - reference.z0_m - reference.v1_m * time;
            expectClose(reference.Ax(z, time), reference.profile_m.value(phase), 1e-13,
                        "equal-media potential");
            expectClose(reference.By(z, time), reference.profile_m.derivative(phase), 1e-13,
                        "equal-media magnetic field");
        }
    }
}  // namespace

int main() {
    try {
        checkProfile();
        checkExperiment(1, 2.13);
        checkExperiment(2.13, 1);
        checkEqualMedia();
        std::cout << "E16/E17 continuum reference checks passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
