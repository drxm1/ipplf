/** @file test_analytic_reference.cpp
 * @brief Check the TE/TM continuum fields, face conditions and initial history tables.
 */
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "experiments/E18E19_oblique_incidence_on_planar_dielectric_interface/E18E19InitialSampler.h"
#include "experiments/E18E19_oblique_incidence_on_planar_dielectric_interface/E18E19Setup.h"

namespace {
    using Reference = chdr::experiments::e18e19::E18E19AnalyticContinuumReference;

    /// Keep failures active in optimized builds and retain the name of the violated relation.
    void expectClose(double actual, double expected, double tolerance, const std::string& name) {
        if (!std::isfinite(actual)
            || std::abs(actual - expected) > tolerance * std::max(1.0, std::abs(expected)))
            throw std::runtime_error(name + ": got " + std::to_string(actual) + ", expected "
                                     + std::to_string(expected));
    }

    /// Independent centered derivative of a potential sample in x, y, z or time.
    double derivative(const Reference& reference, double x, double z, double time,
                      unsigned component, unsigned direction) {
        constexpr double DifferenceStep = 1e-6;
        std::array<double, 4> lower{x, 0, z, time}, upper = lower;
        lower[direction] -= DifferenceStep;
        upper[direction] += DifferenceStep;
        const double a = reference.referenceFourPotentialA(lower[0], lower[1], lower[2], lower[3],
                                                           component + 1);
        const double b = reference.referenceFourPotentialA(upper[0], upper[1], upper[2], upper[3],
                                                           component + 1);
        return (b - a) / (2 * DifferenceStep);
    }

    /** @brief Check both continuum field definitions by numerical differentiation.
     * @details The identities are [chew2014generalizedgauge, Eqs. (5)--(6)],
     * with zero scalar potential and no dependence on y.
     */
    void checkDerivatives(const Reference& reference, double z, double time) {
        constexpr double X = 0.137;
        for (unsigned component = 0; component < 3; ++component)
            expectClose(reference.referenceE(X, 0, z, time, component),
                        -derivative(reference, X, z, time, component, 3), 3e-8,
                        "electric derivative");
        const std::array<double, 3> curl{
            -derivative(reference, X, z, time, 1, 2),
            derivative(reference, X, z, time, 0, 2) - derivative(reference, X, z, time, 2, 0),
            derivative(reference, X, z, time, 1, 0)};
        for (unsigned component = 0; component < 3; ++component)
            expectClose(reference.referenceB(X, 0, z, time, component), curl[component], 3e-8,
                        "magnetic curl");
    }

    /** @brief Check tangential A, E and B at the face.
     * @details These are the source-free interface conditions of
     * [chew2014generalizedgauge, Eqs. (23)--(25)].
     */
    void checkTangentialConditions(const Reference& reference, double time) {
        constexpr double Face = Reference::INTERFACE_Z, X = 0.137;
        const double left = std::nextafter(Face, -std::numeric_limits<double>::infinity());
        for (unsigned component : {0u, 1u}) {
            expectClose(reference.referenceFourPotentialA(X, 0, left, time, component + 1),
                        reference.referenceFourPotentialA(X, 0, Face, time, component + 1), 2e-12,
                        "tangential A");
            expectClose(reference.referenceE(X, 0, left, time, component),
                        reference.referenceE(X, 0, Face, time, component), 2e-12, "tangential E");
            expectClose(reference.referenceB(X, 0, left, time, component),
                        reference.referenceB(X, 0, Face, time, component), 2e-12, "tangential B");
        }
    }

    /** @brief Check the normal displacement and the specified one-sided face values.
     * @details The normal epsilon-A condition is [chew2014generalizedgauge,
     * Eq. (26)]. E19 assigns A_z on the face to vacuum during initialization,
     * while the field formula assigns E_z to dielectric.
     */
    void checkNormalConditions(const Reference& reference, double time) {
        constexpr double Face = Reference::INTERFACE_Z, X = 0.137;
        const double left  = std::nextafter(Face, -std::numeric_limits<double>::infinity());
        const double right = std::nextafter(Face, std::numeric_limits<double>::infinity());
        expectClose(reference.referenceE(X, 0, left, time, 2),
                    Reference::DIELECTRIC_EPSILON * reference.referenceE(X, 0, Face, time, 2),
                    2e-12, "normal D");
        expectClose(
            reference.referenceFourPotentialA(X, 0, Face, time, 3),
            Reference::DIELECTRIC_EPSILON * reference.referenceFourPotentialA(X, 0, right, time, 3),
            2e-12, "normal epsilon A");
    }

    /** @brief Check normalized amplitudes and lossless flux for every spectral term.
     * @details For equal permeability the transmitted fraction is
     * \f$(q'_m/q_m)S_m^2\f$, obtained from the normal plane-wave Poynting
     * flux [ellingson2020ElectromagneticsVol2, Eqs. (3.25), (5.187)] and the
     * potential amplitude in E19-analytic-EB. The same expression uses tau for TE.
     */
    void checkSpectrum(const Reference& reference) {
        expectClose(reference.modes().size(), 327, 0, "mode count");
        double sum = 0;
        for (const auto& mode : reference.modes()) {
            sum += mode.amplitude_m;
            expectClose(mode.reflection_m * mode.reflection_m
                            + mode.qSecond_m / mode.q_m * mode.transmission_m * mode.transmission_m,
                        1, 5e-15, "lossless mode flux");
        }
        expectClose(sum, 1, 1e-14, "normalized spectrum");
        expectClose(reference.measuredWavenumbers().size(), 12, 0, "probe count");
        if (reference.tm_m)
            expectClose(reference.reflection(reference.measuredWavenumbers().back()), 0, 1e-15,
                        "Brewster zero");
    }

    /// Sample both slots through the device functor at the actual component offsets.
    auto sampleInitialTable(const chdr::experiments::e18e19::E18E19InitialSampler& sampler,
                            chdr::LevelContract contract) {
        constexpr int Cells = 7;
        Kokkos::View<double***> values("initial sampler check", 2, 4, Cells);
        Kokkos::parallel_for(
            "sample both potential histories", 2 * 4 * Cells, KOKKOS_LAMBDA(int index) {
                const int slot = index / (4 * Cells), component = (index / Cells) % 4,
                          cell    = index % Cells;
                const auto offset = contract.positionOffset(
                    component == 0 ? chdr::Quantity::Phi : chdr::Quantity::VectorPotential,
                    component == 0 ? 0 : component - 1);
                const double x = 0.137 + offset[0] * sampler.spacingZ_m;
                const double z =
                    sampler.originZ_m
                    + (sampler.firstGlobalZ_m + cell + 0.5 + offset[2]) * sampler.spacingZ_m;
                values(slot, component, cell) =
                    sampler(x, 0, z, sampler.currentTime_m - slot * sampler.timeStep_m, component);
            });
        return Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, values);
    }

    /// Compare tabulated samples with direct sums across the interface at both history times.
    void checkInitialTable(const Reference& reference, chdr::SpatialDiscretization spatial) {
        const chdr::LevelContract contract{chdr::EBReconstructionConvention::R2, spatial};
        const auto sampler =
            chdr::experiments::e18e19::makeInitialSampler(reference, 7, 55, 0, 0.5, contract, 0.25);
        const auto values = sampleInitialTable(sampler, contract);
        for (int slot = 0; slot < 2; ++slot)
            for (unsigned component = 0; component < 4; ++component)
                for (int cell = 0; cell < 7; ++cell) {
                    const auto offset = contract.positionOffset(
                        component == 0 ? chdr::Quantity::Phi : chdr::Quantity::VectorPotential,
                        component == 0 ? 0 : component - 1);
                    const double x = 0.137 + offset[0] * sampler.spacingZ_m;
                    const double z =
                        sampler.originZ_m
                        + (sampler.firstGlobalZ_m + cell + 0.5 + offset[2]) * sampler.spacingZ_m;
                    expectClose(
                        values(slot, component, cell),
                        reference.referenceFourPotentialA(
                            x, 0, z, sampler.currentTime_m - slot * sampler.timeStep_m, component),
                        2e-14, "initial history table");
                }
    }

    /// Check the normal potential tie using the same mesh construction as the study.
    void checkFaceCoordinates() {
        const chdr::LevelContract contract{chdr::EBReconstructionConvention::R2,
                                           chdr::SpatialDiscretization::Staggered};
        const double offset = contract.positionOffset(chdr::Quantity::VectorPotential, 2)[2];
        for (int cells : {24, 48, 96, 192}) {
            const auto grid    = chdr::experiments::e18e19::makeGrid(cells);
            auto box           = chdr::solver::makeBoxFieldContainer(grid.cells_m, grid.lengths(),
                                                                     grid.splitAxis_m, grid.periodic_m);
            const int faceCell = 28 * cells - 1;
            const double face  = box.getMesh().getOrigin()[2]
                                + (faceCell + 0.5 + offset) * box.getMesh().getMeshSpacing()[2];
            expectClose(face, Reference::INTERFACE_Z, 0, "exact interface coordinate");
        }
    }

    /// Exercise incident, reflected and transmitted parts and each declared potential placement.
    void checkExperiment(bool tm) {
        const Reference reference(tm);
        checkSpectrum(reference);
        checkDerivatives(reference, 23.13, 0);
        checkDerivatives(reference, 26.27, 10);
        checkDerivatives(reference, 29.53, 10);
        for (double time : {5.0, 7.0, 9.0}) {
            checkTangentialConditions(reference, time);
            checkNormalConditions(reference, time);
        }
        checkInitialTable(reference, chdr::SpatialDiscretization::Staggered);
        checkInitialTable(reference, chdr::SpatialDiscretization::Collocated);
    }
}  // namespace

int main(int argc, char** argv) {
    ippl::initialize(argc, argv);
    int status = 0;
    try {
        checkFaceCoordinates();
        checkExperiment(false);
        checkExperiment(true);
        std::cout << "E18/E19 continuum reference and history tables passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        status = 1;
    }
    ippl::finalize();
    return status;
}
