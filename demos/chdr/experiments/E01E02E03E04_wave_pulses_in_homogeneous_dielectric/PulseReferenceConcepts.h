/** @file PulseReferenceConcepts.h
 * @brief Named requirements on the analytic continuum references of the pulse run paths.
 */
#pragma once

#include <catalyst_conduit.hpp>
#include <concepts>
#include <string>

#include "core/solver/BoxSetup.h"

namespace chdr::experiments {
    /** @brief Continuum reference of a source-free pulse in the homogeneous dielectric.
     * @details The four-potential call loads the solver histories (solver::FieldComponentSampler,
     * storage order \f$\phi, A_x, A_y, A_z\f$). `referenceFourPotentialA`, `referenceE` and
     * `referenceB` evaluate one component at a native position and time for the measurements,
     * and `v_m` is the wave speed \f$1/\sqrt{\epsilon_r}\f$ that converts travel distances
     * into target times. Every evaluation must be callable inside a Kokkos kernel.
     */
    template <class T>
    concept ContinuumPulseReference =
        solver::FieldComponentSampler<T>
        && requires(const T& reference, double coordinate, unsigned component) {
               {
                   reference.referenceFourPotentialA(coordinate, coordinate, coordinate, coordinate,
                                                     component)
               } -> std::same_as<double>;
               {
                   reference.referenceE(coordinate, coordinate, coordinate, coordinate, component)
               } -> std::same_as<double>;
               {
                   reference.referenceB(coordinate, coordinate, coordinate, coordinate, component)
               } -> std::same_as<double>;
               { reference.v_m } -> std::convertible_to<double>;
           };

    /** @brief Plane pulse of the E01/E02 run path, for example the E01 sampler with any profile.
     * @details The run path constructs it from the domain length \f$L_z\f$, the initial center
     * \f$z_0\f$, the wave speed and the periodicity in \f$z\f$, and the reference records
     * itself in the run record under the experiment name.
     */
    template <class T>
    concept PlanePulseReference =
        ContinuumPulseReference<T>
        && requires(const T& reference, conduit_cpp::Node node, const std::string& experiment) {
               T{double{} /* Lz */, double{} /* z0 */, double{} /* v */, bool{} /* periodic */};
               reference.record(node, experiment);
           };

    /** @brief Dipole pulse of the E03/E04 run path, for example the E03 sampler with any profile.
     * @details The caller constructs it. The reference records its center, speed and profile.
     */
    template <class T>
    concept DipolePulseReference =
        ContinuumPulseReference<T>
        && requires(const T& reference, conduit_cpp::Node node) { reference.record(node); };
}  // namespace chdr::experiments
