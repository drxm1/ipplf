/** @file TiltedDispersion.h
 * @brief Directional Standard-scheme roots for a specified physical interface normal.
 */
#pragma once
#include <numeric>

#include "apps/planar-interface/PlaneSetup.h"

namespace chdr::apps::planarInterface {
    /// Requested tilt and primitive tangent translation, separate from incidence angle.
    struct TiltedRequest {
        PlaneRequest plane_m;
        double tiltDegrees_m = 0, fitResidualLimit_m = 1e-3;
        int primitiveX_m = 0, primitiveZ_m = 0, probeHalfWidth_m = 2;
        int periodMultiplier_m = 0;
        std::optional<int> tangentMode_m;
    };

    /// One signed normal root and its signed normal group velocity.
    struct TiltedRoot {
        double normalK_m, groupVelocity_m;
    };

    /// Resolved geometry; plane_m supplies only common material, mesh, band and timing data.
    struct TiltedRow {
        PlaneRow plane_m;
        std::array<double, 3> requestedShape_m{};
        double tiltDegrees_m = 0, cosine_m = 1, sine_m = 0, kappa_m = 0;
        double facePosition_m = 0;
        int primitiveX_m = 0, primitiveZ_m = 0, multiplier_m = 0, tangentMode_m = 0;
        int periodRows_m = 0, incidentRows_m = 0, probeCells_m = 12, probeHalfWidth_m = 2;
        double rowSpacing_m = 0, normalPeriod_m = 0, tangentPeriod_m = 0, primitivePeriod_m = 0;
        double dielectricWidth_m = 0, incidentWidth_m = 0, sourceDistance_m = 0;
        double fitResidualLimit_m = 1e-3, floquetMargin_m = 0, retainedMargin_m = 0;
        double discardedFieldWeightFraction_m = 0;
        std::vector<std::array<double, 5>> modes_m;
    };

    /// Cartesian wavevector from normal/tangent coordinates [oskooi2013sources, section 4.5].
    inline std::array<double, 2> tiltedWavevector(const TiltedRow& r, double normalK) {
        return {r.cosine_m * normalK - r.sine_m * r.kappa_m,
                r.sine_m * normalK + r.cosine_m * r.kappa_m};
    }

    /// Spatial symbol squared [fallahi2020mithra20fullwavesimulation, Eq. (3.17)].
    inline double tiltedSpatialSquare(const TiltedRow& r, double normalK) {
        const auto k   = tiltedWavevector(r, normalK);
        const double x = planeSpatialSymbol(k[0], r.plane_m.h_m[0]);
        const double z = planeSpatialSymbol(k[1], r.plane_m.h_m[2]);
        return x * x + z * z;
    }

    /// Half the normal derivative of the squared spatial symbol [Trefethen, section 1].
    inline double tiltedSpatialDerivative(const TiltedRow& r, double normalK) {
        const auto k = tiltedWavevector(r, normalK);
        const auto h = r.plane_m.h_m;
        return r.cosine_m * std::sin(k[0] * h[0]) / h[0] + r.sine_m * std::sin(k[1] * h[2]) / h[2];
    }

    /// Frequency of one incident carrier [fallahi2020mithra20fullwavesimulation, Eq. (3.17)].
    inline double tiltedFrequency(const TiltedRow& r, double normalK) {
        const double argument = r.plane_m.dt_m * std::sqrt(tiltedSpatialSquare(r, normalK)) / 2;
        if (!(argument >= 0 && argument < 1))
            throw std::invalid_argument("Tilted carrier leaves the principal time branch");
        return 2 / r.plane_m.dt_m * std::asin(argument);
    }

    /** @brief Normal interval with both Cartesian phases in the strictly convex central region.
     * Deduction: the second derivative of each squared sine symbol is positive for
     * \f$|k_xh_x|,|k_zh_z|<\pi/2\f$ [Fallahi Eq. (3.17)].
     */
    inline std::array<double, 2> tiltedConvexInterval(const TiltedRow& r) {
        const double x = Kokkos::numbers::pi / (2 * r.plane_m.h_m[0]);
        const double z = Kokkos::numbers::pi / (2 * r.plane_m.h_m[2]);
        const double a = (-z - r.cosine_m * r.kappa_m) / r.sine_m;
        const double b = (z - r.cosine_m * r.kappa_m) / r.sine_m;
        return {std::max((-x + r.sine_m * r.kappa_m) / r.cosine_m, std::min(a, b)),
                std::min((x + r.sine_m * r.kappa_m) / r.cosine_m, std::max(a, b))};
    }

    /// Locate the unique grazing minimum in the strictly convex interval.
    inline double tiltedGrazingRoot(const TiltedRow& r) {
        auto interval = tiltedConvexInterval(r);
        if (!(interval[0] < interval[1]) || tiltedSpatialDerivative(r, interval[0]) >= 0
            || tiltedSpatialDerivative(r, interval[1]) <= 0)
            throw std::invalid_argument("Tilted tangent mode has no qualified central minimum");
        for (int iteration = 0; iteration < 64; ++iteration) {
            const double mid = std::midpoint(interval[0], interval[1]);
            interval[tiltedSpatialDerivative(r, mid) < 0 ? 0 : 1] = mid;
        }
        return std::midpoint(interval[0], interval[1]);
    }

    /// Bisect one monotone side of the convex dispersion surface.
    inline double tiltedNormalRoot(const TiltedRow& r, double target, double grazing,
                                   bool forward) {
        const auto interval = tiltedConvexInterval(r);
        double lo = forward ? grazing : interval[0], hi = forward ? interval[1] : grazing;
        if (tiltedSpatialSquare(r, forward ? hi : lo) <= target)
            throw std::invalid_argument("Tilted root leaves the qualified central branch");
        for (int iteration = 0; iteration < 64; ++iteration) {
            const double mid             = std::midpoint(lo, hi);
            const bool below             = tiltedSpatialSquare(r, mid) < target;
            (below == forward ? lo : hi) = mid;
        }
        return std::midpoint(lo, hi);
    }

    /** @brief Two directional roots; no reflection symmetry is assumed.
     * @details The sufficient bound \f$\sqrt\varepsilon W\max h<\sqrt2\f$ puts
     * every possible Cartesian principal root inside the convex region. Thus this
     * pair exhausts that branch. Velocity differentiates Fallahi Eq. (3.17), following
     * [trefethen1982GroupVelocity, section 1, Eqs. (1.4)-(1.5)].
     */
    inline std::array<TiltedRoot, 2> tiltedRoots(const TiltedRow& r, double omega, double epsilon) {
        if (omega == 0 && r.kappa_m == 0)
            return {{{0, -1 / std::sqrt(epsilon)}, {0, 1 / std::sqrt(epsilon)}}};
        const double dt = r.plane_m.dt_m, W = 2 * std::sin(omega * dt / 2) / dt;
        const double target = epsilon * W * W, minimum = tiltedGrazingRoot(r);
        if (!(omega > 0 && omega * dt < Kokkos::numbers::pi)
            || std::sqrt(target) * std::max(r.plane_m.h_m[0], r.plane_m.h_m[2]) >= std::sqrt(2.0)
            || target <= tiltedSpatialSquare(r, minimum))
            throw std::invalid_argument("Tilted band lacks two qualified propagating roots");
        std::array<TiltedRoot, 2> result;
        for (int direction = 0; direction < 2; ++direction) {
            const double k    = tiltedNormalRoot(r, target, minimum, direction != 0);
            result[direction] = {
                k, tiltedSpatialDerivative(r, k) / (epsilon * W * std::cos(omega * dt / 2))};
        }
        return result;
    }

    /** @brief Sufficient exclusion margin for all nonzero material-period Floquet orders.
     * @details Deduction: bound each principal Cartesian phase by its sine symbol, then
     * apply Cauchy-Schwarz to the tangent projection. Order spacing follows the primitive
     * material translation [oskooi2013sources, section 4.5, Eq. (4.30)].
     */
    inline double tiltedOrderMargin(const TiltedRow& r, double omega) {
        const double K = std::sqrt(r.plane_m.epsilon_m) * 2 * std::sin(omega * r.plane_m.dt_m / 2)
                         / r.plane_m.dt_m;
        const auto ratio = [](double x) {
            x = std::min(1.0, x);
            return x == 0 ? 1 : std::asin(x) / x;
        };
        const double x = r.sine_m * ratio(K * r.plane_m.h_m[0] / 2);
        const double z = r.cosine_m * ratio(K * r.plane_m.h_m[2] / 2);
        return 2 * Kokkos::numbers::pi / r.primitivePeriod_m - std::abs(r.kappa_m)
               - K * std::sqrt(x * x + z * z);
    }
}  // namespace chdr::apps::planarInterface
