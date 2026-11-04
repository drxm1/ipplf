/** @file TiltedPulse.h
 * @brief Periodic tilted material samples and exact lattice-mode potential histories.
 */
#pragma once
#include "apps/planar-interface/PlanePulse.h"
#include "apps/planar-interface/TiltedInput.h"

namespace chdr::apps::planarInterface {
    /** @brief Periodic centre-sampled slab with an exact integer normal-row label.
     * @details Twice the centre row is \f$q(2i+1)+\operatorname{sign}(s)p(2k+1)\f$. Cartesian box
     * translations change it by an integer period, including exact interface ties.
     * This samples the closed dielectric interval \f$[\delta-D,\delta]\f$;
     * it introduces no sub-cell average. The signed offset is applied after integer wrapping.
     */
    struct TiltedMaterialMap {
        int x_m, z_m, sign_m, periodRows_m, incidentRows_m;
        double epsilon_m, rowSpacing_m, offset_m;
        explicit TiltedMaterialMap(const TiltedRow& r)
            : x_m(r.primitiveX_m)
            , z_m(r.primitiveZ_m)
            , sign_m(r.sine_m > 0 ? 1 : -1)
            , periodRows_m(r.periodRows_m)
            , incidentRows_m(r.incidentRows_m)
            , epsilon_m(r.plane_m.epsilon_m)
            , rowSpacing_m(r.rowSpacing_m)
            , offset_m(r.facePosition_m) {}

        KOKKOS_INLINE_FUNCTION bool contains(int i, int k) const {
            const std::int64_t row = std::int64_t(z_m) * (2 * std::int64_t(i) + 1)
                                     + std::int64_t(sign_m * x_m) * (2 * std::int64_t(k) + 1);
            const std::int64_t period = 2 * std::int64_t(periodRows_m);
            const auto wrapped        = (row % period + period) % period;
            if (offset_m == 0)
                return wrapped == 0 || wrapped >= 2 * std::int64_t(incidentRows_m);
            const double length = periodRows_m * rowSpacing_m;
            double position     = 0.5 * wrapped * rowSpacing_m - offset_m;
            position -= Kokkos::floor(position / length) * length;
            return position == 0 || position >= incidentRows_m * rowSpacing_m;
        }

        KOKKOS_INLINE_FUNCTION double operator()(material::CellIndex index, material::Position,
                                                 bool& marked) const {
            marked = contains(index[0], index[2]);
            return marked ? epsilon_m : 1;
        }
    };

    /// Cartesian box used by the common manager and field writer.
    inline simulation::RunBox tiltedRunBox(const TiltedRow& r) {
        return {{r.plane_m.cellsX_m, 2, r.plane_m.cellsZ_m}, r.plane_m.h_m};
    }

    /// The anchor lies on \f$\mathbf n\cdot\mathbf x=\delta\f$ with its chosen tangent origin.
    inline std::array<double, 3> tiltedAnchor(const TiltedRow& r) {
        const auto length = tiltedRunBox(r).lengths();
        return {(r.sine_m > 0 ? -0.5 : 0.5) * length[0] + r.cosine_m * r.facePosition_m, 0,
                0.5 * length[2] + r.sine_m * r.facePosition_m};
    }

    /// Record the integer Cartesian winding without narrowing an unrepresentably distant anchor.
    inline int tiltedLift(const TiltedRow& r, unsigned axis) {
        const double winding = std::floor(tiltedAnchor(r)[axis] / tiltedRunBox(r).lengths()[axis]);
        if (!std::isfinite(winding) || std::abs(winding) > std::numeric_limits<int>::max() / 4)
            throw std::invalid_argument(
                "Tilted requested face needs an unsupported Cartesian winding");
        return static_cast<int>(winding);
    }

    /// Normal/tangent coordinates use the lifted actual component position relative to the anchor.
    inline std::array<double, 2> tiltedCoordinates(const TiltedRow& r,
                                                   std::array<double, 3> point) {
        const auto anchor = tiltedAnchor(r);
        point[0] += tiltedLift(r, 0) * r.plane_m.cellsX_m * r.plane_m.h_m[0];
        point[2] += tiltedLift(r, 2) * r.plane_m.cellsZ_m * r.plane_m.h_m[2];
        const double x = point[0] - anchor[0], z = point[2] - anchor[2];
        return {r.cosine_m * x + r.sine_m * z, -r.sine_m * x + r.cosine_m * z};
    }

    /// Select a small Cartesian box around a requested normal probe distance.
    inline diagnostics::ProbeBox tiltedProbeBox(const TiltedRow& r, bool incident) {
        const auto length   = tiltedRunBox(r).lengths();
        const auto anchor   = tiltedAnchor(r);
        const double d      = (incident ? 1 : -1) * r.probeCells_m * r.plane_m.h_m[0];
        const double offset = r.plane_m.polarization_m == Polarization::TM ? 1 : 0.5;
        const int x         = static_cast<int>(std::lround(
            (anchor[0] - tiltedLift(r, 0) * length[0] + r.cosine_m * d) / r.plane_m.h_m[0]
            - offset));
        const int z         = static_cast<int>(std::lround(
            (anchor[2] - tiltedLift(r, 2) * length[2] + r.sine_m * d) / r.plane_m.h_m[2] - offset));
        const int w         = r.probeHalfWidth_m;
        return {{x - w, 0, z - w}, {x + w, 0, z + w}};
    }

    /// Require homogeneous stencils around every recorded component, before allocating fields.
    inline void validateTiltedProbe(const TiltedRow& r, bool incident) {
        const auto box = tiltedProbeBox(r, incident);
        const TiltedMaterialMap map(r);
        if (box.first_m[0] < 1 || box.first_m[2] < 1 || box.last_m[0] + 2 >= r.plane_m.cellsX_m
            || box.last_m[2] + 2 >= r.plane_m.cellsZ_m)
            throw std::invalid_argument("Tilted probe patch crosses a Cartesian box edge");
        for (int i = box.first_m[0] - 1; i <= box.last_m[0] + 1; ++i)
            for (int k = box.first_m[2] - 1; k <= box.last_m[2] + 1; ++k)
                if (map.contains(i, k) == incident)
                    throw std::invalid_argument(
                        "Tilted probe patch intersects a sampled material interface");
    }

    /** @brief Device mode sum with a material-weighted discrete-curl TM start.
     * @details TE uses \f$A_y=\psi\f$. TM uses \f$A_x=-D_z\psi/\varepsilon_x\f$
     * and \f$A_z=D_x\psi/\varepsilon_z\f$,
     * with harmonic face coefficients. Its weighted divergence cancels by commuting
     * differences (deduction from [ryu2016potentialfdtd, Appendix A]; the generalized
     * gauge is [chew2014generalizedgauge, Eqs. (14), (15), (19)]).
     */
    struct TiltedPulse {
        Kokkos::View<double* [4]> modes_m;
        TiltedMaterialMap material_m;
        double hx_m, hz_m, sourceX_m, sourceZ_m;
        bool tm_m;
        explicit TiltedPulse(const TiltedRow& r)
            : modes_m("Tilted pulse modes", r.modes_m.size())
            , material_m(r)
            , hx_m(r.plane_m.h_m[0])
            , hz_m(r.plane_m.h_m[2])
            , sourceX_m(tiltedAnchor(r)[0] + r.cosine_m * r.sourceDistance_m)
            , sourceZ_m(tiltedAnchor(r)[2] + r.sine_m * r.sourceDistance_m)
            , tm_m(r.plane_m.polarization_m == Polarization::TM) {
            auto host = Kokkos::create_mirror_view(modes_m);
            for (size_t m = 0; m < r.modes_m.size(); ++m)
                for (int c = 0; c < 4; ++c)
                    host(m, c) = r.modes_m[m][c];
            Kokkos::deep_copy(modes_m, host);
        }

        KOKKOS_INLINE_FUNCTION double psi(double x, double z, double t) const {
            double value = 0;
            for (size_t m = 0; m < modes_m.extent(0); ++m)
                value += modes_m(m, 3)
                         * Kokkos::cos(modes_m(m, 0) * (x - sourceX_m)
                                       + modes_m(m, 1) * (z - sourceZ_m) - modes_m(m, 2) * t);
            return value;
        }

        /// Centre lookup agrees with MaterialOwner despite rounded component-coordinate arithmetic.
        KOKKOS_INLINE_FUNCTION double epsilon(double x, double z) const {
            return material_m.contains(static_cast<int>(Kokkos::floor(x / hx_m)),
                                       static_cast<int>(Kokkos::floor(z / hz_m)))
                       ? material_m.epsilon_m
                       : 1;
        }

        KOKKOS_INLINE_FUNCTION double operator()(double x, double, double z, double t,
                                                 unsigned c) const {
            if (!tm_m)
                return c == 2 ? psi(x, z, t) : 0;
            if (c != 1 && c != 3)
                return 0;
            const double dx = c == 1 ? hx_m / 2 : 0, dz = c == 3 ? hz_m / 2 : 0;
            const double inverse =
                0.5 * (1 / epsilon(x - dx, z - dz) + 1 / epsilon(x + dx, z + dz));
            const double curl = c == 1 ? (psi(x, z - hz_m / 2, t) - psi(x, z + hz_m / 2, t)) / hz_m
                                       : (psi(x + hx_m / 2, z, t) - psi(x - hx_m / 2, z, t)) / hx_m;
            return inverse * curl;
        }
    };

    /// Stable sinc for the finite geometric time-transform sum.
    inline double tiltedSinc(double x) {
        return std::abs(x) < 1e-5 ? 1 - x * x / 6 + x * x * x * x / 120 : std::sin(x) / x;
    }

    /** @brief Exact dt-weighted transform of a complex exponential at half-step times.
     * Deduction: sum the finite geometric series for \f$t=(k-1/2)\Delta t\f$ in the positive-kernel
     * transform [oskooi2013sources, section 4.3, Eq. (4.6)]. No time trace or repeated
     * device-to-host field copy is required for analytic incident normalization.
     */
    inline std::complex<double> tiltedTimeSum(double difference, double dt, int samples) {
        const double a = difference * dt / 2;
        return dt * samples * tiltedSinc(samples * a) / tiltedSinc(a)
               * std::exp(std::complex<double>(0, a * (samples - 2)));
    }

    /// Analytic incident channel at one true component position and its actual finite time gate.
    inline std::complex<double> tiltedIncidentTransform(const TiltedRow& r,
                                                        std::array<double, 3> point, double omega,
                                                        int samples) {
        const auto anchor           = tiltedAnchor(r);
        const double x              = point[0] - anchor[0] - r.cosine_m * r.sourceDistance_m;
        const double z              = point[2] - anchor[2] - r.sine_m * r.sourceDistance_m;
        const double dt             = r.plane_m.dt_m;
        std::complex<double> result = 0;
        for (const auto& mode : r.modes_m) {
            const double W      = 2 * std::sin(mode[2] * dt / 2) / dt;
            const auto phase    = std::exp(std::complex<double>(0, mode[0] * x + mode[1] * z));
            const auto positive = phase * tiltedTimeSum(omega - mode[2], dt, samples);
            const auto negative = std::conj(phase) * tiltedTimeSum(omega + mode[2], dt, samples);
            result += mode[3]
                      * (r.plane_m.polarization_m == Polarization::TE
                             ? std::complex<double>(0, W / 2) * (positive - negative)
                             : W * W * std::cos(mode[2] * dt / 2) * (positive + negative) / 2.0);
        }
        return result;
    }
}  // namespace chdr::apps::planarInterface
