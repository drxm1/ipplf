/** @file E04Monitor.h
 * @brief E04.meas02 sign-split Poynting power through the fixed interior cube.
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <mpi.h>
#include <stdexcept>

#include "core/base/FieldContainer.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03Plane.h"

namespace chdr::experiments::e04 {
    namespace monitorDetail {
        /// One face's owned quadrature range, expressed in global cell indices.
        struct FaceRange {
            int normalCell_m;
            Kokkos::Array<int, 2> lower_m, upper_m;
        };

        /** @brief Assign each midpoint to the cell immediately below its face.
         * @details The cube faces are one quarter and three quarters across the box.
         * Tangential midpoint indices cover its middle half. The owning cell supplies
         * native bracket samples through at most one halo. On z-only slabs this clips
         * x/y faces in z and assigns each z face to exactly one slab.
         */
        template <class Layout>
        FaceRange faceRange(const Layout& layout, unsigned face) {
            const auto& global    = layout.getDomain();
            const auto local      = layout.getLocalNDIndex();
            const unsigned normal = face / 2, a = (normal + 1) % 3, b = (normal + 2) % 3;
            const int quarter    = global[normal].length() / 4;
            const int normalCell = global[normal].first() + ((face % 2) ? 3 : 1) * quarter - 1;
            if (normalCell < local[normal].first() || normalCell > local[normal].last())
                return {normalCell, {0, 0}, {0, 0}};
            return {normalCell,
                    {std::max(global[a].first() + quarter, local[a].first()),
                     std::max(global[b].first() + quarter, local[b].first())},
                    {std::min(global[a].first() + 3 * quarter, local[a].last() + 1),
                     std::min(global[b].first() + 3 * quarter, local[b].last() + 1)}};
        }

        /** @brief Interpolate one native component to a monitor midpoint.
         * @details The face is half a cell above its normal anchor and lies at tangent
         * cell centres. LevelContract supplies the native E/B offsets
         * [ryu2016potentialfdtd, Appendix A, Eqs. (19)--(24)].
         */
        template <class View>
        KOKKOS_INLINE_FUNCTION double componentAtFace(View view, Kokkos::Array<int, 3> anchor,
                                                      LevelContract contract, Quantity quantity,
                                                      unsigned component, unsigned normal) {
            const auto offset = contract.positionOffset(quantity, component);
            Kokkos::Array<double, 3> delta;
            for (unsigned d = 0; d < 3; ++d)
                delta[d] = (d == normal ? 0.5 : 0.0) - offset[d];
            return e03::planeDetail::trilinear<double>(anchor, delta, [&](auto index) {
                return view(index[0], index[1], index[2])[component];
            });
        }

        /** @brief Return outward flux density after native components are colocated.
         * @details With unit permeability, the flux is
         * \f$(\mathbf E\times\mathbf B)\cdot\mathbf n\f$
         * [ellingson2020ElectromagneticsVol2, section 3.2]. Both inputs have the same
         * R2 half-time. Only the four tangential components enter this face's flux.
         */
        template <class View>
        KOKKOS_INLINE_FUNCTION double fluxAtFace(View electric, View magnetic,
                                                 Kokkos::Array<int, 3> anchor,
                                                 LevelContract contract, unsigned face) {
            const unsigned normal = face / 2, a = (normal + 1) % 3, b = (normal + 2) % 3;
            const double ea =
                componentAtFace(electric, anchor, contract, Quantity::Electric, a, normal);
            const double eb =
                componentAtFace(electric, anchor, contract, Quantity::Electric, b, normal);
            const double ba =
                componentAtFace(magnetic, anchor, contract, Quantity::Magnetic, a, normal);
            const double bb =
                componentAtFace(magnetic, anchor, contract, Quantity::Magnetic, b, normal);
            return (face % 2 ? 1.0 : -1.0) * (ea * bb - eb * ba);
        }

        /// Locally integrated inward/outward flux and the number of invalid samples.
        struct FacePower {
            double inward_m = 0, outward_m = 0;
            std::uint64_t nonfinite_m = 0;
        };

        /// Reduce one face on the field execution space without transferring a field volume.
        template <class Solver>
        FacePower localFacePower(Solver& state, unsigned face) {
            const auto& layout = state.electricField().getLayout();
            const auto range   = faceRange(layout, face);
            FacePower local;
            if (range.lower_m[0] >= range.upper_m[0] || range.lower_m[1] >= range.upper_m[1])
                return local;
            const auto first      = layout.getLocalNDIndex().first();
            const auto electric   = state.electricField().getView(),
                       magnetic   = state.magneticField().getView();
            const auto contract   = state.levelContract();
            const int ghost       = state.electricField().getNghost();
            const unsigned normal = face / 2, a = (normal + 1) % 3, b = (normal + 2) % 3;
            using Execution = typename decltype(electric)::execution_space;
            Kokkos::parallel_reduce(
                "E04 interior monitor face",
                Kokkos::MDRangePolicy<Execution, Kokkos::Rank<2>>(range.lower_m, range.upper_m),
                KOKKOS_LAMBDA(int u, int v, double& inward, double& outward,
                              std::uint64_t& invalid) {
                    Kokkos::Array<int, 3> anchor;
                    anchor[normal]    = range.normalCell_m - first[normal] + ghost;
                    anchor[a]         = u - first[a] + ghost;
                    anchor[b]         = v - first[b] + ghost;
                    const double flux = fluxAtFace(electric, magnetic, anchor, contract, face);
                    if (!Kokkos::isfinite(flux)) {
                        ++invalid;
                        return;
                    }
                    inward += Kokkos::max(-flux, 0.0);
                    outward += Kokkos::max(flux, 0.0);
                },
                local.inward_m, local.outward_m, local.nonfinite_m);
            return local;
        }

        /// Require the E04 cube and a finite positive reference-energy normalization.
        template <class Solver>
        void requireMonitorGrid(Solver& state, double speed, double initialReferenceEnergy) {
            e03::e03MeasurementDetail::requireE03Contract(state.levelContract());
            const auto& domain = state.A_n.getLayout().getDomain();
            const auto h       = state.A_n.get_mesh().getMeshSpacing();
            const auto n       = domain[0].length();
            if (!std::isfinite(speed) || speed <= 0 || !std::isfinite(initialReferenceEnergy)
                || initialReferenceEnergy <= 0 || n % 4 != 0)
                throw std::invalid_argument("Invalid E04 monitor grid or normalization");
            for (unsigned d = 0; d < 3; ++d)
                if (domain[d].length() != n || h[d] != h[0] || std::abs(n * h[d] - 12.0) > 1e-12)
                    throw std::invalid_argument(
                        "E04 monitor requires a uniform 12 by 12 by 12 cube");
        }

        /// Combine face powers and invalid counts before applying the fixed normalization.
        inline std::array<double, 12> globalPowers(std::array<double, 13> values, double scale,
                                                   MPI_Comm communicator) {
            MPI_Allreduce(MPI_IN_PLACE, values.data(), values.size(), MPI_DOUBLE, MPI_SUM,
                          communicator);
            if (values[12] != 0)
                throw std::runtime_error("E04 monitor contains nonfinite field data");
            std::array<double, 12> result;
            for (unsigned i = 0; i < result.size(); ++i) {
                result[i] = values[i] * scale;
                if (!std::isfinite(result[i]))
                    throw std::runtime_error("E04 monitor power is nonfinite");
            }
            return result;
        }
    }  // namespace monitorDetail

    /** @brief E04.meas02 normalized inward/outward power on all six monitor faces.
     * @return Twelve values in x-minus, x-plus, y-minus, y-plus, z-minus, z-plus
     * face order, with inward then outward for each face. All ranks receive the result.
     * @details Exchange reconstructed E and R2 B halos before colocation. The midpoint
     * quadrature sums \f$h^2\max(\mp j_n,0)\f$ over uniquely owned points, then
     * multiplies by \f$w/(vU_{0,h})\f$ with \f$w=1\f$ as defined by E04.typ.
     * The scalar histories can subsequently be integrated in \f$\eta=vt/w\f$.
     */
    template <class Solver>
    std::array<double, 12> measureE04Monitor(Solver& state, double speed,
                                             double initialReferenceEnergy) {
        monitorDetail::requireMonitorGrid(state, speed, initialReferenceEnergy);
        state.electricField().fillHalo();
        state.magneticField().fillHalo();
        std::array<double, 13> values{};
        for (unsigned face = 0; face < 6; ++face) {
            const auto local     = monitorDetail::localFacePower(state, face);
            values[2 * face]     = local.inward_m;
            values[2 * face + 1] = local.outward_m;
            values[12] += local.nonfinite_m;
        }
        const double h = state.A_n.get_mesh().getMeshSpacing()[0];
        return monitorDetail::globalPowers(values, h * h / speed / initialReferenceEnergy,
                                           state.A_n.getLayout().comm);
    }
}  // namespace chdr::experiments::e04
