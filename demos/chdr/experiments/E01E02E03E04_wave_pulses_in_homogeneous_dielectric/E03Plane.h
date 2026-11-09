/** @file E03Plane.h
 * @brief Colocated E03 central sections of the fields and native continuum errors.
 */
#pragma once

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <vector>

#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03Measurements.h"

namespace chdr::experiments::e03 {
    namespace planeDetail {
        /// Tensor-product linear interpolation from the eight native bracketing points.
        template <class Value, class Sample>
        KOKKOS_INLINE_FUNCTION Value trilinear(Kokkos::Array<int, 3> anchor,
                                               Kokkos::Array<double, 3> delta, Sample sample) {
            for (unsigned d = 0; d < 3; ++d) {
                if (delta[d] < 0) {
                    --anchor[d];
                    delta[d] += 1;
                }
            }
            Value result(0.0);
            for (unsigned corner = 0; corner < 8; ++corner) {
                auto index    = anchor;
                double weight = 1;
                for (unsigned d = 0; d < 3; ++d) {
                    const bool upper = (corner & (1u << d)) != 0;
                    index[d] += upper;
                    weight *= upper ? delta[d] : 1 - delta[d];
                }
                if (weight != 0)
                    result += weight * sample(index);
            }
            return result;
        }

        /** @brief Interpolate the numerical sample and its native error separately.
         * @details A centred z=0 plane is half a cell above its owner anchor. A and E use
         * positive-face offsets and B uses edge offsets [ryu2016potentialfdtd, Appendix A].
         * The native difference is formed before interpolation as prescribed by E03.typ.
         */
        template <class View>
        KOKKOS_INLINE_FUNCTION ippl::Vector<double, 2> planeComponent(
            View view, Kokkos::Array<int, 3> anchor, ippl::Vector<double, 3> centre, double h,
            LevelContract contract, const E03Reference& ref, Quantity quantity, double time,
            unsigned c) {
            const auto offset = contract.positionOffset(quantity, c);
            return trilinear<ippl::Vector<double, 2>>(
                anchor, {-offset[0], -offset[1], 0.5 - offset[2]}, [&](auto index) {
                    auto x = centre;
                    for (unsigned d = 0; d < 3; ++d)
                        x[d] += (index[d] - anchor[d] + offset[d]) * h;
                    const double value = view(index[0], index[1], index[2])[c];
                    const double exact = e03MeasurementDetail::referenceValue(ref, quantity, x[0],
                                                                              x[1], x[2], time, c);
                    return ippl::Vector<double, 2>{value, value - exact};
                });
        }

        /// Use one z slab and transfer only the requested plane, including the seam halo.
        template <class Field, class Plane>
        void samplePlane(Field& field, Plane plane, int columns, LevelContract contract,
                         E03Reference ref, Quantity quantity, double time) {
            field.fillHalo();
            const auto local = field.getLayout().getLocalNDIndex();
            const int n = plane.extent_int(0), centreBelow = n / 2 - 1;
            if (centreBelow < local[2].first() || centreBelow > local[2].last())
                return;
            const int ghost = field.getNghost(), k = centreBelow - local[2].first() + ghost;
            const auto view = field.getView();
            const double h  = field.get_mesh().getMeshSpacing()[0];
            using Exec      = typename Field::view_type::execution_space;
            Kokkos::parallel_for(
                "E03 z=0 section", Kokkos::MDRangePolicy<Exec, Kokkos::Rank<2>>({0, 0}, {n, n}),
                KOKKOS_LAMBDA(int i, int j) {
                    const ippl::Vector<double, 3> centre{(i + 0.5) * h, (j + 0.5) * h,
                                                         (centreBelow + 0.5) * h};
                    for (unsigned c = 0; c < 3; ++c) {
                        const auto pair = planeComponent(view, {i + ghost, j + ghost, k}, centre, h,
                                                         contract, ref, quantity, time, c);
                        plane(i, j, columns + 2 * c)     = pair[0];
                        plane(i, j, columns + 2 * c + 1) = pair[1];
                    }
                });
        }

        /// Write one row for each declared plotting point, in the specification's coordinates.
        inline void writePlaneRows(const std::string& directory, int n, double h, long long step,
                                   double time, const std::vector<double>& samples) {
            std::ofstream file;
            file.exceptions(std::ios::failbit | std::ios::badbit);
            file.open(std::filesystem::path(directory) / "plane-z0.csv", std::ios::app);
            file << std::setprecision(std::numeric_limits<double>::max_digits10);
            if (step == 0)
                file << "step,time,x,y,z,Ex,error_Ex,Ey,error_Ey,Ez,error_Ez,"
                        "Bx,error_Bx,By,error_By,Bz,error_Bz\n";
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    file << step << ',' << time << ',' << (i + 0.5) * h - E03BoxCentre << ','
                         << (j + 0.5) * h - E03BoxCentre << ",0";
                    for (unsigned c = 0; c < 12; ++c)
                        file << ',' << samples[12 * (i * n + j) + c];
                    file << '\n';
                }
            }
            file.close();
        }
    }  // namespace planeDetail

    /** @brief Export E and R2 B on one physical z=0 plane with the z-only MPI layout.
     * @details Each point is owned by the rank holding the lower central z cell. Ghosts
     * supply bracket samples but never contribute as additional plotting points.
     */
    template <class Solver>
    void writeE03Plane(Solver& state, const E03Reference& ref, const std::string& directory,
                       long long step) {
        auto& e      = state.electricField();
        const int n  = e.getLayout().getDomain()[0].length();
        using Device = typename EMFieldsType::view_type::device_type;
        Kokkos::View<double***, Kokkos::LayoutRight, Device> plane("E03 field section", n, n, 12);
        const double time = state.labelTime(Quantity::Electric, step);
        planeDetail::samplePlane(e, plane, 0, state.levelContract(), ref, Quantity::Electric, time);
        planeDetail::samplePlane(state.magneticField(), plane, 6, state.levelContract(), ref,
                                 Quantity::Magnetic, time);
        Kokkos::fence();
        const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, plane);
        const auto comm = e.getLayout().comm;
        int rank;
        MPI_Comm_rank(comm, &rank);
        std::vector<double> samples(rank == 0 ? plane.size() : 0);
        MPI_Reduce(host.data(), samples.data(), 12 * n * n, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (rank == 0)
            planeDetail::writePlaneRows(directory, n, e.get_mesh().getMeshSpacing()[0], step, time,
                                        samples);
    }
}  // namespace chdr::experiments::e03
