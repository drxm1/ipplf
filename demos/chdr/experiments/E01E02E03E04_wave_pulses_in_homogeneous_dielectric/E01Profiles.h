/** @file E01Profiles.h
 * @brief Numerical and continuum longitudinal profiles at E01's error observations.
 */
#pragma once

#include <Kokkos_Core.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <mpi.h>
#include <stdexcept>
#include <string>
#include <vector>

#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01Measurements.h"

namespace chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric {
    namespace e01ProfileDetail {
        /// Fixed global cell indices on the transverse mesh, which E01 does not refine.
        inline constexpr int LineCellX = 7, LineCellY = 7;

        /** @brief Extract the owned part of the line without copying a full field to the host.
         * @details Nonowners contribute zeros. Each global z index has exactly one owner;
         * ghosts are excluded, including at the periodic endpoint. The temporary view uses
         * the field's execution/memory spaces and stores one scalar per global z cell.
         */
        template <class Field>
        auto ownedLine(Field& field, unsigned component) {
            const auto local  = field.getLayout().getLocalNDIndex();
            const auto global = field.getLayout().getDomain();
            const auto view   = field.getView();
            const int ghost = field.getNghost(), firstZ = local[2].first() - global[2].first();
            using ExecutionSpace = typename Field::view_type::execution_space;
            Kokkos::View<double*, typename Field::view_type::device_type> line(
                "E01 profile samples", global[2].length());
            if (local[0].first() <= LineCellX && LineCellX <= local[0].last()
                && local[1].first() <= LineCellY && LineCellY <= local[1].last()) {
                const int i = LineCellX - local[0].first() + ghost;
                const int j = LineCellY - local[1].first() + ghost;
                Kokkos::parallel_for(
                    "Extract E01 profile",
                    Kokkos::RangePolicy<ExecutionSpace>(0, local[2].length()),
                    KOKKOS_LAMBDA(int k) { line(firstZ + k) = view(i, j, k + ghost)[component]; });
            }
            Kokkos::fence();
            return Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, line);
        }

        /// Assemble the unique owned samples in global z order on rank zero.
        template <class Field>
        std::vector<double> gatherLine(Field& field, unsigned component) {
            const auto local = ownedLine(field, component);
            const int count  = field.getLayout().getDomain()[2].length();
            std::vector<double> global(count);
            MPI_Reduce(local.data(), global.data(), count, MPI_DOUBLE, MPI_SUM, 0,
                       field.getLayout().comm);
            return global;
        }

        /// Cartesian component of the nonzero E01 reference channel.
        inline unsigned component(Quantity quantity) {
            return quantity == Quantity::Magnetic ? 1 : 0;
        }

        /// CSV channel labels distinguish the three profiles without changing their units.
        inline const char* name(Quantity quantity) {
            if (quantity == Quantity::VectorPotential)
                return "Ax";
            return quantity == Quantity::Electric ? "Ex" : "By";
        }
    }  // namespace e01ProfileDetail

    /** @brief Write profiles.csv on the fixed native line through global cells (7,7).
     * @details E01's fixed transverse mesh places Ax, Ex and By at
     * \f$x=1/2,\ y=15/32\f$. Their longitudinal positions follow LevelContract
     * [ryu2016potentialfdtd, Appendix A, Eqs. (19)--(24)]: By is half a cell beyond
     * Ax/Ex. Its last owned point is at the periodic endpoint, not a duplicate ghost.
     * The caller supplies R2 magneticField() and the same times as the error measurements.
     * @pre A nonempty output directory has already been prepared collectively. An empty
     * directory disables profiles, preserving the error-only CSV on standard output.
     */
    class E01ProfileCsv {
    public:
        E01ProfileCsv(const std::string& directory, const E01Reference& reference,
                      LevelContract contract, MPI_Comm communicator)
            : reference_m(reference)
            , contract_m(contract)
            , enabled_m(!directory.empty()) {
            MPI_Comm_rank(communicator, &rank_m);
            if (!enabled_m || rank_m != 0)
                return;
            file_m.exceptions(std::ios::failbit | std::ios::badbit);
            file_m.imbue(std::locale::classic());
            file_m.open(std::filesystem::path(directory) / "profiles.csv");
            file_m
                << std::setprecision(std::numeric_limits<double>::max_digits10)
                << "component,checkpoint,target_time,step,actual_time,k,x,y,z,numerical,analytic\n";
        }

        /** @brief Collect one component at the same state used for its global error.
         * @pre Every rank calls with the same quantity and observation labels; quantity is
         * VectorPotential, Electric or Magnetic and field is the corresponding solver output.
         */
        template <class Field>
        void write(Field& field, Quantity quantity, const char* checkpoint, double target,
                   long long step, double time) {
            if (!enabled_m)
                return;
            const unsigned component = e01ProfileDetail::component(quantity);
            const unsigned slot      = component + (quantity == Quantity::VectorPotential ? 1 : 0);
            const auto values        = e01ProfileDetail::gatherLine(field, slot);
            if (rank_m == 0) {
                writeRows(field, quantity, component, checkpoint, target, step, time, values);
                file_m.flush();
            }
        }

        /// Close the profile file before the completed run record is published.
        void finish() {
            if (enabled_m && rank_m == 0)
                file_m.close();
        }

    private:
        /// Write native coordinates and evaluate the continuum reference at each numerical sample.
        template <class Field>
        void writeRows(Field& field, Quantity quantity, unsigned component, const char* checkpoint,
                       double target, long long step, double time,
                       const std::vector<double>& values) {
            const auto h = field.get_mesh().getMeshSpacing(), origin = field.get_mesh().getOrigin();
            const auto offset = contract_m.positionOffset(quantity, component);
            const double x    = origin[0] + (e01ProfileDetail::LineCellX + 0.5 + offset[0]) * h[0];
            const double y    = origin[1] + (e01ProfileDetail::LineCellY + 0.5 + offset[1]) * h[1];
            const int firstZ  = field.getLayout().getDomain()[2].first();
            for (std::size_t n = 0; n < values.size(); ++n) {
                const int k        = firstZ + static_cast<int>(n);
                const double z     = origin[2] + (k + 0.5 + offset[2]) * h[2];
                const double exact = e01MeasurementDetail::referenceValue(reference_m, quantity, x,
                                                                          y, z, time, component);
                if (!std::isfinite(values[n]) || !std::isfinite(exact))
                    throw std::runtime_error("E01 profile contains a nonfinite sample");
                file_m << e01ProfileDetail::name(quantity) << ',' << checkpoint << ',' << target
                       << ',' << step << ',' << time << ',' << k << ',' << x << ',' << y << ',' << z
                       << ',' << values[n] << ',' << exact << '\n';
            }
        }

        E01Reference reference_m;
        LevelContract contract_m;
        bool enabled_m;
        int rank_m = 0;
        std::ofstream file_m;
    };
}  // namespace chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric
