/** @file E02Diagnostics.h
 * @brief Initial-normalized residual fields and the native electric probe of E02.
 */
#pragma once

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>

#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01E02FieldSnapshots.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01E02Measurements.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01E02Profiles.h"

namespace chdr::experiments::e01e02 {
    /** @brief Save the two E02 diagnostics without changing the numerical fields.
     * @details E02.typ defines the final/initial discrete norms and the probe at the first
     * electric sample right of z=4w. Norms include all components. The probe transfers
     * only one scalar from the owning device/rank on every step; no field volume is copied.
     * All ranks must enter observe and finish in the same order. The common run output
     * prepares the directory before construction, and its MPI failure handler handles I/O errors.
     */
    class E02Diagnostics {
    public:
        E02Diagnostics(const std::string& directory, PulseSolver& state,
                       const simulation::ChdrCartesianGridSettings& grid)
            : communicator_m(state.A_n.getLayout().comm)
            , initialSquared_m{squaredFields(state)}
            , initialTime_m(state.labelTime(Quantity::Electric, 0)) {
            MPI_Comm_rank(communicator_m, &rank_m);
            for (double value : initialSquared_m)
                if (value <= 0)
                    throw std::runtime_error("E02 initial field norm must be positive");
            locateProbe(state, grid);
            open(probe_m, directory, "probe.csv", "step,actual_time,i,j,k,x,y,z,Ex\n");
            open(residuals_m, directory, "residuals.csv",
                 "measurement,checkpoint,target_time,step,actual_time,initial_time,"
                 "initial_l2,current_l2,relative_residual\n");
        }

        /// Keep the complete signed trace, then measure a final or explicitly partial endpoint.
        void observe(long long step, PulseSolver& state, int finalStep, int protocolSteps,
                     double finalTarget) {
            writeProbe(step, state);
            if (step == finalStep) {
                const bool complete = finalStep == protocolSteps;
                writeResiduals(step, state, complete ? "final" : "partial_final",
                               complete ? finalTarget : state.labelTime(Quantity::Electric, step));
            }
        }

        /// Flush diagnostics before publishing the completed simulation record.
        void finish() {
            if (rank_m == 0) {
                probe_m.close();
                residuals_m.close();
            }
        }

        /// Identify sampling and fixed numerical normalizations for portable analysis.
        void record(output::RunRecordDescription& description) const {
            auto probe    = description.effective_m["probe"];
            probe["file"] = "probe.csv";
            probe["global_cell"].set(
                std::vector<int>{globalCell_m[0], globalCell_m[1], globalCell_m[2]});
            probe["position"].set(std::vector<double>{position_m[0], position_m[1], position_m[2]});
            probe["quantity"]     = "Ex";
            probe["every"]        = 1;
            auto norms            = description.effective_m["residuals"];
            norms["file"]         = "residuals.csv";
            norms["initial_time"] = initialTime_m;
            norms["initial_l2"].set(std::vector<double>{std::sqrt(initialSquared_m[0]),
                                                        std::sqrt(initialSquared_m[1])});
            norms["definition"] = "owned-vector L2 norm relative to its initialized numerical norm";
            norms["magnetic_convention"] = "R2 averaged magnetic field";
        }

    private:
        /// E and magneticField() share the R2 half-step time but have independent native norms.
        static std::array<double, 2> squaredFields(PulseSolver& state) {
            return {pulseSquaredNorm(state.electricField()),
                    pulseSquaredNorm(state.magneticField())};
        }

        /// Set identical native coordinates on all ranks and require exactly one owning rank.
        void locateProbe(PulseSolver& state, const simulation::ChdrCartesianGridSettings& grid) {
            globalCell_m      = {pulseProfileDetail::lineCell(grid.cells_m[0]),
                                 pulseProfileDetail::lineCell(grid.cells_m[1]),
                                 static_cast<int>(std::lround(4.0 / grid.spacing_m[2]))};
            const auto domain = state.electricField().getLayout().getLocalNDIndex();
            const auto offset = state.levelContract().positionOffset(Quantity::Electric, 0);
            const auto origin = state.electricField().get_mesh().getOrigin();
            ownsProbe_m       = true;
            for (unsigned d = 0; d < 3; ++d) {
                ownsProbe_m &=
                    domain[d].first() <= globalCell_m[d] && globalCell_m[d] <= domain[d].last();
                localCell_m[d] =
                    globalCell_m[d] - domain[d].first() + state.electricField().getNghost();
                position_m[d] = origin[d] + (globalCell_m[d] + 0.5 + offset[d]) * grid.spacing_m[d];
            }
            int owners = ownsProbe_m;
            MPI_Allreduce(MPI_IN_PLACE, &owners, 1, MPI_INT, MPI_SUM, communicator_m);
            if (owners != 1)
                throw std::runtime_error("E02 probe must have exactly one owner");
        }

        /// Open one rank-zero CSV using round-trip precision and deterministic numeric formatting.
        void open(std::ofstream& stream, const std::string& directory, const char* filename,
                  const char* header) {
            if (rank_m != 0)
                return;
            stream.exceptions(std::ios::failbit | std::ios::badbit);
            stream.imbue(std::locale::classic());
            stream.open(std::filesystem::path(directory) / filename);
            stream << std::setprecision(std::numeric_limits<double>::max_digits10) << header;
        }

        /// Read one native Ex scalar, on its owner only, and reduce it to the output rank.
        double probeValue(PulseSolver& state) const {
            double local = 0, global = 0;
            if (ownsProbe_m) {
                const auto view      = state.electricField().getView();
                const auto cell      = localCell_m;
                using ExecutionSpace = typename EMFieldsType::view_type::execution_space;
                Kokkos::parallel_reduce(
                    "E02 native Ex probe", Kokkos::RangePolicy<ExecutionSpace>(0, 1),
                    KOKKOS_LAMBDA(int, double& value) {
                        value += view(cell[0], cell[1], cell[2])[0];
                    },
                    local);
            }
            MPI_Reduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, 0, communicator_m);
            return global;
        }

        /// Save actual coordinates and the solver-supplied time without fitting or interpolation.
        void writeProbe(long long step, PulseSolver& state) {
            const double value = probeValue(state);
            if (rank_m != 0)
                return;
            if (!std::isfinite(value))
                throw std::runtime_error("E02 probe contains a nonfinite field value");
            probe_m << step << ',' << state.labelTime(Quantity::Electric, step);
            for (unsigned d = 0; d < 3; ++d)
                probe_m << ',' << globalCell_m[d];
            for (unsigned d = 0; d < 3; ++d)
                probe_m << ',' << position_m[d];
            probe_m << ',' << value << '\n';
        }

        /** @brief Measure residual fields against zero with the fixed initial denominators.
         * @details The discrete norm follows [salari2000CodeVerificationMms, Eq. (15)].
         * These ratios include numerical tails and initialization artifacts as well as
         * reflection; they are not absorbed-energy fractions or isolated reflection coefficients.
         */
        void writeResiduals(long long step, PulseSolver& state, const char* checkpoint,
                            double target) {
            const auto squared = squaredFields(state);
            if (rank_m != 0)
                return;
            const std::array names{"E_L2", "B_L2"};
            const std::array quantities{Quantity::Electric, Quantity::Magnetic};
            for (unsigned c = 0; c < names.size(); ++c)
                residuals_m << names[c] << ',' << checkpoint << ',' << target << ',' << step << ','
                            << state.labelTime(quantities[c], step) << ',' << initialTime_m << ','
                            << std::sqrt(initialSquared_m[c]) << ',' << std::sqrt(squared[c]) << ','
                            << std::sqrt(squared[c] / initialSquared_m[c]) << '\n';
        }

        MPI_Comm communicator_m;
        std::array<double, 2> initialSquared_m;
        double initialTime_m;
        int rank_m       = 0;
        bool ownsProbe_m = false;
        ippl::Vector<int, 3> globalCell_m{}, localCell_m{};
        ippl::Vector<double, 3> position_m{};
        std::ofstream probe_m, residuals_m;
    };
}  // namespace chdr::experiments::e01e02
