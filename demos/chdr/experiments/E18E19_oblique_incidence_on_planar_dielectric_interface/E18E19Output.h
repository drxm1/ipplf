/** @file E18E19Output.h
 * @brief Native field profiles, Fourier measurements and snapshots for E18 and E19.
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <mpi.h>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "core/input/RunInput.h"
#include "core/output/FieldOutput.h"
#include "core/output/RunRecord.h"
#include "core/simulation/SimulationConfig.h"
#include "core/solver/ChdrSolver.h"
#include "experiments/E18E19_oblique_incidence_on_planar_dielectric_interface/E18E19AnalyticContinuumReference.h"

namespace chdr::experiments::e18e19 {
    namespace outputDetail {
        /// Broadcast one rank-zero failure before any rank enters another observation.
        inline void broadcastFailure(std::string failure, MPI_Comm communicator) {
            int length = static_cast<int>(failure.size());
            MPI_Bcast(&length, 1, MPI_INT, 0, communicator);
            failure.resize(length);
            if (length > 0) {
                MPI_Bcast(failure.data(), length, MPI_CHAR, 0, communicator);
                throw std::runtime_error("E18/E19 output failed: " + failure);
            }
        }

        /** @brief Perform one rank-zero file operation and propagate any failure.
         * @details Every rank reaches the broadcast before another observation can start.
         * This follows the collective file publication used by RunRecord and FieldOutput.
         */
        template <class Action>
        void onRoot(MPI_Comm communicator, Action action) {
            int rank = 0;
            MPI_Comm_rank(communicator, &rank);
            std::string failure;
            if (rank == 0) {
                try {
                    action();
                } catch (const std::exception& error) {
                    failure = error.what();
                }
            }
            broadcastFailure(failure, communicator);
        }

        /// Open one table with enough decimal digits to recover every stored double.
        inline void openCsv(std::ofstream& file, const std::string& directory, const char* name,
                            const char* header) {
            file.exceptions(std::ios::failbit | std::ios::badbit);
            file.imbue(std::locale::classic());
            file.open(std::filesystem::path(directory) / name);
            file << std::setprecision(std::numeric_limits<double>::max_digits10) << header << '\n';
        }

        /** @brief Copy the owned part of the longitudinal line at global i=j=0.
         * @details Each rank allocates one scalar per global z cell in the field memory
         * space. A kernel writes its owned samples. The zero-initialized remaining entries
         * allow one MPI sum to assemble the line. Only this line moves to host memory.
         */
        template <class Field>
        auto ownedLine(Field& field, unsigned component) {
            const auto local     = field.getLayout().getLocalNDIndex();
            const auto global    = field.getLayout().getDomain();
            const auto view      = field.getView();
            const int ghost      = field.getNghost();
            const int firstZ     = local[2].first() - global[2].first();
            using ExecutionSpace = typename Field::view_type::execution_space;
            Kokkos::View<double*, typename Field::view_type::device_type> line(
                "E18/E19 profile samples", global[2].length());
            if (local[0].first() <= 0 && local[0].last() >= 0 && local[1].first() <= 0
                && local[1].last() >= 0) {
                const int i = ghost - local[0].first(), j = ghost - local[1].first();
                Kokkos::parallel_for(
                    "E18/E19 extract profile",
                    Kokkos::RangePolicy<ExecutionSpace>(0, local[2].length()),
                    KOKKOS_LAMBDA(int k) { line(firstZ + k) = view(i, j, k + ghost)[component]; });
            }
            Kokkos::fence();
            return Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, line);
        }

        /// Assemble the unique owned samples in global z order on rank zero.
        template <class Field>
        std::vector<double> gatherLine(Field& field, unsigned component) {
            const auto owned = ownedLine(field, component);
            const int count  = field.getLayout().getDomain()[2].length();
            std::vector<double> values(count);
            MPI_Reduce(owned.data(), values.data(), count, MPI_DOUBLE, MPI_SUM, 0,
                       field.getLayout().comm);
            return values;
        }

        /** @brief Locate a field component on the selected longitudinal line.
         * @details The mesh origin and LevelContract determine all component offsets.
         * These are the native staggered locations of the reconstruction
         * [ryu2016potentialfdtd, Appendix A, Eqs. (19)--(24)].
         */
        template <class Field>
        auto position(Field& field, LevelContract contract, Quantity quantity, unsigned component,
                      int k) {
            const auto spacing = field.get_mesh().getMeshSpacing();
            const auto origin  = field.get_mesh().getOrigin();
            const auto offset  = contract.positionOffset(quantity, component);
            return std::array<double, 3>{origin[0] + (0.5 + offset[0]) * spacing[0],
                                         origin[1] + (0.5 + offset[1]) * spacing[1],
                                         origin[2] + (k + 0.5 + offset[2]) * spacing[2]};
        }

        /// Write one measured Fourier ratio with the reference at precisely the same q.
        inline void reflectionRow(std::ofstream& file,
                                  const E18E19AnalyticContinuumReference& reference, double q,
                                  double incident, double reflected) {
            if (!std::isfinite(incident) || incident <= 0 || !std::isfinite(reflected)
                || reflected < 0 || !std::isfinite(reflected / incident))
                throw std::runtime_error("Invalid Fourier amplitudes in reflection.csv");
            const double angle = std::atan(reference.KX / q) * 180 / Kokkos::numbers::pi;
            file << q << ',' << angle << ',' << incident << ',' << reflected << ','
                 << reflected / incident << ',' << std::abs(reference.reflection(q)) << '\n';
        }
    }  // namespace outputDetail

    /** @brief Save the packet's continuum modes exactly as used by initialization.
     * @details The signed coefficient is Gamma in E18 and R in E19. The normalized
     * amplitudes and frequencies retain the inputs of the energy-weighted reference
     * in E18.typ and E19.typ [ellingson2020ElectromagneticsVol2, Sections 5.6--5.7].
     * @pre All ranks call after preparing the run directory.
     */
    inline void writeModes(const std::string& directory,
                           const E18E19AnalyticContinuumReference& reference,
                           MPI_Comm communicator) {
        outputDetail::onRoot(communicator, [&] {
            std::ofstream file;
            outputDetail::openCsv(file, directory, "modes.csv",
                                  "q,omega,q_second,amplitude,reflection_coefficient");
            for (const auto& mode : reference.modes())
                file << mode.q_m << ',' << mode.omega_m << ',' << mode.qSecond_m << ','
                     << mode.amplitude_m << ',' << mode.reflection_m << '\n';
            file.close();
        });
    }

    /** @brief Save E18.meas01 or E19.meas01 from the two recorded Fourier observations.
     * @details The order of rows matches measuredWavenumbers(), with the Brewster sample
     * last. Input amplitudes are the MPI-global magnitudes returned by the measurement.
     * @pre All ranks call with the same initial and final measurements.
     */
    inline void writeReflection(const std::string& directory,
                                const E18E19AnalyticContinuumReference& reference,
                                const std::array<double, 12>& incident,
                                const std::array<double, 12>& reflected, MPI_Comm communicator) {
        outputDetail::onRoot(communicator, [&] {
            std::ofstream file;
            outputDetail::openCsv(file, directory, "reflection.csv",
                                  "q,angle_degrees,incident_magnitude,reflected_magnitude,"
                                  "reflection_magnitude,reference_magnitude");
            const auto wavenumbers = reference.measuredWavenumbers();
            for (std::size_t index = 0; index < wavenumbers.size(); ++index)
                outputDetail::reflectionRow(file, reference, wavenumbers[index], incident[index],
                                            reflected[index]);
            file.close();
        });
    }

    /** @brief Save the energies and the running maximum at the second face.
     * @details Each row is one shared observation at the native field time.
     * field_at_second_face is the maximum through this observation divided by the
     * initial global electric peak, as defined by E18.meas03.
     * @pre All ranks construct, append and close in the same order.
     */
    class MeasurementCsv {
    public:
        MeasurementCsv(const std::string& directory, MPI_Comm communicator)
            : communicator_m(communicator) {
            outputDetail::onRoot(communicator_m, [&] {
                outputDetail::openCsv(file_m, directory, "measurements.csv",
                                      "step,time,U_side1,U_side2,field_at_second_face");
            });
        }

        /// Save one complete observation and make a root I/O failure visible to every rank.
        void row(long long step, double time, double side1, double side2, double secondFace) {
            outputDetail::onRoot(communicator_m, [&] {
                if (!std::isfinite(time) || !std::isfinite(side1) || !std::isfinite(side2)
                    || !std::isfinite(secondFace))
                    throw std::runtime_error("Nonfinite value in measurements.csv");
                file_m << step << ',' << time << ',' << side1 << ',' << side2 << ',' << secondFace
                       << '\n';
                file_m.flush();
            });
        }

        /// Close all CSV data before the completed run record is published.
        void finish() {
            outputDetail::onRoot(communicator_m, [&] {
                file_m.close();
            });
        }

    private:
        MPI_Comm communicator_m;
        std::ofstream file_m;
    };

    /** @brief Write selected native E and R2 B components along global i=j=0.
     * @details Each analytic sample uses its numerical component's coordinate and
     * actual native time. The reference is evaluated only for the gathered line.
     * Full-domain arrays stay in the field memory space until a snapshot is requested.
     * @pre All ranks construct, write and finish in the same order.
     */
    template <class Reference = E18E19AnalyticContinuumReference>
    class ProfileCsv {
    public:
        ProfileCsv(const std::string& directory, const Reference& reference, LevelContract contract,
                   MPI_Comm communicator)
            : reference_m(reference)
            , contract_m(contract)
            , communicator_m(communicator) {
            outputDetail::onRoot(communicator_m, [&] {
                outputDetail::openCsv(file_m, directory, "profiles.csv",
                                      "component,checkpoint,target_time,step,actual_time,k,x,y,z,"
                                      "numerical,analytic");
            });
        }

        /// Gather one component and retain its native placement and actual observation time.
        template <class Field>
        void write(Field& field, Quantity quantity, unsigned component, const char* checkpoint,
                   double target, long long step, double time) {
            if ((quantity != Quantity::Electric && quantity != Quantity::Magnetic) || component > 2)
                throw std::invalid_argument("Profiles require a Cartesian E or B component");
            const auto values = outputDetail::gatherLine(field, component);
            outputDetail::onRoot(communicator_m, [&] {
                writeRows(field, quantity, component, checkpoint, target, step, time, values);
                file_m.flush();
            });
        }

        /// Close the profile table before publishing a completed run record.
        void finish() {
            outputDetail::onRoot(communicator_m, [&] {
                file_m.close();
            });
        }

    private:
        /// Save every owned longitudinal sample, including the native periodic endpoint.
        template <class Field>
        void writeRows(Field& field, Quantity quantity, unsigned component, const char* checkpoint,
                       double target, long long step, double time,
                       const std::vector<double>& values) {
            const int firstZ = field.getLayout().getDomain()[2].first();
            for (std::size_t index = 0; index < values.size(); ++index) {
                const int k = firstZ + static_cast<int>(index);
                const auto point =
                    outputDetail::position(field, contract_m, quantity, component, k);
                const double exact =
                    quantity == Quantity::Electric
                        ? reference_m.referenceE(point[0], point[1], point[2], time, component)
                        : reference_m.referenceB(point[0], point[1], point[2], time, component);
                if (!std::isfinite(values[index]) || !std::isfinite(exact))
                    throw std::runtime_error("Nonfinite value in profiles.csv");
                file_m << (quantity == Quantity::Electric ? 'E' : 'B') << "xyz"[component] << ','
                       << checkpoint << ',' << target << ',' << step << ',' << time << ',' << k
                       << ',' << point[0] << ',' << point[1] << ',' << point[2] << ','
                       << values[index] << ',' << exact << '\n';
            }
        }

        Reference reference_m;
        LevelContract contract_m;
        MPI_Comm communicator_m;
        std::ofstream file_m;
    };

    /** @brief Save E and R2 B at the five prescribed potential-step targets.
     * @details FieldOutput records each array's native component offsets and time
     * from LevelContract. The caller supplies the five step indices, including zero.
     * @pre All ranks construct, observe and finish in the same order.
     */
    class FieldSnapshots {
        using FieldSolver = solver::ChdrSolver<solver::StandardSolverPeriodic>;

    public:
        FieldSnapshots(const input::ChdrYamlRunInput& request,
                       const simulation::ChdrCartesianGridSettings& grid, FieldSolver& state,
                       std::vector<int> prescribedSteps)
            : steps_m(std::move(prescribedSteps))
            , writer_m(std::filesystem::path(request.outputDirectory_m) / "fields",
                       output::OutputGeometry{grid.cells_m, {0, 0, 0}, grid.spacing_m},
                       state.A_n.getLayout().comm, "fields", true,
                       output::fieldOutputScales(request.units_m)) {}

        /// Write the existing reconstructed arrays when the step is a prescribed checkpoint.
        void observe(long long step, FieldSolver& state) {
            if (std::find(steps_m.begin(), steps_m.end(), step) == steps_m.end())
                return;
            const std::array fields{
                output::VectorFieldSnapshot{"E", state.electricField(), Quantity::Electric, 1.0},
                output::VectorFieldSnapshot{"B_R2", state.magneticField(), Quantity::Magnetic,
                                            1.0}};
            const output::FrameRecord frame{static_cast<int>(step),
                                            state.labelTime(Quantity::VectorPotential, step),
                                            state.getDt()};
            writer_m.writeFrame(frame, state.levelContract(), fields);
        }

        /// Publish the complete ParaView series after its last frame.
        void finish() { writer_m.writeSeries(); }

        /// Record the exact step selection and the native field conventions.
        void record(output::RunRecordDescription& description) const {
            auto node       = description.effective_m["field_snapshots"];
            node["enabled"] = true;
            node["prescribed_steps"].set(steps_m);
            node["arrays"]            = "E,B_R2";
            node["sampling_metadata"] = "derived from solver LevelContract";
            node["series"]            = "fields/fields.pvti.series";
        }

    private:
        std::vector<int> steps_m;
        output::FieldOutput writer_m;
    };
}  // namespace chdr::experiments::e18e19
