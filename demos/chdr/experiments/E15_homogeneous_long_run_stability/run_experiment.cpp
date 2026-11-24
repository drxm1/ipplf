/** @file run_experiment.cpp
 * @brief Source-free O1VS propagation of the discrete travelling mode specified in E15.typ.
 */
#include <filesystem>
#include <fstream>
#include <iostream>

#include "core/input/Application.h"
#include "core/output/RunRecord.h"
#include "experiments/E15_homogeneous_long_run_stability/E15Measurements.h"

namespace chdr::experiments::e15 {
    using Json   = nlohmann::json;
    namespace fs = std::filesystem;

    /// Validate the physical case before allocating distributed fields.
    int resolution(const input::ChdrYamlRunInput& request) {
        input::validateParameters(request, {"cells_per_axis", "save_slices"});
        input::validateObservationNames(request, {});
        if ((request.caseName_m != "E15" && request.caseName_m != "E15-S")
            || request.units_m.system_m != "normalized" || request.duration_m
            || request.outputDirectory_m.empty()
            || !request.parameters_m.has_path("cells_per_axis"))
            throw std::invalid_argument(
                "E15 requires cells_per_axis, normalized units and an output directory");
        const auto value = request.parameters_m["cells_per_axis"];
        if (value.dtype().is_number() && value.number_of_elements() == 1) {
            if (request.caseName_m == "E15-S" && value.to_float64() != 16)
                throw std::invalid_argument("E15-S requires 16 cells per axis");
            for (int allowed : {16, 32, 64})
                if (value.to_float64() == allowed)
                    return allowed;
        }
        throw std::invalid_argument("E15 cells_per_axis must be 16, 32 or 64");
    }

    /// Publish complete JSON only after checking the stream, including its close.
    void writeJson(const fs::path& path, const Json& data) {
        const fs::path partial = path.string() + ".tmp";
        std::ofstream stream;
        stream.exceptions(std::ios::badbit | std::ios::failbit);
        stream.open(partial);
        stream << data << '\n';
        stream.close();
        fs::rename(partial, path);
    }

    /** @brief Store the native x-electric field on the central z layer for the illustration.
     * @details Only the slab owner writes. The matrix rows run along y, columns along x;
     * coordinates refer to the first sample and advance by the cell width.
     */
    void saveSlice(Solver& solver, const fs::path& output, Wave wave, long long step) {
        auto& field       = solver.electricField();
        const auto bounds = field.getLayout().getLocalNDIndex();
        const int layer   = wave.cells_m / 2;
        if (layer < bounds[2].first() || layer > bounds[2].last())
            return;
        const int g     = field.getNghost();
        const auto view = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), field.getView());
        const auto offset = solver.levelContract().positionOffset(Quantity::Electric, 0);
        Json values       = Json::array();
        for (int j = 0; j < wave.cells_m; ++j) {
            std::vector<double> row;
            for (int i = 0; i < wave.cells_m; ++i)
                row.push_back(view(i + g, j + g, layer - bounds[2].first() + g)[0]);
            values.push_back(std::move(row));
        }
        const double h = 1.0 / wave.cells_m;
        writeJson(output / ("slice-" + std::to_string(step) + ".json"),
                  {{"n", step},
                   {"t", solver.labelTime(Quantity::Electric, step)},
                   {"h", h},
                   {"x_first", (0.5 + offset[0]) * h},
                   {"y_first", (0.5 + offset[1]) * h},
                   {"z", (layer + 0.5 + offset[2]) * h},
                   {"Ex", values}});
    }

    /// Invert the solver's sample-time label, rounding ties to the later sample.
    int nearestStep(Solver& solver, double time) {
        return static_cast<int>(
            std::lround((time - solver.labelTime(Quantity::Electric, 0)) / solver.getDt()));
    }

    /// Observe all ranks at the same indices. Only rank zero retains the scalar records.
    Json evolve(Simulation& simulation, Wave wave, const input::ChdrYamlRunInput& request,
                int steps, double period, bool slices) {
        auto& solver = simulation.solver();
        int rank;
        MPI_Comm_rank(solver.A_n.getLayout().comm, &rank);
        const int stride = std::lround(period / wave.dt_m);
        std::array<int, 4> frames{0, nearestStep(solver, period / 4),
                                  nearestStep(solver, period / 2),
                                  nearestStep(solver, 3 * period / 4)};
        Json rows = Json::array();
        simulation.advance(steps, [&](long long step, Solver& state, const auto&) {
            if (slices && std::find(frames.begin(), frames.end(), step) != frames.end())
                saveSlice(state, request.outputDirectory_m, wave, step);
            if (step % stride == 0 || step == steps) {
                auto measurement = observe(state, wave, step);
                if (rank == 0)
                    rows.push_back(std::move(measurement));
            }
            if (rank == 0 && step > 0 && step % 10000 == 0)
                std::cout << request.caseName_m << " step " << step << '/' << steps << std::endl;
        });
        return rows;
    }

    /// Save the same E15.meas01–03 history and run record for either choice of time step.
    void writeResults(const input::ChdrYamlRunInput& request,
                      const simulation::ChdrCartesianGridSettings& grid, Solver& solver,
                      const Wave& wave, int steps, Json rows) {
        const MPI_Comm comm = solver.A_n.getLayout().comm;
        int rank;
        MPI_Comm_rank(comm, &rank);
        if (rank == 0)
            writeJson(fs::path(request.outputDirectory_m) / "results.json",
                      {{"N", wave.cells_m},
                       {"period", TwoPi / wave.omega_m},
                       {"kappa_E", wave.kappaE_m},
                       {"kappa_B", wave.kappaB_m},
                       {"phase_lag_rate", wave.omega_m - wave.discreteOmega_m},
                       {"observations", std::move(rows)}});
        auto record =
            output::describeRun("experiment-e15", grid, solver, steps, "periodic travelling wave");
        record.effective_m["material/relative_permittivity"] = Epsilon;
        record.effective_m["measurements/file"]              = "results.json";
        output::writeRunRecord(request.outputDirectory_m, request, record, comm);
    }

    /// Default to 1000 continuum periods for E15 and one million steps for E15-S.
    int runSteps(const input::ChdrYamlRunInput& request, Solver& solver, const Wave& wave) {
        const int refinementSteps = nearestStep(solver, 1000 * (TwoPi / wave.omega_m));
        const int steps =
            request.steps_m.value_or(request.caseName_m == "E15-S" ? 1000000 : refinementSteps);
        if (steps < 0 || steps > (wave.cells_m == 16 ? 1000000 : refinementSteps))
            throw std::invalid_argument("E15 step count outside the specified interval");
        return steps;
    }

    /** @brief Launch the prescribed diagonal mode on a periodic cubic mesh.
     * @details For E15-S only, select the homogeneous dielectric endpoint before loading
     * histories. The checked O1VS step uses the vacuum bound. See E15-S.typ and
     * Fallahi (2020), Eqs. (3.17) and (3.20), with the dielectric wave speed.
     */
    void run(const input::ChdrYamlRunInput& request) {
        const int n    = resolution(request);
        const double h = 1.0 / n;
        const simulation::ChdrCartesianGridSettings grid{
            {n, n, n}, {h, h, h}, 0, {true, true, true}, 2};
        const auto medium = KOKKOS_LAMBDA(material::CellIndex, material::Position, bool& marked) {
            marked = false;
            return Epsilon;
        };
        Simulation simulation(grid, medium, {.timeStep_m = h / 2});
        auto& solver = simulation.solver();
        if (request.caseName_m == "E15-S")
            solver.setTimeStepUnchecked(h * std::sqrt(Epsilon / 3.0));
        const Wave wave(n, solver.getDt());
        const double period = TwoPi / wave.omega_m;
        const int steps     = runSteps(request, solver, wave);
        const bool slices   = input::booleanParameter(request, "save_slices", false);
        output::prepareRunDirectory(request.outputDirectory_m, solver.A_n.getLayout().comm);
        simulation.load(wave);
        writeResults(request, grid, solver, wave, steps,
                     evolve(simulation, wave, request, steps, period, slices));
    }
}  // namespace chdr::experiments::e15

int main(int argc, char** argv) {
    return chdr::input::runApplication(argc, argv, chdr::experiments::e15::run);
}
