/** @file interface-movies.cpp
 * @brief Open-box annulus cases using the common input, observation and recording interfaces.
 */
#include <iomanip>
#include <sstream>
#include <string>
#include <type_traits>

#include "apps/interface-movies/MovieCases.h"
#include "apps/interface-movies/MovieFaceSolver.h"
#include "core/input/Application.h"
#include "core/input/RunInput.h"
#include "core/output/RunRecord.h"
#include "core/simulation/ChDRSimulationManager.h"
#include "opt-in/output/FieldRunOutput.h"

namespace chdr::apps::interfaceMovies {
    /// Use the same face speeds for solver configuration and recorded case settings.
    ippl::Vector<double, 6> movieFaceSpeeds(const MovieCase& selected) {
        const double speed = selected.diagonal_m ? 1 / std::sqrt(selected.epsilon_m) : 1;
        return {1, speed, 1, 1, speed, 1};
    }

    /// Retain the field renderer's mesh, periodicity and boundary-speed metadata.
    std::string movieBoxMetadata(const MovieCase& selected) {
        std::ostringstream text;
        text << std::setprecision(17) << " plane=zx cells=";
        for (unsigned d = 0; d < 3; ++d)
            text << (d ? "," : "") << selected.box_m.cells_m[d];
        text << " h=";
        for (unsigned d = 0; d < 3; ++d)
            text << (d ? "," : "") << selected.box_m.spacing_m[d];
        text << " periodic=";
        for (unsigned d = 0; d < 3; ++d)
            text << (d ? "," : "") << selected.box_m.periodic_m[d];
        text << " split_axis=" << selected.box_m.splitAxis_m << " face_speeds=";
        const auto speeds = movieFaceSpeeds(selected);
        for (unsigned f = 0; f < 6; ++f)
            text << (f ? "," : "") << speeds[f];
        return text.str() + " edge_speed=1";
    }

    /// Record the effective initial-data parameters in the field metadata.
    std::string moviePulseMetadata(const MovieCase& selected) {
        const auto& p = selected.start_m;
        std::ostringstream text;
        text << std::setprecision(17) << " polarization=" << (selected.tm_m ? "TM" : "TE")
             << " epsilon=" << selected.epsilon_m << " centre=" << p.centreX_m << ',' << p.centreZ_m
             << " radius=" << p.radius_m << " wavelength=" << p.wavelength_m
             << " support_half_width=" << 1.5 * p.wavelength_m << " centre_guard=" << p.guard_m
             << " amplitude=" << p.amplitude_m;
        if (selected.diagonal_m)
            text << " material=x>=z";
        else
            text << " rectangle=" << selected.region_m.lowerX_m << ',' << selected.region_m.upperX_m
                 << ',' << selected.region_m.lowerZ_m << ',' << selected.region_m.upperZ_m;
        return text.str();
    }

    /// Keep the existing CSV metadata while the common run record stores structured settings.
    std::string movieMetadata(const MovieCase& selected, const simulation::SolverSettings& settings,
                              const output::FieldRunSettings& observations) {
        std::ostringstream text;
        text << std::setprecision(17) << "case=" << selected.name_m << movieBoxMetadata(selected)
             << " steps=" << selected.steps_m << " every=" << observations.fields_m.every_m
             << " S=" << settings.safetyFactor_m << moviePulseMetadata(selected);
        return text.str();
    }

    /// Preserve the solver-specific face override and the existing TE/TM history construction.
    template <class Solver>
    void loadMovie(simulation::ChDRSimulationManager<Solver>& run, const MovieCase& selected) {
        if constexpr (std::is_same_v<Solver, MovieFaceSolver>)
            run.solver().setMovieFaceSpeeds(movieFaceSpeeds(selected));
        checkMovieStart(selected, run.solver().labelTime(Quantity::VectorPotential, 0));
        if (selected.tm_m)
            run.load(movieTMStart(selected));
        else
            run.load(selected.start_m);
    }

    /// Save the completed solver state and the same schedules used by the output callback.
    template <class Solver>
    void recordMovie(const input::RunInput& request, const MovieCase& selected,
                     simulation::ChDRSimulationManager<Solver>& run, const std::string& metadata,
                     const output::FieldRunSettings& observations) {
        auto description = output::describeRun("interface-movies", selected.box_m, run.solver(),
                                               selected.steps_m, metadata);
        output::recordFieldRunOutput(description, request, observations);
        output::writeRunRecord(request.outputDirectory_m, request, description,
                               run.solver().A_n.getLayout().comm);
    }

    /// Execute one existing case through the shared lifecycle and write its completed run record.
    template <class Solver, class Region>
    void runMovieCase(const input::RunInput& request, MovieCase selected, Region region) {
        const simulation::SolverSettings settings;
        simulation::ChDRSimulationManager<Solver> run(
            selected.box_m, material::TwoMaterialMap{region, selected.epsilon_m, 1.0}, settings);
        selected.steps_m = input::resolvedSteps(request, run.solver().getDt(), selected.steps_m);
        loadMovie(run, selected);
        const auto observations = output::fieldRunSettings(request, selected.frameEvery_m);
        const auto metadata     = movieMetadata(selected, settings, observations);
        output::FieldRunOutput writer(request.outputDirectory_m, selected.box_m, run.solver(),
                                      run.material(), metadata, request.units_m);
        run.advance(selected.steps_m, [&](long long k, auto& state, const auto& owner) {
            writer.observe(k, state, owner, selected.steps_m, observations);
        });
        writer.finish();
        recordMovie(request, selected, run, metadata, observations);
    }

    /// Resolve the selected preset; YAML and compatibility arguments enter the same run path.
    void runConfigured(const input::RunInput& request) {
        input::validateParameters(request, {"smoke"});
        input::validateObservationNames(request, {"fields", "diagnostics"});
        if (request.outputDirectory_m.empty())
            throw std::invalid_argument("A new output directory is required");
        auto selected =
            movieCase(request.caseName_m, input::booleanParameter(request, "smoke", false));
        if (selected.diagonal_m)
            runMovieCase<MovieFaceSolver>(request, selected, DiagonalXZ{});
        else
            runMovieCase<solver::ChdrSolver<solver::StandardSolverPeriodic>>(request, selected,
                                                                             selected.region_m);
    }
}  // namespace chdr::apps::interfaceMovies

int main(int argc, char** argv) {
    return chdr::input::runApplication(argc, argv, chdr::apps::interfaceMovies::runConfigured);
}
