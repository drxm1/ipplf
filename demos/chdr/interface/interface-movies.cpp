/** @file interface-movies.cpp
 * @brief Open-box TE annulus movies using the shared ChDR simulation manager.
 */
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <type_traits>

#include "../output/FieldRunOutput.h"
#include "../simulation/ChDRSimulationManager.h"
#include "MovieCases.h"
#include "MovieFaceSolver.h"
#include "TMAnnulusStart.h"

namespace chdr::interface {
    /// Movie selection and output location; the optional step count permits extending a run.
    struct MovieOptions {
        std::string name_m, directory_m;
        bool smoke_m = false;
        int steps_m  = -1;
    };

    /// Require a complete, nonnegative integer instead of accepting a numeric prefix.
    int parseSteps(const std::string& value) {
        std::size_t used = 0;
        const int result = std::stoi(value, &used);
        if (used != value.size() || result < 0)
            throw std::invalid_argument("--steps requires a nonnegative integer");
        return result;
    }

    /// Parse the few driver options; material and pulse parameters live in the case table.
    MovieOptions movieOptions(int argc, char** argv) {
        if (argc < 2)
            throw std::invalid_argument("Use CASE --output NEW_DIRECTORY [--smoke] [--steps N]");
        MovieOptions result{argv[1], {}};
        for (int i = 2; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--smoke")
                result.smoke_m = true;
            else if (arg == "--output" && i + 1 < argc)
                result.directory_m = argv[++i];
            else if (arg == "--steps" && i + 1 < argc)
                result.steps_m = parseSteps(argv[++i]);
            else
                throw std::invalid_argument(
                    "Use CASE --output NEW_DIRECTORY [--smoke] [--steps N]");
        }
        if (result.directory_m.empty())
            throw std::invalid_argument("A new output directory is required");
        return result;
    }

    /// Use the same face speeds for solver configuration and the recorded experiment input.
    ippl::Vector<double, 6> movieFaceSpeeds(const MovieCase& input) {
        const double speed = input.diagonal_m ? 1 / std::sqrt(input.epsilon_m) : 1;
        return {1, speed, 1, 1, speed, 1};
    }

    /// Record mesh and boundary values from the case rather than assuming the default box.
    std::string movieBoxMetadata(const MovieCase& input) {
        std::ostringstream text;
        text << std::setprecision(17) << " plane=zx cells=";
        for (unsigned d = 0; d < 3; ++d)
            text << (d ? "," : "") << input.box_m.cells_m[d];
        text << " h=";
        for (unsigned d = 0; d < 3; ++d)
            text << (d ? "," : "") << input.box_m.spacing_m[d];
        text << " periodic=";
        for (unsigned d = 0; d < 3; ++d)
            text << (d ? "," : "") << input.box_m.periodic_m[d];
        text << " split_axis=" << input.box_m.splitAxis_m << " face_speeds=";
        const auto speeds = movieFaceSpeeds(input);
        for (unsigned f = 0; f < 6; ++f)
            text << (f ? "," : "") << speeds[f];
        return text.str() + " edge_speed=1";
    }

    /// Preserve the case geometry, pulse scales and sampling cadence with the output data.
    std::string movieMetadata(const MovieCase& input, const simulation::SolverSettings& settings) {
        const auto& p = input.start_m;
        const auto& r = input.region_m;
        std::ostringstream text;
        text << std::setprecision(17) << "case=" << input.name_m << movieBoxMetadata(input)
             << " steps=" << input.steps_m << " every=" << input.frameEvery_m
             << " polarization=" << (input.tm_m ? "TM" : "TE") << " S=" << settings.safetyFactor_m
             << " epsilon=" << input.epsilon_m << " centre=" << p.centreX_m << ',' << p.centreZ_m
             << " radius=" << p.radius_m << " wavelength=" << p.wavelength_m
             << " support_half_width=" << 1.5 * p.wavelength_m << " centre_guard=" << p.guard_m
             << " amplitude=" << p.amplitude_m;
        if (input.diagonal_m)
            text << " material=x>=z";
        else
            text << " rectangle=" << r.lowerX_m << ',' << r.upperX_m << ',' << r.lowerZ_m << ','
                 << r.upperZ_m;
        return text.str();
    }

    /// Load and record a case with its chosen solver and device-copyable material predicate.
    template <class Solver, class Region>
    void runMovieCase(const MovieOptions& options, const MovieCase& input, Region region) {
        const simulation::SolverSettings settings;
        simulation::ChDRSimulationManager<Solver> run(
            input.box_m, material::TwoMaterialMap{region, input.epsilon_m, 1.0}, settings);
        if constexpr (std::is_same_v<Solver, MovieFaceSolver>)
            run.solver().setMovieFaceSpeeds(movieFaceSpeeds(input));
        checkMovieStart(input, run.solver().labelTime(solver::Quantity::VectorPotential, 0));
        if (input.tm_m)
            run.load(movieTMStart(input));
        else
            run.load(input.start_m);
        chdr::output::FieldRunOutput output(options.directory_m, input.box_m, run.solver(),
                                            run.material(), movieMetadata(input, settings));
        run.advance(input.steps_m, [&](long long k, auto& state, const auto& owner) {
            if (k % input.frameEvery_m == 0 || k == input.steps_m)
                output.write(k, state, owner);
        });
        output.finish();
    }

    /// Keep the experimental face override confined to the diagonal demonstration.
    void runMovie(const MovieOptions& options) {
        auto input = movieCase(options.name_m, options.smoke_m);
        if (options.steps_m >= 0)
            input.steps_m = options.steps_m;
        if (input.diagonal_m)
            runMovieCase<MovieFaceSolver>(options, input, DiagonalXZ{});
        else
            runMovieCase<solver::ChdrSolver<solver::StandardSolverPeriodic>>(options, input,
                                                                             input.region_m);
    }
}  // namespace chdr::interface

int main(int argc, char** argv) {
    ippl::initialize(argc, argv);
    int result = 0;
    try {
        chdr::interface::runMovie(chdr::interface::movieOptions(argc, argv));
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        result = 1;
    }
    ippl::finalize();
    return result;
}
