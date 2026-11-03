/** @file test_movie_start.cpp
 * @brief Annulus placement, polarization and finite initial energy in an open material box.
 * @details The tests check the starts of the eight movie cases of interface-movies
 * (MovieCases.h) and the per-face Mur speeds of MovieFaceSolver:
 * - StartPlacement: all eight starts, full and smoke size, pass checkMovieStart(); a start
 *   centre moved onto the dielectric, or near the diagonal staircase at (22.1, 41.9),
 *   throws std::invalid_argument.
 * - UnitFaceSpeeds: after 48 steps, unit face speeds give the potential \f$A_n\f$ exactly
 *   equal to that of the library manager (maximum difference zero); the speeds
 *   (0.7, 0.7, 1, 1, 0.7, 0.7) differ by more than \f$10^{-12}\f$.
 * - InitialFields: the smoke-size half-space TE start has zero potential components 0, 1
 *   and 3 and zero \f$A_y\f$ inside the dielectric, in both histories; its aligned energy
 *   is at least \f$1 - 10^{-12}\f$ times the discrete vector invariant.
 * - TMInitialConstraint: both histories of the smoke-size half-space TM start have a
 *   weighted divergence below \f$10^{-13}\f$ in every cell.
 * Reference: self-consistency, and exact numerical equality with the library manager.
 *
 * CTest: chdr.interface.movies.unit.serial and chdr.interface.movies.unit.mpi2; the
 * two-rank entry does not compare with the one-rank result.
 */
#include <gtest/gtest.h>

#include "apps/interface-movies/MovieCases.h"
#include "apps/interface-movies/MovieFaceSolver.h"
#include "core/simulation/ChDRSimulationManager.h"
#include "opt-in/diagnostics/AlignedFieldEnergy.h"
#include "opt-in/diagnostics/FieldEnergy.h"
#include "tests/core/solver/CoreTestSupport.h"

namespace {
    using namespace chdr;

    /// Detect unwanted potential components or any initial potential inside the dielectric.
    double unexpectedPotential(const AJFourFieldsType& field,
                               const material::MaterialOwner& owner) {
        const auto a       = field.getView();
        const auto epsilon = owner.epsilon();
        return solver::coreTest::globalMaximum(
            field, solver::coreTest::localMaximum(
                       field, KOKKOS_LAMBDA(int i, int j, int k) {
                           const auto value = a(i, j, k);
                           return Kokkos::abs(value[0]) + Kokkos::abs(value[1])
                                  + Kokkos::abs(value[3])
                                  + (epsilon(i, j, k) == 1 ? 0 : Kokkos::abs(value[2]));
                       }));
    }

    /// Eight cases at full and smoke size: the start passes checkMovieStart() and is zero at
    /// its centre; a centre moved onto the dielectric face throws std::invalid_argument.
    TEST(Movies, StartPlacement) {
        for (const auto name : {"half-space", "brick", "slab", "diagonal", "half-space-tm",
                                "brick-tm", "slab-tm", "diagonal-tm"}) {
            for (bool smoke : {false, true}) {
                auto input = apps::interfaceMovies::movieCase(name, smoke);
                EXPECT_NO_THROW(apps::interfaceMovies::checkMovieStart(input, 0));
                EXPECT_EQ(input.start_m(input.start_m.centreX_m, 0, input.start_m.centreZ_m, 0, 2),
                          0);
                input.start_m.centreX_m =
                    input.diagonal_m ? input.start_m.centreZ_m : input.region_m.lowerX_m;
                EXPECT_THROW(apps::interfaceMovies::checkMovieStart(input, 0),
                             std::invalid_argument);
            }
        }
        auto nearStaircase              = apps::interfaceMovies::movieCase("diagonal-tm", true);
        nearStaircase.start_m.centreX_m = 22.1;
        nearStaircase.start_m.centreZ_m = 41.9;
        EXPECT_THROW(apps::interfaceMovies::checkMovieStart(nearStaircase, 0),
                     std::invalid_argument);
    }

    /// Compare owned potential components across all participating ranks.
    double potentialDifference(const AJFourFieldsType& field, const AJFourFieldsType& other) {
        const auto a = field.getView(), b = other.getView();
        const double difference = solver::coreTest::localMaximum(
            field, KOKKOS_LAMBDA(int i, int j, int k) {
                double error = 0;
                for (unsigned q = 0; q < 4; ++q)
                    error = Kokkos::max(error, Kokkos::abs(a(i, j, k)[q] - b(i, j, k)[q]));
                return error;
            });
        return solver::coreTest::globalMaximum(field, difference);
    }

    /// Unit speeds reproduce the library path; non-unit speeds prove the override is reached.
    /// Smoke-size half-space case, 48 steps: difference exactly 0, and above 1e-12.
    TEST(Movies, UnitFaceSpeeds) {
        const auto input = apps::interfaceMovies::halfSpaceMovie(true);
        const auto map   = material::TwoMaterialMap{input.region_m, input.epsilon_m, 1.0};
        simulation::ChDRSimulationManager<> reference(input.box_m, map);
        simulation::ChDRSimulationManager<apps::interfaceMovies::MovieFaceSolver> candidate(
            input.box_m, map),
            control(input.box_m, map);
        control.solver().setMovieFaceSpeeds({0.7, 0.7, 1, 1, 0.7, 0.7});
        reference.load(input.start_m);
        candidate.load(input.start_m);
        control.load(input.start_m);
        reference.run(48);
        candidate.run(48);
        control.run(48);
        EXPECT_EQ(potentialDifference(reference.solver().A_n, candidate.solver().A_n), 0);
        EXPECT_GT(potentialDifference(reference.solver().A_n, control.solver().A_n), 1e-12);
    }

    /** @brief Load both histories and compare the two energy definitions for this TE state.
     * The aligned magnetic square exceeds the cross-level product by a nonnegative square
     * (deduction from [sakkaplangkul2021YeeDispersiveMedia, section 4.1, Eq. (56)]).
     * This checks initialization and diagnostic definitions, not absorption accuracy.
     */
    TEST(Movies, InitialFields) {
        const auto input = apps::interfaceMovies::halfSpaceMovie(true);
        simulation::ChDRSimulationManager<> run(
            input.box_m, material::TwoMaterialMap{input.region_m, input.epsilon_m, 1.0});
        run.load(input.start_m);
        EXPECT_EQ(unexpectedPotential(run.solver().A_n, run.material()), 0);
        EXPECT_EQ(unexpectedPotential(run.solver().A_nm1, run.material()), 0);
        const double aligned = diagnostics::alignedFieldEnergy(run.solver(), run.material());
        const double discrete =
            diagnostics::potentialEnergy(run.solver(), run.material()).vector_m.total();
        EXPECT_TRUE(std::isfinite(aligned));
        EXPECT_GT(discrete, 0);
        EXPECT_GE(aligned, discrete * (1 - 1e-12));
    }

    /// Both TM histories satisfy the weighted divergence constraint before boundary interaction:
    /// smoke-size half-space-tm case, maximum residual below 1e-13, aligned energy positive.
    TEST(Movies, TMInitialConstraint) {
        const auto input = apps::interfaceMovies::movieCase("half-space-tm", true);
        simulation::ChDRSimulationManager<> run(
            input.box_m, material::TwoMaterialMap{input.region_m, input.epsilon_m, 1.0});
        run.load(apps::interfaceMovies::movieTMStart(input));
        for (const auto* field : {&run.solver().A_n, &run.solver().A_nm1}) {
            const solver::O1VSOperator op(*field, run.material(), 1);
            const double residual = solver::coreTest::localMaximum(
                *field, KOKKOS_LAMBDA(int i, int j, int k) { return Kokkos::abs(op.g(i, j, k)); });
            EXPECT_LT(solver::coreTest::globalMaximum(*field, residual), 1e-13);
        }
        EXPECT_GT(diagnostics::alignedFieldEnergy(run.solver(), run.material()), 0);
    }

}  // namespace
