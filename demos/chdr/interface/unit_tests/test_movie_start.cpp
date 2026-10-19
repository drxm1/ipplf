/** @file test_movie_start.cpp
 * @brief Annulus placement, polarization and finite initial energy in an open material box.
 */
#include <gtest/gtest.h>

#include "../../simulation/ChDRSimulationManager.h"
#include "../MovieCases.h"
#include "../MovieFaceSolver.h"
#include "AlignedFieldEnergy.h"
#include "CoreTestSupport.h"
#include "FieldEnergy.h"

namespace {
    using namespace chdr;

    /// Detect unwanted potential components or any initial potential inside the dielectric.
    double unexpectedPotential(const solver::AJFourFieldsType& field,
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

    TEST(Movies, StartPlacement) {
        for (const auto name : {"half-space", "brick", "slab", "diagonal"}) {
            for (bool smoke : {false, true}) {
                auto input = interface::movieCase(name, smoke);
                EXPECT_NO_THROW(interface::checkMovieStart(input, 0));
                EXPECT_EQ(input.start_m(input.start_m.centreX_m, 0, input.start_m.centreZ_m, 0, 2),
                          0);
                input.start_m.centreX_m =
                    input.diagonal_m ? input.start_m.centreZ_m : input.region_m.lowerX_m;
                EXPECT_THROW(interface::checkMovieStart(input, 0), std::invalid_argument);
            }
        }
    }

    /// Unit face speeds must reproduce the existing boundary path exactly, including edges.
    TEST(Movies, UnitFaceSpeeds) {
        const auto input = interface::halfSpaceMovie(true);
        const auto map   = material::TwoMaterialMap{input.region_m, input.epsilon_m, 1.0};
        simulation::ChDRSimulationManager<> reference(input.box_m, map);
        simulation::ChDRSimulationManager<interface::MovieFaceSolver> candidate(input.box_m, map);
        reference.load(input.start_m);
        candidate.load(input.start_m);
        reference.run(48);
        candidate.run(48);
        const auto a            = reference.solver().A_n.getView();
        const auto b            = candidate.solver().A_n.getView();
        const auto& field       = reference.solver().A_n;
        const double difference = solver::coreTest::localMaximum(
            field, KOKKOS_LAMBDA(int i, int j, int k) {
                double error = 0;
                for (unsigned q = 0; q < 4; ++q)
                    error = Kokkos::max(error, Kokkos::abs(a(i, j, k)[q] - b(i, j, k)[q]));
                return error;
            });
        EXPECT_EQ(solver::coreTest::globalMaximum(field, difference), 0);
    }

    /** @brief Load both histories and compare the two energy definitions for this TE state.
     * The aligned magnetic square exceeds the cross-level product by a nonnegative square
     * (deduction from [sakkaplangkul2021YeeDispersiveMedia, section 4.1, Eq. (56)]).
     * This checks initialization and diagnostic definitions, not absorption accuracy.
     */
    TEST(Movies, InitialFields) {
        const auto input = interface::halfSpaceMovie(true);
        simulation::ChDRSimulationManager<> run(
            input.box_m, material::TwoMaterialMap{input.region_m, input.epsilon_m, 1.0});
        run.load(input.start_m);
        EXPECT_EQ(unexpectedPotential(run.solver().A_n, run.material()), 0);
        EXPECT_EQ(unexpectedPotential(run.solver().A_nm1, run.material()), 0);
        const double aligned = solver::alignedFieldEnergy(run.solver(), run.material());
        const double discrete =
            solver::potentialEnergy(run.solver(), run.material()).vector_m.total();
        EXPECT_TRUE(std::isfinite(aligned));
        EXPECT_GT(discrete, 0);
        EXPECT_GE(aligned, discrete * (1 - 1e-12));
    }
}  // namespace
