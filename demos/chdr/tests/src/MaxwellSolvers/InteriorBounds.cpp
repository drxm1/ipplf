/** @file InteriorBounds.cpp
 * @brief Cells written by the interior stencils of StandardFDTDSolver and NonStandardFDTDSolver,
 * and the constructor policy of the boundary modes periodic, absorbing and configured.
 * @details The box has 6 by 5 by 4 cells of spacing (1, 1.25, 0.5) and splits along \f$z\f$.
 * With unit source, zero histories and no boundary update, one step must write exactly
 * \f$\Delta t^2\f$ into every component of the cells with oracle mask 0 (murFaceMask() in
 * BoundaryOracles.h) and leave the sentinel \f$-17\f$ in all other owned cells, for all 8 axis
 * patterns. Two static_assert checks fix the constructor overloads. The other tests check that
 * the solver keeps its copy of the description through initialize(), that the fixed modes
 * keep their all-face meaning, and that invalid descriptions and two-ghost storage throw
 * IpplException before any halo exchange.
 *
 * Reference: exact lattice value; refusal. CTest entries: chdr.src.InteriorBounds.serial and
 * chdr.src.InteriorBounds.mpi2. Only the mpi2 entry checks the refusal of a 6 by 5 by 3 box
 * whose \f$z\f$ split leaves one rank one cell.
 */
#include <gtest/gtest.h>
#include <type_traits>

#include "BoundaryOracles.h"
#include "HaloOrderTest.h"

namespace {
    /** @brief Enforce the boundary policy at overload resolution, not inside constructor bodies.
     * Both fixed modes accept only the three fields. configured requires the
     * fourth argument, including when its description happens to be all periodic.
     */
    template <template <ippl::fdtd_bc> class Solver>
    constexpr bool constructorContracts() {
        using P = fdtd_test::Potential&;
        using F = fdtd_test::Field&;
        using B = const ippl::FDTDBoundary&;
        return std::is_constructible_v<Solver<ippl::periodic>, P, F, F>
               && std::is_constructible_v<Solver<ippl::absorbing>, P, F, F>
               && !std::is_constructible_v<Solver<ippl::configured>, P, F, F>
               && !std::is_constructible_v<Solver<ippl::periodic>, P, F, F, B>
               && !std::is_constructible_v<Solver<ippl::absorbing>, P, F, F, B>
               && std::is_constructible_v<Solver<ippl::configured>, P, F, F, B>;
    }
    static_assert(constructorContracts<fdtd_test::Standard>());
    static_assert(constructorContracts<fdtd_test::NonStandard>());

    /** @brief Expose the interior write set without a subsequent boundary update.
     * A zero history with unit source writes exactly dt squared per component.
     * Untouched samples keep a negative sentinel. This checks which cells the
     * actual stencil writes independently of InteriorBounds::contains itself.
     */
    template <class Base>
    struct InteriorOnly : Base {
        using Base::Base;
        void finishPotentialStep() override {}
    };

    template <class Base>
    void checkInterior(unsigned pattern) {
        const std::array<bool, 3> periodic{bool(pattern & 1), bool(pattern & 2), bool(pattern & 4)};
        fdtd_test::Fields fields({6, 5, 4}, {1, 1.25, 0.5}, periodic);
        InteriorOnly<Base> solver(fields.source, fields.electric, fields.magnetic,
                                  ippl::FDTDBoundary::fromLayout(fields.layout));
        fields.source = 1;
        solver.A_np1  = -17;
        solver.step();
        fdtd_test::owned(solver.A_np1, [&](int i, int j, int k, auto value) {
            const unsigned mask = open_boundary_oracle::murFaceMask(
                {i + 1, j + 1, k + 1}, {6, 5, 4}, {periodic[0], periodic[1], periodic[2]});
            const double expected = mask == 0 ? solver.getDt() * solver.getDt() : -17;
            for (unsigned q = 0; q < 4; ++q) {
                EXPECT_EQ(value[q], expected);
            }
        });
    }
}  // namespace

/// For all 8 axis patterns and both stencils, one step writes \f$\Delta t^2\f$ exactly into the
/// cells with oracle mask 0 and leaves the sentinel -17 in every face, edge and corner cell.
TEST(InteriorBounds, BothStencilsWriteOnlyInterior) {
    for (unsigned pattern = 0; pattern < 8; ++pattern) {
        SCOPED_TRACE(pattern);
        checkInterior<fdtd_test::Standard<ippl::configured>>(pattern);
        checkInterior<fdtd_test::NonStandard<ippl::configured>>(pattern);
    }
}

/** @brief Construction copies the description; resetting potentials preserves that copy.
 * Changing the caller's copy must not change either the kinds or non-default speeds.
 * Every stored potential must recover its field conditions after initialize().
 */
template <class Solver>
void checkBoundaryPersistence() {
    fdtd_test::Fields fields({6, 5, 4}, {1, 1.25, 0.5}, std::array<bool, 3>{false, true, true});
    auto description    = ippl::FDTDBoundary::fromLayout(fields.layout, 0.8);
    const auto expected = description;
    Solver solver(fields.source, fields.electric, fields.magnetic, description);
    description.kind.fill(ippl::FaceKind::periodic);
    description.speed.fill(1);
    solver.initialize();
    EXPECT_EQ(solver.boundary().kind, expected.kind);
    EXPECT_EQ(solver.boundary().speed, expected.speed);
    for (auto* field : {&solver.A_n, &solver.A_nm1, &solver.A_np1}) {
        EXPECT_EQ(field->getFieldBC()[0]->getBCType(), ippl::NO_FACE);
        EXPECT_EQ(field->getFieldBC()[2]->getBCType(), ippl::PERIODIC_FACE);
    }
}

/// Both stencils, open x and Mur speed 0.8: kinds and speeds survive initialize(), and every
/// potential level keeps NO_FACE on the lower x face and PERIODIC_FACE on the lower y face.
TEST(InteriorBounds, BoundarySelectionSurvivesInitialize) {
    checkBoundaryPersistence<fdtd_test::Standard<ippl::configured>>();
    checkBoundaryPersistence<fdtd_test::NonStandard<ippl::configured>>();
}

/** @brief Fixed modes retain their all-face meaning; absorbing never infers periodic faces.
 * The periodic-on-nonperiodic-layout path is an existing supported field-BC use.
 * Any periodically wrapping axis conflicts with the fixed all-Mur description.
 */
template <template <ippl::fdtd_bc> class Solver>
void checkFixedModes() {
    fdtd_test::Fields open({6, 5, 4}, {1, 1.25, 0.5}, false);
    Solver<ippl::periodic> periodic(open.source, open.electric, open.magnetic);
    EXPECT_EQ(periodic.boundary().murMask(), 0u);
    periodic.prepareHistories();
    Solver<ippl::absorbing> absorbing(open.source, open.electric, open.magnetic);
    EXPECT_EQ(absorbing.boundary().murMask(), 63u);
    for (const auto axes :
         {std::array<bool, 3>{true, true, true}, std::array<bool, 3>{false, true, true}}) {
        fdtd_test::Fields wrapped({6, 5, 4}, {1, 1.25, 0.5}, axes);
        EXPECT_THROW((Solver<ippl::absorbing>(wrapped.source, wrapped.electric, wrapped.magnetic)),
                     IpplException);
    }
}

/// Both stencils: periodic has Mur mask 0 and absorbing has mask 63 on an open layout, and
/// absorbing throws IpplException on a layout with a periodic axis.
TEST(InteriorBounds, FixedModesDoNotDeriveFaceKinds) {
    checkFixedModes<fdtd_test::Standard>();
    checkFixedModes<fdtd_test::NonStandard>();
}

/** Solver entry points must reject invalid storage on every rank before halo communication.
 * The three-cell z split is thin on only one rank. Inconsistent descriptions
 * fail during construction; malformed histories must not enter fillHalo.
 * The thin split runs only on two ranks.
 */
TEST(InteriorBounds, SolverRefusesBeforeCommunication) {
    using Solver = fdtd_test::Standard<ippl::absorbing>;
    if (ippl::Comm->size() == 2) {
        fdtd_test::Fields thin({6, 5, 3}, {1, 1, 1}, false);
        EXPECT_THROW((Solver(thin.source, thin.electric, thin.magnetic)), IpplException);
    }
    fdtd_test::Fields fields({6, 5, 4}, {1, 1, 1}, false);
    Solver solver(fields.source, fields.electric, fields.magnetic);
    using Configured = fdtd_test::Standard<ippl::configured>;
    auto invalid     = solver.boundary();
    invalid.kind[0]  = ippl::FaceKind::periodic;
    EXPECT_THROW((Configured(fields.source, fields.electric, fields.magnetic, invalid)),
                 IpplException);
    solver.A_n.updateLayout(fields.layout, 2);
    EXPECT_THROW(solver.prepareHistories(), IpplException);
    fields.electric.updateLayout(fields.layout, 2);
    EXPECT_THROW((Solver(fields.source, fields.electric, fields.magnetic)), IpplException);
}

/** @brief Explicit all-periodic descriptions still receive strict storage validation.
 * Unlike the fixed periodic compatibility mode, configured cannot bypass
 * validation merely because its Mur mask is zero.
 */
TEST(InteriorBounds, ConfiguredPeriodicValidatesStorage) {
    using Solver = fdtd_test::Standard<ippl::configured>;
    fdtd_test::Fields fields({6, 5, 4}, {1, 1.25, 0.5}, true);
    const auto description = ippl::FDTDBoundary::fromLayout(fields.layout);
    Solver solver(fields.source, fields.electric, fields.magnetic, description);
    solver.prepareHistories();
    solver.A_n.updateLayout(fields.layout, 2);
    EXPECT_THROW(solver.prepareHistories(), IpplException);
    fields.magnetic.updateLayout(fields.layout, 2);
    EXPECT_THROW((Solver(fields.source, fields.electric, fields.magnetic, description)),
                 IpplException);
}

int main(int argc, char* argv[]) {
    ippl::initialize(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    ippl::finalize();
    return result;
}
