/** @file test_setup.cpp
 * @brief Independent initialization, local-domain and debug halo-guard tests.
 * @details The sine/history reference is our substitution in Fallahi, MITHRA 2.0,
 * arXiv:2009.13645v1, Eqs. (3.6)-(3.15), printed pp.12-13. The shifted grid is
 * a coordinate/timestep fixture, not a new propagation benchmark. Required
 * stencil faces follow StandardFDTDSolver::step (StandardFDTDSolver.hpp:65-69) and
 * FDTDSolverBase::evaluate_EB (FDTDSolverBase.hpp:113-118).
 * Domain rejection follows HaloCells::exchangeBoundaries (HaloCells.hpp:122-128).
 * Debug guards are tested only
 * when the same production VacuumSetup.cpp is compiled without NDEBUG.
 */
#include "Ippl.h"

#include <Kokkos_MathematicalConstants.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <memory>
#include <stdexcept>

#include "VacuumSetup.h"

namespace {
    using namespace chdr::vacuum;
    using PotentialField        = SourceField_t<double, Dim>;
    constexpr double WaveNumber = 2.0 * Kokkos::numbers::pi_v<double>;

    /// @brief Construct a periodic test layout with caller-selected cells and split axes.
    Fields unitSpacedGrid(const std::array<int, Dim>& cells, std::array<bool, Dim> decompose) {
        ippl::Vector<double, Dim> lower(0.0), spacing(1.0), upper;
        ippl::NDIndex<Dim> domain;
        for (unsigned d = 0; d < Dim; ++d) {
            domain[d] = ippl::Index(cells[d]);
            upper[d]  = cells[d];
        }
        return Fields(spacing, lower, upper, decompose, domain, lower, true);
    }

    /// @brief Demonstrate that domain validation does not assume the demo grid or z split.
    TEST(VacuumDomain, AcceptsDifferentSizesAndSplitDirections) {
        auto fields = unitSpacedGrid({4 * ippl::Comm->size(), 3, 5}, {true, false, false});
        EXPECT_NO_THROW(checkLocalDomain(fields.getFL()));
    }

    /// @brief Reject each possible thin axis, including thin local z partitions.
    TEST(VacuumDomain, RejectsOneOwnedCellOnEveryAxis) {
        for (unsigned axis = 0; axis < Dim; ++axis) {
            std::array<int, Dim> cells{4, 4, 2 * ippl::Comm->size()};
            cells[axis] = axis == 2 ? ippl::Comm->size() : 1;
            auto fields = unitSpacedGrid(cells, {false, false, true});
            EXPECT_THROW(checkLocalDomain(fields.getFL()), std::runtime_error);
        }
    }

    /// @brief Every rank rejects when at least one partition is too thin; MPI2 mixes validity.
    TEST(VacuumDomain, RejectsWhenOnlySomeRanksAreInvalid) {
        auto fields = unitSpacedGrid({4, 4, 2 * ippl::Comm->size() - 1}, {false, false, true});
        EXPECT_THROW(checkLocalDomain(fields.getFL()), std::runtime_error);
    }

    /** @brief Shift coordinates and choose hx<hz to catch a hard-coded axial timestep.
     * @details StandardFDTDSolver::initialize (StandardFDTDSolver.hpp:91-95) uses
     * min(spacing)/2. This fixture retains
     * one z wavelength and tests coordinate/time initialization, not propagation accuracy.
     */
    Fields shiftedGrid() {
        ippl::Vector<double, Dim> lower(-0.25, 0.5, 0.375), upper,
            spacing(1.0 / 256.0, 1.0 / 128.0, 1.0 / 64.0);
        const std::array<int, Dim> cells{8, 4, 64};
        ippl::NDIndex<Dim> domain;
        for (unsigned d = 0; d < Dim; ++d) {
            domain[d] = ippl::Index(cells[d]);
            upper[d]  = lower[d] + cells[d] * spacing[d];
        }
        return Fields(spacing, lower, upper, {false, false, true}, domain, lower, true);
    }

    /// @brief Select owned cells or one face halo; exclude halo edges and corners.
    bool faceOrOwned(const ippl::NDIndex<Dim>& local, int i, int j, int k) {
        return int(i < local[0].first() || i > local[0].last())
                   + int(j < local[1].first() || j > local[1].last())
                   + int(k < local[2].first() || k > local[2].last())
               <= 1;
    }

    /** @brief Independent periodic potential at a global z index and specified time.
     * @pre The global index domain starts at zero, as in both initialization fixtures.
     * @details Wrap the index, then use the actual origin and cell-centered spacing:
     * \f$A_x(z,t)=\sin(2\pi(z-t))/(2\pi)\f$. Our reference follows Fallahi,
     * MITHRA 2.0, Eqs. (3.6)-(3.9); it does not call referencePotential.
     */
    double expectedAx(PotentialField& field, int globalZ, double time) {
        const int count   = field.getLayout().getDomain()[2].length();
        const int wrapped = (globalZ % count + count) % count;
        const double z    = field.get_mesh().getOrigin()[2]
                         + (wrapped + 0.5) * field.get_mesh().getMeshSpacing()[2];
        return std::sin(WaveNumber * (z - time)) / WaveNumber;
    }

    /// @brief Maximum error normalized by potential amplitude 1/k; reject nonfinite values.
    double componentError(const ippl::Vector<double, 4>& value, double expected) {
        double error = 0.0;
        for (unsigned c = 0; c < 4; ++c) {
            if (!std::isfinite(value[c]))
                return std::numeric_limits<double>::infinity();
            error = std::max(error, WaveNumber * std::abs(value[c] - (c == 1 ? expected : 0.0)));
        }
        return error;
    }

    /** @brief Inspect owned cells and one face layer against the independent history.
     * @details Edges/corners are excluded: StandardFDTDSolver::step and
     * FDTDSolverBase::evaluate_EB consume faces only (StandardFDTDSolver.hpp:65-69;
     * FDTDSolverBase.hpp:113-118). Returns a rank-local maximum.
     */
    double historyError(PotentialField& field, double time) {
        const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, field.getView());
        const auto& local = field.getLayout().getLocalNDIndex();
        const int ghost   = field.getNghost();
        double error      = 0.0;
        for (int i = local[0].first() - 1; i <= local[0].last() + 1; ++i)
            for (int j = local[1].first() - 1; j <= local[1].last() + 1; ++j)
                for (int k = local[2].first() - 1; k <= local[2].last() + 1; ++k) {
                    if (!faceOrOwned(local, i, j, k))
                        continue;
                    const auto value =
                        host(i - local[0].first() + ghost, j - local[1].first() + ghost,
                             k - local[2].first() + ghost);
                    error = std::max(error, componentError(value, expectedAx(field, k, time)));
                }
        return error;
    }

    /// @brief Check every scratch component, including all allocated halo cells, is zero.
    double scratchError(PotentialField& field) {
        const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, field.getView());
        double error    = 0.0;
        for (std::size_t i = 0; i < host.extent(0); ++i)
            for (std::size_t j = 0; j < host.extent(1); ++j)
                for (std::size_t k = 0; k < host.extent(2); ++k)
                    error = std::max(error, componentError(host(i, j, k), 0.0));
        return error;
    }

    /** @brief Poison histories, initialize and compare owned cells plus both histories' faces.
     * @pre All layout ranks participate. Each rank receives global errors before assertions.
     * @details We sample our exact wave at 0/-dt to supply the two histories needed by
     * Fallahi, MITHRA 2.0, Eqs. (3.14)-(3.15); scratch must stay exactly zero.
     * Current-history faces feed the stencils; previous-history faces are checked
     * as an additional storage invariant, although those stencils do not consume them.
     */
    void expectInitialization(Fields& fields) {
        checkLocalDomain(fields.getFL());
        initializeFieldStorage(fields);
        VacuumSolver solver(fields.getJ(), fields.getE(), fields.getB());
        solver.A_n   = std::numeric_limits<double>::quiet_NaN();
        solver.A_nm1 = std::numeric_limits<double>::quiet_NaN();
        initializePotentials(solver);
        double errors[3] = {historyError(solver.A_n, 0.0),
                            historyError(solver.A_nm1, -solver.getDt()),
                            scratchError(solver.A_np1)};
        MPI_Allreduce(MPI_IN_PLACE, errors, 3, MPI_DOUBLE, MPI_MAX, fields.getFL().comm);
        EXPECT_LE(errors[0], 1e-10);
        EXPECT_LE(errors[1], 1e-10);
        EXPECT_EQ(errors[2], 0.0);
    }

    /// @brief Check the actual coarse demo, including inter-rank and periodic face values.
    TEST(VacuumInitialization, FillsOwnedCellsAndBothHistoryFaceHalos) {
        auto fields = makePeriodicFieldContainer();
        expectInitialization(fields);
    }

    /// @brief Check shifted coordinates, unequal spacings and the constructor-selected timestep.
    TEST(VacuumInitialization, UsesOriginSpacingAndActualSolverTimestep) {
        auto fields = shiftedGrid();
        expectInitialization(fields);
    }

#ifdef NDEBUG
    /// @brief Fail visibly if the test target was built without the intended debug checks.
    TEST(VacuumHaloGuards, RequiresDebugGuardsInThisTestTarget) {
        FAIL() << "Compile the test target and its VacuumSetup.cpp without NDEBUG";
    }
#else
    /** @brief Fresh periodic storage per guard test; fields outlive the solver.
     * @details Invalid configurations are introduced consistently on all ranks so local
     * guards throw before halo communication. This does not test rank-divergent misuse.
     */
    class VacuumHaloGuards : public ::testing::Test {
    protected:
        void SetUp() override {
            initializeFieldStorage(fields_m);
            solver_m =
                std::make_unique<VacuumSolver>(fields_m.getJ(), fields_m.getE(), fields_m.getB());
        }
        Fields fields_m = makePeriodicFieldContainer();
        std::unique_ptr<VacuumSolver> solver_m;
    };

    /// @brief Accept ordinary one-cell periodic storage and synchronize its faces.
    TEST_F(VacuumHaloGuards, AcceptsValidPeriodicStorage) {
        EXPECT_NO_THROW(preparePotentialHalo(solver_m->A_n));
    }

    /// @brief Reject a field before any view has been allocated.
    TEST_F(VacuumHaloGuards, RejectsUnallocatedStorage) {
        PotentialField empty;
        EXPECT_THROW(preparePotentialHalo(empty), std::logic_error);
    }

    /// @brief Reject two-cell halos, outside this helper's one-cell contract.
    TEST_F(VacuumHaloGuards, RejectsWrongHaloWidth) {
        PotentialField potential;
        potential.initialize(fields_m.getMesh(), fields_m.getFL(), 2);
        potential.setFieldBC(solver_m->A_n.getFieldBC());
        EXPECT_THROW(preparePotentialHalo(potential), std::logic_error);
    }

    /// @brief Reject a layout whose periodic flag disagrees with the helper contract.
    TEST_F(VacuumHaloGuards, RejectsNonperiodicLayout) {
        fields_m.getFL().isAllPeriodic_m = false;
        EXPECT_THROW(preparePotentialHalo(solver_m->A_n), std::logic_error);
    }

    /// @brief Reject a missing face boundary object.
    TEST_F(VacuumHaloGuards, RejectsMissingFace) {
        solver_m->A_n.getFieldBC()[0].reset();
        EXPECT_THROW(preparePotentialHalo(solver_m->A_n), std::logic_error);
    }

    /// @brief Reject a face boundary object with a nonperiodic type.
    TEST_F(VacuumHaloGuards, RejectsNonperiodicFace) {
        solver_m->A_n.getFieldBC()[0] = std::make_shared<ippl::NoBcFace<PotentialField> >(0);
        EXPECT_THROW(preparePotentialHalo(solver_m->A_n), std::logic_error);
    }

    /// @brief Reject a periodic boundary object assigned to the wrong face.
    TEST_F(VacuumHaloGuards, RejectsWrongFaceIndex) {
        solver_m->A_n.getFieldBC()[0] = std::make_shared<ippl::PeriodicFace<PotentialField> >(1);
        EXPECT_THROW(preparePotentialHalo(solver_m->A_n), std::logic_error);
    }
#endif
}  // namespace
