/** @file SerialReference.h
 * @brief Compare a distributed wave run with a rank-zero MPI_COMM_SELF execution.
 * @details The retained output tag "A9 communicator-self witness" identifies
 * records of the maximum per-cell and reduction-relative differences.
 */
#pragma once

#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <iostream>

#include "CoreTestSupport.h"
#include "OwnedStatistics.h"
#include "RankRecords.h"
#include "core/base/SpatialDimension.h"
#include "core/solver/BoxSetup.h"
#include "core/solver/ChdrSolver.h"

static_assert(
    chdr::solver::SpatialDim == 3,
    "SerialReference.h requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver::test {
    struct SerialSnapshot {
        std::vector<double> values_m;
        std::array<double, 3> reductions_m{};
    };
    using SerialState = std::array<SerialSnapshot, 4>;

    /** @brief Capture both per-cell values and genuinely distributed reductions.
     * @details Reducing the gathered vector would hide a change in MPI reduction
     * order. The maximum and the two positive sums therefore run on the field
     * before gatherOwned transfers its values in global index order.
     */
    template <class FieldType>
    SerialSnapshot serialSnapshot(const FieldType& field) {
        SerialSnapshot result;
        result.reductions_m = ownedStatistics(field).values();
        result.values_m     = gatherOwned(field);
        return result;
    }

    /// A zero reference scale requires exact equality; nonfinite values always fail.
    inline double serialRelative(double difference, double scale) {
        const double infinity = std::numeric_limits<double>::infinity();
        if (!std::isfinite(difference) || !std::isfinite(scale))
            return infinity;
        return scale == 0.0 ? (difference == 0.0 ? 0.0 : infinity) : difference / scale;
    }

    /** @brief Compare field entries and their independent global statistics.
     * @details Every entry uses the serial field's global absolute maximum as
     * scale. Each statistic uses its own serial value; serialRelative() handles
     * exact-zero scales and rejects nonfinite data without adding an epsilon floor.
     */
    inline void compareSerialSnapshots(const SerialSnapshot& distributed,
                                       const SerialSnapshot& serial,
                                       std::array<double, 2>& errors) {
        if (distributed.values_m.size() != serial.values_m.size())
            throw std::runtime_error("Serial reference has a different global field size");
        for (std::size_t i = 0; i < serial.values_m.size(); ++i)
            errors[0] = std::max(
                errors[0], serialRelative(std::abs(distributed.values_m[i] - serial.values_m[i]),
                                          serial.reductions_m[0]));
        for (unsigned r = 0; r < 3; ++r)
            errors[1] = std::max(errors[1], serialRelative(std::abs(distributed.reductions_m[r]
                                                                    - serial.reductions_m[r]),
                                                           serial.reductions_m[r]));
    }

    /// Capture both potential histories followed by the electric and magnetic outputs.
    template <class SolverType>
    SerialState serialState(SolverType& solver) {
        return {serialSnapshot(solver.A_n), serialSnapshot(solver.A_nm1),
                serialSnapshot(solver.electricField()), serialSnapshot(solver.magneticField())};
    }

    /** @brief Own a self-communicator mesh and keep its address stable for its fields.
     * @details FELFieldContainer constructs MPI_COMM_WORLD layouts, so this small
     * test owner initializes each object explicitly. It is created only on global
     * rank zero: HaloCells currently consults the global rank for one local index.
     * No global communicator or global tag state is replaced.
     */
    struct SerialFields {
        Mesh_t<SpatialDim> mesh_m;
        ippl::FieldLayout<SpatialDim> layout_m;
        AJFourFieldsType source_m;
        EMFieldsType electric_m, magnetic_m;

        explicit SerialFields(AJFourFieldsType& distributed)
            : mesh_m(distributed.getLayout().getDomain(), distributed.get_mesh().getMeshSpacing(),
                     distributed.get_mesh().getOrigin())
            , layout_m(MPI_COMM_SELF, distributed.getLayout().getDomain(), {false, false, false},
                       std::array<bool, 3>{distributed.getLayout().isPeriodic(0),
                                           distributed.getLayout().isPeriodic(1),
                                           distributed.getLayout().isPeriodic(2)}) {
            source_m.initialize(mesh_m, layout_m);
            electric_m.initialize(mesh_m, layout_m);
            magnetic_m.initialize(mesh_m, layout_m);
            source_m = 0.0;
        }
    };

    /** @brief Replay the distributed run with the same solver on a single-rank layout.
     * @details Copy the distributed material's cached uniform permittivity and
     * actual step; neither is re-estimated from a mean or a new safety factor.
     * The owner-free vacuum path remains owner-free. This reference isolates
     * decomposition effects; it shares the production update and is not an
     * independent reference for the underlying equations.
     */
    template <class SolverType, class AnalyticWave>
    SerialState evolveSerialReference(SerialFields& fields, SolverType& distributed,
                                      const AnalyticWave& wave, int steps,
                                      SpatialDiscretization spatialDiscretization) {
        using Base  = std::conditional_t<SolverType::stencil() == Stencil::Standard,
                                         StandardSolverPeriodic, NonStandardSolverPeriodic>;
        auto handle = coreTest::makeTestSolver<Base>(
            fields.source_m, fields.electric_m, fields.magnetic_m,
            distributed.ebReconstructionConvention(), spatialDiscretization,
            distributed.hasMaterial() ? distributed.homogeneousPermittivity() : 0);
        auto& serial = *handle;
        serial.setTimeStepUnchecked(distributed.getDt());
        loadHistory(serial, wave);
        serial.evaluateInitialFields();
        for (int step = 0; step < steps; ++step) {
            writeSource(serial, wave.sources(), step);
            serial.solve();
        }
        return serialState(serial);
    }

    /// Execute and compare on global rank zero while other ranks wait for the final broadcast.
    template <class SolverType, class AnalyticWave>
    std::array<double, 2> compareReferenceAtRoot(SolverType& solver, const AnalyticWave& wave,
                                                 int steps,
                                                 SpatialDiscretization spatialDiscretization,
                                                 const SerialState& before) {
        if (ippl::Comm->rank() != 0)
            throw std::runtime_error("Serial reference requires global rank zero");
        SerialFields fields(solver.A_n);
        const auto after =
            evolveSerialReference(fields, solver, wave, steps, spatialDiscretization);
        std::array<double, 2> errors{};
        for (unsigned field = 0; field < 4; ++field)
            compareSerialSnapshots(before[field], after[field], errors);
        const auto domain = solver.A_n.getLayout().getDomain();
        std::cout << "A9 communicator-self witness cells=" << domain[0].length() << 'x'
                  << domain[1].length() << 'x' << domain[2].length()
                  << " reading=" << int(spatialDiscretization) << " steps=" << steps
                  << " max_cell_relative=" << errors[0] << " max_reduction_relative=" << errors[1]
                  << '\n';
        return errors;
    }

    /** @brief Compare a completed wave with the same solver type on one rank.
     * @pre Collective on the distributed field communicator; its rank zero is global zero.
     * @details Both executions solve the same potential equations
     * [fallahi2020mithra20fullwavesimulation, Eqs. (3.6)-(3.9)] with identical
     * initial data, sources, time step and component positions. Repeating the
     * calculation on MPI_COMM_SELF isolates changes caused by domain decomposition.
     * The assertions require per-cell relative differences at most
     * \f$10^{-13}\f$ and statistic-relative differences at most \f$10^{-12}\f$,
     * using compareSerialSnapshots() and ownedStatistics() for their definitions.
     * Setting assertChecks to false still computes, prints and broadcasts both
     * errors, allowing the caller to retain a diagnostic record without assertions.
     * @tparam SolverType The same concrete solver type used for both decompositions.
     * @tparam AnalyticWave Analytic potential and source samples accepted by loadHistory.
     * @param solver Distributed state after the requested number of steps.
     * @param wave Initial potential histories and sources for the serial replay.
     * @param steps Number of solves to repeat from the labelled initial histories.
     * @param spatialDiscretization Component placement for the serial reconstruction.
     * @param assertChecks Whether to enforce the stated bounds in addition to recording errors.
     */
    template <class SolverType, class AnalyticWave>
    void verifySerialReference(SolverType& solver, const AnalyticWave& wave, int steps,
                               SpatialDiscretization spatialDiscretization,
                               bool assertChecks = true) {
        const auto comm = solver.A_n.getLayout().comm;
        int rank;
        MPI_Comm_rank(comm, &rank);
        const auto before = serialState(solver);
        std::array<double, 2> errors{};
        if (rank == 0)
            errors = compareReferenceAtRoot(solver, wave, steps, spatialDiscretization, before);
        MPI_Bcast(errors.data(), int(errors.size()), MPI_DOUBLE, 0, comm);
        if (assertChecks) {
            EXPECT_LE(errors[0], 1e-13);
            EXPECT_LE(errors[1], 1e-12);
        }
    }
}  // namespace chdr::solver::test
