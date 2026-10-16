/** @file ReflectionRun.h
 * @brief Shared thirty-step trajectories and transverse projections for reflection tests.
 */
#pragma once
#include "OpenRunSupport.h"
#include "TwoWaveState.h"

namespace chdr::solver::test::reflection {
    /** @brief Selected potential and field in global x,y,z order, on rank zero only.
     * @details TE selects Ay and Ey; TM0 selects Az and averaged By. The two real
     * trajectories supply the real and imaginary parts without mixing histories.
     */
    struct Result {
        std::vector<Complex> potential_m, field_m;
    };

    /// Append one real quadrature from a gathered component-interleaved field.
    inline void combine(std::vector<Complex>& result, const std::vector<double>& data,
                        unsigned components, unsigned component, int quadrature) {
        if (quadrature == 0)
            result.resize(data.size() / components);
        for (std::size_t i = 0; i < result.size(); ++i)
            result[i] += (quadrature == 0 ? Complex(1, 0) : Complex(0, 1))
                         * data[components * i + component];
    }

    /** @brief Evolve both exact quadratures with normal or deliberately defective finishing.
     * @details Control 1 changes only Mur speed by 1e-6; control 2 restores the
     * original halo ordering. The physical step and component time labels come
     * from the shell. Initial histories are prepared once, then solve() advances
     * thirty times. No observation exchanges potentials or repairs control halos.
     * Mesh/layout owners must outlive this routine; source is identically zero.
     */
    inline Result run(const Case& row, SpatialDiscretization reading, bool magnetic, int control,
                      AJFourFieldsType& source, EMFieldsType& electric, EMFieldsType& fieldB) {
        Result result;
        source = 0;
        for (int quadrature = 0; quadrature < 2; ++quadrature) {
            auto solver =
                makeOpenTestSolver(source, electric, fieldB, row.epsilon_m == 1 ? 0 : row.epsilon_m,
                                   EBReconstructionConvention::R2, reading);
            solver->setTimeStepUnchecked(row.dt_m);
            if (control == 1)
                solver->setKernelMurSpeedForTest(1 + 1e-6);
            if (control == 2)
                solver->useShippedOrder();
            loadHistory<StandardSolverPeriodic>(*solver,
                                                TwoWaveState(row, reading, magnetic, quadrature));
            for (int step = 0; step < 30; ++step)
                solver->solve();
            combine(result.potential_m, gatherOwned(solver->A_n), 4, magnetic ? 3 : 2, quadrature);
            combine(result.field_m, gatherOwned(magnetic ? fieldB : electric), 3, 1, quadrature);
        }
        return result;
    }

    /** @brief Remove the transverse carrier at each selected output position.
     * @details The normalized spatial DFT uses every transverse cell. It retains
     * temporal and normal phases, so fit() refers its coefficient to the first
     * supplied layer. Component offsets follow LevelContract; no interpolation
     * between staggered components is introduced. The wave-decomposition model
     * is described by Jones and Parrott (1988), Analysis, equations (1)-(3).
     */
    inline std::vector<Complex> project(const std::vector<Complex>& values, const Case& row,
                                        SpatialDiscretization reading, bool magnetic, bool field,
                                        int first, int count) {
        const LevelContract contract{EBReconstructionConvention::R2, reading};
        const auto offset =
            contract.positionOffset(field ? (magnetic ? Quantity::Magnetic : Quantity::Electric)
                                          : Quantity::VectorPotential,
                                    field ? 1 : (magnetic ? 2 : 1));
        const auto kt = row.tangent();
        std::vector<Complex> result(count);
        for (int i = first; i < first + count; ++i)
            for (int j = 0; j < row.ny_m; ++j)
                for (int k = 0; k < row.nz_m; ++k) {
                    const double phase = kt[1] * (j + .5 + offset[1]) * row.h_m[1]
                                         + kt[2] * (k + .5 + offset[2]) * row.h_m[2];
                    result[i - first] +=
                        values[(i * row.ny_m + j) * row.nz_m + k] * std::exp(Complex(0, -phase));
                }
        for (auto& value : result)
            value /= row.ny_m * row.nz_m;
        return result;
    }

    void verifyRanks(const Case&, SpatialDiscretization, bool magnetic, AJFourFieldsType& source,
                     const Result& distributed);
    void checkFieldReflection(const Case&, SpatialDiscretization, bool magnetic, const Result&);
}  // namespace chdr::solver::test::reflection
