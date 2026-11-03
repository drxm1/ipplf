/** @file ShippedOrderSolver.h
 * @brief Negative control of the halo-order tests: the original finishing order of the
 * potential step. HaloOrder.cpp, HaloOrderPulse.cpp, HaloOrderMixed.cpp and
 * OpenBoundaryRegression.cpp use it.
 */
#ifndef IPPL_TEST_SHIPPED_ORDER_SOLVER_H
#define IPPL_TEST_SHIPPED_ORDER_SOLVER_H

/** @brief Restore the original potential-step order for regression controls.
 * @tparam Solver A StandardFDTDSolver or NonStandardFDTDSolver specialization.
 * @details finishPotentialStep() exchanges the halo of the new level first and applies the
 * boundary update second, so the ghost copies of Mur-updated cells keep their values from
 * before the update. Only the final ordering changes. The production interior and Mur
 * kernels remain in use, so rank differences isolate the stale boundary halo.
 * Initial histories must be exchanged before either control or repaired runs.
 */
template <class Solver>
class ShippedOrderSolver : public Solver {
public:
    using Solver::Solver;

protected:
    void finishPotentialStep() override {
        this->A_np1.fillHalo();
        this->applyBCs();
    }
};

#endif
