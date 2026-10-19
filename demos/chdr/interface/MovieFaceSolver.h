/** @file MovieFaceSolver.h
 * @brief Demo boundary override for a dielectric that reaches two open faces.
 */
#pragma once
#include <cmath>
#include <stdexcept>

#include "../solver/ChdrSolver.h"

namespace chdr::interface {
    /** @brief Supply six face speeds to the existing Mur dispatcher for field movies.
     * @details This demo retains the dispatcher's unit-speed edge and corner updates.
     * It is not a validated heterogeneous absorbing boundary. The stored boundary()
     * descriptor still describes the unit-speed setup; movie metadata records the actual
     * face speeds separately. The physical timestep and interior update are unchanged.
     * See Mur (1981), section IV, Eq. (16), for the face approximation.
     */
    class MovieFaceSolver : public solver::ChdrSolver<solver::StandardSolverPeriodic> {
        using Base = solver::ChdrSolver<solver::StandardSolverPeriodic>;

    public:
        using Base::Base;

        /// Select positive subluminal face speeds before advancing; periodic entries are ignored.
        void setMovieFaceSpeeds(const ippl::Vector<double, 6>& speeds) {
            if (this->completedSolves() != 0)
                throw std::invalid_argument("Set movie face speeds before advancing");
            for (unsigned f = 0; f < 6; ++f)
                if (!std::isfinite(speeds[f]) || speeds[f] <= 0 || speeds[f] > 1)
                    throw std::invalid_argument("Movie face speeds must be in (0,1]");
            speeds_m = speeds;
        }

    protected:
        /// Apply Mur at the physical step, then publish owned samples and close open ghosts.
        void finishPotentialStep() override {
            const auto mask = this->boundary().murMask();
            if (!mask) {
                Base::finishPotentialStep();
                return;
            }
            ippl::Vector<uint32_t, 3> size = this->nr_m;
            size += 2 * this->A_n.getNghost();
            second_order_mur_boundary_conditions boundary;
            boundary.apply(this->A_n, this->A_nm1, this->A_np1, this->dt, size,
                           this->layout_mp->getLocalNDIndex(), mask, speeds_m);
            this->A_np1.fillHalo();
            ippl::fillOpenFaceGhosts(this->A_np1, mask);
        }

    private:
        ippl::Vector<double, 6> speeds_m = {1, 1, 1, 1, 1, 1};
    };
}  // namespace chdr::interface
