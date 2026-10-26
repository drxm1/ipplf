/** @file OpenRunSupport.h
 * @brief Shared source-free open-boundary test setup and deliberate failure controls.
 */
#pragma once
#include "CoreTestSupport.h"
#include "SerialReference.h"
#include "tests/src/MaxwellSolvers/ShippedOrderSolver.h"

namespace chdr::solver::test {
    /** @brief Test-only switches for the old finish order and a mismatched Mur speed.
     * @details Normal operation calls the production finish unchanged. The speed override
     * deliberately bypasses the stability and intersection guards while retaining finite,
     * positive speeds. It supplies negative controls, never a supported solver configuration.
     * Only step completion sees the override; loaded histories retain a validated descriptor.
     */
    template <class Solver>
    class OpenTestSolver : public ShippedOrderSolver<Solver> {
    public:
        using ShippedOrderSolver<Solver>::ShippedOrderSolver;
        void useShippedOrder(bool enabled = true) { shipped_m = enabled; }
        void setKernelMurSpeedForTest(double speed) {
            if (!std::isfinite(speed) || speed <= 0)
                throw IpplException("OpenTestSolver", "Control speed must be finite and positive");
            speed_m = speed;
        }

    protected:
        /** @brief Apply a uniform control speed to faces, edges and corners through the step.
         * @details Normalized Mur weights depend only on \f$v\Delta t\f$: see
         * murFaceWeights, murEdgeWeights and murCornerWeights in BoundaryOracles.h;
         * Mur (1981), section IV, Eq. (16), and Bamberger, Joly and Roberts (1987),
         * section II, Eqs. (2.10)-(2.11). Scaling the boundary step and using unit
         * coefficient speeds therefore changes all three operators, including intersections.
         * The interior has already advanced. Pass a local step and speed copy directly
         * to the Mur dispatcher; neither the solver step nor its const descriptor changes.
         * The shipped-order control omits the final exchange and ghost closure.
         */
        void finishPotentialStep() override {
            if (speed_m > 0 && this->boundary().murMask() != 0)
                finishWithControlSpeed();
            else if (shipped_m)
                ShippedOrderSolver<Solver>::finishPotentialStep();
            else
                Solver::finishPotentialStep();
        }

    private:
        /** @brief Retain the selected exchange order around the boundary-only control.
         * Normal completion publishes final owned samples and closes physical ghosts.
         * The shipped-order control deliberately leaves pre-Mur halo values in place.
         */
        void finishWithControlSpeed() {
            if (shipped_m)
                this->A_np1.fillHalo();
            applyControlBoundary();
            if (!shipped_m) {
                this->A_np1.fillHalo();
                ippl::fillOpenFaceGhosts(this->A_np1, this->boundary().murMask());
            }
        }

        /** @brief Perturb all Mur primitives without mutating production configuration.
         * Unit face speeds with the scaled local step perturb faces, edges and corners
         * alike, because each Mur weight depends only on the product of speed and step;
         * this includes the edge and corner primitives whose coefficient speed is fixed.
         * Periodic faces retain their copied entries and are excluded by the face mask.
         */
        void applyControlBoundary() {
            const auto& boundary = this->boundary();
            ippl::Vector<typename Solver::scalar, 6> speeds;
            for (unsigned face = 0; face < 6; ++face)
                speeds[face] =
                    boundary.kind[face] == ippl::FaceKind::mur ? 1. : boundary.speed[face];
            ippl::Vector<uint32_t, Solver::Dim> trueNr = this->nr_m;
            trueNr += 2 * this->A_n.getNghost();
            const typename Solver::scalar controlDt = this->dt * speed_m;
            second_order_mur_boundary_conditions bcs;
            bcs.apply(this->A_n, this->A_nm1, this->A_np1, controlDt, trueNr,
                      this->layout_mp->getLocalNDIndex(), boundary.murMask(), speeds);
        }
        bool shipped_m = false;
        double speed_m = 0;
    };

    using OpenStandardSolver = OpenTestSolver<ChdrSolver<StandardSolverPeriodic>>;

    /** @brief Build an ordinary shell or a test control on caller-owned field storage.
     * @details Zero epsilon selects owner-free vacuum. A positive value constructs a uniform
     * material; the shell retains the validated scalar value, so the temporary owner can die.
     * Source, output fields, mesh and layout must outlive the returned solver.
     */
    template <class Solver = OpenStandardSolver>
    auto makeOpenTestSolver(
        AJFourFieldsType& source, EMFieldsType& electric, EMFieldsType& magnetic,
        double epsilon                            = 0,
        EBReconstructionConvention reconstruction = EBReconstructionConvention::R1,
        SpatialDiscretization spatial             = SpatialDiscretization::Collocated) {
        if (epsilon == 0)
            return std::make_unique<Solver>(source, electric, magnetic, reconstruction, spatial);
        material::MaterialOwner owner(
            source.get_mesh(), source.getLayout(),
            KOKKOS_LAMBDA(ippl::Vector<int, 3>, ippl::Vector<double, 3>, bool& marked) {
                marked = false;
                return epsilon;
            });
        return std::make_unique<Solver>(source, electric, magnetic, reconstruction, spatial, owner);
    }

    /// Advance only potentials: these checks do not read reconstructed fields or their labels.
    template <class Solver>
    void advancePotentials(Solver& solver) {
        solver.step();
        solver.timeShift();
    }

    /// Absolute owned potential maximum, with nonfinite samples preserved as failures.
    inline double potentialMaximum(const AJFourFieldsType& field) {
        const auto view    = field.getView();
        const double local = coreTest::localMaximum(
            field, KOKKOS_LAMBDA(int i, int j, int k) {
                double maximum = 0;
                for (unsigned q = 0; q < 4; ++q) {
                    const double value = Kokkos::abs(view(i, j, k)[q]);
                    if (!Kokkos::isfinite(value))
                        return std::numeric_limits<double>::infinity();
                    maximum = Kokkos::max(maximum, value);
                }
                return maximum;
            });
        return coreTest::globalMaximum(field, local);
    }
}  // namespace chdr::solver::test
