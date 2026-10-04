#pragma once

#include <Kokkos_Core.hpp>

#include <Kokkos_MathematicalConstants.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
#include <cmath>
#include <limits>
#include <mpi.h>

/** @file VacuumDiagnostics.h
 * @brief Owned-cell measurements for the normalized periodic transverse wave.
 * @details Source-free, three-dimensional, double-precision Standard solver only.
 * The declared mode has wavelength one and Ax amplitude 1/(2*pi). References
 * are independent of referencePotential() and the solver reconstruction.
 * Physics: Fallahi, MITHRA 2.0, arXiv:2009.13645v1, "Wave Equation",
 * Eqs. (3.6)-(3.9); discretization/dispersion: Eqs. (3.11)-(3.17).
 * Code contract: StandardFDTDSolver::step (StandardFDTDSolver.hpp:65-75) and
 * FDTDSolverBase::solve/timeShift/evaluate_EB (FDTDSolverBase.hpp:23-26,60-63,109-127).
 * These references justify our derived modal oracle below.
 * No field is changed, copied to the host, or halo-exchanged by these checks.
 */
namespace chdr::vacuum {

    /// @brief Globally reduced diagnostics; finite_m=false invalidates all comparisons.
    struct Diagnostics {
        int step_m                 = 0;
        double dt_m                = 0.0;
        double timeA_m             = 0.0;
        double timeE_m             = 0.0;
        double timeB_m             = 0.0;
        double relativeE_m         = 0.0;
        double relativeB_m         = 0.0;
        double discreteError_m     = 0.0;
        double zeroError_m         = 0.0;
        double sourceError_m       = 0.0;
        double phaseReal_m         = 0.0;
        double phaseImag_m         = 0.0;
        unsigned long long count_m = 0;
        bool finite_m              = false;
    };

    namespace diagnostic_detail {
        inline constexpr double WaveNumber = 2.0 * Kokkos::numbers::pi_v<double>;

        /// @brief Real and imaginary parts of one complex temporal Fourier coefficient.
        struct Mode {
            double real_m = 0.0;
            double imag_m = 0.0;
        };

        /// @brief Host-prepared scalar state captured by the owned-cell reduction.
        struct Parameters {
            int step_m          = 0;
            int firstZ_m        = 0;
            int ghost_m         = 0;
            double originZ_m    = 0.0;
            double spacingZ_m   = 0.0;
            double dt_m         = 0.0;
            double timeA_m      = 0.0;
            double timeE_m      = 0.0;
            double curlFactor_m = 0.0;
            Mode current_m;
            Mode previous_m;
            bool valid_m = false;
        };

        /** @brief Our closed-form recurrence solution with continuum initial histories.
         * @param omega Discrete angular advance per solver step, on the resolved branch.
         * @param phase Continuum angular advance per solver step.
         * @param step Temporal index, including -1 for the initial older history.
         * @details u(0)=1, u(-1)=exp(i*phase); derived from Fallahi (3.11)-(3.17)
         * and StandardFDTDSolver.hpp:65-69. Both temporal branches are retained.
         */
        inline Mode modeAt(double omega, double phase, int step) {
            const double sine     = std::sin(step * omega);
            const double betaReal = (std::cos(omega) - std::cos(phase)) / std::sin(omega);
            const double betaImag = -std::sin(phase) / std::sin(omega);
            return {std::cos(step * omega) + betaReal * sine, betaImag * sine};
        }

        /// @brief Reject invalid temporal/mesh inputs before evaluating the modal formula.
        inline bool validParameters(const Parameters& p) {
            const double q = WaveNumber * p.spacingZ_m;
            return p.step_m >= 0 && p.ghost_m > 0 && std::isfinite(p.originZ_m)
                   && std::isfinite(p.spacingZ_m) && p.spacingZ_m > 0.0 && std::isfinite(p.dt_m)
                   && p.dt_m > 0.0 && std::isfinite(p.timeA_m) && std::isfinite(p.timeE_m)
                   && q > 0.0 && q < Kokkos::numbers::pi_v<double>;
        }

        /// @brief Complete the positive-frequency mode parameters; no fitted frequency.
        inline void setModes(Parameters& p) {
            const double q        = WaveNumber * p.spacingZ_m;
            const double argument = (p.dt_m / p.spacingZ_m) * std::sin(q / 2.0);
            if (!(argument > 0.0 && argument < 1.0))
                return;
            const double omega = 2.0 * std::asin(argument);
            p.current_m        = modeAt(omega, WaveNumber * p.dt_m, p.step_m);
            p.previous_m       = modeAt(omega, WaveNumber * p.dt_m, p.step_m - 1);
            p.curlFactor_m     = std::sin(q) / q;
            p.valid_m = std::isfinite(p.current_m.real_m) && std::isfinite(p.current_m.imag_m)
                        && std::isfinite(p.previous_m.real_m) && std::isfinite(p.previous_m.imag_m);
        }

        /// @brief Use the actual solver timestep and the potential field's mesh/layout.
        template <class Solver>
        Parameters makeParameters(Solver& solver, int step) {
            const auto& potential = solver.A_n;
            const auto& mesh      = potential.get_mesh();
            Parameters p;
            p.step_m     = step;
            p.firstZ_m   = potential.getLayout().getLocalNDIndex()[2].first();
            p.ghost_m    = potential.getNghost();
            p.originZ_m  = mesh.getOrigin()[2];
            p.spacingZ_m = mesh.getMeshSpacing()[2];
            p.dt_m       = solver.getDt();
            p.timeA_m    = step * p.dt_m;
            p.timeE_m    = p.timeA_m - p.dt_m / 2.0;
            if (validParameters(p))
                setModes(p);
            return p;
        }

        /// @brief Ensure the common owned-cell policy can safely index the supplied field.
        template <class Reference, class Field>
        bool sameStorage(const Reference& reference, const Field& field) {
            if (reference.getNghost() != field.getNghost())
                return false;
            for (unsigned d = 0; d < 3; ++d) {
                if (reference.getView().extent(d) != field.getView().extent(d))
                    return false;
            }
            return &reference.getLayout() == &field.getLayout()
                   && &reference.get_mesh() == &field.get_mesh();
        }

        /// @brief Verify storage compatibility before entering a shared-index reduction.
        template <class Solver, class Fields>
        bool sameStorage(Solver& solver, Fields& fields, const Parameters& p) {
            return p.valid_m && sameStorage(solver.A_n, solver.A_nm1)
                   && sameStorage(solver.A_n, solver.A_np1)
                   && sameStorage(solver.A_n, fields.getE())
                   && sameStorage(solver.A_n, fields.getB())
                   && sameStorage(solver.A_n, fields.getJ());
        }

        /// @brief Global owned-cell count from the global layout; zero flags invalid size.
        template <class Layout>
        unsigned long long expectedCount(const Layout& layout) {
            unsigned long long count = 1;
            for (unsigned d = 0; d < 3; ++d) {
                const auto length = layout.getDomain()[d].length();
                if (length == 0 || length > std::numeric_limits<unsigned long long>::max() / count)
                    return 0;
                count *= length;
            }
            return count;
        }

        /// @brief Kokkos handles alias existing field storage; no field values are deep-copied.
        template <class PotentialView, class FieldView>
        struct Views {
            PotentialView current_m;
            PotentialView previous_m;
            PotentialView next_m;
            PotentialView source_m;
            FieldView electric_m;
            FieldView magnetic_m;
        };

        /// @brief Capture field handles after checking their shared mesh/layout/extents.
        template <class Solver, class Fields>
        auto makeViews(Solver& solver, Fields& fields) {
            return Views{solver.A_n.getView(),    solver.A_nm1.getView(),  solver.A_np1.getView(),
                         fields.getJ().getView(), fields.getE().getView(), fields.getB().getView()};
        }

        /// @brief A single cell's 22 stored components, including both other histories.
        struct Cell {
            Kokkos::Array<double, 4> current_m;
            Kokkos::Array<double, 4> previous_m;
            Kokkos::Array<double, 4> next_m;
            Kokkos::Array<double, 4> source_m;
            Kokkos::Array<double, 3> electric_m;
            Kokkos::Array<double, 3> magnetic_m;
        };

        /// @brief Load one owned cell without reading any halo.
        template <class ViewSet>
        KOKKOS_INLINE_FUNCTION Cell loadCell(const ViewSet& views, int i, int j, int k) {
            Cell cell;
            for (int c = 0; c < 4; ++c) {
                cell.current_m[c]  = views.current_m(i, j, k)[c];
                cell.previous_m[c] = views.previous_m(i, j, k)[c];
                cell.next_m[c]     = views.next_m(i, j, k)[c];
                cell.source_m[c]   = views.source_m(i, j, k)[c];
            }
            for (int c = 0; c < 3; ++c) {
                cell.electric_m[c] = views.electric_m(i, j, k)[c];
                cell.magnetic_m[c] = views.magnetic_m(i, j, k)[c];
            }
            return cell;
        }

        /// @brief Reject nonfinite input before it can disappear inside a maximum.
        KOKKOS_INLINE_FUNCTION bool finiteCell(const Cell& cell) {
            for (int c = 0; c < 4; ++c) {
                if (!Kokkos::isfinite(cell.current_m[c]) || !Kokkos::isfinite(cell.previous_m[c])
                    || !Kokkos::isfinite(cell.next_m[c]) || !Kokkos::isfinite(cell.source_m[c]))
                    return false;
            }
            for (int c = 0; c < 3; ++c) {
                if (!Kokkos::isfinite(cell.electric_m[c]) || !Kokkos::isfinite(cell.magnetic_m[c]))
                    return false;
            }
            return true;
        }

        /// @brief Per-cell analytic references and spatial Fourier factors.
        struct Reference {
            double current_m;
            double previous_m;
            double next_m;
            double electric_m;
            double magnetic_m;
            double continuumE_m;
            double continuumB_m;
            double sine_m;
            double cosine_m;
        };

        /// @brief Temporal-difference part of the independent raw-E oracle; Fallahi (3.9),(3.14).
        KOKKOS_INLINE_FUNCTION double electricReference(const Parameters& p, double sine,
                                                        double cosine) {
            return -(sine * (p.current_m.real_m - p.previous_m.real_m)
                     + cosine * (p.current_m.imag_m - p.previous_m.imag_m))
                   / (WaveNumber * p.dt_m);
        }

        /** @brief Apply independent difference symbols to the closed-form potential mode.
         * @details Our derivation from Fallahi (3.8)-(3.17), matched to
         * FDTDSolverBase.hpp:109-127. Raw E is at timeA-dt/2; raw B is at timeA.
         * A_np1 is zero scratch at step 0 and a current-state duplicate after solve().
         */
        KOKKOS_INLINE_FUNCTION Reference referenceAt(const Parameters& p, double z) {
            const double sine   = Kokkos::sin(WaveNumber * z);
            const double cosine = Kokkos::cos(WaveNumber * z);
            const double current =
                (sine * p.current_m.real_m + cosine * p.current_m.imag_m) / WaveNumber;
            const double previous =
                (sine * p.previous_m.real_m + cosine * p.previous_m.imag_m) / WaveNumber;
            const double electric = electricReference(p, sine, cosine);
            const double magnetic =
                p.curlFactor_m * (cosine * p.current_m.real_m - sine * p.current_m.imag_m);
            return {current,
                    previous,
                    p.step_m == 0 ? 0.0 : current,
                    electric,
                    magnetic,
                    Kokkos::cos(WaveNumber * (z - p.timeE_m)),
                    Kokkos::cos(WaveNumber * (z - p.timeA_m)),
                    sine,
                    cosine};
        }

        /// @brief Check expected values too, before accumulating error norms.
        KOKKOS_INLINE_FUNCTION bool finiteReference(const Reference& ref) {
            return Kokkos::isfinite(ref.current_m) && Kokkos::isfinite(ref.previous_m)
                   && Kokkos::isfinite(ref.next_m) && Kokkos::isfinite(ref.electric_m)
                   && Kokkos::isfinite(ref.magnetic_m) && Kokkos::isfinite(ref.continuumE_m)
                   && Kokkos::isfinite(ref.continuumB_m) && Kokkos::isfinite(ref.sine_m)
                   && Kokkos::isfinite(ref.cosine_m);
        }

        /// @brief Local/global reduction state; sums and maxima have separate MPI operations.
        struct Reduction {
            double errorE_m              = 0.0;
            double normE_m               = 0.0;
            double errorB_m              = 0.0;
            double normB_m               = 0.0;
            double phaseReal_m           = 0.0;
            double phaseImag_m           = 0.0;
            double discreteError_m       = 0.0;
            double zeroError_m           = 0.0;
            double sourceError_m         = 0.0;
            unsigned long long count_m   = 0;
            unsigned long long invalid_m = 0;
        };

        /// @brief Maximum unused potential/E/B component; potential scale is 1/k.
        KOKKOS_INLINE_FUNCTION double inactiveError(const Cell& cell) {
            double error = 0.0;
            for (int c = 0; c < 4; ++c) {
                if (c != 1) {
                    error = Kokkos::max(error, WaveNumber * Kokkos::abs(cell.current_m[c]));
                    error = Kokkos::max(error, WaveNumber * Kokkos::abs(cell.previous_m[c]));
                    error = Kokkos::max(error, WaveNumber * Kokkos::abs(cell.next_m[c]));
                }
            }
            for (int c = 0; c < 3; ++c) {
                if (c != 0)
                    error = Kokkos::max(error, Kokkos::abs(cell.electric_m[c]));
                if (c != 1)
                    error = Kokkos::max(error, Kokkos::abs(cell.magnetic_m[c]));
            }
            return error;
        }

        /// @brief Absolute source error against the identically zero prescribed source.
        KOKKOS_INLINE_FUNCTION double sourceError(const Cell& cell) {
            double error = 0.0;
            for (int c = 0; c < 4; ++c)
                error = Kokkos::max(error, Kokkos::abs(cell.source_m[c]));
            return error;
        }

        /// @brief Check all history Ax values and both nonzero raw field components.
        KOKKOS_INLINE_FUNCTION double discreteError(const Cell& cell, const Reference& ref) {
            double error = WaveNumber * Kokkos::abs(cell.current_m[1] - ref.current_m);
            error =
                Kokkos::max(error, WaveNumber * Kokkos::abs(cell.previous_m[1] - ref.previous_m));
            error = Kokkos::max(error, WaveNumber * Kokkos::abs(cell.next_m[1] - ref.next_m));
            error = Kokkos::max(error, Kokkos::abs(cell.electric_m[0] - ref.electric_m));
            return Kokkos::max(error, Kokkos::abs(cell.magnetic_m[1] - ref.magnetic_m));
        }

        /// @brief Accumulate continuum norms and the signed exp(-ikz) Fourier coefficient.
        KOKKOS_INLINE_FUNCTION void accumulate(const Cell& cell, const Reference& ref,
                                               Reduction& out) {
            const double errorE = cell.electric_m[0] - ref.continuumE_m;
            const double errorB = cell.magnetic_m[1] - ref.continuumB_m;
            out.errorE_m += errorE * errorE;
            out.normE_m += ref.continuumE_m * ref.continuumE_m;
            out.errorB_m += errorB * errorB;
            out.normB_m += ref.continuumB_m * ref.continuumB_m;
            out.phaseReal_m += cell.current_m[1] * ref.cosine_m;
            out.phaseImag_m -= cell.current_m[1] * ref.sine_m;
            const double inactive   = inactiveError(cell);
            const double source     = sourceError(cell);
            const double fieldError = Kokkos::max(discreteError(cell, ref), inactive);
            out.discreteError_m     = Kokkos::max(out.discreteError_m, fieldError);
            out.sourceError_m       = Kokkos::max(out.sourceError_m, source);
            out.zeroError_m         = Kokkos::max(out.zeroError_m, Kokkos::max(inactive, source));
        }

        /// @brief One fused owned-cell Kokkos reduction in the field's execution space.
        template <class ViewSet>
        struct MeasureKernel {
            using value_type = Reduction;
            ViewSet views_m;
            Parameters parameters_m;

            /// @brief Identity for this mixed sum/maximum reduction.
            KOKKOS_INLINE_FUNCTION void init(value_type& value) const { value = Reduction{}; }

            /// @brief Combine partial results without treating maxima as sums.
            KOKKOS_INLINE_FUNCTION void join(value_type& out, const value_type& in) const {
                out.errorE_m += in.errorE_m;
                out.normE_m += in.normE_m;
                out.errorB_m += in.errorB_m;
                out.normB_m += in.normB_m;
                out.phaseReal_m += in.phaseReal_m;
                out.phaseImag_m += in.phaseImag_m;
                out.discreteError_m = Kokkos::max(out.discreteError_m, in.discreteError_m);
                out.zeroError_m     = Kokkos::max(out.zeroError_m, in.zeroError_m);
                out.sourceError_m   = Kokkos::max(out.sourceError_m, in.sourceError_m);
                out.count_m += in.count_m;
                out.invalid_m += in.invalid_m;
            }

            /// @brief Count every owned cell, including invalid cells; skip their arithmetic.
            KOKKOS_INLINE_FUNCTION void operator()(int i, int j, int k, value_type& out) const {
                ++out.count_m;
                const Cell cell = loadCell(views_m, i, j, k);
                const double z  = parameters_m.originZ_m
                                 + (k + parameters_m.firstZ_m - parameters_m.ghost_m + 0.5)
                                       * parameters_m.spacingZ_m;
                const Reference ref = referenceAt(parameters_m, z);
                if (!finiteCell(cell) || !finiteReference(ref)) {
                    ++out.invalid_m;
                    return;
                }
                accumulate(cell, ref, out);
            }
        };

        /// @brief Communicate scalar summaries only; all layout ranks call in this order.
        inline Reduction globalReduce(const Reduction& local, MPI_Comm communicator) {
            double sums[6]   = {local.errorE_m, local.normE_m,     local.errorB_m,
                                local.normB_m,  local.phaseReal_m, local.phaseImag_m};
            double maxima[3] = {local.discreteError_m, local.zeroError_m, local.sourceError_m};
            unsigned long long counts[2] = {local.count_m, local.invalid_m};
            MPI_Allreduce(MPI_IN_PLACE, sums, 6, MPI_DOUBLE, MPI_SUM, communicator);
            MPI_Allreduce(MPI_IN_PLACE, maxima, 3, MPI_DOUBLE, MPI_MAX, communicator);
            MPI_Allreduce(MPI_IN_PLACE, counts, 2, MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);
            return {sums[0],   sums[1],   sums[2],   sums[3],   sums[4],  sums[5],
                    maxima[0], maxima[1], maxima[2], counts[0], counts[1]};
        }

        /// @brief Undefined/nonfinite or reference-zero norms are explicit failures.
        inline double relativeError(double error, double norm) {
            if (!std::isfinite(error) || error < 0.0 || !std::isfinite(norm) || norm <= 0.0)
                return std::numeric_limits<double>::infinity();
            return std::sqrt(error / norm);
        }

        /// @brief Carry actual raw-array timestamps even if later validation fails.
        inline Diagnostics baseDiagnostics(const Parameters& p) {
            Diagnostics result;
            result.step_m  = p.step_m;
            result.dt_m    = p.dt_m;
            result.timeA_m = p.timeA_m;
            result.timeE_m = p.timeE_m;
            result.timeB_m = p.timeA_m;
            return result;
        }

        /// @brief Normalize global sums and expose scalar validity, without pass/fail policy.
        inline Diagnostics finish(const Parameters& p, const Reduction& global,
                                  unsigned long long expected) {
            Diagnostics result     = baseDiagnostics(p);
            result.relativeE_m     = relativeError(global.errorE_m, global.normE_m);
            result.relativeB_m     = relativeError(global.errorB_m, global.normB_m);
            result.discreteError_m = global.discreteError_m;
            result.zeroError_m     = global.zeroError_m;
            result.sourceError_m   = global.sourceError_m;
            result.count_m         = global.count_m;
            const double factor    = global.count_m ? 2.0 * WaveNumber / global.count_m : 0.0;
            result.phaseReal_m     = factor * global.phaseReal_m;
            result.phaseImag_m     = factor * global.phaseImag_m;
            result.finite_m =
                global.invalid_m == 0 && expected > 0 && global.count_m == expected
                && std::isfinite(result.relativeE_m) && std::isfinite(result.relativeB_m)
                && std::isfinite(result.discreteError_m) && std::isfinite(result.zeroError_m)
                && std::isfinite(result.sourceError_m) && std::isfinite(result.phaseReal_m)
                && std::isfinite(result.phaseImag_m);
            return result;
        }
    }  // namespace diagnostic_detail

    /** @brief Measure an initial reconstructed state or a completed solve() result.
     * @param solver Existing Standard solver; no field or history is modified.
     * @param fields Allocated collocated E/B/J, sharing the solver's mesh/layout.
     * @param step Number of completed solves, zero immediately after initial evaluate_EB().
     * @return Global diagnostics, identical on all ranks for valid shared setup.
     * @pre All layout ranks call with the same step after reconstruction; the domain
     * spans one z wavelength and uses zero-based, unit-stride cell indices.
     * @details Potential errors use scale 1/k, field errors unit amplitude; reference
     * zero channels use absolute errors. See file documentation for library/code anchors.
     */
    template <class Solver, class Fields>
    Diagnostics measure(Solver& solver, Fields& fields, int step) {
        const auto parameters = diagnostic_detail::makeParameters(solver, step);
        int valid             = diagnostic_detail::sameStorage(solver, fields, parameters);
        MPI_Allreduce(MPI_IN_PLACE, &valid, 1, MPI_INT, MPI_MIN, fields.getFL().comm);
        if (!valid)
            return diagnostic_detail::baseDiagnostics(parameters);
        const auto views = diagnostic_detail::makeViews(solver, fields);
        diagnostic_detail::Reduction local;
        Kokkos::parallel_reduce(
            "Vacuum owned-cell diagnostics", solver.A_n.getFieldRangePolicy(),
            diagnostic_detail::MeasureKernel<decltype(views)>{views, parameters}, local);
        Kokkos::fence();
        const auto global = diagnostic_detail::globalReduce(local, fields.getFL().comm);
        return diagnostic_detail::finish(parameters, global,
                                         diagnostic_detail::expectedCount(solver.A_n.getLayout()));
    }

    /// @brief Phase result; step zero has no speed/increment and hasSpeed_m=false.
    struct PhaseSample {
        double amplitude_m = 0.0;
        double increment_m = 0.0;
        double speed_m     = 0.0;
        bool hasSpeed_m    = false;
        bool valid_m       = false;
    };

    /** @brief Track signed successive Fourier phases and cumulative mean phase speed.
     * @details The normalized exp(-ikz) coefficient is -i*exp(-ikt) for the continuum
     * wave; this is our Fourier substitution in Fallahi (3.6)-(3.9), not a fitted
     * reference. Positive-z propagation decreases phase. This class measures phase;
     * the caller owns amplitude, direction and speed acceptance thresholds.
     * Each instance belongs to one uninterrupted run starting at step zero.
     */
    class PhaseTracker {
    public:
        /// @brief Consume exactly the next globally reduced sample; never invoke MPI.
        PhaseSample update(const Diagnostics& sample) {
            PhaseSample result;
            result.amplitude_m = std::hypot(sample.phaseReal_m, sample.phaseImag_m);
            if (!validInput(sample, result.amplitude_m))
                return result;
            if (!initialized_m)
                return begin(sample, result);
            const double currentPhase = std::atan2(sample.phaseImag_m, sample.phaseReal_m);
            result.increment_m =
                std::remainder(currentPhase - previousPhase_m, 2.0 * Kokkos::numbers::pi_v<double>);
            accumulated_m += result.increment_m;
            result.speed_m    = -accumulated_m / (diagnostic_detail::WaveNumber * sample.timeA_m);
            result.hasSpeed_m = true;
            result.valid_m    = std::isfinite(result.increment_m) && std::isfinite(result.speed_m);
            remember(sample);
            return result;
        }

    private:
        /// @brief Reject invalid, missing, repeated or inconsistent-time samples.
        bool validInput(const Diagnostics& sample, double amplitude) const {
            if (!sample.finite_m || !std::isfinite(amplitude) || amplitude <= 0.0
                || !std::isfinite(sample.dt_m) || sample.dt_m <= 0.0
                || !std::isfinite(sample.timeA_m))
                return false;
            if (!initialized_m)
                return sample.step_m == 0 && sample.timeA_m == 0.0;
            return sample.step_m == previousStep_m + 1 && sample.dt_m == dt_m
                   && sample.timeA_m == sample.step_m * sample.dt_m;
        }

        /// @brief Establish the initial coefficient; no division by elapsed time.
        PhaseSample begin(const Diagnostics& sample, PhaseSample result) {
            initialized_m  = true;
            dt_m           = sample.dt_m;
            result.valid_m = true;
            remember(sample);
            return result;
        }

        /// @brief Keep only the previous phase and step, not full field snapshots.
        void remember(const Diagnostics& sample) {
            previousStep_m  = sample.step_m;
            previousPhase_m = std::atan2(sample.phaseImag_m, sample.phaseReal_m);
        }

        bool initialized_m     = false;
        int previousStep_m     = -1;
        double previousPhase_m = 0.0;
        double dt_m            = 0.0;
        double accumulated_m   = 0.0;
    };
}  // namespace chdr::vacuum
