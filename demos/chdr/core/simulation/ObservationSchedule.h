/** @file ObservationSchedule.h
 * @brief Sampling intervals in completed-solve indices, independent of any measured quantity.
 */
#pragma once

#include <cstdint>
#include <stdexcept>

namespace chdr::simulation {
    /** @brief Select observations from an inclusive interval of completed solver states.
     * @details State zero is the initialized state. Cadence is measured from first_m;
     * last_m=-1 leaves the upper end open. includeFinal_m additionally selects the run's
     * final state if it lies in the interval, even when it is not on the cadence.
     * All ranks use the same schedule because a measurement may perform MPI collectives.
     * Fourier accumulators that require consecutive samples must retain every_m=1.
     */
    struct ObservationSchedule {
        int every_m          = 1;
        std::int64_t first_m = 0, last_m = -1;
        bool includeFinal_m = false;

        /// Reject invalid intervals before entering any measurement callback.
        void validate() const {
            if (every_m <= 0 || first_m < 0 || last_m < -1 || (last_m >= 0 && last_m < first_m))
                throw std::invalid_argument("Invalid observation cadence or state interval");
        }

        /// Return whether this completed state is selected; an endpoint is never sampled twice.
        bool selects(std::int64_t step, std::int64_t finalStep) const {
            if (step < first_m || (last_m >= 0 && step > last_m))
                return false;
            return (step - first_m) % every_m == 0 || (includeFinal_m && step == finalStep);
        }
    };
}  // namespace chdr::simulation
