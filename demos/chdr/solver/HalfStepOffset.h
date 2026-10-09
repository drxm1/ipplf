/** @file HalfStepOffset.h
 * @brief Distinguish relative half-step counts from absolute step indices.
 */
#pragma once

namespace chdr::solver {
    /** @brief Signed sample-time offset measured in half timesteps.
     * @details Explicit construction and count() make conversion to and from
     * integer arithmetic deliberate. The type records the unit, not which
     * quantity or reconstruction the offset belongs to. It owns no field data
     * and imposes no restriction on the range of the stored integer.
     */
    class HalfStepOffset {
    public:
        /// Construct from a signed number of half timesteps, not a step index.
        explicit constexpr HalfStepOffset(int count)
            : count_m(count) {}

        /// Expose the half-step count for explicit time conversion or comparison.
        constexpr int count() const { return count_m; }

    private:
        int count_m;
    };
}  // namespace chdr::solver
