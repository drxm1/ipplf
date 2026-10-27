/** @file test_mur_layer.cpp
 * @brief Reject unsupported observations on the outermost Mur layer.
 */
#include <gtest/gtest.h>

#include "core/solver/BoxSetup.h"
#include "opt-in/diagnostics/FourierProbe.h"
#include "opt-in/diagnostics/Monitors.h"
#include "tests/core/solver/CoreTestSupport.h"

namespace chdr::diagnostics::test {
    /** @brief Permit potential sampling at Mur while excluding reconstructed outer-layer fields.
     * The same layer rule applies to collocated and staggered readings. A box
     * touching only a periodic face remains valid. Monitors reject a layout with a
     * nonperiodic axis before allocating work arrays with periodic derivative assumptions.
     */
    TEST(MurLayer, ProbeAndMonitorRefusals) {
        auto fields = solver::makeBoxFieldContainer({8, 6, 10}, {8, 6, 10}, 2, {false, true, true});
        initializeFieldStorage(fields);
        for (auto reading : {SpatialDiscretization::Collocated, SpatialDiscretization::Staggered}) {
            const LevelContract levels{EBReconstructionConvention::R2, reading};
            for (auto quantity : {Quantity::Electric, Quantity::Magnetic, Quantity::RawMagnetic}) {
                EXPECT_THROW((FourierProbe(fields.getE(), levels, 0.5, {{0, 0, 0}, {0, 0, 0}},
                                           {{quantity, 0}}, {0.1}, 0, 4)),
                             std::invalid_argument);
                EXPECT_THROW((FourierProbe(fields.getE(), levels, 0.5, {{7, 0, 0}, {7, 0, 0}},
                                           {{quantity, 0}}, {0.1}, 0, 4)),
                             std::invalid_argument);
                EXPECT_NO_THROW((FourierProbe(fields.getE(), levels, 0.5, {{1, 0, 0}, {1, 0, 0}},
                                              {{quantity, 0}}, {0.1}, 0, 4)));
            }
            EXPECT_NO_THROW((FourierProbe(fields.getE(), levels, 0.5, {{0, 0, 0}, {0, 0, 0}},
                                          {{Quantity::Phi, 0}}, {0.1}, 0, 4)));
        }
        EXPECT_THROW(Monitors(fields.getJ()), std::invalid_argument);
    }
}  // namespace chdr::diagnostics::test
