/** @file E04.h
 * @brief Entry point for the all-face Mur spherical pulse experiment.
 */
#pragma once

#include "core/input/RunInput.h"

namespace chdr::experiments::e04 {
    /// Compare periodic propagation and matched-speed Mur exit as specified in E04.typ.
    void runE04(const input::ChdrYamlRunInput& request);
}  // namespace chdr::experiments::e04
