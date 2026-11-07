#include "core/input/Application.h"

namespace chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric {

    void runE01(const chdr::input::ChdrYamlRunInput& input) {}
}  // namespace chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric

int main(int argc, char** argv) {
    return chdr::input::runApplication(
        argc, argv, chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric::runE01);
}
