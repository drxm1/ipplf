/** @file planar-interface.cpp
 * @brief Periodic dielectric packet demonstrations and a plane-wave reflection measurement.
 */
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>

#include "InterfaceRun.h"
#include "PacketStart.h"
#include "PlanePulse.h"
#include "opt-in/output/FieldRunOutput.h"

namespace chdr::interface {
    /// Record the existing packet case in code axis order, with the x-z display plane explicit.
    std::string packetMetadata(const std::string& name, const RunBox& box,
                               const SolverSettings& settings, const material::PlanarSlab& slab,
                               const PacketStart& packet, int steps) {
        std::ostringstream text;
        text << std::setprecision(17) << "case=" << name << " plane=zx cells=" << box.cells_m[0]
             << ',' << box.cells_m[1] << ',' << box.cells_m[2] << " steps=" << steps
             << " h=" << box.spacing_m[0] << ',' << box.spacing_m[1] << ',' << box.spacing_m[2]
             << " S=" << settings.safetyFactor_m << " epsilon=2.13 slab=" << slab.lower_m << ','
             << slab.upper_m << " centre=" << packet.centreX_m << ',' << packet.centreZ_m
             << " widths=" << packet.width_m << ',' << packet.transverseWidth_m
             << " k=" << packet.waveNumber_m << " direction=" << packet.directionX_m << ','
             << packet.directionZ_m << " TM=" << packet.tm_m;
        return text.str();
    }

    /** @brief Run one named, reproducible input; --smoke selects a small four-step box. */
    int demonstration(int argc, char** argv) {
        if (argc < 2)
            throw std::invalid_argument(
                "A case name is required: plane, half-space, slab, brewster-te, brewster-tm");
        const std::string name = argv[1];
        bool smoke             = false;
        std::string directory;
        for (int i = 2; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--smoke")
                smoke = true;
            else if (arg == "--output" && i + 1 < argc)
                directory = argv[++i];
            else
                throw std::invalid_argument("Use CASE [--smoke] [--output NEW_DIRECTORY]");
        }
        if (name == "plane") {
            if (smoke || !directory.empty())
                throw std::invalid_argument("plane accepts no --smoke or --output option");
            const NormalPlane p;
            const auto result = measurePlane(p, false);
            if (ippl::Comm->rank() == 0) {
                const double n = std::sqrt(p.epsilon_m), r = (1 - n) / (1 + n);
                const double kt = 2 * std::asin(n * std::sin(p.omega0_m * p.dt_m / 2) / p.dt_m);
                // Sine matching at the two adjacent cell centres gives the exact TE lattice
                // coefficients [yee1966MaxwellIsotropicMedia, Eqs. (5)-(6); deduction].
                std::cout << std::setprecision(17) << "r=" << result.r_m[0]
                          << " t=" << result.t_m[0] << " Fresnel=" << r << "," << 1 + r
                          << " lattice="
                          << std::sin((p.kx0_m - kt) / 2) / std::sin((p.kx0_m + kt) / 2) << ","
                          << std::sin(p.kx0_m) / std::sin((p.kx0_m + kt) / 2)
                          << " energy drift=" << result.drift_m << '\n';
            }
            return 0;
        }
        if (name != "half-space" && name != "slab" && name != "brewster-te"
            && name != "brewster-tm")
            throw std::invalid_argument("Cases: plane, half-space, slab, brewster-te, brewster-tm");
        const bool brewster = name.starts_with("brewster");
        const int n = smoke ? 32 : 512, steps = smoke ? 4 : brewster ? 1000 : 560;
        const double scale = n / 512.0;
        const material::PlanarSlab slab{0, name == "slab" ? 156 * scale : 0, 230 * scale};
        const RunBox box{{n, 2, n}, {1, 1, 1}};
        const SolverSettings settings;
        InterfaceRun run(box, material::TwoMaterialMap{slab, 2.13, 1.0}, settings);
        PacketStart packet;
        packet.slab_m            = slab;
        packet.centreX_m         = 371 * scale;
        packet.centreZ_m         = n / (brewster ? 4.0 : 2.0);
        packet.length_m          = n;
        packet.width_m           = 29 * scale;
        packet.transverseWidth_m = 45 * scale;
        packet.tm_m              = name == "brewster-tm";
        if (brewster) {
            packet.directionX_m = -1 / std::sqrt(3.13);
            packet.directionZ_m = std::sqrt(2.13 / 3.13);
        }
        run.load(packet);
        std::unique_ptr<chdr::output::FieldRunOutput> output;
        if (!directory.empty())
            output = std::make_unique<chdr::output::FieldRunOutput>(
                directory, box, run.solver(), run.material(),
                packetMetadata(name, box, settings, slab, packet, steps));
        run.advance(steps, [&](long long k, auto& state, const auto& owner) {
            if (output && (k % 10 == 0 || k == steps))
                output->write(k, state, owner);
        });
        if (output)
            output->finish();
        return 0;
    }
}  // namespace chdr::interface
int main(int argc, char** argv) {
    ippl::initialize(argc, argv);
    int result = 0;
    try {
        result = chdr::interface::demonstration(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        result = 1;
    }
    ippl::finalize();
    return result;
}
