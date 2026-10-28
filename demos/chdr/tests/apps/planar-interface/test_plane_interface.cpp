/** @file test_plane_interface.cpp
 * @brief Compare coherent pulse amplitudes with independently derived lattice and Fresnel values.
 */
#include <gtest/gtest.h>

#include "InterfaceOracles.h"
#include "apps/planar-interface/PlanePulse.h"

namespace chdr::apps::planarInterface {
    /** @brief Independently check the invariant's normalization with raw electromagnetic fields.
     * With \f$\Phi=0\f$ and zero gauge residual, the vector invariant equals Yee's electric
     * energy plus the magnetic cross-level product [sakkaplangkul2021YeeDispersiveMedia,
     * section 4.1, Theorem 4.1 and Remark 4.2]. No energy-operator helper enters this sum.
     */
    void checkYee(InterfaceSolver& s, const material::MaterialOwner& owner) {
        EMFieldsType oldB, newB;
        oldB.initialize(s.A_n.get_mesh(), s.A_n.getLayout());
        newB.initialize(s.A_n.get_mesh(), s.A_n.getLayout());
        solver::rawStaggeredMagnetic(s.A_nm1, oldB);
        s.computeRawMagneticField(newB);
        const auto e = s.electricField().getView(), a = oldB.getView(), b = newB.getView();
        const auto fx = owner.derive(material::HarmonicFace<0>{}),
                   fy = owner.derive(material::HarmonicFace<1>{}),
                   fz = owner.derive(material::HarmonicFace<2>{});
        double local = 0, total = 0;
        Kokkos::parallel_reduce(
            "Independent Yee energy", s.A_n.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k, double& sum) {
                const auto E = e(i, j, k), B = a(i, j, k), C = b(i, j, k);
                sum += fx(i, j, k) * E[0] * E[0] + fy(i, j, k) * E[1] * E[1]
                       + fz(i, j, k) * E[2] * E[2] + B[0] * C[0] + B[1] * C[1] + B[2] * C[2];
            },
            local);
        MPI_Allreduce(&local, &total, 1, MPI_DOUBLE, MPI_SUM, s.A_n.getLayout().comm);
        const auto h          = s.A_n.get_mesh().getMeshSpacing();
        const double residual = 0.5 * h[0] * h[1] * h[2] * total
                                    / diagnostics::potentialEnergy(s, owner).vector_m.total()
                                - 1;
        if (ippl::Comm->rank() == 0)
            std::cout << std::setprecision(17) << "Y/W_A-1=" << residual << '\n';
        EXPECT_LE(std::abs(residual), 1e-12);
    }
    /** @brief Assert incident-normalized complex errors at the material face, and energy drift.
     * The oracle defines the frozen mesh, spectrum, separate gates and exact expected values.
     * No fit or reference evolution enters the measurement. Scalar energy is zero for these
     * starts; the vector invariant supplies the positive denominator.
     */
    void checkPlane(const PlaneRun& p, int cadence, double continuumTolerance, int seam = 0) {
        const NormalPlane input{.cellsX_m             = p.cellsX,
                                .cellsY_m             = p.cellsY,
                                .cellsZ_m             = p.cellsZ,
                                .dielectricCells_m    = p.dielectricCells,
                                .h_m                  = p.h,
                                .epsilon_m            = p.epsilon,
                                .courant_m            = p.courant,
                                .ky_m                 = p.ky,
                                .dt_m                 = p.dt,
                                .omega0_m             = p.omega0,
                                .kx0_m                = p.kx0,
                                .sigmaK_m             = p.sigmaK,
                                .pulseCentre_m        = p.pulseCentre,
                                .incidentProbe_m      = p.incidentProbe,
                                .otherProbe_m         = p.otherProbe,
                                .incidentProbeIndex_m = p.incidentProbeIndex,
                                .otherProbeIndex_m    = p.otherProbeIndex,
                                .gateIncident_m       = p.gateIncident,
                                .gateOther_m          = p.gateOther,
                                .steps_m              = p.steps,
                                .incidentFromVacuum_m = p.incidentFromVacuum};
        const auto initialCheck = [&](auto& state, const auto& owner) {
            if (p.polarization != Polarization::TE)
                return;
            if (seam && ippl::Comm->rank() == 0) {
                const auto x = state.A_n.getLayout().getLocalNDIndex()[0];
                EXPECT_EQ(x.first(), 0);
                EXPECT_EQ(x.last() + 1, seam);
            }
            checkYee(state, owner);
        };
        const auto result =
            measurePlane(input, p.polarization == Polarization::TM, cadence, seam, initialCheck);
        ASSERT_EQ(result.modes_m, p.boxModes);
        EXPECT_LE(result.drift_m, 1e-12);
        if (ippl::Comm->rank() != 0)
            return;
        const double in  = p.incidentFromVacuum ? 1 : p.epsilon,
                     out = p.incidentFromVacuum ? p.epsilon : 1;
        const auto lattice =
            latticeCoefficients(p.polarization, in, out, p.h, p.dt, p.ky, 0, p.omega0);
        const auto fresnel = fresnelCoefficients(p.polarization, in, out, p.ky, 0, p.omega0);
        std::cout << std::setprecision(17) << "r=" << result.r_m[0] << " t=" << result.t_m[0]
                  << " lattice errors=" << std::abs(result.r_m[0] - lattice.r) << ","
                  << std::abs(result.t_m[0] - lattice.t)
                  << " Fresnel errors=" << std::abs(result.r_m[0] - fresnel.r) << ","
                  << std::abs(result.t_m[0] - fresnel.t) << " energy drift=" << result.drift_m
                  << " W_Phi=" << result.scalarEnergy_m << " S_abs/W0=" << result.absoluteRatio_m
                  << '\n';
        for (size_t f = 1; f < result.omega_m.size(); ++f) {
            const auto expected =
                latticeCoefficients(p.polarization, in, out, p.h, p.dt, p.ky, 0, result.omega_m[f]);
            std::cout << "band omega=" << result.omega_m[f] << " r=" << result.r_m[f]
                      << " t=" << result.t_m[f] << " lattice=" << expected.r << ',' << expected.t
                      << '\n';
        }
        EXPECT_LE(std::abs(result.r_m[0] - lattice.r), 1e-11);
        EXPECT_LE(std::abs(result.t_m[0] - lattice.t), 1e-11);
        EXPECT_LE(std::abs(result.r_m[0] - fresnel.r), continuumTolerance);
        if (p.polarization == Polarization::TE || !p.incidentFromVacuum) {
            EXPECT_LE(std::abs(result.t_m[0] - fresnel.t), continuumTolerance);
        }
    }
    TEST(PlaneInterface, Normal) {
        checkPlane(normalIncidenceRun(), 1, 5e-3);
    }
    /// The same normal pulse with its rank boundary on the material face, at the same thresholds.
    TEST(PlaneInterface, NormalSeam) {
        if (ippl::Comm->size() != 2)
            GTEST_SKIP() << "The interface seam requires two ranks";
        const auto p = normalIncidenceRun();
        checkPlane(p, 1, 5e-3, p.dielectricCells);
    }
    TEST(PlaneInterface, Brewster) {
        checkPlane(brewsterIncidenceRun(), 10, 1e-2);
    }
    // Registered disabled: reserve this complete configuration until final acceptance.
    TEST(PlaneInterface, DISABLED_Reserved) {
        const auto p = obliqueDielectricRun();
        ASSERT_EQ(ippl::Comm->size(), 2);
        checkPlane(p, 10, 1e-2, p.dielectricCells);
    }
}  // namespace chdr::apps::planarInterface
