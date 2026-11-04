/** @file test_plane_interface.cpp
 * @brief Compare coherent pulse amplitudes with independently derived lattice and Fresnel values.
 * @details Each test runs a plane pulse in a periodic box with measurePlane() of
 * PlanePulse.h. It compares the band-centre reflection r and transmission t at the
 * material face with the exact lattice coefficients and the continuum Fresnel
 * coefficients of InterfaceOracles.h. The recorded run parameters are in PlanePresets.h:
 * - Normal: TE at normal incidence from vacuum onto \f$\varepsilon_r = 2.13\f$, 812 by 2
 *   by 2 cubic cells, dielectric cells 0 to 210 (face at x = 211), Courant factor 0.9,
 *   59 box modes, 1180 solves, one rank.
 * - NormalSeam: the same run on two ranks; rank 0 must own x = 0 to 210, so the rank seam
 *   lies on the material face.
 * - Brewster: TM from vacuum at the Brewster angle onto \f$\varepsilon_r = 2.13\f$, 3581
 *   by 2 by 18 cells of size (1, 2, 2), 1484 dielectric cells, 88 box modes, 6761 solves,
 *   invariant sampled every 10th solve, one rank.
 * - ObliqueDielectricSeam: TM from \f$\varepsilon_r = 3.4\f$ at 18 degrees, 1828 by 2 by 52
 *   cells of size (1, 0.85, 1.25), 949 dielectric cells, rank seam on the material face,
 *   78 box modes, 3789 solves, two ranks.
 *
 * Tolerances: r and t within \f$10^{-11}\f$ of the lattice values; r within
 * \f$5 \times 10^{-3}\f$ (Normal, NormalSeam) or \f$10^{-2}\f$ of Fresnel, and t also,
 * except for TM incidence from vacuum; relative drift of the vector invariant at most
 * \f$10^{-12}\f$; for TE, the invariant equals Yee's energy at the start within
 * \f$10^{-12}\f$ (checkYee). Reference: exact lattice value and continuum analytic
 * (Fresnel).
 *
 * CTest: chdr.interface.normal (Normal), chdr.interface.normal.mpi2 (NormalSeam),
 * chdr.interface.long.brewster (Brewster) and chdr.interface.long.oblique-dielectric.mpi2
 * (ObliqueDielectricSeam). RowInput checks the angle/period input policies before allocation;
 * ResolvedOblique checks three frequencies of a new 30-degree TM row against the lattice
 * reference and an enlarged box/gate calculation within \f$10^{-7}\f$.
 */
#include <gtest/gtest.h>

#include "InterfaceOracles.h"
#include "apps/planar-interface/PlaneInput.h"
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
    /// Compare a measured frequency with the independent aligned lattice reference.
    void checkLatticeBin(const PlaneRow& p, const PlaneResult& result, size_t bin,
                         double tolerance) {
        const auto media    = planeMedia(p);
        const auto expected = latticeCoefficients(p.polarization_m, media[0], media[1], p.h_m,
                                                  p.dt_m, 0, p.kz_m, result.omega_m[bin]);
        std::cout << std::setprecision(17) << "omega=" << result.omega_m[bin]
                  << " r=" << result.r_m[bin] << " t=" << result.t_m[bin]
                  << " lattice errors=" << std::abs(result.r_m[bin] - expected.r) << ','
                  << std::abs(result.t_m[bin] - expected.t) << '\n';
        EXPECT_LE(std::abs(result.r_m[bin] - expected.r), tolerance);
        EXPECT_LE(std::abs(result.t_m[bin] - expected.t), tolerance);
    }

    /// Preserve the recorded continuum comparisons in the fixed Cartesian y basis.
    void checkContinuum(const PlaneRow& p, const PlaneResult& result, double tolerance) {
        const auto media = planeMedia(p);
        const auto expected =
            fresnelCoefficients(p.polarization_m, media[0], media[1], 0, p.kz_m, p.omega0_m);
        EXPECT_LE(std::abs(result.r_m[0] - expected.r), tolerance);
        if (p.polarization_m == Polarization::TE || !p.incidentFromVacuum_m) {
            EXPECT_LE(std::abs(result.t_m[0] - expected.t), tolerance);
        }
    }

    /** @brief Measure the shared resolved row and compare with independent expected values.
     * The preset supplies recorded numerical inputs, not a second copy of measurement logic.
     * The reference formulas remain independent of the resolver and field evolution.
     */
    void checkPlane(const PlaneRow& p, int cadence, double tolerance, size_t modes, int seam = 0) {
        const auto initialCheck = [&](auto& state, const auto& owner) {
            if (seam && ippl::Comm->rank() == 0) {
                const auto x = state.A_n.getLayout().getLocalNDIndex()[0];
                EXPECT_EQ(x.first(), 0);
                EXPECT_EQ(x.last() + 1, seam);
            }
            if (p.polarization_m == Polarization::TE)
                checkYee(state, owner);
        };
        const auto result = measurePlane(p, cadence, seam, initialCheck);
        ASSERT_EQ(result.modes_m, modes);
        EXPECT_LE(result.drift_m, 1e-12);
        if (ippl::Comm->rank() != 0)
            return;
        checkLatticeBin(p, result, 0, 1e-11);
        checkContinuum(p, result, tolerance);
    }
    /// TE pulse at normal incidence onto \f$\varepsilon_r = 2.13\f$, one rank: lattice r and t
    /// within 1e-11, Fresnel within 5e-3, drift at most 1e-12, Yee energy within 1e-12.
    TEST(PlaneInterface, Normal) {
        checkPlane(normalPlaneRow(), 1, 5e-3, 59);
    }
    /// The same normal pulse with its rank boundary on the material face, at the same thresholds.
    TEST(PlaneInterface, NormalSeam) {
        if (ippl::Comm->size() != 2)
            GTEST_SKIP() << "The interface seam requires two ranks";
        const auto p = normalPlaneRow();
        checkPlane(p, 1, 5e-3, 59, p.dielectricCells_m);
    }
    /// TM pulse at the Brewster angle onto \f$\varepsilon_r = 2.13\f$, one rank: lattice r and
    /// t within 1e-11, \f$\lvert r - r_F\rvert \le 10^{-2}\f$ against the Fresnel value
    /// \f$r_F\f$, which is zero at the Brewster angle, drift at most 1e-12; t is not compared
    /// with Fresnel.
    TEST(PlaneInterface, Brewster) {
        checkPlane(brewsterPlaneRow(), 10, 1e-2, 88);
    }
    /// TM pulse from \f$\varepsilon_r = 3.4\f$ at 18 degrees on two ranks, with the rank seam
    /// on the material face: lattice r and t within 1e-11, Fresnel within 1e-2, drift at most
    /// 1e-12.
    TEST(PlaneInterface, ObliqueDielectricSeam) {
        const auto p = dielectricPlaneRow();
        ASSERT_EQ(ippl::Comm->size(), 2);
        checkPlane(p, 10, 1e-2, 78, p.dielectricCells_m);
    }

    /// Selecting a compatible period and actual resolution reproduces the same central mode.
    void checkEquivalentPeriod(const PlaneRow& p) {
        PlaneRequest request;
        request.polarization_m       = p.polarization_m;
        request.cellsPerWavelength_m = p.cellsPerWavelength_m;
        request.transversePeriod_m   = p.cellsZ_m * p.h_m[2];
        const auto equivalent        = resolvePlaneRow(request);
        EXPECT_NEAR(equivalent.omega0_m, p.omega0_m, 1e-14);
        EXPECT_NEAR(equivalent.actualAngleDegrees_m, p.actualAngleDegrees_m, 1e-12);
    }

    /// Request parsing preserves the continuum angle and enforces mutually exclusive inputs.
    TEST(PlaneInterface, RowInput) {
        input::RunInput input;
        input.parameters_m["angle_degrees"]        = 30.0;
        input.parameters_m["cells_per_wavelength"] = 8.0;
        input.parameters_m["polarization"]         = "TM";
        const auto p                               = planeRow(input);
        EXPECT_NEAR(p.actualAngleDegrees_m, 30, 1e-12);
        EXPECT_NEAR(p.kz_m * p.cellsZ_m * p.h_m[2], 2 * Kokkos::numbers::pi, 1e-12);
        checkEquivalentPeriod(p);
        input.parameters_m["transverse_period"] = 24.0;
        EXPECT_THROW(planeRow(input), std::invalid_argument);
        PlaneRequest invalid;
        invalid.angleDegrees_m = 90;
        EXPECT_THROW(resolvePlaneRow(invalid), std::invalid_argument);
        invalid.angleDegrees_m   = 30;
        invalid.cutoffFraction_m = 0;
        EXPECT_THROW(resolvePlaneRow(invalid), std::invalid_argument);
    }

    /// Compare both layouts against independent lattice coefficients across the measured band.
    void checkExpandedPlane(const PlaneRow& p, const PlaneRow& expanded) {
        const auto first = measurePlane(p, 10), second = measurePlane(expanded, 10);
        if (ippl::Comm->rank() != 0)
            return;
        for (size_t bin = 0; bin < first.omega_m.size(); ++bin) {
            checkLatticeBin(p, first, bin, 1e-7);
            checkLatticeBin(expanded, second, bin, 1e-7);
            EXPECT_LE(std::abs(first.r_m[bin] - second.r_m[bin]), 1e-7);
            EXPECT_LE(std::abs(first.t_m[bin] - second.t_m[bin]), 1e-7);
        }
    }

    /// Denser speed sampling and enlarged layers/gates qualify one non-preset oblique TM row.
    TEST(PlaneInterface, ResolvedOblique) {
        PlaneRequest request;
        request.polarization_m       = Polarization::TM;
        request.angleDegrees_m       = 30;
        request.cellsPerWavelength_m = 8;
        const auto p                 = resolvePlaneRow(request);
        request.timingSamples_m      = 4097;
        const auto dense             = resolvePlaneRow(request);
        EXPECT_LE(std::abs(p.gateIncident_m - dense.gateIncident_m), 1);
        EXPECT_LE(std::abs(p.gateOther_m - dense.gateOther_m), 1);
        request.separationFactor_m = 1.2;
        request.gatePaddingSteps_m = 16;
        const auto expanded        = resolvePlaneRow(request);
        checkExpandedPlane(p, expanded);
    }
}  // namespace chdr::apps::planarInterface
