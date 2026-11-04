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
#include "apps/planar-interface/TiltedMeasurement.h"
#include "tests/core/solver/support/CoreTestSupport.h"

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
        ASSERT_TRUE(result.valid(bin)) << result.reason_m[bin];
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

    /// A fixed physical wavelength and face fit the planned three-resolution normal ladder.
    void checkPhysicalRequests() {
        input::RunInput input;
        input.parameters_m["cell_shape_y"] = 2.0;
        input.parameters_m["cell_shape_z"] = 2.0;
        for (int n : {20, 40, 80})
            for (double fraction : {0.0, 0.25, 0.5}) {
                input.parameters_m["cells_per_wavelength"] = n;
                input.parameters_m["cell_size"]            = 20.0 / n;
                input.parameters_m["face_position"]        = 500 + fraction * 20 / n;
                const auto p                               = planeRow(input);
                EXPECT_NEAR(p.omega0_m, 2 * std::numbers::pi / (20 * std::sqrt(2.13)), 1e-15);
                EXPECT_EQ(p.dielectricCells_m, 25 * n + (fraction == 0.5));
                EXPECT_EQ(p.facePosition_m, 500 + fraction * p.h_m[0]);
            }
        input.parameters_m["face_position"] = 1.0;
        EXPECT_THROW(planeRow(input), std::invalid_argument);
        input.parameters_m["cell_size"] = 0.0;
        EXPECT_THROW(planeRow(input), std::invalid_argument);
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
        checkPhysicalRequests();
    }

    /** @brief Check sampled cells and the harmonic weights actually used by the TM initializer.
     * The cell predicate includes upper-face ties. Harmonic coefficients follow the existing
     * MaterialOwner rule [eymard2000FiniteVolume, section 7.1, Example 7.1].
     */
    void checkPlaneMaterial(InterfaceSolver& state, const material::MaterialOwner& owner,
                            const PlaneRow& p) {
        const PlanePulse pulse(p);
        const auto coords =
            material::detail::coordinates(state.A_n.get_mesh(), state.A_n.getLayout());
        const auto epsilon = owner.epsilon(), face = owner.derive(material::HarmonicFace<0>{});
        double local = 0, error = 0;
        Kokkos::parallel_reduce(
            "Plane sampled material", state.A_n.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k, double& worst) {
                const double x = coords.centre(coords.global(i, j, k))[0];
                const double a = pulse.epsilon(x), b = pulse.epsilon(x + coords.spacing_m[0]);
                worst = Kokkos::fmax(worst, Kokkos::abs(a - epsilon(i, j, k)));
                worst = Kokkos::fmax(worst, Kokkos::abs(2 / (1 / a + 1 / b) - face(i, j, k)));
            },
            Kokkos::Max<double>(local));
        MPI_Allreduce(&local, &error, 1, MPI_DOUBLE, MPI_MAX, state.A_n.getLayout().comm);
        EXPECT_LE(error, 1e-14);
        EXPECT_EQ(owner.statistics().marked_m,
                  std::uint64_t(p.dielectricCells_m) * p.cellsY_m * p.cellsZ_m);
    }

    /** @brief Move independent represented-face coefficients to the requested phase origin.
     * Deduction from the incident/reflected/transmitted propagation phases
     * [schneider2010UnderstandingFdtd, section 5.8.2, Eqs. (5.96)-(5.99)].
     */
    LatticeCoefficients shiftedPlaneReference(const PlaneRow& p, double omega) {
        const auto media = planeMedia(p);
        auto reference = latticeCoefficients(p.polarization_m, media[0], media[1], p.h_m, p.dt_m, 0,
                                             p.kz_m, omega);
        const double shift = p.facePosition_m - p.dielectricCells_m * p.h_m[0];
        const double sign  = p.incidentFromVacuum_m ? -1 : 1;
        reference.k1 *= sign;
        reference.k2 *= sign;
        reference.r *= std::exp(Complex(0, -2 * shift) * reference.k1);
        reference.t *= std::exp(Complex(0, shift) * (reference.k2 - reference.k1));
        return reference;
    }

    /// Synthetic fields at the probes must recover the given requested-face amplitudes.
    void checkSyntheticPlanePhase(const PlaneRow& p, double omega,
                                  const LatticeCoefficients& reference) {
        const Complex incoming(0.7, -0.4);
        const double xi        = p.incidentProbe_m - p.facePosition_m;
        const double xt        = p.otherProbe_m - p.facePosition_m;
        const Complex incident = incoming * std::exp(Complex(0, xi) * reference.k1);
        const Complex total =
            incident + incoming * reference.r * std::exp(Complex(0, -xi) * reference.k1);
        const Complex transmitted =
            incoming * reference.t * std::exp(Complex(0, xt) * reference.k2);
        const auto result = planeFaceCoefficients(p, omega, incident, total, transmitted);
        EXPECT_LE(std::abs(result[0] - reference.r), 1e-14);
        EXPECT_LE(std::abs(result[1] - reference.t), 1e-14);
    }

    /// A measured shifted face must match the independent aligned lattice result and phase shift.
    void checkShiftedPlane(const PlaneRow& p) {
        const auto result = measurePlane(p, 10, 0, [&](auto& state, const auto& owner) {
            checkPlaneMaterial(state, owner, p);
        });
        EXPECT_LE(result.drift_m, 1e-12);
        if (ippl::Comm->rank() != 0)
            return;
        for (size_t bin = 0; bin < result.omega_m.size(); ++bin) {
            ASSERT_TRUE(result.valid(bin)) << result.reason_m[bin];
            const auto reference = shiftedPlaneReference(p, result.omega_m[bin]);
            checkSyntheticPlanePhase(p, result.omega_m[bin], reference);
            EXPECT_LE(std::abs(result.r_m[bin] - reference.r), 1e-10);
            EXPECT_LE(std::abs(result.t_m[bin] - reference.t), 1e-10);
            std::cout << "shift/h=" << (p.facePosition_m / p.h_m[0] - p.dielectricCells_m)
                      << " lattice errors=" << std::abs(result.r_m[bin] - reference.r) << ','
                      << std::abs(result.t_m[bin] - reference.t) << '\n';
        }
    }

    /// Quarter-cell shifts retain the map; the inclusive half-cell tie moves it by one cell.
    TEST(PlaneInterface, ShiftedFace) {
        PlaneRequest request;
        request.cellsPerWavelength_m = 8;
        request.h_m                  = {1, 2, 2};
        for (const auto polarization : {Polarization::TE, Polarization::TM})
            for (double fraction : {0.0, 0.25, 0.5}) {
                request.polarization_m = polarization;
                request.facePosition_m = 200 + fraction;
                const auto p           = resolvePlaneRow(request);
                EXPECT_EQ(p.dielectricCells_m, 200 + (fraction == 0.5));
                checkShiftedPlane(p);
            }
        request.incidentFromVacuum_m = false;
        request.facePosition_m       = 500.5;
        checkShiftedPlane(resolvePlaneRow(request));
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

    /// Scaling every length and time together preserves dimensionless coefficients and cells.
    TEST(PlaneInterface, ScaledCell) {
        PlaneRequest request;
        request.cellsPerWavelength_m = 8;
        request.h_m                  = {1, 2, 2};
        request.facePosition_m       = 200;
        const auto original          = resolvePlaneRow(request);
        request.h_m                  = {0.5, 1, 1};
        request.facePosition_m       = 100;
        const auto scaled            = resolvePlaneRow(request);
        EXPECT_EQ(scaled.cellsX_m, original.cellsX_m);
        EXPECT_EQ(scaled.gateIncident_m, original.gateIncident_m);
        EXPECT_DOUBLE_EQ(scaled.omega0_m, 2 * original.omega0_m);
        EXPECT_DOUBLE_EQ(scaled.dt_m, original.dt_m / 2);
        checkExpandedPlane(original, scaled);
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
    /// Changing only the observation grid must not change the default pulse, box or gate.
    TEST(PlaneInterface, BandInput) {
        input::RunInput input;
        input.parameters_m["frequency_bins"] = 41;
        const auto p = planeRow(input), original = normalPlaneRow();
        EXPECT_EQ(p.cellsX_m, original.cellsX_m);
        EXPECT_EQ(p.pulseCentre_m, original.pulseCentre_m);
        EXPECT_EQ(p.gateIncident_m, original.gateIncident_m);
        EXPECT_EQ(p.omega_m.size(), 41);
        input.parameters_m["frequency_bins"] = 4;
        EXPECT_THROW(planeRow(input), std::invalid_argument);
    }

    /// A coarse near-grazing row can propagate discretely while lacking a continuum incident mode.
    TEST(PlaneInterface, BandDomain) {
        PlaneRequest request;
        request.angleDegrees_m       = 85;
        request.cellsPerWavelength_m = 8;
        const auto p                 = resolvePlaneRow(request);
        EXPECT_TRUE(planePropagationReason(p, p.omega_m[0]).empty());
        EXPECT_TRUE(planeLatticePropagationReason(p, p.omega_m[1]).empty());
        EXPECT_LT(p.omega_m[1], std::abs(p.kz_m));
        EXPECT_EQ(planePropagationReason(p, p.omega_m[1]), "continuum_incident_cutoff");
    }

    /// Weak and nonpropagating bins preserve their raw data but expose no plausible coefficient.
    void checkInvalidPlaneBins(const PlaneResult& result, size_t firstInvalid) {
        for (size_t bin = firstInvalid; bin < result.omega_m.size(); ++bin) {
            EXPECT_FALSE(result.valid(bin));
            EXPECT_TRUE(finitePlaneSpectrum(result.incident_m[bin]));
            EXPECT_TRUE(finitePlaneSpectrum(result.total_m[bin]));
            EXPECT_TRUE(finitePlaneSpectrum(result.transmitted_m[bin]));
            EXPECT_TRUE(std::isnan(result.r_m[bin].real()));
            EXPECT_TRUE(std::isnan(result.t_m[bin].real()));
        }
        EXPECT_EQ(result.reason_m[firstInvalid], "weak_incident");
        EXPECT_EQ(result.reason_m[firstInvalid + 1], "weak_incident");
        EXPECT_EQ(result.reason_m[firstInvalid + 2], "continuum_incident_cutoff");
    }

    /// A populated band is checked independently, including explicit unusable sampled frequencies.
    TEST(PlaneInterface, BandSpectrum) {
        PlaneRequest request;
        request.polarization_m       = Polarization::TM;
        request.angleDegrees_m       = 30;
        request.cellsPerWavelength_m = 8;
        request.frequencyBins_m      = 41;
        auto p                       = resolvePlaneRow(request);
        p.omega_m.push_back(planeFrequency(p, p.kx0_m - 7.5 * p.sigmaK_m));
        p.omega_m.push_back(planeFrequency(p, p.kx0_m + 7.5 * p.sigmaK_m));
        p.omega_m.push_back(0.1 * p.omega0_m);
        const auto result = measurePlane(p, 10);
        if (ippl::Comm->rank() != 0)
            return;
        for (size_t bin = 0; bin < 41; ++bin)
            checkLatticeBin(p, result, bin, 1e-8);
        checkInvalidPlaneBins(result, 41);
    }

    /// Every retained carrier is exactly Cartesian-periodic and has the requested tangent mode.
    void checkTiltedCarriers(const TiltedRow& r) {
        const auto& p = r.plane_m;
        EXPECT_GT(r.floquetMargin_m, 0);
        EXPECT_GT(r.retainedMargin_m, 0);
        for (const auto& mode : r.modes_m) {
            const double x = mode[0] * p.cellsX_m * p.h_m[0] / (2 * std::numbers::pi);
            const double z = mode[1] * p.cellsZ_m * p.h_m[2] / (2 * std::numbers::pi);
            EXPECT_NEAR(x, std::round(x), 1e-12);
            EXPECT_NEAR(z, std::round(z), 1e-12);
            EXPECT_NEAR(-r.sine_m * mode[0] + r.cosine_m * mode[1], r.kappa_m, 1e-14);
            EXPECT_NEAR(tiltedRoots(r, mode[2], 1)[0].normalK_m, mode[4], 1e-12);
        }
    }

    /// A compatible enlarged period preserves physical spectral controls exactly.
    void checkPinnedTiltedControl(TiltedRequest input, const TiltedRow& baseline) {
        input.primitiveX_m               = baseline.primitiveX_m;
        input.primitiveZ_m               = baseline.primitiveZ_m;
        input.periodMultiplier_m         = 2 * baseline.multiplier_m;
        input.tangentMode_m              = 2 * baseline.tangentMode_m;
        input.plane_m.separationFactor_m = 1.2;
        input.plane_m.probeCells_m += 4;
        const auto control = resolveTiltedRow(input);
        EXPECT_DOUBLE_EQ(control.kappa_m, baseline.kappa_m);
        EXPECT_DOUBLE_EQ(control.plane_m.omega0_m, baseline.plane_m.omega0_m);
        EXPECT_EQ(control.plane_m.omega_m, baseline.plane_m.omega_m);
        EXPECT_GT(control.dielectricWidth_m, baseline.dielectricWidth_m);
        EXPECT_GT(control.incidentWidth_m, baseline.incidentWidth_m);
        checkTiltedCarriers(control);
    }

    /// Two non-diagonal normals and an oblique incidence retain distinct directional roots.
    TEST(PlaneInterface, TiltedLayout) {
        TiltedRequest input;
        input.plane_m.h_m             = {1, 2, 2};
        input.plane_m.probeCells_m    = 12;
        input.plane_m.timingSamples_m = 65;
        for (double alpha : {17.0, 33.0}) {
            input.tiltDegrees_m = alpha;
            const auto r        = resolveTiltedRow(input);
            checkTiltedCarriers(r);
            std::cout << "tilt=" << alpha << " cells=" << r.plane_m.cellsX_m << 'x'
                      << r.plane_m.cellsZ_m << " steps=" << r.plane_m.steps_m << '\n';
        }
        input.plane_m.angleDegrees_m = 30;
        const auto r                 = resolveTiltedRow(input);
        checkTiltedCarriers(r);
        const auto roots = tiltedRoots(r, r.plane_m.omega0_m, 1);
        EXPECT_GT(std::abs(roots[0].normalK_m + roots[1].normalK_m), 1e-5);
        EXPECT_LT(roots[0].groupVelocity_m, 0);
        EXPECT_GT(roots[1].groupVelocity_m, 0);
        checkPinnedTiltedControl(input, r);
    }

    /// Both loaded histories satisfy the existing material-weighted divergence operator.
    void checkTiltedConstraint(InterfaceSolver& state, const material::MaterialOwner& owner) {
        for (const auto* field : {&state.A_n, &state.A_nm1}) {
            const solver::O1VSOperator op(*field, owner, 1);
            const double local = solver::coreTest::localMaximum(
                *field, KOKKOS_LAMBDA(int i, int j, int k) { return Kokkos::abs(op.g(i, j, k)); });
            EXPECT_LT(solver::coreTest::globalMaximum(*field, local), 1e-10);
        }
    }

    /// Compare every measured point with the exact finite-gate incident lattice mode sum.
    void checkTiltedHomogeneous(const TiltedRow& r, const TiltedResult& result) {
        if (ippl::Comm->rank() != 0)
            return;
        double error = 0, scale = 0;
        for (bool incident : {true, false}) {
            const auto& probe = incident ? result.incidentProbe_m : result.transmittedProbe_m;
            const int gate    = incident ? r.plane_m.gateIncident_m : r.plane_m.gateOther_m;
            for (size_t bin = 0; bin < r.plane_m.omega_m.size(); ++bin)
                for (size_t point = 0; point < probe.coordinates_m.size(); ++point) {
                    const auto expected = tiltedIncidentTransform(r, probe.coordinates_m[point],
                                                                  r.plane_m.omega_m[bin], gate);
                    error = std::max(error, std::abs(probe.value(point, 0, bin) - expected));
                    scale = std::max(scale, std::abs(expected));
                }
        }
        std::cout << "homogeneous relative field-spectrum error=" << error / scale << '\n';
        EXPECT_LT(error / scale, 1e-10);
        EXPECT_LT(result.drift_m, 1e-11);
    }

    /// Sampled material and both weighted-curl histories retain periodic geometry.
    void checkTiltedMaterialStart(TiltedRequest input) {
        input.plane_m.epsilon_m      = 2.13;
        input.plane_m.polarization_m = Polarization::TM;
        input.plane_m.facePosition_m = input.tiltDegrees_m > 0 ? 0.37 : -0.37;
        const auto r                 = resolveTiltedRow(input);
        const TiltedMaterialMap map(r);
        for (int i : {0, r.plane_m.cellsX_m / 2, r.plane_m.cellsX_m - 1})
            for (int k : {0, r.plane_m.cellsZ_m / 2, r.plane_m.cellsZ_m - 1}) {
                EXPECT_EQ(map.contains(i, k), map.contains(i + r.plane_m.cellsX_m, k));
                EXPECT_EQ(map.contains(i, k), map.contains(i, k + r.plane_m.cellsZ_m));
            }
        auto run = tiltedRun(r);
        run.load(TiltedPulse(r));
        checkTiltedConstraint(run.solver(), run.material());
    }

    /// The negative-incidence row also crosses both Cartesian winding directions.
    TiltedRow homogeneousTiltedRow(TiltedRequest& input) {
        input.plane_m.angleDegrees_m = input.tiltDegrees_m > 0 ? 0 : -20;
        auto r                       = resolveTiltedRow(input);
        if (input.tiltDegrees_m < 0) {
            input.plane_m.facePosition_m = -2 * r.normalPeriod_m - 0.37;
            r                            = resolveTiltedRow(input);
            EXPECT_NE(tiltedLift(r, 0), 0);
            EXPECT_NE(tiltedLift(r, 2), 0);
        }
        const auto anchor = tiltedAnchor(r);
        EXPECT_NEAR(r.cosine_m * anchor[0] + r.sine_m * anchor[2], r.facePosition_m, 1e-10);
        return r;
    }

    /// Exercise both polarization histories, signed angles/offsets and true component/time labels.
    TEST(PlaneInterface, TiltedHomogeneous) {
        TiltedRequest input;
        input.plane_m.h_m             = {1, 2, 2};
        input.plane_m.epsilon_m       = 1;
        input.plane_m.tailSigmas_m    = 4;
        input.plane_m.probeCells_m    = 12;
        input.plane_m.timingSamples_m = 33;
        input.probeHalfWidth_m        = 1;
        for (double alpha : {17.0, -33.0}) {
            input.tiltDegrees_m          = alpha;
            input.plane_m.facePosition_m = alpha > 0 ? 0.37 : -0.37;
            input.plane_m.polarization_m = alpha > 0 ? Polarization::TE : Polarization::TM;
            const auto r                 = homogeneousTiltedRow(input);
            auto run                     = tiltedRun(r);
            const auto result            = measureTiltedRun(run, r, 10, {}, checkTiltedConstraint);
            checkTiltedHomogeneous(r, result);
            checkTiltedMaterialStart(input);
        }
    }
}  // namespace chdr::apps::planarInterface
