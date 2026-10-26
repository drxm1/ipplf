# ChDR potential solver

`ChdrSolver<Base>` in `core/solver/ChdrSolver.h` extends IPPL's Standard and
NonStandard potential solvers in normalized units. It supplies field reconstruction,
time labels and the optional timestep setter. Without a material owner, its
update remains the inherited vacuum step. A bound material selects its update
through the host formulation interface.

The residual monitors, the Fourier probe and the field energies are opt-in
diagnostics under `demos/chdr/opt-in/diagnostics`, documented in the README there.
`ChdrSolver` is under `demos/chdr/core/solver`; shared boundary operators are in
`src/MaxwellSolvers`.

`ChdrSolver` supports **both spatial readings**; applications select their
sampling convention explicitly. The supported common-time contract is
R2; R2Prime is its relabelled equivalence experiment, and R1 is the inherited
regression/control. Homogeneous dielectric materials and mixed periodic/Mur
boundaries are supported. Particles and material interfaces require further work.

The constructor requires an explicit reconstruction. Select `EBReconstructionConvention::R2`
for the common-time contract. Both collocated and staggered routes are retained
as selectable implementations for later material and sourced applications; the
choice for a particular result will follow its benchmarks. `ChdrSolver`
supports Standard and NonStandard in the collocated reading, and Standard
in the staggered reading. Record the selected `SpatialDiscretization` and stencil with every result.

## Homogeneous material

Construct `chdr::material::MaterialOwner` on the fields' mesh and layout, then pass
it after the explicit reconstruction and spatial reading in the solver constructor.
`FormulationParameters` selects `MaterialFormulation::Homogeneous` and `alpha_m = 1` by
default. The guard requires `alpha_m = 1`, the Standard stencil and a finite constant
relative permittivity \f$\varepsilon\ge1\f$. This is the value sampled into the
owner's `epsilon()` field, not a derived coefficient. Homogeneity is checked
over owned samples across all ranks; halo copies do not contribute statistics.
The guard reads cached rank-identical statistics before
the base constructor enters collective operations. Nonuniform maps are permanently
outside this formulation's validity. Future interface formulations provide their
own guards, updates and bounds through `MaterialFormulation.h`.

The four potential components are \f$\psi=(\phi,\mathbf A)\f$, and the corresponding
physical source remains \f$(\rho,\mathbf J)\f$. For unit permeability,
\f$\partial_t^2\psi=\varepsilon^{-1}\nabla^2\psi+
(\rho/\varepsilon^2,\mathbf J/\varepsilon)\f$
[chew2014generalizedgauge, section 2, Eqs. (11)–(13)]. `HomogeneousUpdate` copies
\f$(\rho/\varepsilon,\mathbf J)\f$ into its own source field and calls the inherited
kernel with \f$\Delta t/\sqrt\varepsilon\f$ (deduction by substitution).
`ChdrSolver` restores the physical timestep and source binding on every exit path;
all time labels and reconstructions use the physical step. Unit permittivity
therefore reproduces the vacuum path bit for bit. The workspace stores four doubles
per allocated cell. Each step launches one device copy/scaling kernel over owned
cells, with no additional MPI exchange.

`maxStableTimeStep()` returns the Standard vacuum bound multiplied by
\f$\sqrt\varepsilon\f$ [fallahi2020mithra20fullwavesimulation, Eq. (3.20), with
medium speed]. `setSafetyFactor()` uses that bound. Construction retains the
inherited default timestep; select the desired safety factor before loading
histories. `Monitors` refuses a material solver unless its permittivity and its
gauge parameter are one, because its residual identities are vacuum identities
(`opt-in/diagnostics/README.md`). Medium convergence tests evaluate their own
gauge and Gauss norms.

Mesh, layout and source/output owners must retain their addresses for the
`ChdrSolver` instance's lifetime. The material owner is borrowed only during construction; the update
retains its scalar value and private workspace. Record `formulationName()`,
`formulationParameters().alpha_m` and `homogeneousPermittivity()` beside stencil,
spatial reading and reconstruction. Medium tests cover source scaling, phase,
stability and continuum convergence using independent expected values in
`tests/core/solver/MediumOracles.h`. Rank comparisons repeat the same update on one and two ranks.

## Periodic and Mur boundaries

`BoxSetup.h` constructs boxes with independent periodic axes. Its default is
periodic in every direction. These flags select the layout's exchange topology,
not an electromagnetic boundary condition. `ChdrSolver` selects Mur conditions
on the two global faces normal to each nonperiodic axis. The following example
uses Mur faces normal to x, periodic y and z axes, and decomposition along z:

```cpp
using namespace chdr::solver;
auto fields = makeBoxFieldContainer({40, 4, 48}, {40., 4., 48.}, 2,
                                   {false, true, true});
initializeFieldStorage(fields);
ChdrSolver<StandardSolverPeriodic> solver(fields.getJ(), fields.getE(), fields.getB(),
    EBReconstructionConvention::R2, SpatialDiscretization::Staggered);
solver.setSafetyFactor(0.9);
```

`ChdrSolver` supports Mur faces with Standard only. Periodic NonStandard remains
available with the collocated reading. Source and output fields must refer to the
same mesh and layout objects, with one ghost layer. Every rank needs at least two
owned cells per axis; each nonperiodic axis needs at least three global cells.
Global indices start at zero.

`ChdrSolver` retains the alias names `StandardSolverPeriodic` and
`NonStandardSolverPeriodic` for existing callers, but both select `ippl::configured`.
`ChdrSolver` derives and validates the face description before constructing its base.
Direct library users construct `StandardFDTDSolver<..., ippl::configured>` or
`NonStandardFDTDSolver<..., ippl::configured>` with
`(source, electric, magnetic, const FDTDBoundary&)`. The `periodic` and `absorbing`
specializations accept only the three fields; `absorbing` selects Mur on every
face and rejects periodically wrapping layouts; it does not infer face kinds from
the layout. `boundary()` is read-only. There is no
`setBoundary()`; construct a new solver to change faces or speeds.

`FDTDBoundary` describes six faces. Face `2*d` is lower, and `2*d+1` is upper.
Periodic faces agree with the layout and occur in opposing pairs. Mur faces
use the second-order condition of Mur [mur1981AbsorbingBoundaryFdtd, section IV,
Eq. (16)]. Edges and corners use the compatibility conditions of
[bamberger1987SecondOrderCorner, section II, Eqs. (2.10)–(2.11)].

`MurSpeed::physical(v)` specifies a speed in the time coordinate of
`ChdrSolver::getDt()`, with vacuum propagation speed normalized to one.
For `Homogeneous`, let \f$n=\sqrt{\varepsilon_r}\f$: the potential update uses
\f$\Delta t_{\mathrm{update}}=\Delta t/n\f$
[chew2014generalizedgauge, section 2, Eqs. (11)–(13); `HomogeneousUpdate.h`].
The Mur coefficients depend on speed times timestep
[mur1981AbsorbingBoundaryFdtd, section IV, Eq. (16)]. `ChdrSolver` therefore
passes \f$nv\f$ to Mur, so that
\f$(nv)\Delta t_{\mathrm{update}}=v\Delta t\f$ (deduction).
Without homogeneous timestep rescaling, the conversion factor is one.

`MurSpeed::matched()` instead passes exactly one to the Mur coefficient
constructors, equal to the propagation speed in the reused potential-update
equations. This corresponds to speed \f$1/n\f$ in the time coordinate of
`getDt()`. The API name does not assert perfect absorption. Pass the policy
after the spatial reading, or after the owner and formulation parameters for
a material solver. Material metadata and owner identity are validated before
reading permittivity; conversion and boundary validation finish before base
construction, not through a later setter.

Direct library callers supply `FDTDBoundary` speeds in the time coordinate of
the potential update; the library does not apply the conversion above. The
speed passed to each Mur face must be finite, positive and below \f$\sqrt{2}\f$.
This continuum guard is not a discrete stability proof
[trefethen1986OneWayWellPosed, Theorem 2, p. 424]. If Mur is selected on two or
three axes, every selected speed must be exactly one: the retained edge/corner
updates use the unit-interior-speed compatibility equations
[bamberger1987SecondOrderCorner, section II, Eqs. (2.10)–(2.11)]. With Mur on
only one axis, its two opposing faces do not intersect, and direct library
callers may supply separate validated speeds. A Mur–periodic intersection does
not impose the unit-speed restriction.

Ghost cells are extra storage around each rank's owned subdomain. `fillHalo()`
supplies inter-rank ghosts and periodic-wrap ghosts. Separate boundary helpers
supply exterior ghosts beyond nonperiodic global faces. The outermost owned
cells at those faces are not ghosts: Mur updates those owned cells.

With Mur faces selected, each step updates interior cells, then the selected
Mur faces, edges and corners. It exchanges the completed potentials before
filling exterior physical ghosts. `fillOpenFaceGhosts()` uses
\f$g=3u_0-3u_1+u_2\f$, where the three samples start at the outermost owned
cell and proceed inward. This gives a second-order one-sided derivative through
the unchanged centred reconstruction [fornberg1988FiniteDifferenceFormulas,
Table 3, p. 704]. Axis sweeps include transverse ghosts to complete edge and
corner values. `loadHistory()` prepares both stored potential levels with the
same halo exchange and exterior-ghost rule.
After loading arrays directly, call `markHistoryLoaded()` and `prepareHistories()`
before reconstruction or stepping. Preparation does not apply Mur to loaded owned samples.

Material ghosts wrap along periodic axes and copy the nearest material sample
beyond nonperiodic global faces. The Fourier probe and the monitors have their
own rules for the outermost Mur layer and for nonperiodic layouts
(`opt-in/diagnostics/README.md`). Keep particle deposits off the Mur layer: its update does not consume
a source, and halo accumulation does not wrap deposits across nonperiodic
global boundaries.

`MurReflection`, `MurFieldReflection` and `MurReflectionRefinement` compare
propagating, grazing and evanescent modes with independent discrete reflection
values. `MurEdgeCornerFamilies` checks local edge/corner compatibility.
`MurStability` checks fixed long trajectories, not stability for arbitrary
initial data; `OpenBeam` compares a finite beam in an all-Mur box with a
reference run using periodic transverse axes and Mur faces normal to x. The
disabled reflection holdout is configured separately.

Tests labelled `record` also print dense three-mesh reflection tables, vacuum
gauge and Gauss residuals, the field error from an unfilled physical ghost, and
the drift from unrestricted initial boundary data. The positive potential
diagnostic in `MurStabilityRecords` is defined by its temporal differences and
spatial links; it is not a conserved energy. `OpenGaugeRecords` distinguishes
the closed physical layer from an interior comparison with the plane-wave
oracle. These records carry their own labels and normalization and impose no
additional acceptance bounds.

## Units and labels

Use normalized units \f$c=\epsilon_0=\mu_0=1\f$, with beam direction z and the
future interface normal x, with the dielectric on \f$x<0\f$. Electron charge is
negative in FEL-style source input. For conversion back to SI, \f$c\f$, \f$\epsilon_0\f$
and \f$\mu_0\f$ in the formulas below denote their SI vacuum values. Given a
length scale \f$L_0\f$ and potential scale \f$\Phi_0\f$, the conversions are
\f$\mathbf x_{\rm SI}=L_0\mathbf x\f$,
\f$t_{\rm SI}=L_0t/c\f$,
\f$\phi_{\rm SI}=\Phi_0\phi\f$,
\f$\mathbf A_{\rm SI}=\Phi_0\mathbf A/c\f$,
\f$\mathbf E_{\rm SI}=\Phi_0\mathbf E/L_0\f$,
\f$\mathbf B_{\rm SI}=\Phi_0\mathbf B/(cL_0)\f$,
\f$\rho_{\rm SI}=\epsilon_0\Phi_0\rho/L_0^2\f$,
\f$\mathbf J_{\rm SI}=\epsilon_0c\Phi_0\mathbf J/L_0^2\f$.
These are dimensional substitutions into the potential equations of
[Fallahi, section 3.1.1, Eqs. (3.6)–(3.9)][fallahi]. No SI input adapter is
implemented here.

After k completed solves, `A_n` contains slot k and `A_nm1` slot k-1.
Component zero is the scalar potential; components 1–3 are the vector
potential. The source field stores charge density followed by current.
The caller loads both histories and the source at their declared labels.
Sources present before solve k+1 are those of slot k. Initialization resets
history, so initialize before loading it.

The entries below are half-step offsets m relative to k: the time is
\f$(k+m/2)\Delta t\f$. `LevelContract` owns these labels. The R2 labels follow
[Fallahi, sections 3.1.2, 3.1.4, 3.2.2][fallahi] and
[Christlieb et al., Eqs. (25)–(29)][gauge]. R2Prime shifts the scalar-potential
and charge slots by one step while retaining the same reconstructed fields
when its histories are consistent (deduction from the componentwise update).

Bitwise R2/R2Prime equality assumes matching floating-point contraction in both
kernels; changing the rounded operation sequence can change the result
[Jeannerod and Rump, Eqs. (1.1)–(1.2)][roundoff]. The tests retain bitwise assertions.
Only a demonstrated, recorded contraction difference permits the fallback
\f$\max\lvert E_{\mathrm{R2}}-E_{\mathrm{R2Prime}}\rvert
\le 10^{-15}\max\lvert E_{\mathrm{R2}}\rvert\f$.
Both maxima cover owned cells, components and ranks, rather than normalizing each
component separately. Record the architecture and contraction flags with any fallback.

| Quantity | R1 | R2 | R2Prime |
|---|---:|---:|---:|
| Vector potential, current | 0 | 0 | 0 |
| Scalar potential, charge | 0 | +1 | −1 |
| Electric field | mixed: −1 and 0 | −1 | −1 |
| Output magnetic field | 0 | −1 | −1 |
| Raw magnetic field | 0 | 0 | 0 |

R1 electric output has no single time label for a time-dependent scalar
potential. A probe channel for it therefore requires the explicit
`phiIsZero` declaration. The R1 convergence control is compared at the
vector-potential difference midpoint; this deliberately exposes its scalar-term
time offset rather than assigning a single time to the combined field.

`SpatialDiscretization::Collocated` samples everything at cell centres.
`SpatialDiscretization::Staggered` is supported only for R2 and Standard: scalar potential
and charge at centres, vector potential/current/electric components on their
positive faces, magnetic components on the transverse positive edges.
`LevelContract::positionOffset` reports the extra offset from the cell centre
in cell units. In homogeneous vacuum with \f$\epsilon=\mu=\alpha=1\f$,
[Chew, section 2, Eqs. (11)–(13), and the paragraph after Eq. (18)][chew]
give the Lorenz gauge and the component wave equations used here. Component
positions and staggered differences follow [Ryu et al., Appendix A, Eqs. (19)–(24)][ryu];
the arrays and potential update are shared.

## Operators

For collocated reconstruction, the spatial first difference is
\f$D_cf_i=(f_{i+1}-f_{i-1})/(2h)\f$, with symbol
\f$i\sin(kh)/h\f$ [Trefethen, section 5.1, Eqs. (5.1.8)–(5.1.9)][trefethen]. The Standard wave
Laplacian has symbol \f$-\sum_d K_d^2\f$, where
\f$K_d=2\sin(k_dh_d/2)/h_d\f$. The NonStandard symbol is
\f$-[S_z(K_x^2+K_y^2)+K_z^2]\f$, with
\f$S_z=1-4\mathcal A\sin^2(k_zh_z/2)\f$ and IPPL's coefficient
\f$\mathcal A=\tfrac14[1+0.02/((h_z/h_x)^2+(h_z/h_y)^2)]\f$
[Fallahi, Eqs. (3.17), (3.21)–(3.25)][fallahi].


For the staggered reading, the forward gradient maps centre scalars to faces,
and the forward curl maps face vectors to edges. Backward divergence returns
face vectors to centres; magnetic divergence instead uses forward differences
from edges to vertices. The backward curl returns edge vectors to faces for
the corresponding field equation. These are the staggered differences of
[Ryu et al., Appendix A][ryu].

Their gradient/divergence composition is the Standard Laplacian: in each
direction, direct substitution gives
\f$D^-D^+f_i=(f_{i+1}-2f_i+f_{i-1})/h^2\f$. Cross-direction differences commute
on the periodic uniform grid, so the divergence of the curl and the curl of the
gradient vanish to round-off in both readings. These are deductions from the
staggered difference definitions and the lattice identities in
[Ryu et al., Appendix A][ryu] and [Clemens and Weiland, section 2,
Eqs. (18)–(24)][fit]. The collocated Laplacian differs from the composition
of centred first differences, whereas the staggered composition is exactly the
Standard Laplacian. The residual monitors that check these identities on
evolved data, with their normalization and time slots, are documented in
`opt-in/diagnostics/README.md`.

Centred first differences annihilate checkerboard modes, so derivative-based
checks in the collocated reading do not observe every grid-scale potential mode
[Trefethen, section 5.3, p. 203][trefethen].

For Standard, the symbol of the collocated mismatch is
\f$L_h-D_c\cdot D_c=-\sum_d4\sin^4(k_dh_d/2)/h_d^2\f$.
Its leading smooth-mode value is \f$-\sum_dh_d^2k_d^4/4\f$; hence on a cubic
mesh its magnitude relative to \f$k^2\f$ is \f$(kh)^2/4\f$ along an axis,
\f$(kh)^2/8\f$ on a face diagonal and \f$(kh)^2/12\f$ on the body diagonal.
At the all-axis checkerboard the mismatch magnitude is \f$\sum_d4/h_d^2\f$.
These are deductions by subtraction and Taylor expansion of the symbols above
[Trefethen, section 5.1, Eqs. (5.1.8)–(5.1.9)][trefethen]. For NonStandard,
subtract the same centred-gradient symbol from its stated wave symbol.

Kokkos kernels operate on the field execution/memory space. Evolution and
stencil passes exchange halos through IPPL. Do not capture a host solver
`this` pointer inside a device kernel.

Load histories after selecting the time step, then call
`evaluateInitialFields()`. Before each solve, write the source at the labels
returned by `ChdrSolver`. `solve()` allocates no diagnostic storage; the
residual monitors, the Fourier probe and the field energies are separate
objects and functions (`opt-in/diagnostics/README.md`).

## Time step

Without a setter call the inherited time step is retained. The checked
safety factor S sets \f$\Delta t=S\Delta t_{\max}\f$, with
\f$\Delta t_{\max}=(\sum h_d^{-2})^{-1/2}\f$ for Standard and
\f$\Delta t_{\max}=h_z\f$ for admissible NonStandard aspect ratios
[Fallahi, Eqs. (3.19), (3.20), (3.24)][fallahi]. Standard accepts
\f$0<S<1\f$; NonStandard accepts \f$0<S\le1\f$. Setters refuse changes after
histories have been loaded or reconstructed; initialize before changing their
time spacing. The unchecked entry point is for tests and records only.
No filter is added. NonStandard's axial
Nyquist endpoint is a repeated root; its recorded growth follows the
discrete recurrence (deduction from Eq. (3.24)).

## Tests and evidence

The solver tests are in `demos/chdr/tests/core/solver` and build the executable
`test-chdr-solver` with `BUILD_TESTING` and `IPPL_ENABLE_UNIT_TESTS`.
The project Nix configuration also enables CHDR, FEL, FFT and SOLVERS with
SERIAL and OPENMP. Use at most two compiler jobs for local development.

From the IPPL checkout, retain the usual toolchain and dependency options when
configuring. For an existing configured build, enable the relevant modules with:

```sh
cmake -S . -B build -DBUILD_TESTING=ON -DIPPL_ENABLE_UNIT_TESTS=ON \
  -DIPPL_ENABLE_CHDR=ON -DIPPL_ENABLE_FEL=ON \
  -DIPPL_ENABLE_FFT=ON -DIPPL_ENABLE_SOLVERS=ON
cmake --build build --parallel 2 --target test-chdr-solver
ctest --test-dir build --no-tests=error -R '^chdr\.solver\.' \
  -LE 'long|record' --output-on-failure -j1
ctest --test-dir build --no-tests=error -R '^chdr\.solver\.' \
  -L 'long|record' --output-on-failure -j1
```

The last selection runs the longer ladders and supplementary records separately.
Use `-V` to retain numerical records from passing tests. Python3 and an MPI
launcher are required for complete rank comparisons; unavailable comparisons
appear as skipped tests. The tests of the opt-in diagnostics (`chdr.diagnostics.*`)
and of the programs under `apps/` (`chdr.vacuum.*`, `chdr.interface.*`) are
described in their own READMEs.

The CTests separate bounded unit checks, long ladders and records. Tests run
on one and two MPI ranks with one OpenMP thread per rank. Numerical assertions
use global reductions. Rank fixtures compare global-index owned-cell dumps,
including the source/potential seam cases. The rank-invariance field-statistic
checks (`tests/core/solver/OwnedStatistics.h`) use balanced pairwise local sums
before their plain-double MPI sum. This keeps their rounding depth logarithmic
in the cell count; the fixed sum over at most four components adds a bounded
local depth (deduction from [Jeannerod and Rump, p. 7][roundoff]). The evolved
fields are unaffected. Maxima still use MPI_MAX. The same helper serves both
the rank-record writer and the serial-reference comparison. Configure `CHDR_SOLVER_RECORD_ROOT`
to relocate generated records; its default is
`build/demos/chdr/tests/core/solver/runs`. It contains serial and MPI2
subdirectories. CTest assigns each producer's `CHDR_SOLVER_RECORD_DIR`; a direct
executable invocation can set that variable itself. Preserve these directories
with the run log.

| Test suites | Main checks |
|---|---|
| `Manufactured` | Source sign, all components and time labels; wrong-label control |
| `Levels`, `GaugeWave`, `LevelEquivalence` | Discrete gauge cancellation, R1 closed form and bitwise R2/R2Prime equivalence |
| `Convergence`, `LongConvergence`, `StaggeredConvergence`, `LongStaggeredConvergence` | Mixed source-free and analytic sourced waves, compared at each output's label |
| `LongGaugeGrowth` | Gauge norm in the first and last of ten periods |
| `Dispersion`, `RecordDispersion` | Single-root cumulative phase and angular dispersion |
| `TimeStep`, `TimeStepLong`, `TimeStepRecords` | Permitted range, effective safety factor, bounded/growing histories and endpoint records |
| `ComovingT7` | Co-moving mode and bunch-train reconstruction versus exact symbols and leading-order predictions |
| `RankGather`, `RankCases`, `SerialReference`, `StableStatistics` | Owned fields and well-conditioned reductions across rank layouts |
| `ChdrSolverContract` | Construction contract of `ChdrSolver`, and its evolution against IPPL's Standard solver on the vacuum-wave setup |
| `Operators`, `CompactOperators`, `BoxSetup` | Difference symbols, periodic halos and geometry |
| `ChdrStaggered` | Staggered difference composition, component positions and the Standard staggered R2 reconstruction |

The NonStandard source-amplitude record contains two host-side symbol predictions.
The manufactured test checks that the update adds its source without averaging;
the smooth sourced-wave record does not test evanescent near-field response.
That application remains a separate experiment.

The `DISABLED_Holdout` suite is excluded from normal runs. It uses a different
wave vector and anisotropic cells and is run once after development checks
and their evidence are complete. The collocated case checks field/Gauss
convergence, discrete identities and rank agreement; the staggered case is an
additional record. Do not tune against its result. The exact inputs and
numerical assertions are defined in `Convergence.h` and the named test bodies.

The four upstream Maxwell examples remain compile-only. The Standard example's
update, dispersion and rank claims are covered here by `Manufactured`,
`Dispersion` and `RankCases`. The NonStandard convergence example and the two
image demos retain their compile-only registration; this module does not
modify `test/maxwell`.
Whether to convert the shared NonStandard convergence example into an asserted
runtime test remains unresolved (an open decision). The mapping
above records current coverage and does not settle that upstream change.

## Library references

The citation keys used in source comments identify these works:

- `fallahi2020mithra20fullwavesimulation`: Arya Fallahi (2020), [MITHRA 2.0: A Full-Wave Simulation Tool for Free Electron Lasers][fallahi].
- `christlieb2024GaugeConserving`: Andrew J. Christlieb, William A. Sands and Stephen R. White (2025), [A Particle-in-Cell Method for Plasmas with a Generalized Momentum Formulation, Part III: A Family of Gauge Conserving Methods][gauge].
- `kreiss2002SecondOrderWave`: Heinz-Otto Kreiss, N. Anders Petersson and Jacob Yström (2002), [Difference Approximations for the Second Order Wave Equation](https://doi.org/10.1137/S0036142901397435), SIAM J. Numer. Anal. 40, 1940–1967.
- `chew2014generalizedgauge`: Weng Cho Chew (2014), [Vector Potential Electromagnetic Theory with Generalized Gauge for Inhomogeneous Anisotropic Media][chew].
- `ryu2016potentialfdtd`: Christopher J. Ryu, Aiyin Y. Liu, Wei E. I. Sha and Weng Cho Chew (2016), [Finite-Difference Time-Domain Simulation of the Maxwell–Schrödinger System][ryu].
- `trefethen1996FiniteDifferenceSpectral`: Lloyd N. Trefethen (1996), [Finite Difference and Spectral Methods for Ordinary and Partial Differential Equations][trefethen].
- `schneider2010UnderstandingFdtd`: John B. Schneider (2010), [Understanding the Finite-Difference Time-Domain Method][schneider]; equation and printed-page references use the author PDF dated 2026-03-22.
- `clemens2001FiniteIntegration`: Markus Clemens and Thomas Weiland (2001), [Discrete Electromagnetism with the Finite Integration Technique][fit].
- `jeannerod2013ImprovedErrorBounds`: Claude-Pierre Jeannerod and Siegfried M. Rump (2013), [Improved Error Bounds for Inner Products in Floating-Point Arithmetic][roundoff].
- `oskooi2010meep`: Ardavan F. Oskooi et al. (2010), [Meep: A Flexible Free-Software Package for Electromagnetic Simulations by the FDTD Method][meep].
- `jones1988MultiPointImpedance`: Michael G. Jones and Tony L. Parrott (1988), [Evaluation of a Multi-Point Method for Determining Acoustic Impedance](https://ntrs.nasa.gov/citations/19880017529), NASA TM-100637; Analysis, equations (1)–(3), p. 3.
- `heinzel2002SpectrumSpectralDensity`: Gerhard Heinzel, Albrecht Rüdiger and Roland Schilling (2002), [Spectrum and Spectral Density Estimation by the Discrete Fourier Transform (DFT), Including a Comprehensive List of Window Functions and Some New Flat-Top Windows][window].
- `mur1981AbsorbingBoundaryFdtd`: Gerrit Mur (1981), [Absorbing Boundary Conditions for the Finite-Difference Approximation of the Time-Domain Electromagnetic-Field Equations](https://doi.org/10.1109/TEMC.1981.303970).
- `bamberger1987SecondOrderCorner`: Alain Bamberger, Patrick Joly and Jean E. Roberts (1987), [Second Order Absorbing Boundary Conditions for the Wave Equation: A Solution for the Corner Problem](https://inria.hal.science/inria-00075909).
- `fornberg1988FiniteDifferenceFormulas`: Bengt Fornberg (1988), [Generation of Finite Difference Formulas on Arbitrarily Spaced Grids](https://doi.org/10.1090/S0025-5718-1988-0935077-0).
- `trefethen1986OneWayWellPosed`: Lloyd N. Trefethen and Laurence Halpern (1986), [Well-Posedness of One-Way Wave Equations and Absorbing Boundary Conditions](https://doi.org/10.1090/S0025-5718-1986-0856695-2).

[fallahi]: https://arxiv.org/abs/2009.13645
[gauge]: https://arxiv.org/abs/2410.18414
[chew]: https://arxiv.org/abs/1406.4780
[ryu]: https://doi.org/10.1109/JMMCT.2016.2605378
[trefethen]: https://people.maths.ox.ac.uk/trefethen/pdetext.html
[schneider]: https://eecs.wsu.edu/~schneidj/ufdtd/
[fit]: https://doi.org/10.2528/PIER00080103
[roundoff]: https://doi.org/10.1137/120894488
[meep]: https://doi.org/10.1016/j.cpc.2009.11.008
[window]: https://hdl.handle.net/11858/00-001M-0000-0013-557A-5
