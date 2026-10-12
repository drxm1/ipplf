# Periodic ChDR solver audit

This module adds a periodic, normalized-vacuum shell over IPPL's Standard and
NonStandard potential solvers. The inherited `step()` remains the update;
`ChdrSolver<Base>` supplies reconstruction of fields and time labels and the optional
time-step setter.

Five residual monitors and Fourier probe are separate objects that might be useful for measurements along the line.
All work is under `demos/chdr/solver`; the existing vacuum demo and `src/MaxwellSolvers` are reference implementations.

The module provides **both spatial readings**; applications select their
sampling convention explicitly. The supported common-time contract is
R2; R2Prime is its relabelled equivalence experiment, and R1 is the inherited
regression/control. Material, particles and open boundaries belong to later
units.

The constructor requires an explicit reconstruction. Select `EBReconstructionConvention::R2`
for the common-time contract. Both collocated and staggered routes are retained
as selectable implementations for later material and sourced applications; the
choice for a particular result will follow its benchmarks. The present vacuum
module supports Standard and NonStandard in the collocated reading, and Standard
in the staggered reading. Record the selected `SpatialDiscretization` and stencil with every result.

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

## Operators and monitors

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
on the periodic uniform grid, giving the curl/divergence cancellations below.
These are deductions from the staggered difference definitions and the lattice identities
in [Ryu et al., Appendix A][ryu] and [Clemens and Weiland, section 2][fit].

| Monitor or identity | Collocated reading | Staggered Standard reading |
|---|---|---|
| Divergence of magnetic curl | Zero to round-off | Zero to round-off |
| Faraday, raw magnetic fields at adjacent whole steps | Zero to round-off | Zero to round-off |
| Gauge and matching continuity | Same wave operator relates them | Same relationship with staggered differences |
| Gauss defect for R2 | Laplacian/composed-gradient difference plus negative gauge time difference | Negative gauge time difference |
| Laplacian equals gradient/divergence composition | Approximate for smooth resolved modes; exact mismatch given below | Exact Standard composition |

The exact curl/divergence cancellations are deductions from commuting
periodic differences [Clemens and Weiland, section 2, Eqs. (18)–(24)][fit].
The gauge propagation and Gauss relationship are the discrete substitutions
of [Christlieb et al., Lemma 3.3 and Theorem 3.5][gauge]. The collocated
Laplacian differs from the composition of centred first differences, so its
Gauss defect need not vanish. The explicit formulas and time slots are in
the documentation of `Monitors` and `evaluateGauss`.

Centred first differences annihilate checkerboard modes, so the collocated
gradient/curl monitors do not observe every grid-scale potential mode
[Trefethen, section 5.3, p. 203][trefethen]. Dedicated Ampère and null-mode
projection diagnostics are deferred to later source-coupling work.

To preserve zero Gauss defect in sourced staggered runs, initialize compatible
histories with zero matching discrete gauge residual and satisfy the matching
discrete source continuity equation (deduction from
[Christlieb et al., Theorem 3.3 and Theorem 3.5, Eq. (43)][gauge]).

For Standard, the symbol of the collocated mismatch is
\f$L_h-D_c\cdot D_c=-\sum_d4\sin^4(k_dh_d/2)/h_d^2\f$.
Its leading smooth-mode value is \f$-\sum_dh_d^2k_d^4/4\f$; hence on a cubic
mesh its magnitude relative to \f$k^2\f$ is \f$(kh)^2/4\f$ along an axis,
\f$(kh)^2/8\f$ on a face diagonal and \f$(kh)^2/12\f$ on the body diagonal.
At the all-axis checkerboard the mismatch magnitude is \f$\sum_d4/h_d^2\f$.
These are deductions by subtraction and Taylor expansion of the symbols above
[Trefethen, section 5.1, Eqs. (5.1.8)–(5.1.9)][trefethen]. For NonStandard,
subtract the same centred-gradient symbol from its stated wave symbol.

Monitor normalization evaluates the same expression with absolute inputs
and additions in place of subtractions. It reports both the maximum of the
cellwise ratio and the ratio of global maxima. Zero divided by zero is
defined as zero. This avoids normalizing an exact identity by its nearly
zero result. [Jeannerod and Rump, Eq. (1.6), Theorem 4.2][roundoff] motivate the
absolute scale for separately rounded inner products. With fused multiply-add
(FMA) contraction or nested stencil evaluation, the round-off reasoning instead
follows by applying the standard model, Eqs. (1.1)–(1.2), to the composed
expression (deduction, assuming round-to-nearest without underflow or overflow);
Theorem 4.2 does not directly cover those cases. Only owned cells
enter reductions. MPI exchanges plain doubles. The Gauss check retains the
source actually consumed by the solve, so replacing the next source does
not change its interpretation. R2Prime's Gauss report is one solve late.

One ghost layer and an all-periodic layout are required. After external
history writes, prepare both histories' halos. Monitor differences of E/B
and source arrays require their halos; composed differences exchange the
intermediate field before the second pass. Kokkos kernels operate on the
field execution/memory space. Evolution and stencil passes exchange halos
through IPPL; diagnostics additionally use the documented host mirrors,
MPI reductions and gathers. Do not capture a host solver `this` pointer inside
a device kernel.

The rank-invariance field-statistic checks use balanced pairwise local sums before their
plain-double MPI sum. This keeps their rounding depth logarithmic in the
cell count; the fixed sum over at most four components adds a bounded local
depth (deduction from [Jeannerod and Rump, p. 7][roundoff]). The evolved fields
are unaffected. Maxima still use MPI_MAX. The same helper serves both the
rank-record writer and the serial-reference comparison.

Load histories after selecting the time step, then call
`evaluateInitialFields()`. Before each solve, write the source at the labels
returned by the shell. For diagnostics, `Monitors::advance(solver)` performs
the before/after snapshots around the solve. Seed its previous source when
the first continuity report is needed. `solve()` itself remains usable without
allocating monitor work fields.
Monitoring retains five four-component fields and one three-component field,
including halos, in addition to the solver storage. Each evaluation launches
several stencil passes and halo exchanges; it is intended as an opt-in diagnostic.
Continuity is available only for adjacent source snapshots from the same solver
history. A skipped monitoring interval or history reload invalidates that pairing.

## Time step and probe

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

The Fourier probe computes the rectangular-rule time transform with a positive
time-exponential and no forward normalization:
\f$\widehat f=\Delta t\sum_n w_nf_n\exp(i\omega t_n)\f$.
Each channel reads its own time label and spatial offset from the contract.
The gate is explicit and half-open. Frequencies are cycles per step in
\f$(0,1/2)\f$; phases use the integer half-step index and FMA product remainder.
Windows are rectangle and periodic Hann
\f$w_n=[1-\cos(2\pi n/N)]/2\f$
[Oskooi et al., section 5.1, Eq. (6)][meep];
[Heinzel et al., Eq. (17), Appendix B][window].
Accumulators are local per owned cell; output reduces zero-padded real and
imaginary doubles with an ownership count of exactly one. Metadata includes
the gate, time step, frequency, window sums, labels and sampling coordinates.

The full convention uses the spatial forward kernel
\f$\exp(-i\mathbf k\cdot\mathbf x)\f$ and one inverse factor
\f$(2\pi)^{-1}\f$ per transformed variable; fields have phase
\f$\exp(i\mathbf k\cdot\mathbf x-i\omega t)\f$.
The probe transforms time only and applies no spatial area or volume factor.

For evolving fields, construct `SolverFourierProbe` after `loadHistory` and call
its parameterless `sample()` after each solve. The binding reads the solver's
completed-solve count, fields, time step and `LevelContract`; a gate containing
count zero reconstructs the initial E/B automatically. It allocates separate raw
B storage only for a requested raw-B channel. History reload or reinitialization
requires a new binding. The lower-level `FourierProbe::sample(step, ...)` remains
available for independently generated data; `step` must be the completed-solve
count, with initial fields reconstructed before sample zero.

For P global cells, C channels and F frequencies, gather allocates 2PCF doubles
on each rank and another 2PCF doubles on rank zero. Indices and PC coordinates
are also global on each rank; device accumulators use 2P_local CF doubles.
This implementation suits bounded point/plane diagnostics; account for those
global buffers before choosing a large surface or frequency list. Host phase
updates wait for the preceding sample kernel, including on asynchronous backends.

The `DISABLED_Holdout` suite is excluded from normal runs. It uses a different
wave vector and anisotropic cells and is run once after development checks
and their evidence are complete. The collocated case checks field/Gauss
convergence, discrete identities and rank agreement; the staggered case is an
additional record. Do not tune against its result. The exact inputs and
numerical assertions are defined in `Convergence.h` and the named test bodies.

## Library references

The citation keys used in source comments identify these works:

- `fallahi2020mithra20fullwavesimulation`: Arya Fallahi (2020), [MITHRA 2.0: A Full-Wave Simulation Tool for Free Electron Lasers][fallahi].
- `christlieb2024GaugeConserving`: Andrew J. Christlieb, William A. Sands and Stephen R. White (2025), [A Particle-in-Cell Method for Plasmas with a Generalized Momentum Formulation, Part III: A Family of Gauge Conserving Methods][gauge].
- `chew2014generalizedgauge`: Weng Cho Chew (2014), [Vector Potential Electromagnetic Theory with Generalized Gauge for Inhomogeneous Anisotropic Media][chew].
- `ryu2016potentialfdtd`: Christopher J. Ryu, Aiyin Y. Liu, Wei E. I. Sha and Weng Cho Chew (2016), [Finite-Difference Time-Domain Simulation of the Maxwell–Schrödinger System][ryu].
- `trefethen1996FiniteDifferenceSpectral`: Lloyd N. Trefethen (1996), [Finite Difference and Spectral Methods for Ordinary and Partial Differential Equations][trefethen].
- `schneider2010UnderstandingFdtd`: John B. Schneider (2010), [Understanding the Finite-Difference Time-Domain Method][schneider]; equation and printed-page references use the author PDF dated 2026-03-22.
- `clemens2001FiniteIntegration`: Markus Clemens and Thomas Weiland (2001), [Discrete Electromagnetism with the Finite Integration Technique][fit].
- `jeannerod2013ImprovedErrorBounds`: Claude-Pierre Jeannerod and Siegfried M. Rump (2013), [Improved Error Bounds for Inner Products in Floating-Point Arithmetic][roundoff].
- `oskooi2010meep`: Ardavan F. Oskooi et al. (2010), [Meep: A Flexible Free-Software Package for Electromagnetic Simulations by the FDTD Method][meep].
- `heinzel2002SpectrumSpectralDensity`: Gerhard Heinzel, Albrecht Rüdiger and Roland Schilling (2002), [Spectrum and Spectral Density Estimation by the Discrete Fourier Transform (DFT), Including a Comprehensive List of Window Functions and Some New Flat-Top Windows][window].

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
