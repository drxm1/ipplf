# ChDR opt-in diagnostics

The six headers in this directory measure a running `ChdrSolver<Base>` (`core/solver/ChdrSolver.h`).
They change neither its update nor its potential histories. All names are in namespace
`chdr::solver`.

- `Monitors.h` evaluates residuals of the discrete identities of the vacuum potential update.
- `MonitorExpressions.h` supplies the residual arithmetic, its normalization and the reductions.
- `MonitorRecords.h` evaluates the Gauss-law defect and compares it with its discrete formula.
- `FourierProbe.h` accumulates running time transforms of selected field components.
- `FieldEnergy.h` sums the two discrete potential invariants of the O1-VS formulation.
- `AlignedFieldEnergy.h` sums the electromagnetic energy of the time-aligned staggered R2 fields.

The diagnostics are opt-in. No file under `core/` includes them, and `ChdrSolver::solve()`
allocates and evaluates none of them. A program that uses them links the CMake target
`chdr-diagnostics`, an INTERFACE library that links `chdr-core`, and includes each header by its
path from `demos/chdr`, for example `#include "opt-in/diagnostics/Monitors.h"`. Units and labels
follow `core/solver/README.md`, section "Units and labels": \f$c=\epsilon_0=\mu_0=1\f$, and after
\f$k\f$ completed solves `A_n` holds potential slot \f$k\f$ and `A_nm1` slot \f$k-1\f$. The monitors
and the probe take their time pairing, sample times and component positions from the solver's
`LevelContract` (`core/base/LevelContract.h`). `FourierProbe` and the static functions
`Monitors::gauge` and `Monitors::continuity` also accept fields without a solver. Every call that
reduces data is collective: all ranks of the field layout make it in the same order. Kernels run
on the field execution space; only field halos, reduction scalars and the gathered probe output
cross MPI. Outside the tests, `apps/planar-interface/PlanePulse.h` uses the probe and
`potentialEnergy`, and `opt-in/output/FieldRunOutput.h` uses both energies. Only tests use
`Monitors`.

## Residual monitors (`Monitors.h`, `MonitorExpressions.h`)

`Monitors` brackets one `solve()` and reports how closely the evolved data satisfy the discrete
identities of the vacuum potential update. It checks the discretization and the field
reconstruction; it is not a physical observable.

```cpp
using namespace chdr::solver;
Monitors monitors(solver.A_n);                           // all-periodic layout only
monitors.seedPreviousSource(previous);                   // optional: source of the preceding solve
const MonitorReport& report = monitors.advance(solver);  // beforeSolve, solve(), afterSolve
```

Let \f$k\f$ be the completed-solve count after the bracketed solve, \f$\phi[k]\f$ and
\f$\mathbf A[k]\f$ the scalar and vector parts of `A_n`, and \f$\phi[k-1]\f$, \f$\mathbf A[k-1]\f$
those of `A_nm1`. \f$D_h\f$ is the centred difference \f$D_c\f$ in the collocated reading, or the
compact forward/backward pair in the staggered reading [ryu2016potentialfdtd, Appendix A,
Eqs. (19)-(20)]. \f$C_h\f$ is the curl built from the same differences (forward differences in the
staggered reading), with the Cartesian component order of [schneider2010UnderstandingFdtd,
section 2.5, Eq. (2.24), p. 16]. `MonitorReport` holds one `ResidualMetrics` per residual:

- `divMagnetic_m` is \f$D_h\cdot\mathbf B\f$ of the output magnetic field; `divRawMagnetic_m` is
  \f$D_h\cdot(C_h\mathbf A[k])\f$. The staggered reading uses forward edge-to-vertex differences.
- `faraday_m[d]` is \f$[(C_h\mathbf A[k])_d-(C_h\mathbf A[k-1])_d]/\Delta t+(C_h\mathbf E)_d\f$
  with the output E. The vector-potential difference cancels the raw-curl difference, and the
  scalar gradient has zero discrete curl because coordinate differences commute (deduction from
  [fallahi2020mithra20fullwavesimulation, Eqs. (3.8)-(3.9)]).
- `gauge_m` is
  \f[G_w=\frac{\phi[k]-\phi[k-1]}{\Delta t}+D_h\cdot\bigl(w\mathbf A[k]+(1-w)\mathbf A[k-1]\bigr),\f]
  with \f$w\f$ from `LevelContract::gaugeWeight()`: 1 for R2, 1/2 for R1, 0 for R2Prime. The
  vector term then has the time of the scalar difference, \f$t_k\f$, \f$t_{k-1/2}\f$ or
  \f$t_{k-1}\f$ with \f$t_k=k\Delta t\f$ (deduction from [christlieb2024GaugeConserving, Eq. (29)]
  and the time table of `LevelContract`). The staggered divergence is the backward face-to-centre
  difference [ryu2016potentialfdtd, Appendix A, Eq. (23)].
- `continuity_m` is \f$(\rho_n-\rho_{n-1})/\Delta t+D_h\cdot(w\mathbf J_n+(1-w)\mathbf J_{n-1})\f$,
  where \f$n\f$ and \f$n-1\f$ denote the sources consumed by the bracketed and by the preceding
  solve. It is the source counterpart of \f$G_w\f$, with the same time pairing (deduction from
  [christlieb2024GaugeConserving, Lemma 3.3]). `beforeSolve` evaluates it;
  `continuityAvailable_m` states whether the preceding source is known.
- `gaussDefect_m` and `gauss_m` are described in the next section.

For a residual of an output field (B or E), the scale in the normalization below is the absolute
evaluation of that field's reconstruction formula from the potentials, not the stored field.

| Residual | Collocated reading | Staggered reading (R2, Standard) |
|---|---|---|
| `divMagnetic_m`, `divRawMagnetic_m` | Zero to round-off | Zero to round-off |
| `faraday_m`, raw magnetic fields of adjacent slots | Zero to round-off | Zero to round-off |
| `gauge_m` and `continuity_m` | The gauge residual follows the discrete wave equation, driven by the continuity residual | Same relation with staggered differences |
| `gaussDefect_m` for R2 | Laplacian/composed-gradient difference minus the time difference of the auxiliary gauge | Minus the time difference of the auxiliary gauge |
| `gauss_m` | Zero to round-off | Zero to round-off |

The exact curl/divergence cancellations are deductions from commuting periodic differences
[clemens2001FiniteIntegration, section 2, Eqs. (18)-(24)]. The relation between the gauge and
continuity residuals is the discrete form of [christlieb2024GaugeConserving, Lemma 3.3,
Eq. (38)]; both residuals are reported, not required to vanish. The monitors do not detect every
error. Centred first differences annihilate checkerboard modes, so vanishing derivative-based
residuals do not exclude grid-scale potentials [trefethen1996FiniteDifferenceSpectral,
section 5.3, p. 203], and there is no Ampère residual. A common shift of the charge and current
time labels preserves the continuity identity, so the monitors cannot validate source time labels
(`LevelContract::halfStepLabel`); analytic source comparisons must fix them.

**Normalization.** Each residual is a `monitorDetail::Expression` that carries its value \f$R\f$
and its absolute evaluation \f$S\f$: the same expression with coefficient magnitudes and with every
subtraction replaced by the addition of input magnitudes. Nested stencils propagate this scale;
`monitorDetail::Workspace` exchanges the halos of both parts before the outer stencil. A residual
that cancels to round-off is thus compared with the size of the cancelling terms, not with its
own nearly zero value. `ResidualMetrics` reports \f$\max(\lvert R\rvert/S)\f$
(`cell_m`), \f$\max\lvert R\rvert/\max S\f$ (`global_m`), \f$\max\lvert R\rvert\f$ (`maxAbs_m`),
\f$\max S\f$ (`maxScale_m`), and the Euclidean norms of \f$R\f$ (`l2_m`) and \f$S\f$ (`scaleL2_m`)
without a cell-volume factor. Zero divided by zero is zero; a nonzero residual with zero scale, or
a nonfinite residual, gives an infinite ratio. The scale is a deduction from the absolute error
scale of [jeannerod2013ImprovedErrorBounds, section 1, Eq. (1.6)]. It is a cancellation-free
diagnostic scale, not an error bound: Theorem 4.2 there covers inner products with separately
rounded operations, not nested expressions that the compiler contracts to fused multiply-add.
Only owned cells enter the reductions; maxima use `MPI_MAX` and squared norms `MPI_SUM`, on plain
doubles. A maximum is exact in any reduction order, so the four maxima do not depend on the
decomposition; the last bits of `l2_m` and `scaleL2_m` can (deduction).

**Contract.**

- The constructor takes a four-component field on the solver's mesh and layout. It throws
  `std::invalid_argument` for a layout that is not periodic on every axis, before it allocates
  work arrays. The identities rest on periodic index shifts, which commute. Next to a Mur face
  the stencils would read the outermost Mur layer and the exterior ghosts, which hold a one-sided
  extrapolation (`core/solver/README.md`, section "Periodic and Mur boundaries"). A correct
  boundary residual needs more than that closure: a scale rule for the closed values, the
  boundary terms of the closed operators and a separate boundary reduction. The record test
  `OpenGaugeRecords` (`chdr.solver.open.record.gauge`) evaluates gauge and Gauss residuals near
  Mur faces with its own normalization.
- `beforeSolve` throws `std::logic_error` for a material solver whose `homogeneousPermittivity()`
  or `formulationParameters().alpha_m` differs from one; a nonuniform map reports a NaN
  permittivity and is refused too. The identities are those of the vacuum update. In a
  homogeneous medium the gauge condition and Gauss's law are
  \f$\nabla\cdot\mathbf A+\varepsilon\partial_t\phi=0\f$ and
  \f$\nabla\cdot(\varepsilon\mathbf E)=\rho\f$ [chew2014generalizedgauge, section 2, Eqs. (4),
  (11)], so the vacuum residuals of an exact medium solution do not vanish. At
  \f$\varepsilon=1\f$ the generalized gauge
  \f$\nabla\cdot(\varepsilon\mathbf A)=-\chi\partial_t\phi\f$ with \f$\chi=\alpha\varepsilon^2\f$
  [chew2014generalizedgauge, Eqs. (15), (18)] becomes
  \f$\nabla\cdot\mathbf A+\alpha\partial_t\phi=0\f$, and the O1-VS update adds
  \f$-(1-1/\alpha)D_d^+(D^-\cdot\mathbf A)\f$ to the effective source of vector component \f$d\f$
  (deduction from `O1VSOperator::curlCorrection` and `core/formulation/O1VSUpdate.h`). With
  \f$\varepsilon=1\f$ and \f$\alpha=1\f$ these forms reduce to the vacuum forms, so such a
  material solver is accepted; `O1VS.VacuumBits` shows that
  the O1-VS update then reproduces the vacuum potentials and fields bit for bit. The refusal
  occurs before any input is captured; it rejects the diagnostic, not the solver. Medium tests
  evaluate their own gauge and Gauss residuals (`tests/core/solver/WaveNorms.h`).
- Write the source and make the potential halos valid before `beforeSolve`. It copies `A_nm1`,
  the output E and the source; at count zero it first calls `reconstructEB()`. `afterSolve`
  evaluates the step and keeps its source for the next continuity residual; `advance` calls
  `beforeSolve`, `solve()` and `afterSolve`. `evaluate` re-evaluates the retained step from
  `consumedSource()`, so a later write to the live source does not change the result.
- `afterSolve` and `evaluate` throw `std::logic_error` unless the solver is the same object, with
  the same history revision and time step, exactly one solve later. A second `beforeSolve`
  without `afterSolve` also throws. Monitoring can start at any count and can skip solves.
  Continuity needs the source of the immediately preceding solve of the same solver, history
  revision and time step; otherwise `continuityAvailable_m` is false, which is not a measured
  zero. `seedPreviousSource(source)` declares that source for the next step only.
- All readings that `ChdrSolver` accepts are supported: R1, R2 and R2Prime with the Standard or
  NonStandard stencil in the collocated reading, and R2 with Standard in the staggered reading.
  The formulas assume the inherited vacuum step and sources written before `solve()`.

**Cost and output.** A `Monitors` object allocates 23 doubles per cell, halos included: five
four-component fields (the `A_nm1` copy, two source copies, the value and scale arrays of the
workspace) and one three-component copy of E. A monitored step with continuity available makes
four full-field device copies, 14 halo fills and nine residual reductions with two
`MPI_Allreduce` calls each; without continuity it makes two halo fills and one reduction fewer.
No full-field host mirror is made. `afterSolve`, `advance` and `evaluate` return the
`MonitorReport`, identical on every rank; `report()` returns the last one. `Monitors` writes no
file.

## Gauss formula check (`MonitorRecords.h`)

`evaluateGauss(solver, oldestA, consumedSource, previousElectric, workspace)` returns a
`GaussReport`. `Monitors::afterSolve` stores its `defect_m` as `gaussDefect_m` and its
`identity_m` as `gauss_m`.

- `defect_m` is the physical defect \f$D_h\cdot\mathbf E-\rho\f$, with \f$\rho\f$ from the retained
  source of the completed solve. R2 and R1 use the new output E. R2Prime uses the E from before
  the solve, because its forward scalar difference needs the newly evolved scalar slot; its Gauss
  report is therefore one solve late.
- `identity_m` is the defect minus its discrete prediction. For collocated R2 after solve
  \f$k\f$, substitution of the scalar wave update gives
  \f[D_c\cdot\mathbf E[k]-\rho_s[k-1]=L_h\phi[k-1]-D_c\cdot D_c\phi[k-1]
     -\frac{\Gamma_0[k]-\Gamma_0[k-1]}{\Delta t},\f]
  where \f$\rho_s[j]\f$ is the charge of the source that advances slot \f$j\f$ to \f$j+1\f$,
  \f$\Gamma_0[k]=(\phi[k]-\phi[k-1])/\Delta t+D_c\cdot\mathbf A[k]\f$ is the auxiliary backward
  gauge, and \f$L_h\f$ is the wave Laplacian of the solver's stencil (deduction from
  [christlieb2024GaugeConserving, Theorem 3.5] and [fallahi2020mithra20fullwavesimulation,
  Eqs. (3.6)-(3.9)]). R1 takes its scalar gradient from the current slot, so its double gradient
  is \f$D_c\cdot D_c\phi[k]\f$ while the wave-operator term keeps \f$\phi[k-1]\f$; R1 uses
  \f$\Gamma_0\f$, not its reported mean gauge \f$G_{1/2}\f$. R2Prime pairs each scalar difference
  with the preceding vector slot. In the staggered reading the backward divergence of the forward
  gradient equals the Standard Laplacian, so the prediction is
  \f$-(\Gamma_0[k]-\Gamma_0[k-1])/\Delta t\f$ alone, with the backward divergence in
  \f$\Gamma_0\f$ (deduction by substitution from [ryu2016potentialfdtd, Appendix A,
  Eqs. (19)-(23)]).

A small `identity_m` shows agreement with the formula; it does not imply a small `defect_m`. The
collocated defect contains the mismatch \f$L_h-D_c\cdot D_c\f$, whose symbol `core/solver/README.md`
gives; its leading smooth-mode term is of order \f$h^2\f$, and it does not vanish at grid scale.
In the staggered reading the defect equals \f$-(\Gamma_0[k]-\Gamma_0[k-1])/\Delta t\f$ up to
`identity_m`, and \f$\Gamma_0\f$ is the reported R2 gauge residual. The defect therefore stays zero
in a sourced staggered run when the histories start with zero discrete gauge residual and the
sources satisfy the matching discrete continuity equation (deduction from
[christlieb2024GaugeConserving, Theorem 3.3 and Theorem 3.5, Eq. (43)]). The check is part of the
`Monitors` cost: one E halo fill, one gradient fill with two halo exchanges and two reductions.

## Fourier probe (`FourierProbe.h`)

`FourierProbe` accumulates a running time transform at each cell of an inclusive global-index box
(`ProbeBox`: point, line, plane or slab) for selected channels and frequencies, without storing a
time history. `apps/planar-interface/PlanePulse.h` computes reflection and transmission
coefficients from two point probes.

```cpp
const ProbeBox point{{ix, iy, iz}, {ix, iy, iz}};
SolverFourierProbe probe(solver, point, {{Quantity::Electric, 2}}, {nu}, 0, N);
probe.sample();                             // count 0: reconstructs the initial E and B first
for (std::int64_t n = 1; n < N; ++n) {
    solver.solve();
    probe.sample();                         // completed-solve count n
}
const ProbeRecord record = probe.gather();  // collective; values on rank 0 only
```

**Transform.** `LevelContract` gives channel \f$c\f$ the half-step label \f$m_c\f$: the sample
\f$f_c[k]\f$ represents the time \f$(k+m_c/2)\Delta t\f$. With \f$\nu\f$ in cycles per step and
\f$\omega=2\pi\nu/\Delta t\f$, the probe computes
\f[\widehat f_c=\Delta t\sum_{n=0}^{N-1}w_n\,f_c[k_0+n]\,e^{2\pi i\nu(k_0+n+m_c/2)}\f]
over the half-open gate \f$k_0\le k<k_0+N\f$ [oskooi2010meep, section 5.1, Eq. (6);
oskooi2013sources, section 4.3, Eq. (4.6), with its normalization converted to forward factor
one]. There is no division by the sample count or the window sum and no area or volume factor.
`ProbeRecord::CONVENTION` records the full convention: time kernel \f$e^{+i\omega t}\f$, spatial
kernel \f$e^{-i\mathbf k\cdot\mathbf x}\f$, forward factor one, inverse factor \f$1/(2\pi)\f$ per
variable. The probe transforms time only.

The positive time exponential is the convention of both cited sources. For a real field the
transform with the opposite sign is the complex conjugate, so spectra of the two conventions must
not be combined without conjugation. The inverse transforms have the kernels \f$e^{-i\omega t}\f$
and \f$e^{i\mathbf k\cdot\mathbf x}\f$; a real plane wave therefore enters the spectrum at
\f$\omega>0\f$ through its component \f$e^{i(\mathbf k\cdot\mathbf x-\omega t)}\f$ (deduction). The
factor \f$\Delta t\f$ is the rectangle rule: \f$\widehat f_c\f$ is the discrete-time Fourier
transform of the gated samples. For a sampled harmonic signal the sum is a finite geometric
series with a closed form, which `test_probe_full.cpp` uses as its independent reference. The
trapezoid rule would differ only by half weights on the first and last samples, which are
negligible when the gate contains the whole pulse (deduction).

**Windows.** `ProbeWindow::Rectangle` has \f$w_n=1\f$. Use it for a pulse that the gate contains
completely: any other window weights the pulse samples unequally and changes the transform of the
pulse (deduction). `ProbeWindow::Hann` tapers both gate ends and so reduces the spectral leakage
of a signal that does not decay inside the gate [heinzel2002SpectrumSpectralDensity, section 8].
It is the periodic Hann window on gate-local indices, \f$w_n=[1-\cos(2\pi n/N)]/2\f$
[heinzel2002SpectrumSpectralDensity, section 8.1, Eq. (17)]. Its cosine completes one period over
the gate, so it has the cosine-sum form of [heinzel2002SpectrumSpectralDensity, Appendix B,
Eq. (48)]: for a harmonic signal the Hann sum is \f$\tfrac12D(\mu)-\tfrac14[D(\mu+1/N)+D(\mu-1/N)]\f$,
with \f$D\f$ the rectangle sum at frequency offset \f$\mu\f$, and the window sums have the closed
forms below (deduction). The record stores \f$S_1=\sum_nw_n\f$ and \f$S_2=\sum_nw_n^2\f$
[heinzel2002SpectrumSpectralDensity, Eqs. (19), (20)]. For Hann, \f$S_1=N/2\f$ when \f$N\ge2\f$
and \f$S_2=3N/8\f$ when \f$N\ge3\f$; \f$N=1\f$ gives \f$S_1=S_2=0\f$, and \f$N=2\f$ gives
\f$S_2=1\f$.

**Phase.** The host computes each phase from the integer half-step index \f$j=2k+m_c\f$: the
rounded product \f$p=(\nu/2)j\f$, its remainder \f$r=\mathrm{fma}(\nu/2,j,-p)\f$, and the argument
\f$2\pi[(p-\mathrm{nint}(p))+r]\f$, where \f$\mathrm{nint}\f$ is the nearest integer. The phase is
periodic in cycles [oskooi2010meep, Eq. (6)]; the construction is a deduction under the
round-to-nearest model [jeannerod2013ImprovedErrorBounds, section 1, Eq. (1.1)]. Each sample is
phased from its integer index, never from an accumulated time or a phase recurrence, so the
argument error does not grow with the step count (deduction). The gate end \f$k_0+N\f$ is at most
\f$2^{51}-2\f$, so the half-step index is exact in binary64.

**Channels.** A `ProbeChannel` names a `Quantity` and a Cartesian component: `Electric`,
`Magnetic` (output B), `RawMagnetic`, `Phi` (component zero) or `VectorPotential`. The probe does
not interpolate: it reads each channel at its own staggered position and records per-channel
coordinates \f$x_d=x_{0,d}+(i_d+1/2+o_d)h_d\f$, with \f$o_d\f$ from `LevelContract::positionOffset`
[ryu2016potentialfdtd, Appendix A, Eqs. (19)-(24)]. R1 electric output has no single time label
for a time-dependent scalar potential, so an R1 electric channel requires the explicit
`phiIsZero`.

**Contract.**

- The constructor throws `std::invalid_argument` for an invalid gate (nonfinite or nonpositive
  \f$\Delta t\f$, negative first step, no samples, gate end beyond \f$2^{51}-2\f$), an unknown
  window, no frequencies, a frequency outside \f$(0,1/2)\f$ or with nonfinite \f$\omega\f$, no
  channels, a `Charge` or `Current` channel, an invalid component, an R1 electric channel without
  `phiIsZero`, a box outside the global domain or with a last index below its first, nonfinite
  geometry or nonpositive spacing, or an output larger than one MPI count. It also refuses E, B
  and raw-B channels in a box that touches the outermost layer of a nonperiodic axis: that layer
  supports the field reconstruction but is not an interior field sample. Potential channels remain
  valid there.
- `sample(step, electric, magnetic, rawMagnetic, currentA)` ignores steps outside the gate. Inside
  it, call it once per step, in order, on every rank; a duplicate or skipped step throws
  `std::logic_error`. `step` is the completed-solve count; at count zero, reconstruct the initial
  E and B first. All fields share the constructor field's layout and halo width.
- `gather()` is collective and requires the complete gate (`std::logic_error` otherwise). A cell
  without exactly one owner makes every rank throw `std::runtime_error`.

`SolverFourierProbe` binds a probe to one solver and takes \f$\Delta t\f$, the `LevelContract` and
the history revision from it. Construct it after the time step is selected and the histories are
loaded; it throws `std::invalid_argument` for a gate that starts before `completedSolves()`. Call
its parameterless `sample()` at every completed-solve count in the gate, including count zero
before the first solve, where it calls `evaluateInitialFields()`. For a raw-B channel it allocates one
three-component field and calls `computeRawMagneticField()` in each sample; the output B is not
changed. After a history reload, `initialize()` or a time-step change, `sample()` throws
`std::logic_error`. The solver and its mesh and layout must outlive the binding.

**Cost and output.** For \f$P\f$ global box cells (\f$P_{\mathrm{local}}\f$ owned), \f$C\f$ channels
and \f$F\f$ frequencies, each rank stores \f$P\f$ global indices and \f$PC\f$ coordinates on the
host; the device holds \f$2P_{\mathrm{local}}CF\f$ accumulator doubles and a \f$2CF\f$ phase table.
Each sample waits for the preceding sample kernel, evaluates the phases on the host, copies them
to the device and launches one kernel over \f$P_{\mathrm{local}}C\f$ work items with \f$2F\f$
multiply-adds each, without MPI. `gather()` allocates \f$2PCF\f$ doubles on each rank and another
\f$2PCF\f$ on rank zero, and calls `MPI_Reduce` twice and `MPI_Bcast` once; budget these buffers
for a large plane. Each cell's temporal sum comes from its single owner, and the other ranks
contribute zero. Adding zeros is exact, so the gathered value equals the owner's sum for every
rank count (deduction). Real and imaginary parts travel as `MPI_DOUBLE`, because IPPL's type table
maps `Kokkos::complex<double>` to `MPI_CXX_FLOAT_COMPLEX` (`src/Communicate/DataTypes.h`), a type
of half the size. `gather()` returns a `ProbeRecord` whose values are complete on communicator
rank zero only; other ranks receive the metadata with empty `values_m` and `ownership_m`.
`values_m` is ordered cell, channel, frequency, real/imaginary part, with cells in x, y, z order
and z contiguous; `value(cell, channel, frequency)` returns one complex value. The metadata holds
the gate, \f$\Delta t\f$, the window and its sums, the `LevelContract`, `phiIsZero_m`, origin,
spacing, the frequencies in cycles per step and the angular frequencies \f$2\pi\nu/\Delta t\f$,
the channel labels and offsets, global indices, coordinates and ownership counts. The probe
writes no file.

## Potential invariants (`FieldEnergy.h`)

`potentialEnergy(solver, owner)` returns an `EnergyRecord` with the vector invariant \f$W_A\f$
(`vector_m`) and the scalar invariant \f$W_\Phi\f$ (`scalar_m`) of the O1-VS update. With \f$a\f$,
\f$b\f$ the vector parts and \f$\phi\f$, \f$\phi'\f$ the scalar parts of `A_n` and `A_nm1`,
\f[W_A=\frac{\Delta V}{2}\sum_{\mathrm{owned}}\Bigl[\sum_d\varepsilon_d\Bigl(\frac{a_d-b_d}{\Delta t}\Bigr)^2
   +\sum_d(C_ha)_d(C_hb)_d+\frac{g(a)\,g(b)}{\chi}\Bigr],\f]
\f[W_\Phi=\frac{\Delta V}{2}\sum_{\mathrm{owned}}\Bigl[\chi\Bigl(\frac{\phi-\phi'}{\Delta t}\Bigr)^2
   +\sum_d\varepsilon_d(D_d^+\phi)(D_d^+\phi')\Bigr],\f]
with \f$\Delta V=h_xh_yh_z\f$, the harmonic face permittivity \f$\varepsilon_d\f$, the backward
divergence \f$g=G(\varepsilon a)\f$ of the face-weighted vector potential,
\f$\chi=\alpha\varepsilon^2\f$ with the cell permittivity, the forward differences \f$D_d^+\f$ and
the forward curl \f$C_h\f$ [chew2014generalizedgauge, Eqs. (15)-(18); ryu2016potentialfdtd,
Appendix A]. \f$C_ha\f$ and \f$C_hb\f$ are the raw magnetic fields of the two slots. These are the
cross-level invariants of the leapfrog update [kreiss2002SecondOrderWave, section 2, form before
Eq. (2.10)], specialized to the masses and stiffness of Chew's equations
[chew2014generalizedgauge, Eqs. (16)-(18)]. After solve \f$k\f$, \f$W_A\f$ is labelled
\f$(k-1/2)\Delta t\f$ (`vectorTime_m`) and \f$W_\Phi\f$ \f$k\Delta t\f$ (`scalarTime_m`).

Without sources in a periodic box both invariants are constant for any two-level data; no initial
gauge constraint is assumed. Their drift therefore measures the round-off of the update and of
the sum, also for random histories. With Mur faces the same sums are diagnostics, not conserved
quantities: they include the boundary samples but no flux term. \f$W_A\f$ and \f$W_\Phi\f$ are
not in general the field energy: \f$W_A\f$ pairs the magnetic fields of two adjacent slots,
contains the gauge term \f$g(a)g(b)/\chi\f$, and its first term omits the scalar-gradient part of
E. With \f$\phi=0\f$ and zero gauge residual, \f$W_A\f$ equals the Yee electric energy plus the
cross-level magnetic product [sakkaplangkul2021YeeDispersiveMedia, section 4.1, Theorem 4.1,
Remark 4.2, Eq. (56)]; `PlaneInterface.Normal` checks this equality to \f$10^{-12}\f$ against an
independent sum. `apps/planar-interface/PlanePulse.h` uses \f$W_A\f$ as the positive reference of
its drift record.

Each `EnergyParts` splits its sum into `vacuum_m` and `dielectric_m`; `total()` adds them. A term
at the scalar position is dielectric when its cell permittivity differs from one, a face term
when \f$\varepsilon_d\ne1\f$, and an edge term when any of the four adjoining cells is
dielectric; both interfaces of a dielectric slab therefore count as dielectric. `absolute_m` sums
the magnitudes of the elementary terms: a cancellation scale, not an energy.

The rounding-error bound of a plain sum grows with the number of terms times the sum of their
magnitudes [jeannerod2013ImprovedErrorBounds, section 1, Eq. (1.5)]; for many cells it can exceed
the drift that a conservation check must resolve. All sums are therefore compensated: TwoSum keeps
the rounding residual of each addition [ogita2005AccurateSumDot, Algorithm 3.1]. MPI may change
the evaluation order of a reduction and advises gathering the operands when the order matters
[mpi2023Standard41, section 6.9.1, p. 226]. `MPI_Allgather` therefore collects twelve doubles per
rank, and every rank merges them in increasing rank order, so the record is identical on every
rank. Compensation reduces, but does not remove, the dependence on the partition.

`potentialEnergy` throws `std::invalid_argument` unless the solver uses `MaterialFormulation::O1VS`,
which requires the Standard stencil, the staggered reading and R2. `owner` must be the
`MaterialOwner` passed to the solver constructor (not checked), and both potential slots need
valid halos. The O1-VS update derives the face and inverse-square coefficient fields from that
owner at construction, and the owner caches them; with the same owner the call allocates no field.
It launches one device reduction over owned cells, copies no field to the host and makes one
`MPI_Allgather`.

## Aligned field energy (`AlignedFieldEnergy.h`)

`alignedFieldEnergy(solver, owner)` approximates the electromagnetic energy at one time,
\f$\tfrac12\int(\varepsilon\lvert\mathbf E\rvert^2+\lvert\mathbf B\rvert^2)\,dV\f$ in normalized
units with unit permeability [ellingson2020ElectromagneticsVol2, section 3.1, Eqs. (3.20),
(3.22)], by
\f[W_h=\frac{\Delta V}{2}\sum_{\mathrm{owned},d}\varepsilon_dE_d^2
   +\frac{\Delta V}{2}\sum_{\mathrm{owned},d}\overline B_d^{\,2}.\f]
\f$W_h\f$ is the quantity to compare with a field energy; \f$W_A\f$ is the quantity that the
update conserves. Each square is sampled at its own component position: E on positive faces with
the harmonic face permittivity, \f$\overline{\mathbf B}\f$ on edges. These are separate volume
sums, not an energy density at a common point; the component-wise sums follow the discrete inner
products of [sakkaplangkul2021YeeDispersiveMedia, section 3, Eq. (30)]. "Aligned" refers to time
only: after solve \f$k\f$ the staggered R2 output is
\f$\overline{\mathbf B}=(C_h\mathbf A[k]+C_h\mathbf A[k-1])/2\f$ at the electric-field time
\f$(k-1/2)\Delta t\f$ [fallahi2020mithra20fullwavesimulation, Eqs. (3.50)-(3.52)]. Mur boundary
samples receive the full volume weight, without boundary quadrature weights or flux. \f$W_h\f$ is
not asserted to be conserved.

The magnetic part of \f$W_h\f$ exceeds the cross-level product in \f$W_A\f$ by
\f[\frac{\Delta V}{8}\sum_{\mathrm{owned},d}\bigl[(C_ha)_d-(C_hb)_d\bigr]^2\ge0\f]
(deduction by expanding the square of the mean; the cross-level form is that of
[sakkaplangkul2021YeeDispersiveMedia, section 4.1, Eq. (56)]). For a resolved field of angular
frequency \f$\omega\f$ the difference is of relative order \f$(\omega\Delta t)^2\f$ (deduction by
Taylor expansion). `Movies.InitialFields` checks \f$W_h\ge W_A\f$ for a TE start state with
\f$\phi=0\f$.

The function throws `std::invalid_argument` unless the solver uses the staggered reading and R2;
it does not require O1-VS. The fields must be reconstructed for the current state, and `owner`
must be the solver's material (not checked). The first request of a face-permittivity field from
`owner` is collective and allocates it; later requests use the cache. The call launches one
plain-sum reduction and one `MPI_Allreduce` and returns the same value on every rank. Its
summation order is not fixed, so the last bits can depend on the decomposition (deduction from
[mpi2023Standard41, section 6.9.1]).

`FieldRunOutput::write` (`opt-in/output/FieldRunOutput.h`) calls both energy functions, and rank
zero appends one row to `energy.csv`. The first line is `# ` and the case metadata; the header is
`step,time,vacuum,dielectric,total,scalar,absolute,field_energy`. `time` is `vectorTime_m`, the
electric-field time \f$(k-1/2)\Delta t\f$ of \f$W_A\f$ and \f$W_h\f$. `vacuum`, `dielectric`,
`total` and `absolute` belong to \f$W_A\f$, `scalar` is \f$W_\Phi\f$ at \f$k\Delta t\f$, and
`field_energy` is \f$W_h\f$. No column includes the flux through open faces. A nonfinite
\f$W_A\f$ or \f$W_h\f$ stops the run with `std::runtime_error`.

## Tests

`tests/opt-in/diagnostics` builds the executable `test-chdr-diagnostics` when `IPPL_ENABLE_CHDR`,
`BUILD_TESTING` and `IPPL_ENABLE_UNIT_TESTS` are on. It also compiles
`apps/vacuum-wave/VacuumSetup.cpp` and includes support headers from `tests/core/solver`. Configure
the build as in `core/solver/README.md`, then:

```sh
cmake --build build --parallel 2 --target test-chdr-diagnostics
ctest --test-dir build --no-tests=error -R '^chdr\.diagnostics\.' --output-on-failure -j1
```

Each suite has a `.serial` entry and, with an MPI launcher, a `.mpi2` entry on two ranks. All
entries use one OpenMP thread per rank, a 60 s timeout and the labels `chdr;diagnostics;unit`,
plus `mpi` for two ranks and `probe;mpi` for `chdr.diagnostics.probe.ranks`. The suites print
their measured maxima; `ctest -V` shows them.

| CTest entry | Test file and filter | Main checks |
|---|---|---|
| `chdr.diagnostics.energy` | `test_interface_energy.cpp`, `InterfaceEnergy.*` (serial: `NormalSplit`) | O1-VS dielectric slab, \f$\varepsilon=2.13\f$, random histories, 50 source-free steps: relative drift of \f$W_A\f$ and \f$W_\Phi\f$ at most \f$10^{-12}\f$; against a run on `MPI_COMM_SELF`, bitwise equal histories and energy parts within \f$10^{-14}\f$ of the initial total; decomposition normal and tangential to the slab |
| `chdr.diagnostics.monitors` | `test_monitors.cpp`, `ChdrMonitors.*` | Collocated random and Nyquist data, R1, R2 and R2Prime, Standard and NonStandard: magnetic divergence and Faraday at most \f$10^{-12}\f$; a planted E error is detected; exact zeros for zero data; planted continuity values; continuity unavailable after gaps and reloads |
| `chdr.diagnostics.monitor-records` | `test_monitor_records.cpp`, `ChdrMonitorRecords.*` | Collocated `gauss_m` at most \f$10^{-12}\f$ for each reconstruction and stencil; result unchanged after a live-source overwrite |
| `chdr.diagnostics.probe` | `test_probe.cpp`, `FourierProbe.*` | Refusals, gate order, FMA phase against an integer reference, Hann weights and sums, collective ownership failure (two ranks only) |
| `chdr.diagnostics.open-refusals` | `test_open_probe.cpp`, `OpenShell.ProbeLayerAndMonitorRefusals` | In both readings, field probes refused on the outermost Mur layer and potentials accepted there; `Monitors` refused on an open layout |
| `chdr.diagnostics.probe_wiring` | `test_probe_full.cpp`, four `FourierProbeT8` tests | Label-shift factor \f$e^{i\pi\nu}\f$; channel labels and coordinates for each supported reading and reconstruction; staggered plane against the closed form; `SolverFourierProbe` on an evolving mode in five solver configurations |
| `chdr.diagnostics.probe_gates` | `test_probe_full.cpp`, `FourierProbeT8.FrozenGatesAndNaNGhosts` | Two planes (rectangle and Hann) and six points over 100000 steps with NaN in every halo: each value within \f$10^{-12}R\f$ of a long-double closed form; a wrong half-step label differs by at least \f$10^{-6}R\f$ |
| `chdr.diagnostics.probe.ranks` | `compare_probe_ranks.py` | Serial against two-rank output of `FrozenGatesAndNaNGhosts` |

\f$R=\Delta t\,S_1\sum_s\lvert a_s\rvert\f$ bounds the transform of the test signal with tone
amplitudes \f$a_s\f$ (deduction from [oskooi2010meep, Eq. (6)] and
[heinzel2002SpectrumSpectralDensity, Eq. (19)]). `chdr.diagnostics.probe.ranks` runs
`FourierProbeT8.FrozenGatesAndNaNGhosts` on one and on two ranks with `CHDR_PROBE_OUTPUT` set.
Rank zero writes `# communicator_size N` and one line per value: probe, cell, channel, frequency,
\f$R\f$, real and imaginary part. The script requires that the producer passes in both runs,
complete coverage (262 cells, 11 channels, 6 frequencies), identical \f$R\f$, finite values and a
difference of at most \f$10^{-14}R\f$; it reports bitwise equality but does not require it.
Without Python3 or an MPI launcher the entry is registered as skipped, with the label
`unavailable`. The cache variable
`CHDR_DIAGNOSTICS_RECORD_ROOT` (default `build/demos/chdr/tests/opt-in/diagnostics/runs`) sets the
record directory; the comparison writes `probe-ranks1.txt`, `probe-ranks2.txt` and their logs to
its subdirectory `probe-ranks`.

Other tests use the diagnostics with the solver:

- `ChdrStaggered` (`chdr.solver.staggered.*`): staggered identities, `gauss_m` and continuity.
- The vacuum convergence suites (`chdr.solver.convergence.*`, `chdr.solver.long.Convergence.*`):
  identities on the last step.
- `RankCases` (`chdr.solver.rank_cases.*`) and `LongGaugeGrowth`
  (`chdr.solver.long.GaugeGrowth.*`): residuals of seam cases, gauge norm over ten periods.
- `MediumGuard` (`chdr.solver.medium.Guard.*`) and `O1VS.Refusals` (`chdr.formulation.o1vs.*`):
  monitor acceptance at \f$\varepsilon=1\f$, \f$\alpha=1\f$ and refusal otherwise.
- `PlaneInterface.Normal` (`chdr.interface.normal`): \f$W_A\f$ against an independent Yee sum.
- `Movies.InitialFields` (`chdr.interface.movies.unit.*`): \f$W_h\ge W_A\f$ for a TE start state.

## Library references

The citation keys used in this README and in the headers of this directory identify these works:

- `fallahi2020mithra20fullwavesimulation`: Arya Fallahi (2020), [MITHRA 2.0: A Full-Wave Simulation Tool for Free Electron Lasers][fallahi].
- `christlieb2024GaugeConserving`: Andrew J. Christlieb, William A. Sands and Stephen R. White (2025), [A Particle-in-Cell Method for Plasmas with a Generalized Momentum Formulation, Part III: A Family of Gauge Conserving Methods][gauge].
- `chew2014generalizedgauge`: Weng Cho Chew (2014), [Vector Potential Electromagnetic Theory with Generalized Gauge for Inhomogeneous Anisotropic Media][chew].
- `ryu2016potentialfdtd`: Christopher J. Ryu, Aiyin Y. Liu, Wei E. I. Sha and Weng Cho Chew (2016), [Finite-Difference Time-Domain Simulation of the Maxwell–Schrödinger System][ryu].
- `clemens2001FiniteIntegration`: Markus Clemens and Thomas Weiland (2001), [Discrete Electromagnetism with the Finite Integration Technique][fit].
- `trefethen1996FiniteDifferenceSpectral`: Lloyd N. Trefethen (1996), [Finite Difference and Spectral Methods for Ordinary and Partial Differential Equations][trefethen].
- `schneider2010UnderstandingFdtd`: John B. Schneider (2010), [Understanding the Finite-Difference Time-Domain Method][schneider]; equation and printed-page references use the author PDF dated 2026-03-22.
- `jeannerod2013ImprovedErrorBounds`: Claude-Pierre Jeannerod and Siegfried M. Rump (2013), [Improved Error Bounds for Inner Products in Floating-Point Arithmetic][roundoff].
- `oskooi2010meep`: Ardavan F. Oskooi et al. (2010), [Meep: A Flexible Free-Software Package for Electromagnetic Simulations by the FDTD Method][meep].
- `oskooi2013sources`: Ardavan Oskooi and Steven G. Johnson (2013), [Electromagnetic Wave Source Conditions](https://arxiv.org/abs/1301.5366), chapter 4 of *Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology*, Artech House.
- `heinzel2002SpectrumSpectralDensity`: Gerhard Heinzel, Albrecht Rüdiger and Roland Schilling (2002), [Spectrum and Spectral Density Estimation by the Discrete Fourier Transform (DFT), Including a Comprehensive List of Window Functions and Some New Flat-Top Windows][window].
- `kreiss2002SecondOrderWave`: Heinz-Otto Kreiss, N. Anders Petersson and Jacob Yström (2002), [Difference Approximations for the Second Order Wave Equation](https://doi.org/10.1137/S0036142901397435), SIAM J. Numer. Anal. 40, 1940–1967.
- `ogita2005AccurateSumDot`: Takeshi Ogita, Siegfried M. Rump and Shin'ichi Oishi (2005), [Accurate Sum and Dot Product](https://doi.org/10.1137/030601818), SIAM J. Sci. Comput. 26, 1955–1988.
- `mpi2023Standard41`: Message Passing Interface Forum (2023), [MPI: A Message-Passing Interface Standard, Version 4.1](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report.pdf).
- `sakkaplangkul2021YeeDispersiveMedia`: Puttha Sakkaplangkul and Vrushali A. Bokil (2021), [Convergence Analysis of Yee-FDTD Schemes for 3D Maxwell's Equations in Linear Dispersive Media](https://doi.org/10.4208/ijnam.v18.n4.p524), Int. J. Numer. Anal. Model. 18, 524–568.
- `ellingson2020ElectromagneticsVol2`: Steven W. Ellingson (2020), [Electromagnetics, Vol. 2](https://doi.org/10.21061/electromagnetics-vol-2), Virginia Tech Publishing.

[fallahi]: https://arxiv.org/abs/2009.13645
[gauge]: https://arxiv.org/abs/2410.18414
[chew]: https://arxiv.org/abs/1406.4780
[ryu]: https://doi.org/10.1109/JMMCT.2016.2605378
[fit]: https://doi.org/10.2528/PIER00080103
[trefethen]: https://people.maths.ox.ac.uk/trefethen/pdetext.html
[schneider]: https://eecs.wsu.edu/~schneidj/ufdtd/
[roundoff]: https://doi.org/10.1137/120894488
[meep]: https://doi.org/10.1016/j.cpc.2009.11.008
[window]: https://hdl.handle.net/11858/00-001M-0000-0013-557A-5
