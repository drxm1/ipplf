#import "@preview/physica:0.9.8": *
#set page(paper: "a4", margin: 22mm, numbering: "1")
#set text(size: 10.5pt)
#set par(justify: true)
#set math.equation(numbering: "(1)")
#show heading: set block(above: 1.3em, below: 0.6em)

#let mur = link("https://doi.org/10.1109/TEMC.1981.303970")[Mur (1981), §§III–V]
#let chew = link("https://arxiv.org/abs/1406.4780")[Chew, §2, Eqs. (1)–(13)]
#let ellingson = link("https://doi.org/10.21061/electromagnetics-vol-2")[Ellingson, §§3.1–3.2]
#let e03 = link("E03.typ")[E03]

= Experiment E04 detailed specification\*

E04 measures what happens when the exact spherical pulses of #e03 reach an all-Mur boundary.
The questions are how much electromagnetic energy returns to an interior region,
how much field energy remains after the outgoing pulse should have left,
and whether the boundary generates Gauss or Lorenz-gauge defects.
These are separate observations.
A decrease in total energy alone does not establish accurate absorption.

We keep E03's two radial polarizations, medium, initialization, cubic geometry and refinement ladder.
Only the boundary conditions and run duration change.
This document specifies measurements and required capabilities;
it does not claim that the all-face dielectric boundary is already implemented or validated.

== Reference and solver configuration

Use the complete TE\_r and TM\_r potentials and fields in #e03,
including the radial near-field terms and the nonzero scalar potential of TM\_r.
The continuum reference is now the single outgoing pulse in unbounded homogeneous material,
restricted to the computational box.
Do not periodically repeat it after boundary contact.
Its exact fields remain known throughout E04 and become exactly zero inside the box after the pulse leaves.

#table(
  columns: (1fr, 2fr), inset: 5pt, align: left + top, stroke: 0.4pt,
  table.header([*Setting*], [*Value*]),
  [Solver], [O1VS, $alpha=1$, Standard stencil, staggered placement, R2],
  [Material and sources], [$epsilon=2.13$, $mu=1$, $rho=0$, $vb(J)=vb(0)$],
  [Units and speed], [$c_"vac"=epsilon_0=mu_0=1$, $v=1/sqrt(epsilon)$],
  [Domain and pulse], [$Omega=[-6w,6w]^3$, $R_0=2w$, $A_0=w=1$],
  [Polarizations], [TE\_r and TM\_r separately, $vb(d)=frac(lr((1,2,3)),sqrt(14))$],
  [Grid], [$N_x=N_y=N_z=96,144,216,324$; $h=12w/N$],
  [Time step], [$h_t=h/2$, selected before loading histories],
  [Boundaries], [Second-order componentwise Mur on all six faces, matched speed $v$],
  [Duration], [$T_4=16w/v approx 23.351w/c_"vac"$],
)

The same homogeneous speed must reach the face, edge and corner updates for all four potential components.
Use the physical speed $v$ with the O1VS timestep;
do not substitute the vacuum speed or adjust the timestep to hide a boundary-speed mismatch.
The detailed discrete closure, including its exterior-ghost rule, belongs in the run record.
Mur's scalar boundary construction is given by #mur;
its application component by component to potentials does not itself prove preservation of the electromagnetic constraints in #chew.
That is one reason for measuring both polarizations and the constraint residuals.

=== Required implementation before the experiment can run

At the time of this specification, `src/MaxwellSolvers/FDTDBoundary.h` rejects non-unit kernel speeds when more than one axis is open.
The edge and corner calls in `src/MaxwellSolvers/AbsorbingBC.h` retain their default unit speed.
Consequently the requested O1VS all-Mur configuration at $epsilon=2.13$ is not yet supported.
The existing matched-speed single-axis path does not establish the all-face case.

Before collecting E04 evidence,
the implementation must support the common dielectric speed consistently on intersecting faces,
retain the existing unit-speed behaviour for its current users,
and check the resulting edge/corner treatment.
This is a prerequisite, not an instruction to bypass validation.
Any change to the boundary closure must be identified in the study record.

The existing generic `opt-in/diagnostics/Monitors.h` also requires a periodic vacuum configuration.
It cannot supply E04's material-aware open-boundary residuals unchanged.
The requirements below define the quantities independently of that interface.
The E03 analytic initializer, its scalar time placement and its continuum refinement study must be checked before interpreting E04.
No solver, monitor or boundary implementation is part of this specification change.

== Timing and comparison with E03

Write $eta=v t/w$.
The exact pulse occupies
$
  (1+eta)w<r<(3+eta)w.
$ <e04-support>
The following event times are geometric deductions from this support and the cube;
they are not statements about the exact numerical domain of dependence.

#table(
  columns: (2.2fr, 1fr), inset: 5pt, stroke: 0.4pt,
  table.header([Event], [$eta$]),
  [E03 final comparison, before boundary contact], [2],
  [Front first reaches a face centre], [3],
  [Front first reaches an outer edge midpoint], [$6sqrt(2)-3 approx 5.485$],
  [Front first reaches an outer corner], [$6sqrt(3)-3 approx 7.392$],
  [Rear has passed every point of the outer box], [$6sqrt(3)-1 approx 9.392$],
  [Late-energy targets], [10, 12, 16],
)
The finite-mesh Mur layers occupy component-dependent native positions near the nominal faces.
Save their actual locations and the actual sample times;
the geometric times above are continuum guides, not exact first-contact indices for every component.

Through $eta=2$, repeat all E03 observations at the identical indices and native positions.
Compare both runs with the continuum and directly with each other.
For example, save
$
  d_(E,2)=norm(vb(E)_"Mur"-vb(E)_"periodic")_(2,h)
    /norm(vb(E)_"ref")_(2,h),
$
and the corresponding magnetic and potential differences at their own times.
This is a control for changes introduced by the boundary implementation or initialization.
Expect a small, refinement-dependent difference before continuum contact;
do not require exact equality, since discrete dispersive tails can reach the boundary earlier.
Investigate a difference comparable to the continuum error before treating the later result as an absorption measurement.

After contact, E03's small periodic box is no longer a free-space reference:
its pulse wraps around.
Use the unbounded analytic solution and the fixed normalizations below.
A larger periodic box is an optional diagnostic only if it is large enough to exclude return during the chosen observation window,
and uses identical local spacing, timestep and initialization.
Establish independence from its boundary by the discrete domain of dependence or a further box-enlargement comparison;
continuum travel time alone is insufficient.
It is not a mandatory additional refinement ladder.

== Fixed field and energy normalizations

Use E03's componentwise native sample positions and R2 half-time for E and B.
Let $t_0=-h_t/2$ be the initial reconstructed-field time.
Define the reference energy quadrature
$
  U_(0,h)=h^3/2 sum_(bold(i) in cal(M),c)
  (epsilon abs(E_(c,bold(i))^"ref" (t_0))^2+abs(B_(c,bold(i))^"ref" (t_0))^2).
$ <e04-energy-normalization>
This is a sum over owned native samples, with no ghosts or duplicated MPI contributions.
Each component is squared at its own location; a common vector position is not required for this volume quadrature.
Use the same rule for both cases and every mesh,
record $U_(0,h)$ and check that it approaches a common value under refinement.
The two continuum cases have equal energy by E03's duality relation.
Their finite-grid quadratures need not agree exactly.

#block(breakable: false)[
@e04-energy-normalization is a discrete approximation to
$
  U=1/2 integral_Omega (epsilon norm(vb(E))^2+norm(vb(B))^2) dif V,
$
]
the electromagnetic energy in a lossless nondispersive medium (#ellingson, Eqs. (3.20), (3.22)).
The fixed initial scale prevents a vanishing denominator after the reference pulse has left.
Also save $(U_h^"num" (t_0))/U_(0,h)$ to show the initial reconstruction error.
Do not silently renormalize each later observation to the numerical energy still in the box.

== E04.meas01: continuum discrepancy and energy left in the box

At each field observation, define an energy-weighted field error
$
  e_"field" (t)=sqrt(
  (h^3/(2U_(0,h))) sum_(bold(i) in cal(M),c)
  (epsilon abs(E_(c,bold(i))^"num" (t)-E_(c,bold(i))^"ref" (t))^2
  +abs(B_(c,bold(i))^"num" (t)-B_(c,bold(i))^"ref" (t))^2)).
$ <e04-field-error>
This is a norm of the field difference.
Subtract fields before squaring; the difference of two total energies is a different observable and contains different interference terms.
The definition uses the energy weights above and the global norms of E03.
At finite resolution it includes initialization and bulk propagation errors as well as boundary effects;
subtracting the analytic field does not algebraically isolate the boundary contribution.
Save the separate electric and magnetic contributions as well.

Also save the numerical electromagnetic energy and its remaining fraction
$
  U_h (t)=h^3/2 sum_(bold(i) in cal(M),c)
  (epsilon abs(E_(c,bold(i))^"num" (t))^2+abs(B_(c,bold(i))^"num" (t))^2),
  quad F_"remain" (t)=(U_h (t))/U_(0,h).
$ <e04-remaining-energy>
Use R2-averaged B, not raw integer-time B.
`opt-in/diagnostics/AlignedFieldEnergy.h` implements this type of native, time-aligned sum.
It gives every owned sample a full volume weight, including boundary samples;
it is a quadrature estimate, not an exactly conserved discrete invariant or a boundary-flux balance.
The potential energies in `FieldEnergy.h` are different quantities.

For $eta>6sqrt(3)-1$, the continuum reference is zero throughout the box,
so $e_"field"^2=F_"remain"$ with the same discrete sample set.
Report $F_"remain"$ at $eta=10,12,16$ and its sampled maximum over $10<=eta<=16$.
For this maximum, use all native R2 samples whose actual times lie in the interval,
and include energies linearly interpolated from bracketing samples to the two endpoints.
Before the exit time, report $e_"field"$ and $F_"remain"$ together.

*Purpose:* quantify departure from free propagation and the field energy that remains after the exact pulse leaves.
The remaining energy includes reflected radiation, numerical tails and constraint-related fields.
It is not a reflected-power coefficient,
and $1-F_"remain"$ is not by itself an absorption efficiency.

== E04.meas02: energy returning through an interior surface

Use the fixed monitor cube
$
  V=[-3w,3w]^3, quad Sigma=partial V.
$
The initial pulse fits inside this cube.
Its outgoing reference field has cleared the entire monitor surface once
$
  eta>3sqrt(3)-1 approx 4.196.
$
The earliest geometrical return from the nearest outer face to the nearest monitor face has $eta=6$:
the leading front travels $3w$ to the outer face and another $3w$ back.
Thus the gate $5<=eta<=16$ begins after the direct reference field has left the monitor,
and before the first geometrical reflected arrival.
This separation is why we use this smaller interior cube rather than a surface immediately next to Mur.

#block(breakable: false)[
At each face quadrature point, interpolate every component of E and R2 B to the same physical point and the same half-time.
Then compute the Poynting flux density
$
  j_n=(vb(E) times vb(B)) dot vb(n)_Sigma,
$
where $vb(n)_Sigma$ is the outward normal of the monitor cube.
]
Here $mu=1$, so $vb(H)=vb(B)$ (#ellingson, §3.2).
Forming this cross product from native array entries at different locations is not an acceptable flux measurement.

Use midpoint quadrature on each face with spacing $h$ in its two tangent directions,
and weights $h^2$.
Each face has $(N/2)^2$ points at tangent coordinates
$-3w+(j+1/2)h$, $j=0,dots,N/2-1$.
The normal coordinate is exactly $plus.minus 3w$.
All four meshes admit this construction.
Use multilinear interpolation from each component's own native grid,
with positions obtained from `LevelContract`.
Internal MPI halos are communication data;
each face point has one owner and contributes once to the reduction.
The monitor is far from exterior ghosts.

Save inward and outward powers separately for every monitor face:
$
  P_"in" (t)=sum_(Sigma,"points") h^2 max(-j_n,0),
  quad P_"out" (t)=sum_(Sigma,"points") h^2 max(j_n,0).
$
The sums may also be combined over all six faces.
Split the sign locally before summation,
so simultaneous incoming and outgoing fluxes on different parts of the surface do not cancel.
The signed net outward power is $P_"out"-P_"in"$.

The primary return measure is
$
  F_"in"=(1/U_(0,h)) integral_(5w/v)^(16w/v) P_"in" (t) dif t.
$ <e04-return-fraction>
Call this the _finite-window inward crossing fraction_.
Also retain the cumulative curve and each face's contribution.
This makes the arrival times and directional differences visible.

Sample powers at every solver step and integrate with the trapezoidal rule.
Retain samples bracketing the exact gate endpoints and interpolate the powers linearly to those endpoints.
The final half-time may therefore require a solver step after $T_4$.
Do not silently integrate from the nearest sample instead of the specified physical endpoint.
Apply the same rule to all meshes.

As a numerical floor, record the maximum inward power and the inward energy over $4.5<=eta<=5.5$,
where the direct continuum field is zero and the geometrical return has not arrived.
This measures dispersive precursors, residual outgoing tails and interpolation effects actually present in the run.
It is not a proof that the numerical field must be zero there.
Do not subtract a fitted floor from @e04-return-fraction;
show it alongside the return signal and examine its refinement trend.

*Purpose:* distinguish energy travelling back toward the interior from energy merely remaining near the outer boundary.
Some oblique reflections miss the monitor, and a field can cross it more than once.
Near fields and constraint defects can also contribute to local flux.
Therefore $F_"in"$ is neither the total reflected fraction of the outer box nor an angle-resolved Mur reflection coefficient.
The fixed gate covers the geometrical return estimates along face normals, edge bisectors and body diagonals.
It does not cover every multiply reflected path;
the time history must show whether a significant late tail remains at its end.

#block(breakable: false)[
== E04.meas03: boundary and interior constraint defects

Use E03's material-aware residuals
$
  Q_h^(n-1/2)=D_h^- dot (epsilon vb(E)^(n-1/2)),
  quad G_h^n=D_h^- dot vb(A)^n+epsilon (phi^n-phi^(n-1))/h_t.
$
Their continuum values are zero (#chew).
The scalar slots remain at $(n+1/2)h_t$ and $(n-1/2)h_t$;
the gauge residual is labelled at $n h_t$.
Save the raw $L^2$ and maximum norms,
the maximum's coordinates,
and the dimensionless $q_2,q_infinity,g_2,g_infinity$ defined in E03 with the fixed initial reference scales.
The definitions apply to a subset by restricting only the numerator's samples;
keep the same full-domain denominators for all regions and times.
]

Separate the spatial regions as follows.

- *Interior:* residual centres in $[-5w,5w]^3$.
  This is a fixed physical crop on every mesh.
- *Boundary region with owned stencils:* remaining residual centres whose complete difference stencil uses only globally owned physical samples.
  A sample across an MPI partition is still a physical sample and belongs here after halo exchange.
- *Exterior-dependent region:* residual centres whose stencil uses at least one exterior ghost value.
  Report these separately, with the actual exterior continuation rule, or mark them unavailable.

Never treat an omitted exterior-dependent residual as zero.
Do not invent analytic exterior values for the numerical field to reduce a boundary defect.
This classification concerns the immediate residual stencil only.
An owned electric sample near an upper face can already depend on exterior scalar ghosts through reconstruction;
an owned stencil does not imply boundary-independent evolution or reconstruction.
Retain that closure provenance in the record.
For each region save its sample count and quadrature volume,
so the changing size of the exterior-dependent set is visible under refinement.
The fixed interior crop distinguishes a defect confined to the boundary from one that propagates into the domain.

Measure the initial residuals and the pre-contact E03 control first.
Direct continuum sampling already gives nonzero discrete residuals for these curved fields;
boundary-generated error must be interpreted against that baseline.
Do not subtract residual norms from each other.
Save representative signed residual maps at $eta=2,4,6,8,10,12,16$ to locate their onset and propagation.
If these maps implicate an edge or corner,
use the component fields and per-face flux histories to guide a focused follow-up;
this one spherical run does not uniquely separate every part of the boundary closure.

*Purpose:* identify loss of charge-free Gauss consistency or Lorenz-gauge compatibility,
especially in TM\_r where the scalar contribution matters.
A reflected Maxwell wave can satisfy both constraints.
Conversely, a large constraint defect can exist even when little energy remains.
Neither diagnostic replaces the field and flux measurements.

== Sampling, output and computational cost

Save the initial state and all E03 pre-contact checkpoints.
For the extended run, use targets
$
  eta in {3,4,5,6,8,10,12,16}.
$
At these targets save field errors, energy, constraints and central-plane component/error sections.
Use E03's nearest-native-time selection separately for A, scalar potential and E/B,
and retain their actual times.
The signed residual maps additionally include the $eta=2$ control.
For the return flux and the late-energy maximum, sample every solver step;
evaluate gauge and Gauss scalar norms at this cadence too, with distinct time labels.
Dense scalar histories are required even when volume snapshots are sparse.

Retain full-domain E and R2 B snapshots at least at $eta=2,6,10,16$ on one documented mesh per case.
On the larger meshes, plane sections plus dense scalar histories are sufficient unless a discrepancy needs a volume inspection.
Raw B may be saved separately with its integer time metadata.
Do not overwrite E03 output or mix TE\_r and TM\_r records.

E04 lasts eight times as long as E03 at the same spacing,
before accounting for the additional reductions and flux interpolation.
Run one coarse case through the complete exit/return interval first,
then estimate memory, runtime and output volume for the full ladder from measurements.
The four-level physical protocol remains unchanged if runs are scheduled in stages.
An unfinished expensive level is a limitation of the evidence, not a reason to substitute anisotropic refinement.

#block(breakable: false)[
For the flux implementation, one smooth-reference check of colocation/quadrature and one small MPI comparison are sufficient initial checks.
On a representative run, compare the every-step integral with the integral from every second sample using the same endpoints.
If the difference is material relative to the mesh trend, investigate temporal resolution before drawing an absorption conclusion.
Preserve inputs, source revision and local modifications, boundary configuration,
rank/thread layout, build, machine, timing and peak memory with all measurements.
]

== Analysis and interpretation

For each polarization and mesh, tabulate:

- the E03 pre-contact control errors and Mur–periodic differences;
- $F_"in"$, its per-face contributions and the pre-return floor;
- $F_"remain"$ at $eta=10,12,16$ and its maximum on that late interval;
- the initial, pre-contact and late maxima of the interior and boundary constraint norms.

Plot the time histories with face, edge and corner encounter times marked.
Show the field/error sections with common physical coordinates and comparable colour scales.
Plot each fixed-target observable against $h$.
Report consecutive differences and any apparent nonzero limit;
do not compute an order of decay to zero for a quantity whose continuum boundary model already reflects.

For orientation, inserting incident and reflected scalar plane waves into Mur's second continuum approximation
(#mur, Eq. (12)) gives the amplitude magnitude
$
  abs(r_2(theta))=((1-cos(theta))/(1+cos(theta)))^2,
$
where $theta$ is the incidence angle from the face normal.
This is a deduction for a single plane face with matched speed.
It predicts zero at normal incidence and nonzero reflection at oblique incidence;
the corresponding planar power ratio is $abs(r_2)^2$.
It is not an exact reference value for E04's broadband curved pulse,
potential reconstruction, or intersecting-face closure.
Consequently mesh refinement can approach a nonzero boundary error even when the discretization is converging.

A resolved return signal with small, decreasing constraint defects supports an interpretation dominated by reflected radiation.
A persistent interior Gauss/gauge defect requires a constraint interpretation as well.
A signal near the pre-return numerical floor is unresolved and should be reported as a bound or sensitivity limit,
not as exact zero reflection.
If the finest meshes do not establish a trend, retain the results and identify the missing refinement or targeted control.
No arbitrary absorption percentage is imposed in advance.

E04 supplies an absorption and constraint study for this pulse, dielectric, box and observation window.
It does not establish a universal Mur reflection coefficient,
all-angle accuracy, a rigorous energy balance or long-time stability.
Those claims require additional experiments selected from the actual findings.

== Sources

The exact pulse, derivation, staggered sampling, norm conventions and full Chew reference are in #e03.
The additional primary sources are held in the project literature library:

- `mur1981AbsorbingBoundaryFdtd`: G. Mur, _Absorbing Boundary Conditions for the Finite-Difference Approximation of the Time-Domain Electromagnetic-Field Equations_, IEEE Transactions on Electromagnetic Compatibility 23 (1981), 377–382. #mur.
- `ellingson2020ElectromagneticsVol2`: S. W. Ellingson, _Electromagnetics, Vol. 2_ (2020). Energy density and Poynting flux: #ellingson.
