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
#let bamberger = link("https://inria.hal.science/inria-00075909")[Bamberger et al., §II, Eqs. (2.10)–(2.11)]

= Experiment E04 detailed specification\*

E04 measures energy returning to an interior region, energy remaining in the box,
and Gauss/gauge defects when the spherical pulse of #e03 meets an all-Mur boundary.
These observations distinguish interior contamination from energy left near the outer boundary.

== Reference and solver configuration

Retain E03's dipole pulse, medium, direct-history initialization,
O1VS/Standard/R2 configuration, cubic mesh and timestep.
The continuum reference is the single outgoing pulse in unbounded material, restricted to the box.
It is not periodically repeated after boundary contact.
The input parameter `profile_sharpness` selects the radial profile, as in #e03.
The study uses the same value as the E03 study, which is 8.

#table(
  columns: (1fr, 2fr), inset: 5pt, align: left + top, stroke: 0.4pt,
  table.header([*Setting*], [*Value*]),
  [Inherited grid], [$Omega=[-6w,6w]^3$, $N=96,144,216,324$, $h_t=h/2$],
  [Boundary change], [Second-order componentwise Mur on all six faces],
  [Boundary speed], [$v=1/sqrt(2.13)$ for faces, edges and corners],
  [Duration], [$T_4=16w/v approx 23.351w/c_"vac"$],
)

The face condition follows #mur, with edge/corner compatibility conditions from #bamberger.
Their componentwise application to potentials does not establish preservation of Maxwell's constraints (#chew).

=== Boundary implementation

`FDTDBoundary.h` validates a common kernel speed on intersecting open faces.
`AbsorbingBC.h` supplies that speed to the edge and corner updates for all four potential components.
O1VS uses the physical timestep and needs kernel speed $v$.
`MurSpeed::matched()` selects one, so it does not select the required dielectric speed here.
A successful single-axis E02 run does not validate this extension.
The periodic-vacuum `Monitors` interface also needs material/open-boundary diagnostics as defined below.
Check E03's initialization and propagation baseline before interpreting E04.

== Timing and comparison with E03

Write $eta=v t/w$.
The exact pulse occupies
$
  (1+eta)w<r<(3+eta)w.
$ <e04-support>
The event times below are geometric deductions from this support and the cubes.
They describe particular continuum paths, not exact numerical arrival indices.

#table(
  columns: (2.2fr, 1fr), inset: 5pt, stroke: 0.4pt,
  table.header([Event], [$eta$]),
  [E03 final comparison, before boundary contact], [2],
  [Front first reaches a face centre], [3],
  [Front first reaches an outer edge midpoint], [$6sqrt(2)-3 approx 5.485$],
  [Front first reaches an outer corner], [$6sqrt(3)-3 approx 7.392$],
  [Rear has passed every point of the outer box], [$6sqrt(3)-1 approx 9.392$],
  [Rear has cleared the monitor $[-3w,3w]^3$], [$3sqrt(3)-1 approx 4.196$],
  [Earliest face-normal return to the monitor], [6],
  [Late-energy targets], [10, 12, 16],
)
Save the actual Mur-layer positions and sample times, which differ between native components.

#block(breakable: false)[
At initialization and the E03 target $eta=2$, compare E and R2 B with E03 at identical native positions and indices.
Reuse its continuum error rows and save the direct differences
$
  d_(E,2)=norm(vb(E)_"Mur"-vb(E)_"periodic")_(2,h)
    /norm(vb(E)_"ref")_(2,h),
$
and $d_(B,2)$ defined analogously.
]
Discrete tails may reach the boundary before continuum contact.
A difference comparable to the continuum error requires investigation before interpreting absorption.

After contact, the small periodic box has wraparound and is no longer a free-space control.
If needed for diagnosis, use a larger box with identical local spacing, timestep and initialization,
and establish boundary independence by the discrete domain of dependence or a further enlargement comparison.

== Fixed field and energy normalizations

Use E03's componentwise native sample positions and R2 half-time for E and B.
Let $t_0=-h_t/2$ be the initial reconstructed-field time.
Define the reference energy quadrature
$
  U_(0,h)=h^3/2 sum_(bold(i) in cal(M),c)
  (epsilon abs(E_(c,bold(i))^"ref" (t_0))^2+abs(B_(c,bold(i))^"ref" (t_0))^2).
$ <e04-energy-normalization>
Count each owned native sample once across MPI ranks, without ghosts.
Components are squared at their own locations.
Record $U_(0,h)$ once and check its refinement trend.

#block(breakable: false)[
@e04-energy-normalization is a discrete approximation to
$
  U=1/2 integral_Omega (epsilon norm(vb(E))^2+norm(vb(B))^2) dif V,
$
]
the electromagnetic energy in a lossless nondispersive medium (#ellingson, Eqs. (3.20), (3.22)).
This fixed scale remains nonzero after exit.
The startup ratio $(U_h^"num" (t_0))/U_(0,h)$ measures energy bias,
which can be small even when the startup field error is substantial.

== E04.meas01: continuum discrepancy and energy left in the box

At each field observation, define an energy-weighted field error
$
  e_"field" (t)=sqrt(
  (h^3/(2U_(0,h))) sum_(bold(i) in cal(M),c)
  (epsilon abs(E_(c,bold(i))^"num" (t)-E_(c,bold(i))^"ref" (t))^2
  +abs(B_(c,bold(i))^"num" (t)-B_(c,bold(i))^"ref" (t))^2)).
$ <e04-field-error>
Subtract fields before squaring.
This error includes initialization, propagation and boundary effects at finite resolution.
Separate electric/magnetic contributions are optional diagnostics.

Define the numerical electromagnetic energy and save its remaining fraction
$
  U_h (t)=h^3/2 sum_(bold(i) in cal(M),c)
  (epsilon abs(E_(c,bold(i))^"num" (t))^2+abs(B_(c,bold(i))^"num" (t))^2),
  quad F_"remain" (t)=(U_h (t))/U_(0,h).
$ <e04-remaining-energy>
Use R2-averaged B.
`opt-in/diagnostics/AlignedFieldEnergy.h` implements this type of native, time-aligned sum.
It gives every owned sample a full volume weight, including boundary samples.
This is a quadrature estimate, not a discrete invariant or boundary-flux balance.

For $eta>6sqrt(3)-1$, the continuum reference is zero throughout the box,
so $e_"field"^2=F_"remain"$ with the same discrete sample set.
Retain one $F_"remain"$ trace, with $U_(0,h)$ stored separately.
Report $e_"field"$ during contact and exit, without duplicating it after the two quantities become algebraically redundant.
The sampling table below defines the late targets and maximum from the same trace.

*Purpose:* quantify departure from free propagation and the field energy that remains after the exact pulse leaves.
The remaining energy includes reflected radiation, numerical tails and constraint-related fields.
It does not by itself give a reflection coefficient or absorption efficiency.

== E04.meas02: energy returning through an interior surface

Use the fixed monitor cube
$
  V=[-3w,3w]^3, quad Sigma=partial V.
$
The initial pulse fits inside this cube.
The timing table places the return gate $5<=eta<=16$ after continuum clearance and before the first face-normal return.
This monitor measures interior contamination.
For a single specular ray reflected at $x=6w$, the image source is at $x=12w$.
Entry at $x=3w$ requires both tangential slopes to have magnitude at most $1/3$,
so the largest admitted incidence angle is $arctan(sqrt(2)/3) approx 25.24 degree$.
This geometric deduction explains why more oblique first reflections can miss the monitor.
It is not an angular bound for the actual curved field or edge/corner returns.

#block(breakable: false)[
At each face quadrature point, interpolate every component of E and R2 B to the same physical point and the same half-time.
Then compute the Poynting flux density
$
  j_n=(vb(E) times vb(B)) dot vb(n)_Sigma,
$
where $vb(n)_Sigma$ is the outward normal of the monitor cube.
]
Here $mu=1$, so $vb(H)=vb(B)$ (#ellingson, §3.2).

Use midpoint quadrature on each face with spacing $h$ in its two tangent directions,
and weights $h^2$.
Each face has $(N/2)^2$ points at tangent coordinates
$-3w+(j+1/2)h$, $j=0,dots,N/2-1$.
The normal coordinate is exactly $plus.minus 3w$.
All four meshes admit this construction.
Use multilinear interpolation from each component's own native grid,
with positions obtained from `LevelContract`.
Synchronize reconstructed E/B halos before interpolation.
Each face point has one MPI owner and contributes once.

Compute inward and outward powers separately on each monitor face:
$
  P_"in" (t)=sum_(Sigma,"points") h^2 max(-j_n,0),
  quad P_"out" (t)=sum_(Sigma,"points") h^2 max(j_n,0).
$
Split the sign locally before summation,
so simultaneous incoming and outgoing fluxes on different parts of the surface do not cancel.
Sum face contributions to obtain the total powers.
Retaining the outward trace also permits inspection of recrossing.

Save one normalized power history per face and sign,
from which totals, cumulative integrals and the return measure are derived:
$
  p_"in,out" (eta)=frac(w P_"in,out" (w eta/v),v U_(0,h)),
  quad F_"in"=integral_5^16 p_"in" (eta) dif eta.
$ <e04-return-fraction>
Call this the _finite-window inward crossing fraction_.
Use the same $p_"in"$ trace in the disjoint precursor interval $4.5<=eta<5$ and the return gate.
The precursor contains numerical tails and interpolation effects after continuum clearance.
Compare power levels directly, or compare mean rates if using intervals of different duration.
An earlier precursor does not bound later return and must not be subtracted as a fitted background.

*Purpose:* distinguish energy travelling back toward the interior from energy merely remaining near the outer boundary.
Some oblique reflections miss the monitor, and a field can cross it more than once.
Near fields and constraint defects can also contribute to local flux.
Therefore $F_"in"$ is neither the total reflected fraction of the outer box nor an angle-resolved Mur reflection coefficient.
The geometric path estimates do not determine the collected fraction.
Inspect the trace for a late tail at the gate's end.

#block(breakable: false)[
== E04.meas03: boundary and interior constraint defects

Use E03's material-aware residuals
$
  Q_h^(n-1/2)=D_h^- dot (epsilon vb(E)^(n-1/2)),
  quad G_h^n=D_h^- dot vb(A)^n+epsilon (phi^n-phi^(n-1))/h_t.
$
Their continuum values are zero (#chew).
The scalar slots remain at $(n+1/2)h_t$ and $(n-1/2)h_t$.
The gauge residual is labelled at $n h_t$.
Save each raw $L^2$ and maximum norm once, with E03's fixed reference scales.
Normalize subset norms with the same full-domain denominators.
The Gauss/gauge bulk identity in E03 does not justify replacing either measurement on a Mur row.
]

Separate the spatial regions as follows.

- *Interior:* residual centres in $[-5w,5w]^3$.
  This is a fixed physical crop on every mesh.
- *Accessible boundary region:* remaining centres whose stencil inputs are supplied by physical samples or the recorded potential closure.
  Samples across MPI partitions remain physical after halo exchange.
- *Unavailable centres:* omit a residual whose stencil needs unavailable exterior samples, and record the count separately for each residual.

Synchronize E halos before $Q_h$ and potential halos before $G_h$.
The quadratic potential extrapolation in `OpenFaceClosure.h` does not define exterior E/B samples.
An omitted residual is unavailable, not zero.
Owned E and B near a boundary may already depend on extrapolated potentials through reconstruction.
Record that provenance and the sample counts, from which quadrature volumes follow.
A shrinking region can give a falling volume norm without a falling local amplitude, so retain the maximum too.

Use the numerical startup residual as the baseline, not a separately sampled analytic-field residual.
Measure the selected states below first.
Add denser samples, maximum locations or signed maps when they are needed to locate an observed defect.

*Purpose:* identify loss of charge-free Gauss consistency or Lorenz-gauge compatibility.
A reflected Maxwell wave can satisfy both constraints.
Conversely, a large constraint defect can exist even when little energy remains.

== Sampling, output and computational cost

#table(
  columns: (1fr, 2.6fr), inset: 5pt, align: left + top, stroke: 0.4pt,
  table.header([Output], [Cadence and interpretation]),
  [Startup and control], [E/B errors and direct E03 differences at state zero and target $eta=2$.],
  [Field discrepancy], [$e_"field"$ at $eta=3,4,5,6,8$, using the nearest R2 half-time.],
  [Remaining energy], [Every R2 half-time. Linearly interpolate this trace to exact $eta=10,12,16$ and report its maximum over $[10,16]$, including interpolated endpoints.],
  [Monitor powers], [Every R2 half-time from a sample at or before $eta=4.5$ through one at or after $eta=16$. Interpolate to exact interval endpoints and use trapezoidal integration in $eta$.],
  [Gauss and gauge], [State zero and targets $eta=2,3,4,6,8,10,12,16$. Use nearest R2 half-times for Q and integer times for G, retaining actual labels.],
  [Spatial output], [Add field/error sections or signed residual maps when the scalar measurements need localization. No routine full-volume or raw-B output.],
)
The late maximum is the maximum of the piecewise-linear energy trace over the stated window.
Sparse residual samples support comparisons at those times, not whole-run maxima.
There is no separate late scalar-potential target or relative potential error after pulse exit.

Advance through $ceil(T_4/h_t+1/2)$ to bracket the final R2 endpoint:
375, 561, 842 and 1262 steps on the four meshes.
Run $N=96,144$ with the diagnostics enabled before selecting further refinement.
Preserve unresolved return traces as well as resolved ones.
Use the pilot to assess whether higher resolution will separate interior return from numerical contamination,
and to estimate memory and runtime.
#block(breakable: false)[
At $N=324$, the run requires about $4.29 times 10^10$ owned cell advances by this schedule.
Keep the physical geometry and duration fixed, and report any incomplete ladder.
]

#block(breakable: false)[
For the flux implementation, one smooth-reference check of colocation/quadrature and one small MPI comparison are sufficient initial checks.
On a representative run, compare the every-step integral with the integral from every second sample using the same endpoints.
If the difference is material relative to the mesh trend, investigate temporal resolution before drawing an absorption conclusion.
Use E03's run-record checklist, adding the Mur speed and potential-ghost closure.
]

== Analysis and interpretation

For each mesh, tabulate:

- the E03 pre-contact control errors and Mur–periodic differences;
- $F_"in"$, its per-face contributions and the precursor/return normalized-power traces;
- $F_"remain"$ at $eta=10,12,16$ and its maximum on that late interval;
- the initial and selected evolved interior/boundary constraint norms.

Plot the time histories with face, edge and corner encounter times marked.
Plot the fixed-target mesh comparisons, combining related quantities rather than creating one plot for every target.
Report consecutive differences and any apparent nonzero limit.
Do not compute an order of decay to zero for a quantity whose continuum boundary model already reflects.

For orientation, inserting incident and reflected scalar plane waves into Mur's second continuum approximation
(#mur, Eq. (12)) gives the amplitude magnitude
$
  abs(r_2(theta))=((1-cos(theta))/(1+cos(theta)))^2,
$
where $theta$ is the incidence angle from the face normal.
This is a deduction for a single plane face with matched speed.
It predicts zero at normal incidence and nonzero reflection at oblique incidence.
The corresponding planar power ratio is $abs(r_2)^2$.
It is not an exact reference value for E04's broadband curved pulse,
potential reconstruction, or intersecting-face closure.
Consequently mesh refinement can approach a nonzero boundary error even when the discretization is converging.

A resolved return signal with small, decreasing constraint defects supports an interpretation dominated by reflected radiation.
A persistent interior Gauss/gauge defect requires a constraint interpretation as well.
A return trace comparable to the precursor trace is unresolved at the measured precursor level,
not bounded by that earlier signal.
If the finest meshes do not establish a trend, retain the results and identify the missing refinement or targeted control.

E04 supplies an absorption and constraint study for this pulse, dielectric, box and observation window.
It does not establish a universal Mur reflection coefficient,
all-angle accuracy, a rigorous energy balance or long-time stability.
Those claims require additional experiments selected from the actual findings.

== Sources

The exact pulse, derivation, staggered sampling, norm conventions and full Chew reference are in #e03.
The additional primary sources are held in the project literature library:

- `mur1981AbsorbingBoundaryFdtd`: G. Mur, _Absorbing Boundary Conditions for the Finite-Difference Approximation of the Time-Domain Electromagnetic-Field Equations_, IEEE Transactions on Electromagnetic Compatibility 23 (1981), 377–382. #mur.
- `ellingson2020ElectromagneticsVol2`: S. W. Ellingson, _Electromagnetics, Vol. 2_ (2020). Energy density and Poynting flux: #ellingson.
- `bamberger1987SecondOrderCorner`: A. Bamberger, P. Joly and J. E. Roberts, _Second Order Absorbing Boundary Conditions for the Wave Equation: A Solution for the Corner Problem_, INRIA RR-0644 (1987). #bamberger.
