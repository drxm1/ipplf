#import "@preview/physica:0.9.8": *
#set math.equation(numbering: "(1)")
#set page(numbering: "1")

= Experiment E02 detailed specs\*

E02 measures the return and residual fields when E01's transverse pulse
leaves through a Mur boundary. Use #link("E01.typ")[E01.typ] for the units,
material, profile, continuum field definitions, initialization, solver
configuration and native sampling conventions.

== Setup

Keep $A_0=w=1$, $z_0=2w$, $epsilon_r=2.13$, O1VS with $alpha=1$,
Standard/staggered/R2 and $h_t=h_z/2$. The E02 settings are:

#table(
  columns: (1fr, 2.7fr),
  inset: 5pt,
  align: left + top,
  stroke: 0.5pt,
  [Domain], [$L_x=L_y=w$, $L_z=6w$],
  [Transverse mesh], [$N_x=N_y=8$, matching the completed E01 study],
  [Longitudinal ladder], [$N_"pulse"=32,64,128,256,512$;
    $h_z=2w/N_"pulse"$ and $N_z=3N_"pulse"$],
  [Boundaries], [Periodic in $x,y$; Mur at the left ($-z$) and right ($+z$) faces],
  [Mur speed], [$v_"Mur"=v=1/sqrt(epsilon_r) approx 0.68519$ on both faces],
)

Run one 3D E02 calculation per level, five runs in total.
There is no additional 3D control or repeat of the long E01 study.

The compact pulse initially occupies $[w,3w]$ and, after $2w$ of travel,
$[3w,5w]$. Its leading edge reaches the geometric right face at $3w/v$;
the pulse has left by $5w/v$. These times follow from E01's support and speed.

=== Boundary convention

Select #raw("MurSpeed::physical(v)") explicitly: the current O1VS
#raw("MurSpeed::matched()") supplies kernel speed one. Use the existing second
Mur approximation and face recurrence
(#link("https://doi.org/10.1109/TEMC.1981.303970")[Mur (1981), Eqs. (12), (16)]).
The continuum normal-incidence pulse satisfies this boundary condition
without reflection when the speeds match; the discrete return is measured.

For $A_x$, the first and last owned samples have the cell-centre
$z$-coordinates $h_z/2$ and $L_z-h_z/2$. Their $x$ offset is as in E01.
The centred Mur recurrence acts between each boundary sample and its
inward neighbour, so use
$
  z_b=L_z-h_z
$
for the right recurrence plane in the arrival estimates. Keep the existing
quadratic ghost extension, including its effect on reconstructed $B_y$.
These conventions are implemented in
#link("../../../../src/MaxwellSolvers/FDTDBoundary.h")[FDTDBoundary.h],
#link("../../../../src/MaxwellSolvers/AbsorbingBC.h")[AbsorbingBC.h] and
#link("../../../../src/MaxwellSolvers/OpenFaceClosure.h")[OpenFaceClosure.h].

== References and duration

For the continuum comparison, use E01's compact profile
$p(z-z_0-v t)$ without periodic wrapping, with the same potential and field
formulas. Initialize the two potential histories at $0,-h_t$ as in E01.

For the numerical baseline, reuse the saved E01 profiles at matching times.
At additional times, evaluate the existing one-dimensional discrete E01
reference offline on its original periodic grid:
$L_z^"E01"=44w$, $N_z^"E01"=22N_"pulse"$, with the same spacing,
timestep and initial histories. Use its values only at E02's sample
positions in $[0,6w]$; the 3D E02 grid remains $N_z=3N_"pulse"$.
This recurrence follows from
#link("https://arxiv.org/pdf/2009.13645v1")[Fallahi (2020), §3.1.2,
Eqs. (3.11)–(3.15)], using the homogeneous material speed from
#link("https://arxiv.org/pdf/1406.4780v1")[Chew (2014), Eq. (13)].

Reconstruct its $E_x$ and R2 $B_y$ with E01's discrete operators; all other
components are zero. Compare the evaluator with archived E01
#raw("travel_2w") profiles on all five levels and record the discrepancy.
Use field values at matching positions and times, not differences between
scalar error measures or interpolation between sparse saved checkpoints.

Stop at the R2 sample nearest
$
  T_"end"=7.5w/v, quad
  n_"end"=floor(T_"end"/h_t+1), quad
  t_"end"=(n_"end"-1/2)h_t.
$
The baseline's periodic wrap cannot return during this interval:
the shortest wraparound path from $[w,3w]$ into $[0,6w]$ is $39w$.
At the recurrence's maximum influence speed $h_z/h_t=2c_"vac"$,
this takes approximately $19.5w/c_"vac"$, beyond
$T_"end" approx 10.95w/c_"vac"$. The few-cell history/reconstruction
offsets leave ample margin. This is a discrete stencil-support deduction.

=== Pre-contact comparison

Save E01.meas01–04 and the longitudinal profiles at initialization,
$w/v$ and $2w/v$, using E01's separate nearest-time rules for potentials
and fields. Evaluate the unwrapped continuum reference at the actual times.

At $2w/v$, compare with the archived 8-by-8 E01 profiles at the same native
samples, restricting their owned indices to $0<=k<3N_"pulse"$.
This retains the $B_y$ endpoint at $z=6w$. E01's archived global errors
cover its full $44w$ box; they are not E02-domain norms.
After the continuum pulse leaves, use the fixed normalizations below
rather than dividing by its vanishing instantaneous norm.

== Additional measurements

=== E02.meas01: Returned electric pulse

Use the native $E_x$ probe at
$
  k_m=4w/h_z, quad z_m=4w+h_z/2,
$
on the transverse line $(i,j)=(3,3)$. Obtain component positions and times
from #raw("LevelContract"). Check transverse uniformity and inactive components
in the full fields before interpreting this line as representative.

At every R2 step, save E02's trace and the discrete E01 baseline:
$
  I_h (t_n)=E_(x,h)^"E01" (z_m,t_n), quad
  D_h (t_n)=E_(x,h)^"E02" (z_m,t_n)-I_h (t_n),
  quad t_n=(n-1/2)h_t.
$
Here $h$ labels the refinement level. The incident and reflected pulse
centres have the continuum arrival estimates
$
  t_I=(z_m-z_0)/v, quad t_R=(2z_b-z_0-z_m)/v.
$
Use the windows
$
  W_I=[t_I-1.25w/v,t_I+1.25w/v], quad
  W_R=[t_R-1.25w/v,t_R+1.25w/v].
$
They approach $[0.75,3.25]w/v$ and $[4.75,7.25]w/v$ under refinement.
The reflected pulse first reaches the left recurrence plane at
$(9w-3h_z)/v$, after the final observation. These support-based estimates
do not imply compact numerical tails.

Let $Q_h [f;W]$ be the trapezoidal integral of $f^2$ over $W$.
Insert the exact gate endpoints by linear interpolation of the signed
trace before squaring. Use every saved sample, without fitted shifts,
amplitude adjustments, tapers or baseline removal. Report
$
  R_Q=(Q_h [D_h;W_R])/(Q_h [I_h;W_I]), quad r_2=sqrt(R_Q),
$
$
  r_"peak"=(max_(t in W_R) abs(D_h (t)))/
            (max_(t in W_I) abs(I_h (t))).
$
The maxima include the interpolated endpoints. Reject empty gates or zero
incident denominators. $R_Q$ measures the integrated squared return;
$r_"peak"$ also exposes a narrow high-amplitude return.

For isolated continuum plane waves, E01's field relation and
$vb(S)=vb(E) times vb(H)$ give $abs(S_z)=E_x^2/v$
(#link("https://doi.org/10.21061/electromagnetics-vol-2")[Ellingson (2020),
§3.2, Eq. (3.25)]). Thus $R_Q$ has a reflected-fluence interpretation when
the pulses are separated and the initialization artifact discussed below
is negligible. It is not an exact finite-grid energy balance.

==== Window check

Recompute with half-width $1.5w/v$ instead of $1.25w/v$ in both windows,
denoting these wider windows by $W^+$. They remain separated and fit
within the saved interval at every level. Define
$
  g[f;W]=(Q_h [f;W^+]-Q_h [f;W])/(Q_h [f;W^+]).
$
Require $g[I_h;W_I]<0.01$, $g[D_h;W_R]<0.01$ and
$abs(r_"peak"^+-r_"peak")/r_"peak"^+<0.01$ for a resolved return.
This checks window capture, not Mur accuracy. If it fails, inspect the
traces before changing the gates consistently across levels. Do not extend
through a left-face interaction; revise the geometry if necessary.
Retain failed checks and report noise-level returns as unresolved.

=== E02.meas02: Residual electric and magnetic fields

At $t_"end"$, define for $F in {E,B}$
$
  d_F=sqrt(
    (sum_(bold(i) in sans(M),c)
      abs(F_(c,bold(i))^"E02" (t_"end")
        -F_(c,bold(i))^"E01" (t_"end"))^2)
    /(sum_(bold(i) in sans(M),c)
      abs(F_(c,bold(i))^"E01" (-h_t/2))^2)
  ).
$
Here $sans(M)$ contains all owned E02 cells, including boundary rows.
Use every vector component at its native offset, R2 $B$, no ghosts and
one contribution per owned sample across MPI ranks. The one-dimensional
baseline is uniform across the transverse mesh. Its initial norm over
$sans(M)$ remains fixed.

This reuses E01's discrete $L^2$ convention
(#link("https://www.osti.gov/servlets/purl/759450")[Salari and Knupp (2000),
§4.1, Eq. (15)]). It detects near-face residue and unwanted components
missed by the probe. The reflected pulse also contributes: $d_E,d_B$ are
residual-field ratios, not separate reflection coefficients or absorbed
fractions. Use the final difference profiles to locate the residuals.

== Interpretation limits

Continuum-history sampling can introduce a small left-going numerical
component because of discrete dispersion
(deduction from #link("https://eecs.wsu.edu/~schneidj/ufdtd/ufdtd.pdf")[Schneider,
§7.4, Eqs. (7.34)–(7.36)]). Its left-Mur return can overlap the intended
right-face return. E01 subtraction retains this contribution, whose size
has not been quantified. Both measurements therefore include effects of
both Mur faces. Investigate it with a one-dimensional boundary calculation
only if it affects the conclusion; no extra 3D control is required.

Gauss residuals remain sanity checks. For transverse $E_x (z,t)$ and
$B_y (z,t)$, the corresponding divergences vanish even with reflection
(deduction from #link("https://arxiv.org/pdf/1406.4780v1")[Chew,
Eqs. (3)–(4)]). A zero residual does not measure absorption.

The conclusions apply to this homogeneous, normally incident pulse and
the stated duration, not to oblique incidence, intersecting open faces,
general constraint compatibility or long-time stability.

== Output and analysis

Save the following at each resolution:

- E01.meas01–04 and $A_x,E_x,B_y$ profiles at initialization, $w/v$ and
  $2w/v$; include continuum and discrete E01 references.
- E02, discrete E01 and difference probe traces at every R2 step.
- Profiles at $4w/v$, $6w/v$ and $7.5w/v$, with both references, using
  E01's quantity-specific nearest-time rules.
- Full E02 snapshots at initialization, $2w/v$ and those three later
  targets, at the nearest field step: $vb(E)$, R2 $vb(B)$ and separately
  labelled raw $vb(B)$. Save initial/final one-dimensional baseline
  profiles for reproducing $d_E,d_B$; no 3D baseline snapshots are needed.
- $R_Q,r_2,r_"peak",d_E,d_B$, both gate results, endpoints, denominators,
  actual times and reference-evaluator discrepancy. Retain E01's run
  metadata and add the resolved Mur speeds, baseline period and evaluator revision.

Plot the pre-contact comparisons, signed incident/return traces with gates,
final electric/magnetic difference profiles including the boundary rows,
and the refinement curves. Show $r_2,r_"peak",d_E,d_B$ against $h_z/w$;
list the squared-amplitude quantity $R_Q$ separately.

For each positive resolved metric $m$, report
$
  p_"obs"=log_2 (m(h_z)/m(h_z/2))
$
(#link("https://www.osti.gov/servlets/purl/759450")[Salari and Knupp, §4.2,
Eq. (25)]). The slope of $R_Q$ is twice that of $r_2$.
Report plateaus, growth and electric/magnetic disagreements rather than
assuming a convergence order from Mur's local truncation error.
