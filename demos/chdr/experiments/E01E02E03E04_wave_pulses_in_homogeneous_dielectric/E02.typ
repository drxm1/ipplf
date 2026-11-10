#import "../../docs/manual/preamble.typ": *

= Experiment E02 detailed specs\* <experiment-e02>

E02 measures the return and residual fields when E01's transverse pulse
leaves through a Mur boundary. Use #link(<experiment-e01>)[E01.typ] for the units,
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
  [Longitudinal ladder], [$N_"pulse"=32,64,128,256,512,1024$;
    $h_z=2w/N_"pulse"$ and $N_z=3N_"pulse"$],
  [Boundaries], [Periodic in $x,y$; Mur at the left ($-z$) and right ($+z$) faces],
  [Mur speed], [$v_"Mur"=v=1/sqrt(epsilon_r) approx 0.68519$ on both faces],
)

Run one 3D E02 calculation per level, six runs in total.
There is no additional 3D control or repeat of the long E01 study.

The compact pulse initially occupies $[w,3w]$ and, after $2w$ of travel,
$[3w,5w]$. Its leading edge reaches the geometric right face at $3w/v$;
the pulse has left by $5w/v$. These times follow from E01's support and speed.

=== Boundary convention

Select #raw("MurSpeed::physical(v)") explicitly: the current O1VS
#raw("MurSpeed::matched()") supplies kernel speed one. Use the existing second
Mur approximation and face recurrence
(Mur (1981) @mur1981AbsorbingBoundaryFdtd[Eqs. (12), (16)]).
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
#link("../../src/MaxwellSolvers/FDTDBoundary.h")[FDTDBoundary.h],
#link("../../src/MaxwellSolvers/AbsorbingBC.h")[AbsorbingBC.h] and
#link("../../src/MaxwellSolvers/OpenFaceClosure.h")[OpenFaceClosure.h].

== Continuum comparison and duration

Use E01's compact profile $p(z-z_0-v t)$ without periodic wrapping,
with the same potential and field formulas. Initialize the two potential
histories at $0,-h_t$ as in E01.

Stop at the R2 sample nearest
$
  T_"end"=7.5w/v, quad
  n_"end"=floor(T_"end"/h_t+1), quad
  t_"end"=(n_"end"-1/2)h_t.
$
The continuum pulse has left by $5w/v$. Its first right-face reflection
reaches the left recurrence plane at $(9w-3h_z)/v$, after the final sample
on all six levels. These times follow from E01's pulse support and the
recurrence planes above; numerical tails need not share that support.

=== Pre-contact comparison

Save E01.meas01–04 and the longitudinal profiles at initialization,
$w/v$ and $2w/v$, using E01's separate nearest-time rules for potentials
and fields. Evaluate the unwrapped continuum reference at the actual times.

At $2w/v$, compare with the archived 8-by-8 E01 profiles at the same native
samples, restricting their owned indices to $0<=k<3N_"pulse"$.
This retains the $B_y$ endpoint at $z=6w$. E01's archived global errors
cover its full box; they are not E02-domain norms. No additional numerical
reference evaluation is required.

== Additional measurements

=== E02.meas01: Residual electric and magnetic fields

Compare the final fields with zero. For $F in {E,B}$, define
$
  d_F=sqrt(
    (sum_(bold(i) in sans(M),c)
      abs(F_(c,bold(i))^"E02" (t_"end"))^2)
    /(sum_(bold(i) in sans(M),c)
      abs(F_(c,bold(i))^"E02" (-h_t/2))^2)
  ).
$
Here $sans(M)$ contains all owned E02 cells, including boundary rows.
Use every vector component at its native offset, R2 $B$, no ghosts and
one contribution per owned sample across MPI ranks. Each denominator is
the corresponding numerical field's initial norm and remains fixed.

This uses E01's discrete $L^2$ convention
(Salari and Knupp (2000)
@salari2000CodeVerificationMms[§4.1, Eq. (15)]). The two ratios measure the total field remaining after
exit, including reflection, numerical tails and initialization artifacts.
They are not isolated reflection coefficients or absorbed-energy fractions.
Use the spatial profiles to distinguish a returning pulse from residue
near a boundary.

=== E02.meas02: Electric probe trace

Record $E_x$ at the native probe
$
  k_m=4w/h_z, quad z_m=4w+h_z/2,
$
on the transverse line $(i,j)=(3,3)$, at every R2 step from initialization
through $n_"end"$. Save the signed value, step, actual position and time
$t_n=(n-1/2)h_t$, using #raw("LevelContract").

The trace records the incident and returning signals for inspection or
later analysis. No window integration, reflection ratio or reference
subtraction is required. Check transverse uniformity and inactive
components in the full fields before interpreting the line as representative.

== Interpretation

Compact support is a property of the continuum pulse. Numerical dispersion
can leave tails after its predicted exit time
(Schneider @schneider2010UnderstandingFdtd[§7.4]).
The residual norms include these tails and effects from both Mur faces.
Report what remains; investigate its origin only if the results require it.

Gauss residuals remain sanity checks. For transverse $E_x (z,t)$ and
$B_y (z,t)$, the corresponding divergences vanish even with reflection
(deduction from Chew
@chew2014generalizedgauge[Eqs. (3)–(4)]). A zero residual does not measure absorption.

The conclusions apply to this homogeneous, normally incident pulse and
the stated duration, not to oblique incidence, intersecting open faces,
general constraint compatibility or long-time stability.

== Output and analysis

Save the following at each resolution:

- E01.meas01–04 and numerical/continuum $A_x,E_x,B_y$ profiles at
  initialization, $w/v$ and $2w/v$.
- The signed E02 probe trace at every R2 step.
- Numerical profiles at $4w/v$, $6w/v$ and $7.5w/v$, using E01's
  quantity-specific nearest-time rules. Show zero as the post-exit reference.
- Full E02 snapshots at initialization, $2w/v$ and those three later
  targets, at the nearest field step: $vb(E)$, R2 $vb(B)$ and separately
  labelled raw $vb(B)$.
- $d_E,d_B$, their initial denominators and actual final time.
  Retain E01's run metadata and add the resolved Mur speed on each face.

Plot the pre-contact comparisons, signed probe trace, post-exit spatial
profiles including the boundary rows, and $d_E,d_B$ against $h_z/w$.
For each positive resolved metric $m$, report
$
  p_"obs"=log_2 (m(h_z)/m(h_z/2))
$
(Salari and Knupp @salari2000CodeVerificationMms[§4.2,
Eq. (25)]). Report plateaus, growth, unresolved values and electric/magnetic
disagreements rather than assuming an order from Mur's local truncation error.
