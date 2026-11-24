#import "../../docs/manual/preamble.typ": *
#import "obliquePacketValues.typ": *

= Experiment E18-S detailed specs <experiment-e18s>

E18-S repeats #link(<experiment-e18>)[E18] with three differences.
+ The time step of E16-S (@e16s-time-step), $h_t = 0.99 thin h_t^"max"$.
+ A second reference for the reflection magnitude:
  the Fourier sums of E18.meas01 applied to the continuum field (E18-S.meas04).
+ A Richardson extrapolation from the two finest levels (E18-S.meas05).

Everything that is not stated otherwise below is the same as in E18:
the packet and its continuum solution, the domain, the refinement levels, the initialization,
the measurements E18.meas01–03 and their analysis.

== Why a second reference <e18s-why>
In E18 the observed order of $abs(r^"num" (q) - abs(Gamma(q)))$ stays close to two at most of the twelve wavenumbers.
At the largest angle it falls on the last refinements.
This is the behaviour of an error with a part that does not decrease with $h$.

The sums of E18.meas01 extend over the vacuum side only, at two fixed times.
The run ends shortly after the slowest plane waves of the measurement have left the face (E18, "Duration").
Applied to the continuum field, the two sums therefore do not give $abs(Gamma(q))$ exactly (deduction):
the window cuts the packets and their tails, and neighbouring wavenumbers contribute to each sum.
This part belongs to the measurement and not to the update, and it does not decrease with $h$.
E18-S.meas04 measures it, so that the error of the update can be separated from it.

== Experiment E18-S time step and duration
The cells are cubes, so $h_t^"max" = h\/sqrt(3)$ and
$
  h_t = 0.99 thin h \/ sqrt(3) approx #str(calc.round(0.99 / calc.sqrt(3), digits: 4)) thin h
$
on every level. The program selects it with `ChdrSolver::setSafetyFactor(0.99)`, as E16-S.
The ratio $h_t\/h$ is the same on all levels.
The duration is $T = 28 thin L_x\/c_"vac"$ as in E18.
$T\/h_t$ is not an integer, so the last step is the nearest-time rule of E01 for the R2 fields,
$
  n_T = floor(T\/h_t + 1),
$
instead of $56 N_x$.
The energies of E18.meas02 are sampled every $N_x$ steps and at $n_T$.

#let e18s-levels = (24, 48, 96, 192)
#let e18s-last-step(cells) = calc.floor(28 * cells * calc.sqrt(3) / 0.99 + 1)
#align(center)[
  #table(
    columns: 3,
    align: center,
    inset: 6pt,
    stroke: 0.5pt,
    table.header([$N_x$], [$n_T$], [$n_T$ of E18]),
    ..e18s-levels.map(cells => ([#cells], [#e18s-last-step(cells)], [#(56 * cells)])).flatten(),
  )
]

The time step does not change the reflection coefficient of a mode (@e16s-what-changes).
It changes the phase speed and the group speed, and therefore where the packets are at $t_(n_T - 1\/2)$.

== Experiment E18-S measured quantities
E18.meas01–03 are used unchanged.

=== E18-S.meas04: Reflection magnitude of the continuum field through the same sums
Let $E_y^"ref"$ be the continuum field of @E18-analytic-EB, with the expression for $z < z_f$.
With the cells $sans(M)_1$, the positions $(x_bold(i), z_bold(i))$ and the wavenumbers $q_j$ of E18.meas01, define
$
  I^"ref" (q) &= sum_(bold(i) in sans(M)_1) E_y^"ref" (x_bold(i), z_bold(i), -h_t\/2) exp(-i(k_x x_bold(i) + q z_bold(i))), \
  P^"ref" (q) &= sum_(bold(i) in sans(M)_1) E_y^"ref" (x_bold(i), z_bold(i), (n_T - 1\/2) h_t) exp(-i(k_x x_bold(i) - q z_bold(i))), \
  r^"win" (q) &= abs(P^"ref" (q)) / abs(I^"ref" (q)).
$ <E18S-window-reference>
The positions and the two times are those of the numerical samples of E18.meas01,
taken from `LevelContract`.
The program evaluates $E_y^"ref"$ with the continuum reference that it already uses for the initialization and for `profiles.csv`.
This costs one evaluation of the #modes.len() plane waves at each sample of side 1 at two times.

With $r^"win"$ the deviation of the measured magnitude splits into two parts,
$
  r^"num" (q) - abs(Gamma(q)) = underbrace(r^"num" (q) - r^"win" (q), "update") + underbrace(r^"win" (q) - abs(Gamma(q)), "measurement").
$ <E18S-split>

*Purpose:* the first part is the error of the update.
Its observed order is the order of O1VS for the reflection at each angle,
free of the part of @e18s-why.
The second part shows how large that part is.
It depends on the window, the duration and the sample positions, and only weakly on $h$ (deduction).

=== E18-S.meas05: Richardson extrapolation
For each $q_j$, and for the energy fractions $rho_S$ of E18.meas02 at the final step,
we form the extrapolation of E16-S.meas04 from the two finest levels,
$
  r^"R" (q_j) = (4 r^"num"_(192) (q_j) - r^"num"_(96) (q_j))/3, quad
  rho_S^"R" = (4 rho_(S,192) - rho_(S,96))/3,
$ <E18S-richardson>
where the subscript is $N_x$ (Slater @slater2021SpatialGridConvergence[section "Richardson Extrapolation"]).
$r^"R"$ estimates the limit $h -> 0$ of $r^"num"$.
This limit is $r^"win"$ and not $abs(Gamma)$, because the sums are the same.
We compare $r^"R"$ with $r^"win"$ of the finest level and with $abs(Gamma)$,
and $rho_1^"R"$ with $macron(R)$ (@E18-energy-fractions).
As a check we also form @E18S-richardson from the levels $48$ and $96$.

The extrapolation of a magnitude assumes that $r^"num" - r^"win"$ has the same sign on both levels
and decreases with order two.
The analysis reports both conditions for each $q_j$.

*Purpose:* the most accurate reflection magnitudes that the four levels allow,
and a test that the error of the update is a pure $h^2$ term at every angle.

== Expected behaviour
- $r^"num"$ of E18-S differs from $r^"num"$ of E18 on the same level by much less than $abs(r^"num" - r^"win")$,
  because the time step does not enter the reflection of a mode (@e16s-what-changes).
  The observed order of $abs(r^"num" - abs(Gamma))$ at the largest angle therefore stays below two, as in E18.
- $r^"win" - abs(Gamma)$ is almost the same on all levels.
- $abs(r^"num" - r^"win")$ decreases with order two at all twelve wavenumbers, including the largest angle.
- $abs(r^"R" - r^"win")$ on the finest level is much smaller than $abs(r^"num" - r^"win")$.

If the order of $abs(r^"num" - r^"win")$ also falls at the largest angle,
the part that does not decrease with $h$ comes from the update and not from the measurement.

== Experiment E18-S saved output
The inputs are `inputs/E18-S-n<N>.yaml` with the case `E18-S` and the four levels of E18.
`reflection.csv` gets three columns:
`reference_incident_magnitude` ($abs(I^"ref")$), `reference_reflected_magnitude` ($abs(P^"ref")$)
and `window_reference_magnitude` ($r^"win"$).
`measurements.csv`, `modes.csv`, `profiles.csv` and `run.json` are written as in E18,
and `run.json` records $h_t$ and $n_T$.
The inputs do not enable the field frames.
E18-S.meas05 is formed by the analysis from `reflection.csv` and `measurements.csv` of the levels $48$, $96$ and $192$.

== Experiment E18-S analysis
For each $q_j$ we report $abs(Gamma)$, $r^"win" - abs(Gamma)$ on the finest level,
$abs(r^"num" - r^"win")$ on each level with its observed order (Salari and Knupp @salari2000CodeVerificationMms[§4.2, Eq. (25)]),
the observed order of $abs(r^"num" - abs(Gamma))$ for comparison with E18,
and $abs(r^"R" - r^"win")$.
We compare $r^"num"$ with E18 on each level,
and we report $rho_1$, $rho_2$, $rho_S^"R"$ and E18.meas03 as in E18.
