#import "../../docs/manual/preamble.typ": *

= Experiment E16-S detailed specs <experiment-e16s>

E16-S repeats #link(<experiment-e16>)[E16] with two differences.
+ The time step is the largest one that the solver accepts with a margin of one percent, instead of $h_t = h_z\/2$.
+ A Richardson extrapolation of the two energy fractions from the two finest levels is added (E16-S.meas04).

Everything that is not stated otherwise below is the same as in E16:
the reference (@E16-analytic-potentials, @E16-analytic-EB, @E16-energy-fractions), the profile, the media,
the domain, the position of the face, the refinement levels, the initialization,
the measurements E16.meas01–03 and their analysis.

== Experiment E16-S time step <e16s-time-step>
Let
$
  h_t^"max" = (sum_d h_d^(-2))^(-1\/2)
$
be the stability bound of the Standard stencil
(Fallahi @fallahi2020mithra20fullwavesimulation[Eq. (3.20)], with $c_"vac" = 1$).
`O1VSUpdate::maxStableTimeStep` returns this bound for O1VS.
E16-S uses
$
  h_t = C h_t^"max", quad C = 0.99,
$ <E16S-time-step>
and selects it with `ChdrSolver::setSafetyFactor(0.99)`.
The program does not compute the bound itself.
The checked setter excludes $C = 1$, where the shortest wave on the mesh grows linearly in time.

The transverse cells of E16 are $h_x = h_y = w\/8$, so
$
  h_t / h_z = C / sqrt(1 + 2 (h_z\/h_x)^2).
$
This ratio is not the same on all levels.
It grows with refinement, because the transverse terms of the bound lose weight.

The duration $T = 3w\/v_1$ and the three observation targets are those of E16.
The nearest-time rule of E01 selects the observed steps, so the last completed step is
$
  n_T = floor(T\/h_t + 1).
$
The energies of E16.meas03 are also sampled every $N_"pulse"\/16$ steps, as in E16.

#let eps-dielectric = 2.13
#let e16s-levels = (128, 256, 512, 1024, 2048)
// h_t/h_z with h_z = 2w/N_pulse and h_x = h_y = w/8.
#let step-ratio(cells) = 0.99 / calc.sqrt(1 + 2 * calc.pow((2 / cells) / (1 / 8), 2))
#let last-step(cells) = calc.floor(3 / (step-ratio(cells) * 2 / cells) + 1)
// value rounded to `digits` decimals, with trailing zeros.
#let decimals(value, digits) = {
  let parts = str(calc.round(value, digits: digits)).split(".")
  let fraction = if parts.len() > 1 { parts.at(1) } else { "" }
  parts.at(0) + "." + fraction + "0" * (digits - fraction.len())
}

#block(breakable: false, width: 100%)[#align(center)[
  #table(
    columns: 6,
    align: center,
    inset: 6pt,
    stroke: 0.5pt,
    table.header(
      [$N_"pulse"$], [$h_t\/h_z$], [$n_T$], [$n_T$ of E16],
      [$b_1$], [$b_2$],
    ),
    ..e16s-levels.map(cells => (
      [#cells], [#decimals(step-ratio(cells), 5)], [#last-step(cells)], [#(3 * cells + 1)],
      [#decimals(1 - calc.pow(step-ratio(cells), 2), 3)],
      [#decimals(1 - calc.pow(step-ratio(cells), 2) / eps-dielectric, 3)],
    )).flatten(),
  )
]]
$b_1$ and $b_2$ are the coefficients of the phase error of @E16S-phase-error in vacuum and in the dielectric.
With $h_t = h_z\/2$, E16 has $b_1 = 0.75$ and $b_2 = #decimals(1 - 0.25 / eps-dielectric, 3)$ on every level.
A run of E16-S takes about half the steps of the corresponding run of E16.

== What the time step changes and what it does not change <e16s-what-changes>

*The phase speed changes.*
The fields of E16 do not vary along $x$ and $y$.
For one wave $e^(i(k z - omega t))$ in a medium with the speed $v$,
the update then gives the dispersion relation
$
  sin(omega h_t \/ 2) = (v h_t)/h_z sin(k h_z \/ 2)
$
(Schneider @schneider2010UnderstandingFdtd[§7.4, Eq. (7.36)]).
We expand both sides to fourth order and obtain (deduction)
$
  omega/(v k) - 1 = - (k h_z)^2/24 thin b + cal(O)((k h_z)^4), quad
  b = 1 - (v h_t \/ h_z)^2.
$ <E16S-phase-error>
The table above gives $b$ for both media.
A larger step lowers $b$ in vacuum from $0.75$ to a few hundredths.
In the dielectric $b$ stays above $1 - 1\/epsilon_2$,
because the vacuum limits the time step on the same mesh.
E16.meas01–02 contain this phase error, so we expect them to decrease on every level.

*The reflection at the face does not change.*
For a mode $e^(-i omega n h_t)$ of the update,
the second difference in time acts as a multiplication by $-Omega^2$ with
$
  Omega = 2/h_t sin(omega h_t \/ 2),
$
which is the left side of the dispersion relation of Fallahi
@fallahi2020mithra20fullwavesimulation[Eq. (3.17)].
The update then requires that the mode solves the spatial equations of the semi-discrete update at the frequency $Omega$ (deduction).
The incident wave fixes $Omega$ through the relation of the first medium.
The reflected and the transmitted wave, and their amplitudes, follow from the spatial equations at this $Omega$.
$h_t$ does not appear in them, so the reflection coefficient of one mode does not depend on $h_t$ (deduction).
The R2 reconstruction multiplies the incident and the reflected part of a mode with the same factor.
We therefore expect $delta_1$ and $delta_2$ of E16.meas03 to agree with E16 on each level, to leading order.

== Experiment E16-S measured quantities
E16.meas01–03 are used unchanged.

=== E16-S.meas04: Richardson extrapolation of the energy fractions
Let $rho_S (N)$ be the energy fraction $rho_S^(n-1\/2)$ of E16.meas03 at the final observation on the level $N = N_"pulse"$,
and let $rho_1^"cont" = Gamma^2$ and $rho_2^"cont" = 1 - Gamma^2$ (@E16-energy-fractions).
From the two finest levels we form
$
  rho_S^"R" = (4 rho_S (2048) - rho_S (1024))/3, quad
  delta_S^"R" = rho_S^"R"/rho_S^"cont" - 1, quad S in {1, 2}.
$ <E16S-richardson>
This is the Richardson extrapolation with the refinement ratio $2$ for a second-order quantity
(Slater @slater2021SpatialGridConvergence[section "Richardson Extrapolation"]).
It assumes $rho_S (N) = rho_S^"cont" + g_S h_z^2 + cal(O)(h_z^3)$ with $g_S$ independent of $h_z$,
which the observed order of $delta_S$ in E16 supports.
The $h_z^2$ term then cancels, and $rho_S^"R"$ has the error $cal(O)(h_z^3)$,
or $cal(O)(h_z^4)$ if the odd powers vanish
@slater2021SpatialGridConvergence[section "Richardson Extrapolation"].

The two finest levels have slightly different ratios $h_t\/h_z$ (table above).
Because the reflection does not depend on $h_t$ (@e16s-what-changes), this does not change $g_S$.
As a check we form @E16S-richardson also from the levels $512$ and $1024$,
and we report the observed order of $abs(delta_S^"R")$ between the two pairs.

*Purpose:* the extrapolation gives the most accurate energy fractions that the five levels allow.
It also tests that the error of $rho_S$ is a pure $h_z^2$ term on the two finest levels.
We do not extrapolate the field errors of E16.meas01–02.
A norm of an error does not converge to a value that the extrapolation could estimate.

== Experiment E16-S saved output
The inputs are `inputs/E16-S-n<N>.yaml` with the case `E16-S` and the five levels of E16.
Each run writes `measurements.csv`, `profiles.csv` and `run.json` as in E16.
`run.json` records $C$ and $h_t$.
The inputs do not enable the field frames and the movie profiles of E16.
E16-S.meas04 is formed by the analysis from `measurements.csv` of the levels $512$, $1024$ and $2048$.

== Experiment E16-S analysis
We report the analysis of E16 for E16-S.
For each level we compare with E16 at the same level:
the ratios $e_(F,S)$ of E16 over $e_(F,S)$ of E16-S for $F in {E, B}$ and $S in {1, 2}$,
and $delta_S$ of E16-S over $delta_S$ of E16.
We report $rho_S^"R"$ and $delta_S^"R"$ of E16-S.meas04 against $delta_S (2048)$,
and the order between the two pairs.
