#import "../../docs/manual/preamble.typ": *

= Experiment E17 detailed specs <experiment-e17>

E17 repeats #link(<experiment-e16>)[E16] with the two media exchanged.
The pulse starts in the dielectric and meets the face to vacuum at normal incidence.
Everything that is not stated otherwise below is the same as in E16:
the reference (@E16-analytic-potentials, @E16-analytic-EB, @E16-energy-fractions), the profile,
the solver configuration, the domain, the position of the face, the refinement levels,
the initialization, the measurements E16.meas01–03, the saved output and the analysis.

== Experiment E17 media and coefficients
$
  epsilon_1 = 2.13, quad epsilon_2 = 1.
$

#let n-dielectric = calc.sqrt(2.13)
#let gamma-e17 = (n-dielectric - 1) / (n-dielectric + 1)
#let rounded(value, digits: 4) = str(calc.round(value, digits: digits))

With $n_1 = sqrt(2.13)$ and $n_2 = 1$, @E16-coefficients gives
$Gamma approx #rounded(gamma-e17)$ and $tau approx #rounded(1 + gamma-e17)$.
Three things differ from E16.
- The reflected pulse keeps the sign of the incident pulse.
- The transmitted profile is stretched by $n_1\/n_2 approx #rounded(n-dielectric)$.
- The transmitted potential is larger than the incident potential, $tau > 1$.
  The energy fractions $Gamma^2$ and $1 - Gamma^2$ are the same as in E16.

== Experiment E17 duration and the second face
The pulse travels with $v_1 = c_"vac"\/n_1$, so the target duration
$
  T = 3w/v_1 = 3 n_1 w\/c_"vac"
$
is longer than in E16 by the factor $n_1$.
With $h_t = h_z\/2$ a run takes about $3 n_1 N_"pulse"$ steps.
At $t = T$ the reflected pulse occupies $[1.5w, 3.5w]$ as in E16.
The transmitted pulse occupies $z_f + (n_1\/n_2) [0.5w, 2.5w]$,
which is about $[#rounded(4 + 0.5 * n-dielectric, digits: 2) w, #rounded(4 + 2.5 * n-dielectric, digits: 2) w]$.
The continuum fields stay exactly zero within $1.5w$ of the second face during the whole run.

== Experiment E17 resolution
The carrier wavelength is $lambda_c$ in the first medium, which is now the dielectric,
and $lambda_c n_1\/n_2$ in vacuum.
The shortest wavelength of E17 is therefore $lambda_c$, with $lambda_c\/h_z = 16, dots, 256$ cells on the five levels.
The frequency of the carrier is lower than in E16 by the factor $n_1$.

== Experiment E17 initialization
At $t = 0$ and at $t = -h_t$ the pulse lies in the dielectric, $0.5w$ away from the face and $1.5w$ away from the second face.
The loaded histories are those of E01 with $epsilon_r = 2.13$ and the profile of E16.
