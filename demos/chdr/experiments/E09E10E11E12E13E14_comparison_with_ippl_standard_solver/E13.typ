#import "../../docs/manual/preamble.typ": *

= Experiment E13 detailed specs <experiment-e13>

E13 repeats #link(<experiment-e07>)[E07] in the dielectric, $epsilon_r = 2.13$, with the four solver configurations of E11 (@E09-configurations, @E11-standard-r2prime).
The O1VS configuration is E07 itself, run again so that the study holds all configurations.
Everything else, including the measurements and the sections, is as in E11.

== Experiment E13 equivalence
In a uniform medium without sources, every potential component of Chew's homogeneous equations satisfies
$nabla^2 psi - epsilon_r partial_t^2 psi = 0$ (Chew @chew2014generalizedgauge[§2, Eqs. (12)–(13)]).
With $tau = h_t\/sqrt(epsilon_r)$ the discrete step in the medium is the vacuum step with the time step $tau$ (`HomogeneousUpdate.h`, a deduction).
The Homogeneous configuration therefore advances the potentials with IPPL's vacuum update and the step $tau$.
The Standard configuration does the same on the collocated layout, and its `evaluate_EB()` divides the difference of $vb(A)$ by the physical step $h_t$.
On a uniform medium O1VS is equivalent to the Standard update with this scaled step (E01–E04 report).
E13 tests this: if O1VS reduces to the scaled update, the O1VS–Homogeneous difference of E11.meas01 is at rounding level.

*Purpose:* show that the O1VS formulation in a homogeneous dielectric equals IPPL's Standard update with the scaled time step,
and compare its accuracy with the collocated Standard solver.
