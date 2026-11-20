#import "../../docs/manual/preamble.typ": *

= Experiment E14 detailed specs <experiment-e14>

E14 repeats #link(<experiment-e08>)[E08] in the dielectric, $epsilon_r = 2.13$, with the four solver configurations of E11 (@E09-configurations, @E11-standard-r2prime).
The O1VS configuration is E08 itself, run again so that the study holds all configurations.
Everything else, including the measurements, is as in E12, with the speed and the duration of E08.

== Experiment E14 boundary
The Mur speed is the physical speed $v = 1\/sqrt(2.13)$ on all faces, edges and corners.
For the Homogeneous formulation the solver multiplies it by $sqrt(epsilon_r)$ before it passes it to the inherited boundary,
because that boundary sees the step $tau = h_t\/sqrt(epsilon_r)$ (`MurSpeed::kernel`).
The runs advance through 375, 561, 842 and 1262 steps, as in E04.

*Purpose:* show that the update equivalence of E13 also holds with Mur faces, edges and corners,
and compare the energy left in the box with the collocated Standard solver.
