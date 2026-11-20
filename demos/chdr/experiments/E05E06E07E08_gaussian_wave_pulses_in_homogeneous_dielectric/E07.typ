#import "../../docs/manual/preamble.typ": *

#let chew = [Chew @chew2014generalizedgauge[§2, Eqs. (5)–(13)]]

= Experiment E07 detailed specs <experiment-e07>

E07 repeats #link(<experiment-e03>)[E03] with the Gaussian profile of #link(<experiment-e05>)[E05].
Everything that is not stated otherwise below is the same as in E03:
the dipole direction $vb(d)$, the potentials and fields (@e03-dipole-potentials, @e03-dipole-fields),
the box $[-6w, 6w]^3$, the solver configuration, the levels $N = 96, 144, 216, 324$, the time step,
the duration $T_3 = 2w\/v$, the initialization including the scalar-potential times,
the measurements E03.meas01–06, the schedule, the spatial output and the analysis.
E07 has no `profile_sharpness` parameter.

== Experiment E07 radial generator
We keep the generator $psi = f(r - R_0 - v t)\/r$ of E03, with $f = w^2 p$ and $p$ from E05 (@E05-profile).
The derivatives of $f$ are those of E05 (@E05-profile-derivatives), multiplied by $w^2$.
The radial derivatives, potentials and fields of E03 do not change.

In E03 the compact support keeps the center outside the pulse.
A Gaussian $f$ is not zero at $r = 0$, so near the center $psi$ behaves like $f(-R_0 - v t)\/r$.
The function $(f(r - R_0 - v t) - f(-r - R_0 - v t))\/r$ has an odd numerator,
so it is smooth at the center and solves the homogeneous wave equation everywhere
(deduction from the homogeneous wave equations in #chew).
Our generator differs from it by the inward-traveling term $f(-r - R_0 - v t)\/r$,
which is singular at the center.
Thus the reference also contains the field of a point source at the center.

No native sample lies at the center.
Every component position in the table of E03 has at least one coordinate at a cell center,
and with an even $N$ no cell center lies on a coordinate plane through the center.
Thus $r >= h\/2$ at every sample.
At these samples and at the startup times $-h_t <= t <= h_t$,
the inward-traveling term is below $10^(-9)$ of the peak values of $vb(E)$, $vb(B)$, $vb(A)$ and $phi$ on all four meshes.
After a travel of $w\/4$ it is below $10^(-13)$.
We neglect it, as we neglect the periodic copies below.

=== Periodic copies
E03 evaluates only the central copy of the periodic reference.
With the Gaussian the other copies are not exactly zero.
At every selected time $t <= T_3 + h_t\/2$, every other copy has $abs(s) >= 1.9w = 7.6 sigma$ everywhere in the box,
where $p$, $sigma abs(p')$ and $sigma^2 abs(p'')$ are below $2 times 10^(-11)$ of their maxima.
So we evaluate only the central copy, as in E03.
At $T_3$ the effective support of E05 reaches the faces,
so the final comparison is before contact only up to the bound of E05.

== Experiment E07 levels and checks
E03's pilot measured memory and time per step on the same four meshes, so E07 runs all four levels directly.

Before the study we differentiate the potential pair numerically at generic interior points
and compare the result with the field evaluator, as in E03.
The check of the zero continuation in E03 does not apply.

== Experiment E07 video output
For videos, every level also saves the $z=0$ section of $vb(E)$ and R2 $vb(B)$ of E03 every 1, 1, 2 and 2 steps on the four meshes,
including the final state.
We do not save full-domain frames, because two-dimensional videos are sufficient.
