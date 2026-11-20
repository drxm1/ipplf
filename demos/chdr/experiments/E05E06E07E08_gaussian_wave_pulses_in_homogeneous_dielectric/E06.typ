#import "../../docs/manual/preamble.typ": *

= Experiment E06 detailed specs <experiment-e06>

E06 repeats #link(<experiment-e02>)[E02] with the Gaussian profile of #link(<experiment-e05>)[E05].
Everything that is not stated otherwise below is the same as in E02:
the domain with $L_z = 6w$, the transverse mesh, the six levels,
the Mur boundaries with $v_"Mur" = v$ on both faces and their convention,
the initialization, the duration $T_"end" = 7.5w\/v$, the measurements E02.meas01–02, the output and the analysis.

== Experiment E06 profile and timing
We use the profile of E05 (@E05-profile) without periodic wrapping.
At $t = 0$ the effective support of E05 is $[0, 4w]$.
So there is no time without contact as in E02:
at $t = 0$, the pulse touches the left face with values below $10^(-12)$ of their maxima (see E05).

The times below follow from the effective support, as those of E02 follow from the support.
The front of the effective support reaches the right face at $2w\/v$ (E02: $3w\/v$).
The pulse center crosses this face at $4w\/v$,
and the effective support has left the box at $6w\/v$ (E02: $5w\/v$).
The first right-face reflection of the front reaches the left recurrence plane at $(8w - 3h_z)\/v$ (E02: $(9w - 3h_z)\/v$),
which is after the final sample on all six levels.
At $t_"end"$ the incident pulse inside the box has $abs(s) >= 3.5w = 14 sigma$,
where its values are below $10^(-38)$ of their maxima.

== Experiment E06 pre-contact comparison
We save E01.meas01–04 and the longitudinal profiles at initialization, $w\/v$ and $2w\/v$, as in E02.
At $2w\/v$ the front of the effective support reaches the right face,
so these continuum comparisons are before contact only up to the bound of E05.

E02 also compared its profiles at $2w\/v$ with those of the periodic E01 run.
There the Mur faces changed the propagation before contact only very little,
by $5 times 10^(-4)$ of the peak at $N_"pulse" = 32$ and by $1.5 times 10^(-8)$ at $N_"pulse" = 1024$.
So E06 does not repeat this comparison with E05,
and the E06 study needs no E05 study as a reference.

== Experiment E06 video output
For videos, the levels 32, 128 and 512 also save the full-domain frames of E02 every $2 N_"pulse"\/32$ steps,
about 175 frames per level at the same physical times.
The frames at the targets of E02 are kept.

== Experiment E06 interpretation
The Gaussian has no exact support,
but at $t_"end"$ the continuum pulse inside the box is below $10^(-38)$ of its maximum.
So, as in E02, the residual norms $d_E$ and $d_B$ measure the numerical tails,
the return from both Mur faces and the initialization artifacts.
