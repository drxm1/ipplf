#import "../../docs/manual/preamble.typ": *

= Experiment E08 detailed specs <experiment-e08>

E08 repeats #link(<experiment-e04>)[E04] with the pulse of #link(<experiment-e07>)[E07].
Everything that is not stated otherwise below is the same as in E04:
the second-order Mur condition with speed $v$ on all faces, edges and corners,
the box, the levels, the time step, the duration $T_4 = 16w\/v$, the periodic control (now E07),
the field and energy normalizations, the measurements E04.meas01–03, the monitor cube, the sampling table and the analysis.

== Experiment E08 timing
Write $eta = v t\/w$ as in E04.
The event times of E04 follow from the support $(1 + eta) w < r < (3 + eta) w$.
For E08 we use the effective support of E05,
$
  eta w < r < (4 + eta) w.
$
This gives the following events. The E04 column is repeated for comparison.

#align(center)[
  #table(
    columns: (2.2fr, 1fr, 1fr), inset: 5pt, stroke: 0.4pt,
    table.header([Event], [E04 $eta$], [E08 $eta$]),
    [E03/E07 final comparison], [2], [2],
    [Front first reaches a face center], [3], [2],
    [Front first reaches an outer edge midpoint], [$6 sqrt(2) - 3 approx 5.485$], [$6 sqrt(2) - 4 approx 4.485$],
    [Front first reaches an outer corner], [$6 sqrt(3) - 3 approx 7.392$], [$6 sqrt(3) - 4 approx 6.392$],
    [Rear has passed every point of the outer box], [$6 sqrt(3) - 1 approx 9.392$], [$6 sqrt(3) approx 10.392$],
    [Rear has cleared the monitor $[-3w, 3w]^3$], [$3 sqrt(3) - 1 approx 4.196$], [$3 sqrt(3) approx 5.196$],
    [Earliest face-normal return to the monitor], [6], [5],
    [Late-energy targets], [10, 12, 16], [10, 12, 16],
  )
]

We keep all observation times, the precursor interval $4.5 <= eta < 5$ and the return gate $5 <= eta <= 16$ of E04,
so that the results of the two experiments can be compared directly.
With the Gaussian, three statements of E04 then hold only up to a bound:

- At $eta = 2$ the effective support reaches the face centers,
  so the comparison with E07 is before contact only up to the bound of E05.
- At $eta = 4.5$ the outgoing pulse has $abs(s) >= 1.3w = 5.2 sigma$ on the monitor surface,
  where $p < 1.4 times 10^(-6)$ and $sigma abs(p')$, $sigma^2 abs(p'')$ are below $4 times 10^(-5)$ of their maxima.
  At $eta = 5$ it has $abs(s) >= 1.8w = 7.2 sigma$, where these values are below $3 times 10^(-10)$.
  So the precursor interval contains a small part of the outgoing continuum pulse, in addition to the numerical tails.
  The face-normal return of the front of the effective support reaches the monitor at the start of the return gate.
- At $eta = 10$ the outgoing pulse still touches the outer corners, with $abs(s) >= 1.6w = 6.4 sigma$,
  where $p$, $sigma abs(p')$ and $sigma^2 abs(p'')$ are below $10^(-7)$ of their maxima.
  From $eta = 6 sqrt(3)$ on, the continuum fields in the box are below the bound of E05.
  We treat $e_"field"^2$ and $F_"remain"$ as redundant from there, as E04 does after the exact pulse has left.

== Experiment E08 levels
E04's pilot measured memory and time per step on the same four meshes, so E08 runs all four levels directly,
with 375, 561, 842 and 1262 steps as in E04.

== Experiment E08 video output
For videos, every level saves the $z=0$ section of $vb(E)$ and R2 $vb(B)$ of E03 every 2, 3, 4 and 7 steps on the four meshes,
including the final state.
This gives about 200 sections per mesh at nearly the same physical times.
We do not save full-domain frames, because two-dimensional videos are sufficient.
