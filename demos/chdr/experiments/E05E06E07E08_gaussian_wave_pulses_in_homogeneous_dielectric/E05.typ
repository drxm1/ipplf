#import "../../docs/manual/preamble.typ": *

= Experiment E05 detailed specs <experiment-e05>

E05 repeats #link(<experiment-e01>)[E01] with a Gaussian potential profile.
Everything that is not stated otherwise below is the same as in E01:
the units, the material, the potentials and fields (@E01-analytic-potentials, @E01-analytic-EB),
the solver configuration, the domain, the transverse mesh, the initialization,
the time steps and comparison times, the measurements E01.meas01–04 and the saved output.

== Experiment E05 profile
We replace the compact profile of E01 by
$
  p(s) = A_0 exp(-s^2\/(2 sigma^2)), quad A_0 = 1, quad sigma = w\/4.
$ <E05-profile>
Differentiation gives
$
  p'(s) = -s/sigma^2 p(s), quad p''(s) = (s^2/sigma^4 - 1/sigma^2) p(s).
$ <E05-profile-derivatives>
The proof in E01 that the potential family solves the homogeneous equations only needs $a in cal(C)^2$,
so it also holds for this profile.

We choose $sigma = w\/4$ because the sharpened profile of E03 with $beta = 8$ is close to $A_0 exp(-8 s^2\/w^2)$ near its center (@e03-profile-sharpness),
which is @E05-profile.

The Gaussian has no compact support, so $w$ is no longer a support half-width.
We keep $w$ as the unit of length of the geometry,
so that the domain, the initial center and the travel distance of E01 do not change.
Where $abs(s) >= 2w = 8 sigma$, the values $p$, $sigma abs(p')$ and $sigma^2 abs(p'')$ are below $10^(-12)$ of their maxima.
We call the interval $abs(s) < 2w$ the _effective support_ of the pulse.
E06–E08 use it where E02–E04 used the support.

=== Periodic extension
We keep the periodic profile $a(s) = sum_(m in ZZ) p(s - m L_z)$ of E01
and its evaluation from the nearest copy, $a(s) = p(r)$.
A Gaussian copy is never exactly zero.
For $abs(r) <= L_z\/2$, every other copy is at a distance of at least $L_z\/2 = 22w = 88 sigma$,
where $p < exp(-3872)$.
This is below the smallest positive double-precision number,
so the nearest-copy evaluation is exact in double precision.
The support test $abs(r) >= w$ of the E01 pseudocode is removed.

=== Independence of the transverse directions
As in E01, $A_x$ depends only on $z$ and $t$.
Thus the fields do not vary in $x$ and $y$,
and the loaded histories satisfy the discrete homogeneous Lorenz gauge as shown for E01.
We keep $N_x = N_y = 8$ on every level.
To check that the results do not depend on the transverse mesh,
we run the coarsest level a second time with $N_x = N_y = 4$.
Its errors and profiles must agree with those of the $8 times 8$ run to rounding.

== Experiment E05 refinement
We use the six levels of E02:

#align(center)[
  #table(
    columns: 3,
    align: center,
    inset: 6pt,
    stroke: 0.5pt,
    table.header([$N_"pulse"$], [$h_z/w$], [$N_z$]),
    [32],  [$1/16$],  [704],
    [64],  [$1/32$],  [1408],
    [128], [$1/64$],  [2816],
    [256], [$1/128$], [5632],
    [512], [$1/256$], [11264],
    [1024],[$1/512$], [22528],
  )
]

$N_"pulse" := 2w\/h_z$ keeps its definition from E01.
For E05 it is the number of cells across $8 sigma$, which is half the width of the effective support.
E01 needed the levels 2048 and 4096 because its profile reached second order only there.
We do not run them for E05.

== Experiment E05 video output
For videos, the levels 32, 128 and 512 save the full-domain frames of E01 every $10 N_"pulse"\/32$ steps instead of every quarter run.
This gives about 190 frames per level at the same physical times.
The other levels keep the quarter-run frames of E01.

== Experiment E05 analysis
For each measurement and each pair of successive levels, we report the observed order
$
  p_"obs" = log_2 (e(h_z)\/e(h_z\/2))
$
(Salari and Knupp @salari2000CodeVerificationMms[§4.2, Eq. (25)]),
at the final time $T$ and after the travel of $2w$.
