// Standalone draft. Compile in this directory with `typst compile E15.typ`;
// it reads no file outside this directory.
// To include E15 in the manual: replace the standalone setup below by
// `#import "../../docs/manual/preamble.typ": *`, delete the bibliography call
// at the end and E15-references.bib, and include this file from
// experiments/experimentDefinitions.typ like E01 to E04.

// --- standalone setup (begin) ---
#import "@preview/physica:0.9.8": *
// Same colour as `proposal` in docs/manual/preamble.typ.
#let proposal(body) = text(fill: rgb("1f4e9c"), body)
#set document(title: [Experiment E15: long-run stability])
#set page(paper: "a4", margin: (x: 24mm, y: 24mm), numbering: "1")
#set text(font: "Libertinus Serif", size: 11pt, lang: "en", fill: rgb("202428"))
#set par(justify: true, leading: 0.6em)
#set heading(numbering: "1.1", supplement: [Section])
#set math.equation(numbering: "(1)")
#show math.equation: set text(font: "New Computer Modern Math")
#show heading: set text(weight: "semibold")
#show raw: set text(font: "DejaVu Sans Mono", size: 0.85em)
// --- standalone setup (end) ---

= Experiment E15 specs <experiment-e15>

#proposal[
All quantities use the normalized solver units $c_"vac" = epsilon_0 = mu_0 = 1$.
]

== Goal and scope

#proposal[
E15 tests whether the O1VS update stays stable and accurate over a long run.
A single plane wave travels through a periodic box for 1000 periods,
and we compare it with the analytic continuum wave at regular times.
The wave leaves through one face and enters through the opposite face without change,
so the continuum reference is known at all times and no boundary closure takes part.

A single error norm mixes two effects that we want to separate.
A change of amplitude must not occur in a stable, non-dissipative scheme.
A phase lag from numerical dispersion is expected, and it grows linearly with time.
E15 therefore measures the amplitude and the phase of the numerical wave relative to the continuum wave,
and everything else that the numerical solution contains.

We use one discrete Fourier mode instead of a compact pulse as in E01.
Under numerical dispersion a pulse spreads,
and after enough transits its continuum error saturates
and no longer shows what happens to amplitude and phase.
A single mode keeps its shape;
the update can only change its amplitude and phase.
Its discrete solution is also known in closed form (@e15-expected).

The potentials of the wave contain a gauge part:
a longitudinal part of $vb(A)$ and a nonzero scalar potential $phi$.
Their contributions to $vb(E)$ cancel, so $vb(E)$ and $vb(B)$ are those of a transverse wave.
The run therefore also uses the scalar-potential update,
the weighted-divergence term of the vector-potential update,
the gradient term of $vb(E)$
and the time placement of $phi$ relative to $vb(A)$.

E15 uses the core solver and the opt-in diagnostics only, and no code of other experiments.
Every measured sample uses its component position and time label from `LevelContract`.
Out of scope are open boundaries, material interfaces, sources and charges.
]

== Analytic continuum reference

#proposal[
The medium is nonmagnetic, lossless, nondispersive and homogeneous with $epsilon_r = 2.13$ on the whole domain.
There is no free charge or current.
The domain is the cube $[0,L]^3$ with $L = 1$, periodic in all three directions.

We use the travelling plane wave
$
  vb(A)(vb(x),t) = (A_0 vu(p) + a_parallel vu(k)) cos theta, quad
  phi(vb(x),t) = v a_parallel cos theta, quad
  theta = vb(k) dot vb(x) - omega t,
$ <e15-wave>
with amplitudes $A_0 = a_parallel = 1$ and
$
  vb(k) = k_c (1,1,1), quad k_c = (2 pi)/L, quad
  vu(k) = 1/sqrt(3) (1,1,1), quad
  vu(p) = 1/sqrt(2) (1,-1,0), quad
  omega = v abs(vb(k)), quad v = 1/sqrt(epsilon_r).
$
The wave travels along the space diagonal, so all three mesh directions take part.
It has one wavelength per box length along each axis and is therefore $L$-periodic in $x$, $y$ and $z$.
Its wavelength is $lambda = 2 pi \/ abs(vb(k)) = L\/sqrt(3)$
and its period $T_p = 2 pi \/ omega = sqrt(epsilon_r) L \/ sqrt(3) approx 0.8426 L$.

The transverse part $A_0 vu(p) cos theta$ carries the electromagnetic wave.
The gauge part $a_parallel vu(k) cos theta$ and $phi$ follow from it by the gauge transformation
$vb(A) -> vb(A) + grad Lambda$, $phi -> phi - partial_t Lambda$
with $Lambda = (a_parallel \/ abs(vb(k))) sin theta$,
which leaves $vb(E)$ and $vb(B)$ unchanged (deduction from Chew @chew2014generalizedgauge[§1, Eqs. (5)–(6)]).
With $vb(E) = -partial_t vb(A) - grad phi$ and $vb(B) = curl vb(A)$ we obtain
$
  vb(E) = -A_0 omega vu(p) sin theta, quad
  vb(B) = -A_0 (vb(k) times vu(p)) sin theta,
$ <e15-EB>
with $vb(k) times vu(p) = k_c (1,1,-2) \/ sqrt(2)$.
So $vb(E)$ has an $x$ and a $y$ component, and $vb(B)$ has all three components.

Since $vb(k) dot vu(p) = 0$, only the gauge part contributes to the divergence:
$div vb(A) = -a_parallel abs(vb(k)) sin theta$.
Since $epsilon_r v omega = abs(vb(k))$, we have $epsilon_r partial_t phi = a_parallel abs(vb(k)) sin theta$,
so $div vb(A) = -epsilon_r partial_t phi$ and the homogeneous Lorenz gauge holds
(Chew @chew2014generalizedgauge[Eq. (11)]).
Since $epsilon_r omega^2 = abs(vb(k))^2$, every component of $vb(A)$ and $phi$ satisfies the homogeneous wave equation
(Chew @chew2014generalizedgauge[§2, Eqs. (12)–(13)]).
@e15-wave is therefore a source-free Maxwell solution in Lorenz gauge.

For the projection in @e15-measurements we also use the complex fields
$
  vb(cal(E)) = i A_0 omega vu(p) e^(i theta), quad
  vb(cal(B)) = i A_0 (vb(k) times vu(p)) e^(i theta),
$
whose real parts are the fields of @e15-EB.
]

== The wave on the mesh and its dispersion relation <e15-dispersion>

#proposal[
The mesh has cubic cells, $h_x = h_y = h_z = h$.
In the staggered placement, $A_d$ lives on the positive $d$-face of cell $bold(i)$,
at $vb(x)_bold(i) + (h\/2) bold(e)_d$, and $phi$ at the cell centre $vb(x)_bold(i)$
(Ryu et al. @ryu2016potentialfdtd[Appendix A, Eqs. (19)–(23)]; `LevelContract::positionOffset`).
By the harmonic representation of Schneider @schneider2010UnderstandingFdtd[§7.3],
a staggered difference along $d$ acts on a sampled mode $e^(i vb(k) dot vb(x))$ as multiplication by $i tilde(k)_d$, with
$
  tilde(k)_d = 2/h sin((k_d h)/2).
$
Since all components of $vb(k)$ are equal, $tilde(vb(k))$ is parallel to $vu(k)$ (deduction).
This is the reason for cubic cells and for equal components of $vb(k)$.
It has two consequences (deduction):

- The sampled transverse part has zero discrete divergence, because $tilde(vb(k)) dot vu(p) = 0$.
- The sampled gauge part is a discrete gradient, because it points along $tilde(vb(k))$.
  Its discrete curl therefore vanishes.

For $alpha = 1$ and the homogeneous medium, the O1VS update of each part reduces to the scalar leapfrog update,
with the spatial operator acting as multiplication by $-abs(tilde(vb(k)))^2$:
the staggered curl–curl on the transverse part,
the gradient of the weighted divergence on the gauge part,
and the weighted Laplacian on $phi$ (deduction from the potential equations of Chew @chew2014generalizedgauge[Eqs. (15)–(18)] with $chi = alpha epsilon_r^2$).
Substitution of $e^(i (vb(k) dot vb(x) - omega_h t))$ into the leapfrog update gives one discrete angular frequency $omega_h$ for all three parts,
$
  epsilon_r (2/h_t)^2 sin^2((omega_h h_t)/2) = abs(tilde(vb(k)))^2 = 3 (2/h)^2 sin^2((k_c h)/2).
$ <e15-dispersion-relation>
This is the three-dimensional form of the FDTD dispersion relation
(Schneider @schneider2010UnderstandingFdtd[§7.4, Eq. (7.34)] in one dimension; the three-dimensional form follows by the same substitution, deduction).
It gives $omega_h < omega$ for the meshes below.

In R2, the electric field after $n$ steps is
$
  vb(E)^(n-1/2) = -(vb(A)^n - vb(A)^(n-1))/h_t - grad_h phi^(n-1/2),
$
where $phi^(n-1/2)$ is the previous scalar slot, at $t_(n-1/2)$ (`LevelContract`).
For the gauge part, at the face positions and with $psi = vb(k) dot vb(x) - omega_h t_(n-1/2)$,
the time-difference term is $-a_parallel tilde(omega) vu(k) sin psi$ with $tilde(omega) = (2\/h_t) sin(omega_h h_t\/2)$,
and the gradient term is $v a_parallel tilde(vb(k)) sin psi$.
By @e15-dispersion-relation, $tilde(omega) = v abs(tilde(vb(k)))$, and $tilde(vb(k)) = abs(tilde(vb(k))) vu(k)$,
so the two terms cancel exactly (deduction).
In the same way, the discrete Lorenz-gauge residual vanishes.
The discrete $vb(E)$ and $vb(B)$ therefore come from the transverse part alone,
and the discrete Gauss residual of $epsilon_r vb(E)$ vanishes, in exact arithmetic.

The cancellation needs $phi$ at $t_(n-1/2)$.
If the previous scalar slot held $phi$ at $t_n$ instead,
$vb(E)$ would contain a longitudinal part, and the residual $r$ of E15.meas03 would be
$sqrt(2) (a_parallel \/ A_0) (tilde(omega)\/omega) sin(omega_h h_t\/4) approx omega h_t \/ (2 sqrt(2))$,
about $8 times 10^(-2)$ at $N = 16$ (deduction).
]

== Mesh, time step, duration and ranks

#proposal[
We use $N in {16, 32, 64}$ cells per axis and $h = L\/N$.
$N$ is also the number of cells per wavelength along each axis.
The time step is $h_t = h\/2$ on every level.
The leapfrog update is stable when $max_j abs(lambda_j) h_t^2 < 4$
(Kreiss et al. @kreiss2002SecondOrderWave[§2, Eq. (2.8)]),
which for this homogeneous medium gives $h_t < h sqrt(epsilon_r \/ 3)$ (deduction).
The vacuum bound $h_t < h\/sqrt(3)$ (Fallahi @fallahi2020mithra20fullwavesimulation[Eq. (3.20)]) is stricter, and $h_t = h\/2$ satisfies both.

The run lasts $P = 1000$ periods, $T_"end" = P T_p approx 842.6$.
As in E01, we stop at the R2 sample nearest to $T_"end"$:
$
  n_"end" = floor(T_"end"\/h_t + 1), quad t_"end" = (n_"end" - 1/2) h_t.
$

#align(center)[
  #table(
    columns: 4,
    align: center,
    inset: 6pt,
    stroke: 0.5pt,
    table.header([$N$], [$h_t$], [$n_"end"$], [cell updates]),
    [16], [$1\/32$], [26964], [$1.1 times 10^8$],
    [32], [$1\/64$], [53928], [$1.8 times 10^9$],
    [64], [$1\/128$], [107855], [$2.8 times 10^10$],
  )
]

We use eight MPI ranks with decomposition along $z$, the one axis that the core decomposes,
and four OpenMP threads per rank.
At $N = 16$ each rank owns two cell layers.
From the rates of the E03 runs we estimate about ten minutes for all three levels on the PC.
]

== Initialization

#proposal[
We load the two vector-potential levels and the two scalar-potential slots of the discrete wave:
$
  (A_d^n)_(bold(i) + 1/2 bold(e)_d)
  = (A_0 p_d + a_parallel hat(k)_d) cos(vb(k) dot vb(x)_(d,bold(i)) - omega_h n h_t), quad n in {0, -1},
$
$
  phi_bold(i)^m = v a_parallel cos(vb(k) dot vb(x)_bold(i) - omega_h m h_t), quad m in {1/2, -1/2}.
$
Here $vb(x)_(d,bold(i)) = vb(x)_bold(i) + (h\/2) bold(e)_d$ is the face position of $A_d$,
and $omega_h$ comes from @e15-dispersion-relation.
In R2, the current scalar slot holds $phi$ at $t_(1/2)$ and the previous slot at $t_(-1/2)$;
the loader takes these times from `LevelContract`.
$vb(A)^0$ equals the continuum sample of @e15-wave.
The other three samples differ from the continuum samples at their times by at most $abs(omega - omega_h) h_t$ relative to their amplitudes.

This is initialization option B of `../E01E02E03E04_wave_pulses_in_homogeneous_dielectric/experiment_definition_draft.md`, not option A of E01.
The reason is as follows.
Continuum samples at $t = 0$ and $t = -h_t$ start a small backward wave beside the forward wave,
of relative amplitude about $(omega - omega_h)\/(2 omega)$
(deduction; about $5 times 10^(-4)$ at $N = 32$).
The measured amplitude would then oscillate with angular frequency about $2 omega$,
and a slow change of amplitude below that level would be hidden.
With option B the discrete solution is exactly the forward wave,
so the amplitude is constant to rounding error and a slow change becomes visible.
]

== Expected behaviour <e15-expected>

#proposal[
With this initialization the discrete solution is
$vb(A)^n = (A_0 vu(p) + a_parallel vu(k)) cos(vb(k) dot vb(x) - omega_h n h_t)$ at the face positions
and $phi^m = v a_parallel cos(vb(k) dot vb(x) - omega_h m h_t)$ at the cell centres,
for all $n$ and all scalar times $m$, in exact arithmetic (deduction from @e15-dispersion-relation).
The gauge part does not contribute to the R2 fields (@e15-dispersion).
The R2 fields at $t_(n-1/2) = (n - 1/2) h_t$ then differ from the continuum fields @e15-EB in two ways only (deduction):

+ The amplitude is multiplied by a constant factor.
  It is $kappa_E = 2 sin(omega_h h_t\/2) \/ (omega h_t)$ for $vb(E)$, from the time difference of $vb(A)$,
  and $kappa_B = (sin(k_c h\/2)\/(k_c h\/2)) cos(omega_h h_t\/2)$ for $vb(B)$, from the staggered curl and the R2 average.
+ The phase lags by $delta(t) = (omega - omega_h) t$, which grows linearly with time.

Below the stability bound, leapfrog conserves a positive discrete energy
(Kreiss et al. @kreiss2002SecondOrderWave[§2, Eqs. (2.8)–(2.10)]).
Rounding errors made in each step therefore stay bounded while they propagate,
and their accumulated size after $n$ steps is at most of order $n epsilon_"mach"$ relative to the wave (deduction).

#align(center)[
  #table(
    columns: 6,
    align: center,
    inset: 6pt,
    stroke: 0.5pt,
    table.header([$N$], [$(omega - omega_h)\/omega$], [$kappa_E - 1$], [$kappa_B - 1$], [$delta(t_"end")$ in rad], [$e_F (t_"end")$]),
    [16], [$4.18 times 10^(-3)$], [$-6.41 times 10^(-3)$], [$-1.31 times 10^(-2)$], [26.27], [1.07],
    [32], [$1.04 times 10^(-3)$], [$-1.61 times 10^(-3)$], [$-3.30 times 10^(-3)$], [6.55], [0.26],
    [64], [$2.60 times 10^(-4)$], [$-4.02 times 10^(-4)$], [$-8.25 times 10^(-4)$], [1.64], [1.46],
  )
]
Here $e_F$ is the continuum field error of @e15-derived-error.
The last column shows why a final-time continuum error is not a useful long-run measure:
it does not decrease with $N$, because the phase lag passes multiples of $2 pi$ at different times.
The phase-lag rate $(omega - omega_h)\/omega$ decreases by a factor of four per refinement, as expected for a second-order scheme.
]

== Measured quantities <e15-measurements>

#proposal[
All three measurements use the R2 fields $vb(E)$ and $vb(B)$ after $n$ completed steps,
with time label $t_(n-1/2) = (n - 1/2) h_t$ (`LevelContract::halfStepLabel`).
$E_d$ lives on the positive $d$-face, $vb(x)_bold(i) + (h\/2) bold(e)_d$,
and $B_d$ on the edge $vb(x)_bold(i) + (h\/2)(bold(e)_p + bold(e)_q)$ with $p, q != d$
(`LevelContract::positionOffset`).
The complex reference fields $vb(cal(E)), vb(cal(B))$ are evaluated at the same positions and at the same time.
Let $sans(M)$ be the set of all owned cell indices.
All sums below run over $bold(i) in sans(M)$ and $d in {x,y,z}$,
with every sample at its own position.

=== E15.meas01 and E15.meas02: amplitude ratio and phase lag of $vb(E)$ and $vb(B)$

We define the complex mode coefficients
$
  c_E (t_(n-1/2)) = (2 sum E_d^"num" overline(cal(E)_d)) / (sum abs(cal(E)_d)^2), quad
  c_B (t_(n-1/2)) = (2 sum B_d^"num" overline(cal(B)_d)) / (sum abs(cal(B)_d)^2),
$
and from them the amplitude ratios and phase lags
$
  rho_E = abs(c_E), quad delta_E = arg c_E, quad
  rho_B = abs(c_B), quad delta_B = arg c_B,
$
where $delta_E$ and $delta_B$ are unwrapped across consecutive observations.
E15.meas01 is $c_E$ and E15.meas02 is $c_B$.
If the numerical $vb(E)$ is the continuum $vb(E)$ with amplitude factor $rho$ and phase lag $delta$, then $c_E = rho e^(i delta)$,
and likewise for $vb(B)$.
On the periodic mesh the real and imaginary parts of $vb(cal(E))$ are orthogonal and equally large,
because the sum of $e^(2 i vb(k) dot vb(x))$ over all samples is zero (deduction).
Therefore $op("Re")(c_E vb(cal(E)))$ is the least-squares fit of the continuum wave to the numerical $vb(E)$,
and likewise for $vb(B)$.
We fit $vb(E)$ and $vb(B)$ separately because their constant amplitude factors differ (@e15-expected).
One common coefficient would leave their difference, about $3 times 10^(-3)$ at $N = 16$, in the residual of E15.meas03.

*Purpose:* $rho_E$ and $rho_B$ show whether the wave grows or decays during the run.
The expected values are the constants $kappa_E$ and $kappa_B$ of @e15-expected.
A trend indicates instability or dissipation.
$delta_E$ and $delta_B$ measure the accumulated numerical dispersion.
Both are expected to equal $(omega - omega_h) t_(n-1/2)$.

=== E15.meas03: residual outside the wave

$
  r(t_(n-1/2)) = sqrt(
    (sum [epsilon_r (E_d^"num" - op("Re")(c_E cal(E)_d))^2 + (B_d^"num" - op("Re")(c_B cal(B)_d))^2])
    / (sum [epsilon_r (op("Re") cal(E)_d)^2 + (op("Re") cal(B)_d)^2])
  ).
$
The weights $epsilon_r$ and $1$ are those of the electromagnetic energy density
(Ellingson @ellingson2020ElectromagneticsVol2[§3.1, Eqs. (3.20), (3.22)]),
as in `opt-in/diagnostics/AlignedFieldEnergy.h`.
Since $epsilon_r omega^2 = abs(vb(k))^2$, the electric and magnetic parts of the denominator are equal.

*Purpose:* measure everything in the numerical solution that is not the original wave:
accumulated rounding error, a mode that grows from rounding error, a longitudinal (charge-like) field,
or a defect at rank boundaries.
A longitudinal field also appears when the gauge part of $vb(A)$ and $phi$ stop cancelling in $vb(E)$,
for example through a wrong time level of $phi$ (@e15-dispersion)
or a defect in the scalar update or in the weighted-divergence term.
The expected value is at rounding level, at most of order $n epsilon_"mach"$.
A longitudinal field is orthogonal to the transverse reference,
because summation by parts moves the discrete gradient onto the zero discrete divergence of $vb(cal(E))$ (deduction).
It contributes fully to $r$.
We therefore do not measure the Gauss residual separately.

=== Derived continuum field error <e15-derived-error>

The relative continuum error with the same weights,
$
  e_F = sqrt(
    (sum [epsilon_r (E_d^"num" - op("Re") cal(E)_d)^2 + (B_d^"num" - op("Re") cal(B)_d)^2])
    / (sum [epsilon_r (op("Re") cal(E)_d)^2 + (op("Re") cal(B)_d)^2])
  ),
$
follows from the three measurements,
$
  e_F^2 = 1/2 abs(1 - c_E)^2 + 1/2 abs(1 - c_B)^2 + r^2,
$
by the orthogonality above and the equal electric and magnetic parts (deduction).
The analysis computes it; the run does not measure it separately.
It is the discrete $L^2$ norm of Salari and Knupp @salari2000CodeVerificationMms[§4.1, Eq. (15)] with energy weights.

=== Schedule and saved output

#table(
  columns: (1fr, 4fr),
  inset: 6pt,
  align: left + top,
  stroke: 0.5pt,
  [*Scalar table*],
  [
    E15.meas01–03 at $n = 0$, at every $N_"obs" = floor(T_p\/h_t + 1/2)$ steps (about one period), and at $n_"end"$.
    Each row records $N$, $n$ and the actual time label $t_(n-1/2)$.
    The observation at $n = 0$ has the time label $-h_t\/2$.
  ],
  [*Full-domain fields*],
  [
    Numerical $vb(E)$ and $vb(B)$ (R2), $vb(A)$ (time label $t_n$) and the previous scalar slot $phi$ (time label $t_(n-1/2)$),
    at $n = 0$ and $n = n_"end"$,
    with their component positions and time labels from `LevelContract`.
    They localize the content of $r$ if $r$ grows, and the potentials show which update is at fault.
  ],
)
]

== Agreement and disagreement

#proposal[
We expect, on every level:
- $rho_E$ and $rho_B$ constant over the run and equal to $kappa_E$ and $kappa_B$ to rounding level;
- $delta_E$ and $delta_B$ on the line $(omega - omega_h) t$ to rounding level;
- $r$ at rounding level, growing at most linearly in $n$.

A trend in $rho_E$ or $rho_B$, a departure of $delta_E$ or $delta_B$ from the line, or a faster growth of $r$
indicates an instability or a defect in the implementation.
We then use the saved fields to locate it.
]

== Limits <e15-limits>

#proposal[
E15 tests one mode, one direction, one polarization and one Courant number.
A mode that is absent at the start grows only from rounding error.
If such a mode grows by a factor $1 + gamma$ per step, $r$ shows it once $e^(gamma n_"end")$ is clearly above one,
so E15 detects growth rates per step above about $5\/n_"end"$ (deduction; about $5 times 10^(-5)$ at $N = 64$).
Random initial data would excite all modes at once, but it has no continuum reference.

The gauge part of $vb(A)$ and $phi$ evolve in separate updates.
A defect in either one breaks their cancellation in $vb(E)$ and appears in $r$,
but $r$ does not show which update is at fault; the saved potentials do.
E15 uses $alpha = 1$ only.
For $alpha != 1$ the gauge part travels at $v\/sqrt(alpha)$, slower or faster than the physical wave,
and $vb(E)$ stays transverse (deduction from the generalized gauge of Chew @chew2014generalizedgauge[Eqs. (15)–(18)]).
It makes no statement about Mur faces, material interfaces or sources.
]

// --- standalone bibliography (delete when included in the manual) ---
#bibliography("E15-references.bib", style: "ieee", title: [References])
