#import "../../docs/manual/preamble.typ": *

= Experiment E16 detailed specs <experiment-e16>

E16 sends a transverse plane pulse from vacuum at normal incidence onto a planar face of the dielectric.
At the face the pulse splits into a reflected pulse and a transmitted pulse.
The main result is the energy on each side of the face after the split,
compared with the continuum energy fractions, and its convergence.
We also compare $vb(E)$ and $vb(B)$ with the continuum solution separately on each side of the face.
The pulse travels only a short distance, so that little numerical dispersion enters.

Everything that is not stated otherwise below is the same as in #link(<experiment-e01>)[E01]:
the units, the notation for the staggered samples and their time levels,
the solver configuration, the transverse domain and mesh,
the initialization method and the nearest-time rule for the observations.

== Experiment E16 analytic continuum reference

The two media are nonmagnetic, lossless and nondispersive.
There is no free charge or current.
The face is the plane $z = z_f$, and
$
  epsilon_r (z) = cases(
    epsilon_1 & "if" z < z_f,
    epsilon_2 & "if" z >= z_f
  ), quad
  n_j := sqrt(epsilon_j), quad v_j := c_"vac"\/n_j, quad j in {1,2}.
$ <E16-media>
E16 uses $epsilon_1 = 1$ (vacuum) and $epsilon_2 = 2.13$.
We keep $epsilon_1$ and $epsilon_2$ general in the formulas, because E17 exchanges them.

We define three phases,
$
  s_"i" := z - z_0 - v_1 t, quad
  s_"r" := 2 z_f - z - z_0 - v_1 t, quad
  s_"t" := z_f - z_0 - v_1 t + n_2/n_1 (z - z_f).
$ <E16-phases>
The phase $s_"r"$ is the mirror image of $s_"i"$ at the face.
In $s_"t"$ the factor $n_2\/n_1$ compresses the profile in the second medium when $n_2 > n_1$.
At $z = z_f$ the three phases are equal.

The potentials are
$
  vb(A)(z,t) = vb(hat(x)) cases(
    a(s_"i") + Gamma a(s_"r") & "if" z < z_f,
    tau a(s_"t") & "if" z >= z_f
  ), quad phi = 0,
$ <E16-analytic-potentials>
with
$
  Gamma = (n_1 - n_2)/(n_1 + n_2), quad tau = 1 + Gamma = (2 n_1)/(n_1 + n_2).
$ <E16-coefficients>
These are the reflection and transmission coefficients of the electric field at normal incidence
(Ellingson @ellingson2020ElectromagneticsVol2[§5.1, Eqs. (5.16), (5.17), (5.27)]).
They do not depend on the frequency.
Therefore the reflected and the transmitted pulse keep the profile of the incident pulse.

As in E01, $vb(E) = -partial_t vb(A) - grad phi$ and $vb(B) = curl vb(A)$
@chew2014generalizedgauge[§2, Eqs. (5)–(6)] give
$
  vb(E)(z,t) = vb(hat(x)) v_1 cases(
    a'(s_"i") + Gamma a'(s_"r") & "if" z < z_f,
    tau a'(s_"t") & "if" z >= z_f
  ), quad
  vb(B)(z,t) = vb(hat(y)) cases(
    a'(s_"i") - Gamma a'(s_"r") & "if" z < z_f,
    (n_2\/n_1) tau a'(s_"t") & "if" z >= z_f
  ).
$ <E16-analytic-EB>

#let n-dielectric = calc.sqrt(2.13)
#let gamma-e16 = (1 - n-dielectric) / (1 + n-dielectric)
#let rounded(value, digits: 4) = str(calc.round(value, digits: digits))

For E16, $n_1 = 1$ and $n_2 = sqrt(2.13)$ give
$Gamma approx #rounded(gamma-e16)$ and $tau approx #rounded(1 + gamma-e16)$.
The reflected pulse has the opposite sign of the incident pulse.

=== Proof: The potentials solve the equations of the two media
We use Chew's equations for an inhomogeneous medium with $mu_r = 1$.

*Gauge.*
Neither $A_x$ nor $epsilon_r$ depends on $x$, so $div (epsilon_r vb(A)) = partial_x (epsilon_r A_x) = 0$.
With $phi = 0$ this satisfies the generalized gauge $div (epsilon_r vb(A)) = -chi partial_t phi$ @chew2014generalizedgauge[Eq. (15)].

*Wave equation in each medium.*
Without the gauge term and without current,
the vector equation @chew2014generalizedgauge[Eq. (17)] reduces to
$
  partial_z^2 A_x - epsilon_j partial_t^2 A_x = 0 quad "in medium" j.
$
For $a in cal(C)^2$, the terms $a(s_"i")$ and $a(s_"r")$ give $a'' - epsilon_1 v_1^2 a'' = 0$.
The term $a(s_"t")$ gives $(n_2\/n_1)^2 a'' - epsilon_2 v_1^2 a'' = 0$.

*Conditions at the face.*
The tangential components of $vb(A)$ and of $mu_r^(-1) curl vb(A)$ are continuous @chew2014generalizedgauge[§3, Eqs. (23)–(24)].
Here these are $A_x$ and $partial_z A_x$.
Since the three phases are equal at $z = z_f$, the two conditions are
$
  1 + Gamma = tau, quad 1 - Gamma = n_2/n_1 tau,
$
and @E16-coefficients is their solution.

*Before the pulse reaches the face.*
The profile below vanishes for $abs(s) >= w$.
While $v_1 t < z_f - z_0 - w$, we have $s_"r" > w$ for all $z < z_f$ and $s_"t" > w$ for all $z >= z_f$.
Only the incident term is then nonzero, which is the pulse of E01 in the first medium.

=== Energy of the three pulses
The energy density is $u = (epsilon_r norm(vb(E))^2 + norm(vb(B))^2)\/2$
@ellingson2020ElectromagneticsVol2[§3.1, Eqs. (3.20), (3.22)].
With @E16-analytic-EB the incident pulse has $u = a'(s_"i")^2$,
the reflected pulse $u = Gamma^2 a'(s_"r")^2$
and the transmitted pulse $u = (n_2\/n_1)^2 tau^2 a'(s_"t")^2$.
We integrate over $z$ and use $dif s_"t" = (n_2\/n_1) dif z$.
After the two pulses have separated, which is the case for $v_1 t > z_f - z_0 + w$,
the energies on the two sides are
$
  U_1 = Gamma^2 U_"i", quad
  U_2 = n_2/n_1 tau^2 U_"i" = (1 - Gamma^2) U_"i",
$ <E16-energy-fractions>
where $U_"i"$ is the energy of the incident pulse.
These are the fractions of the reflected and the transmitted power density
@ellingson2020ElectromagneticsVol2[§5.1, Eqs. (5.22), (5.24)].
For E16, $Gamma^2 approx #rounded(gamma-e16 * gamma-e16)$ and $1 - Gamma^2 approx #rounded(1 - gamma-e16 * gamma-e16)$.

=== Specification of the profile
We use the compact envelope $p_beta$ of E03 (@e03-profile-sharpness) with $beta = 8$, $A_0 = 1$ and $w = 1$,
and we multiply it by a carrier:
$
  a(s) = p_beta (s) cos(k_c s), quad k_c = (2 pi)/lambda_c, quad lambda_c = w/4.
$ <E16-profile>
Differentiation gives
$
  a'(s) = p_beta '(s) cos(k_c s) - k_c p_beta (s) sin(k_c s),
$ <E16-profile-derivative>
with $p_beta '$ as in E03.
Both vanish for $abs(s) >= w$.
The input parameter `carrier_wavelength` sets $lambda_c$.

We choose this profile for three reasons.
- The carrier gives the pulse a wavelength.
  It is $lambda_c$ in the first medium and $lambda_c n_1\/n_2$ in the second medium.
  We can then state the resolution in cells per wavelength.
- The envelope with $beta = 8$ is close to a Gaussian with $sigma = w\/4$ near its center.
  The flanks of the E01 profile ($beta = 1$) were too steep for the coarse levels.
- The support is compact.
  The pulse is exactly zero at the second face of the periodic box (see below).

Unlike E01, we do not extend the profile periodically.

#import "@preview/cetz:0.3.4": canvas, draw

// Envelope sharpness and carrier wavenumber of the figure, with w = 1.
#let figure-beta = 8
#let figure-kc = 8 * calc.pi

// The pair (a, a'/k_c) at u = s/w.
#let profile-values(u) = {
  if calc.abs(u) >= 1 {
    (0, 0)
  } else {
    let q = 1 - u * u
    let envelope = calc.exp(figure-beta * (1 - 1 / q))
    let slope = -2 * figure-beta * u * envelope / (q * q)
    let phase = figure-kc * u
    (envelope * calc.cos(phase),
      slope * calc.cos(phase) / figure-kc - envelope * calc.sin(phase))
  }
}

#let profile-axes(position, ordinate) = {
  draw.rect((0, 0), (5.2, 3), stroke: 0.35pt + luma(65%))
  draw.line(position(-1.35, 0), position(1.35, 0), stroke: 0.35pt + luma(75%))
  for u in (-1, 0, 1) {
    let x = position(u, 0).at(0)
    draw.line((x, 0), (x, -0.08), stroke: 0.35pt + luma(65%))
    draw.content((x, -0.16), anchor: "north", [$#u$])
  }
  for value in (-1, 0, 1) {
    let y = position(0, value).at(1)
    draw.line((0, y), (-0.08, y), stroke: 0.35pt + luma(65%))
    draw.content((-0.16, y), anchor: "east", [$#value$])
  }
  draw.content((0, 3.16), anchor: "south-west", ordinate)
  draw.content((2.6, -0.54), anchor: "north", [$s\/w$])
}

#let profile-panel(derivative: false) = {
  let ordinate = if derivative { [$a'(s)\/(k_c A_0)$] } else { [$a(s)\/A_0$] }
  let position(u, value) = (5.2 * (u + 1.35) / 2.7, 3 * (value + 1.15) / 2.3)
  let points = range(1001).map(i => {
    let u = -1.35 + 2.7 * i / 1000
    position(u, profile-values(u).at(if derivative { 1 } else { 0 }))
  })
  set text(size: 9pt)
  canvas(length: 10mm, {
    profile-axes(position, ordinate)
    draw.line(..points, stroke: 0.8pt)
  })
}

#figure(
  grid(columns: (1fr, 1fr), column-gutter: 8pt, align: center,
    profile-panel(), profile-panel(derivative: true)),
  caption: [Profile of E16 and its derivative, in dimensionless form.
    The fields $E_x$ and $B_y$ are proportional to the derivative.],
) <fig-e16-profile>

== Experiment E16 solver configuration

The configuration is that of E01:
O1VS with $alpha = 1$, Standard stencil with staggered placement, R2,
periodic boundaries in all directions and $h_t = h_z\/2$.
Only the material differs.
The solver samples @E16-media at the cell centers
(`TwoMaterialMap` with `PlanarSlab` in `core/material`),
with $mu_r = 1$, $rho = 0$ and $vb(J) = vb(0)$.

== Experiment E16 computational domain and refinement

The transverse domain and mesh are those of E01,
$L_x = L_y = w = 1$ and $N_x = N_y = 8$.
Along the propagation direction we choose
$
  z in [0, L_z], quad L_z = 10w, quad z_f = 4w, quad z_0 = 2.5w.
$
The first medium occupies $0 <= z < z_f$ and the second medium $z_f <= z < L_z$.

=== The second face of the periodic box
The box is periodic along $z$.
It therefore contains a second face at $z = 0$, which is the same plane as $z = L_z$.
The reference @E16-analytic-potentials describes one face only.
It is the solution in the box as long as the fields are zero at the second face.
The target duration is
$
  T = 3w/v_1.
$
The pulse occupies $[1.5w, 3.5w]$ at $t = 0$.
Its front reaches the face at $0.5w\/v_1$, and the two pulses are separated after $2.5w\/v_1$.
At $t = T$ the reflected pulse occupies $[1.5w, 3.5w]$ again,
and the transmitted pulse occupies $z_f + (n_1\/n_2) [0.5w, 2.5w]$,
which is about $[#rounded(4 + 0.5 / n-dielectric, digits: 2) w, #rounded(4 + 2.5 / n-dielectric, digits: 2) w]$ for E16.
During the whole run the continuum fields are exactly zero within $1.5w$ of the second face.

=== The travel distance
The pulse starts as close to the face as its support allows with a margin of $0.5w$,
and the run ends $0.5w\/v_1$ after the two pulses have separated.
Each part of the pulse then travels $3w$ in the first medium, which is #(3 * 4) wavelengths $lambda_c$.
We keep this distance short because the numerical dispersion shifts the phase of the pulse in proportion to the travelled distance
(Schneider @schneider2010UnderstandingFdtd[§7.4]).
E01 measured this growth without a face.
The energy on each side does not depend on the phase.
It is therefore the main result, and the field errors are secondary.

=== The position of the face on the mesh
As in E01, $N_"pulse" := 2w\/h_z$ is the number of cells across the support of the pulse.
Then $N_z = 5 N_"pulse"$, and the cell centers are $z_k = (k + 1\/2) h_z$ for $k = 0, dots, N_z - 1$.
The face lies at $z_f = 2 N_"pulse" h_z$.
This is the cell face between the cells $k = 2 N_"pulse" - 1$ and $k = 2 N_"pulse"$ on every level.
All cells with $z_k < z_f$ get $epsilon_1$, and all others get $epsilon_2$.
The face of the sampled material is therefore at $z_f$ on every level.
E16 does not test a face inside a cell.

=== Refinement levels
We refine $h_z$ and $h_t$ and keep the transverse mesh fixed, as in E01.

#let cells-in-dielectric(cells) = rounded(cells / n-dielectric, digits: 1)
#align(center)[
  #table(
    columns: 5,
    align: center,
    inset: 6pt,
    stroke: 0.5pt,
    table.header(
      [$N_"pulse"$], [$h_z\/w$], [$N_z$],
      [$lambda_c\/h_z$], [$lambda_c\/(n_2 h_z)$],
    ),
    [128],  [$1/64$],   [640],   [16],  [#cells-in-dielectric(16)],
    [256],  [$1/128$],  [1280],  [32],  [#cells-in-dielectric(32)],
    [512],  [$1/256$],  [2560],  [64],  [#cells-in-dielectric(64)],
    [1024], [$1/512$],  [5120],  [128], [#cells-in-dielectric(128)],
    [2048], [$1/1024$], [10240], [256], [#cells-in-dielectric(256)],
  )
]
The last two columns are the cells per wavelength in vacuum and in the dielectric.
The physical domain, the profile and the target duration are the same on all levels.
We use eight MPI ranks along $z$ and four OpenMP threads per rank.

=== Time steps and comparison times
We set $h_t = h_z\/2$ on every level.
With $v_1 = c_"vac"$ this gives $T\/h_t = 3 N_"pulse"$ steps.
The observation targets are
$
  tau in {0, T/2, T}.
$
For each target we select the nearest sample of $vb(A)$ and of the R2 fields with the rule of E01,
and we evaluate the reference at the actual time of each sample.
At $tau = 0$ the pulse has not reached the face.
At $T\/2$ its center is at the face.
At $T$ the reflected and the transmitted pulse are separated.
The energies of E16.meas03 are also sampled every $N_"pulse"\/16$ steps, which gives their history during the split.

== Experiment E16 initialization
We sample the continuum potential history as in E01,
$
  vb(A)^0 eq.triple_"stagg" vb(A)_"exact" (0), quad vb(A)^(-1) eq.triple_"stagg" vb(A)_"exact" (-h_t),
$
with @E16-analytic-potentials and $phi = 0$ in both slots.
At these two times only the incident term is nonzero.
The pulse lies in the first medium, $0.5w$ away from the face and $1.5w$ away from the second face.
The loaded histories are therefore those of E01 in a uniform medium with $epsilon_1$.
The discrete gauge holds in exact arithmetic for the same reason as in E01:
neither $A_x$ nor its permittivity varies along $x$.

== What E16 exercises in the update
The vector potential has only the component $A_x$, which is tangential to the face and does not vary along $x$.
The weighted divergence of the gauge @chew2014generalizedgauge[Eq. (15)] is then zero.
Its discrete form in `O1VSOperator` is a backward difference along $x$ and is zero as well.
The gauge term of the vector equation and the scalar potential stay zero during the run.
The permittivity of $A_x$ is the harmonic mean of two cells that are neighbors along $x$ (`HarmonicFace`).
This is the cell value $epsilon_r (z_k)$.
The update of $A_x$ is then the discrete wave equation of one component with the cell permittivity.

At normal incidence there is no plane of incidence, so the cases TE and TM coincide.
For the update E16 is of the TE kind, because the electric field is tangential to the face.
E16 therefore tests the update of a tangential component across the face.
An update that only scales the wave equation of each component by the cell permittivity gives the same result for this field,
so E16 cannot distinguish O1VS from it.
E16 does not test the component of $vb(A)$ normal to the face,
the harmonic mean of two different permittivities, the gauge term or the scalar potential.
These parts need the TM case at oblique incidence, with the electric field in the plane of incidence.

== Experiment E16 measured quantities and their purpose

Let $sans(M)_1$ be the set of the cells with $z_k < z_f$ and $sans(M)_2$ the set of the other cells.
A sample belongs to the side of the cell index $bold(i)$ that stores it.
The samples of $B_y$ at $z_k + h_z\/2 = z_f$ lie on the face and belong to side 1.
Both expressions for $B_y$ in @E16-analytic-EB agree there.
The numerical samples, the reference samples, their positions and their times are defined as in E01.meas01 and E01.meas02.

=== E16.meas01: Relative electric-field error on each side in the discrete $L^2$ norm
For $S in {1,2}$ and $n$ completed steps, define the two sums
$
  N_(E,S)^(n-1/2) = sum_(bold(i) in sans(M)_S,c)
    abs(E_(c,bold(i)+1/2 bold(e)_c)^("num",n-1/2) - E_(c,bold(i)+1/2 bold(e)_c)^("ref",n-1/2))^2, quad
  D_(E,S)^(n-1/2) = sum_(bold(i) in sans(M)_S,c)
    abs(E_(c,bold(i)+1/2 bold(e)_c)^("ref",n-1/2))^2,
$
and from them
$
  e_(E,S)^(n-1/2) = sqrt(N_(E,S)^(n-1/2) / D_(E,S)^(n-1/2)), quad
  e_(E)^(n-1/2) = sqrt((N_(E,1)^(n-1/2) + N_(E,2)^(n-1/2)) / (D_(E,1)^(n-1/2) + D_(E,2)^(n-1/2))).
$
The global error $e_E$ is E01.meas01 with the reference of E16.
We do not form $e_(E,2)$ while $D_(E,2) = 0$, which is the case before the pulse reaches the face.

*Purpose:* measure the reflected and the transmitted pulse separately.
After the separation, side 1 contains only the reflected pulse and side 2 only the transmitted pulse.
Each error is then relative to the norm of its own pulse.
The reflected pulse carries only the fraction $Gamma^2$ of the energy,
so the global error would hide an error of the reflected pulse.
A reflected pulse with the wrong sign gives $e_(E,1) approx 2$.

=== E16.meas02: Relative magnetic-field error on each side in the discrete $L^2$ norm
Define $N_(B,S)$, $D_(B,S)$, $e_(B,S)$ and $e_B$ as in E16.meas01,
with the R2 magnetic samples $B_(c,bold(i)+bold(delta)_c^B)^(n-1/2)$ of E01.meas02.

*Purpose:* the reflected magnetic field has the sign opposite to the reflected electric field (@E16-analytic-EB),
and the transmitted magnetic field has the additional factor $n_2\/n_1$.
The magnetic error checks both.

=== E16.meas03: Field energy on each side
Let $epsilon_bold(i)$ be the relative permittivity of cell $bold(i)$. Define
$
  U_S^(n-1/2) = (h_x h_y h_z)/2 sum_(bold(i) in sans(M)_S,c)
    (epsilon_bold(i) abs(E_(c,bold(i)+1/2 bold(e)_c)^("num",n-1/2))^2
     + abs(B_(c,bold(i)+bold(delta)_c^B)^("num",n-1/2))^2), quad
  rho_S^(n-1/2) = U_S^(n-1/2) / (U_1^(-1/2) + U_2^(-1/2)).
$
The denominator is the energy at the initial observation.
For $E_x$ and $E_y$ the cell value $epsilon_bold(i)$ is the permittivity of the component.
For $E_z$ it differs from the face value only in the layer at the face, and the reference of $E_z$ is zero.

*Purpose:* this is the main result of E16. It compares the energy split with @E16-energy-fractions.
After the separation the continuum values are $rho_1 = Gamma^2$ and $rho_2 = 1 - Gamma^2$.
We report the relative deviations at the final observation,
$
  delta_1 = rho_1/Gamma^2 - 1, quad delta_2 = rho_2/(1 - Gamma^2) - 1.
$ <E16-energy-deviations>
The sum $rho_1 + rho_2 - 1$ is the relative change of the sampled energy.
This sampled energy is not the conserved energy of the leapfrog update,
so the sum can differ from one by a discretization error.

== Experiment E16 saved output

Each run writes the following files, and `run.json` with its configuration.

#table(
  columns: (1fr, 4fr),
  inset: 6pt,
  align: left + top,
  stroke: 0.5pt,
  [`measurements.csv`],
  [
    One row per quantity and observation,
    with the target time, the selected step index and the actual time.
    The quantities are the eight sums `E_error_side1`, `E_reference_side1`, `E_error_side2`, `E_reference_side2`
    and the same four names with `B` ($N_(F,S)$ and $D_(F,S)$ of E16.meas01–02),
    and the two energies `U_side1` and `U_side2` (E16.meas03).
    The two energies are also written every $N_"pulse"\/16$ steps.
    The analysis forms the ratios $e_(F,S)$, $e_F$ and $rho_S$ from these values.
  ],
  [`profiles.csv`],
  [
    At the same observations,
    the numerical and the analytic profiles of $A_x$, $E_x$ and R2-averaged $B_y$ along one line parallel to $z$,
    as in E01.
  ],
  [`fields/`],
  [
    The full-domain arrays `E`, `B_R2` and `B_raw` at the three observations.
  ],
)

== Experiment E16 analysis
The main result is $rho_1$ and $rho_2$ at $T$ on each level, with their relative deviations $delta_1$ and $delta_2$ (@E16-energy-deviations).
For each pair of successive levels we report the observed order
$
  p_"obs" = log_2 (e(h_z)\/e(h_z\/2))
$
(Salari and Knupp @salari2000CodeVerificationMms[§4.2, Eq. (25)])
of $abs(delta_1)$ and $abs(delta_2)$.
We plot the history of $rho_1$ and $rho_2$ during the split.

As secondary results we report $e_(E,1)$, $e_(E,2)$, $e_(B,1)$ and $e_(B,2)$ at $T$ with their observed orders,
and we plot the profiles at the three observations.
These errors contain the phase error of the numerical dispersion over the travel distance as well as the error from the face.
