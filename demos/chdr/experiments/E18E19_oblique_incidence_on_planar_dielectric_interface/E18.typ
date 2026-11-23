#import "../../docs/manual/preamble.typ": *
#import "obliquePacketValues.typ": *

= Experiment E18 detailed specs <experiment-e18>

E18 measures the reflection coefficient of the planar face between vacuum and the dielectric as a function of the angle of incidence.
The electric field is perpendicular to the plane of incidence (the TE case).
A short wave packet comes from vacuum.
It contains plane waves with angles of incidence between about $35 degree$ and $60 degree$.
We measure how much of each plane wave is reflected and compare this with the Fresnel coefficient at its angle.

The experiment does not reuse the measurements of E16.
It takes from #link(<experiment-e16>)[E16] only the two media (@E16-media) with $epsilon_1 = 1$ and $epsilon_2 = 2.13$,
the solver configuration, and the rule that the face lies on a cell face.

== Why the packet is a sum of plane waves with one tangential wavenumber
The plane of incidence is the $x z$ plane, and the box is periodic along $x$ with the period $L_x$.
A field fits into the box only if its tangential wavenumber is a multiple of $2 pi\/L_x$. We use
$
  k_x = (2 pi)/L_x.
$
A plane pulse with a compact profile is not periodic along $x$ when it travels obliquely.
A single plane wave with this $k_x$ is periodic, but it fills the whole box and never leaves the face.
We therefore build a packet from plane waves that share $k_x$ and differ in their wavenumber $q$ normal to the face.
Each plane wave has its own angle of incidence $theta = arctan(k_x\/q)$.
We make the spectrum broad on purpose, so that one run covers a range of angles.

== Experiment E18 analytic continuum solution

=== The plane waves
For $m = -M, dots, M$ let
$
  q_m = q_0 + m Delta q, quad
  omega_m = v_1 sqrt(k_x^2 + q_m^2), quad
  q'_m = sqrt(epsilon_2 omega_m^2\/c_"vac"^2 - k_x^2).
$ <E18-modes>
Here $q_m$ and $q'_m$ are the normal wavenumbers in the first and in the second medium.
The tangential wavenumber and the frequency are the same on both sides,
which is the law of refraction @ellingson2020ElectromagneticsVol2[§5.6, Eqs. (5.122), (5.146)].
The three phases of plane wave $m$ are
$
  psi_("i",m) &= k_x x + q_m (z - z_c) - omega_m t, \
  psi_("r",m) &= k_x x - q_m (z - z_f) + q_m (z_f - z_c) - omega_m t, \
  psi_("t",m) &= k_x x + q'_m (z - z_f) + q_m (z_f - z_c) - omega_m t.
$ <E18-phases>
They are equal at the face $z = z_f$.
At $t = 0$ all incident phases agree on the plane $z = z_c$, which is the center of the packet.

=== Potentials and fields
$
  vb(A)(x,z,t) = vb(hat(y)) sum_(m=-M)^M a_m cases(
    cos psi_("i",m) + Gamma_m cos psi_("r",m) & "if" z < z_f,
    tau_m cos psi_("t",m) & "if" z >= z_f
  ), quad phi = 0,
$ <E18-analytic-potentials>
$
  Gamma_m = Gamma(q_m), quad tau_m = 1 + Gamma_m, quad
  Gamma(q) = (q - q')/(q + q').
$ <E18-coefficients>
With the wave impedance $eta_j = 1\/n_j$ and $cos psi^"i" = q\/(n_1 omega)$, $cos psi^"t" = q'\/(n_2 omega)$,
these are the TE coefficients of Ellingson @ellingson2020ElectromagneticsVol2[§5.6, Eqs. (5.142), (5.147)].
The fields follow from $vb(E) = -partial_t vb(A)$ and $vb(B) = curl vb(A)$:
$
  E_y &= -sum_m a_m omega_m cases(
    sin psi_("i",m) + Gamma_m sin psi_("r",m),
    tau_m sin psi_("t",m)
  ), \
  B_x &= sum_m a_m cases(
    q_m (sin psi_("i",m) - Gamma_m sin psi_("r",m)),
    q'_m tau_m sin psi_("t",m)
  ), quad
  B_z = -k_x sum_m a_m cases(
    sin psi_("i",m) + Gamma_m sin psi_("r",m),
    tau_m sin psi_("t",m)
  ).
$ <E18-analytic-EB>
In each pair the upper line holds for $z < z_f$ and the lower line for $z >= z_f$.
All other components are zero.

=== Proof: Each plane wave solves the equations of the two media
*Gauge.* Neither $A_y$ nor $epsilon_r$ depends on $y$, so $div (epsilon_r vb(A)) = 0$,
and the generalized gauge @chew2014generalizedgauge[Eq. (15)] holds with $phi = 0$.

*Wave equation.* The vector equation @chew2014generalizedgauge[Eq. (17)] reduces to
$(partial_x^2 + partial_z^2) A_y - epsilon_j partial_t^2 A_y = 0$ in medium $j$.
Each term satisfies it because $k_x^2 + q_m^2 = epsilon_1 omega_m^2$ and $k_x^2 + q'^2_m = epsilon_2 omega_m^2$ (@E18-modes, $c_"vac" = 1$).

*Conditions at the face.* $A_y$ and $partial_z A_y$ are continuous @chew2014generalizedgauge[§3, Eqs. (23)–(24)].
With equal phases at the face this gives $1 + Gamma_m = tau_m$ and $q_m (1 - Gamma_m) = q'_m tau_m$,
and @E18-coefficients is the solution.

The sum over $m$ is a solution because the equations are linear.

=== The packet
The amplitudes follow the profile $p_beta$ of E03 (@e03-profile-sharpness), used in the wavenumber:
$
  a_m = A_0 (hat(p)(u_m)) / (sum_(m') hat(p)(u_(m'))), quad
  u_m = (q_m - q_0)/Q, quad
  hat(p)(u) = cases(
    exp(beta(1 - 1/(1 - u^2))) & "if" abs(u) < 1,
    0 & "if" abs(u) >= 1.
  )
$ <E18-amplitudes>
#table(
  columns: (1fr, 2.4fr),
  inset: 6pt,
  align: left + top,
  stroke: 0.5pt,
  table.header([*Parameter*], [*Value*]),
  [Central wavenumber], [$q_0 = k_x$, a central angle of incidence of $45 degree$],
  [Half-width of the spectrum], [$Q = 0.85 thin k_x$],
  [Sharpness], [$beta = 8$],
  [Spacing], [$Delta q = 2 pi\/Lambda$ with $Lambda = 4 L_z$],
  [Number of plane waves], [$M = Q\/Delta q = #half-count$, so $2M + 1 = #modes.len()$ plane waves],
  [Amplitude], [$A_0 = 1$],
)
The spectrum is exactly zero outside $q_0 plus.minus Q$, so the packet contains no plane wave with $q <= 0.15 thin k_x$.
Such plane waves would travel almost parallel to the face and would not leave it during the run.
With $beta = 8$ the spectrum is close to a Gaussian with the standard deviation $Q\/4$.
At $t = 0$ the packet is a few periods $L_x$ long and centered at $z_c$.
Each plane wave moves along $z$ with $v_1 cos theta$ (deduction from @E18-modes, $partial omega\/partial q = v_1 q\/sqrt(k_x^2 + q^2)$),
so the packet becomes longer during the run.

The packet is not compact in space, and the sum repeats along $z$ with the period $Lambda$.
The parameters are chosen so that the field stays below $10^(-5)$ of its peak at the second face during the run.
E18.meas03 records the actual value.

=== The reflection coefficient against the angle
The reference of the main measurement is the magnitude of the Fresnel coefficient $abs(Gamma(q))$ of @E18-coefficients,
with $omega = v_1 sqrt(k_x^2 + q^2)$ and $q' = sqrt(epsilon_2 omega^2 - k_x^2)$, at the angle $theta = arctan(k_x\/q)$.
At the central angle $Gamma = #rounded(gamma-te(plane-wave(q-center)))$.

=== Energy fractions
Two plane waves with different $q_m$ are orthogonal over the period $Lambda$.
The energy of the incident packet is therefore the sum of the energies of its plane waves, which are proportional to $a_m^2 omega_m^2$
(energy density from @ellingson2020ElectromagneticsVol2[§3.1, Eqs. (3.20), (3.22)]).
Plane wave $m$ reflects the fraction $Gamma_m^2$ of its power @ellingson2020ElectromagneticsVol2[§5.6, Example 5.7].
After the packets have separated, the energies on the two sides are
$
  U_1 = macron(R) thin U_"i", quad U_2 = (1 - macron(R)) U_"i", quad
  macron(R) = (sum_m a_m^2 omega_m^2 Gamma_m^2)/(sum_m a_m^2 omega_m^2).
$ <E18-energy-fractions>
The second relation is the conservation of energy in lossless media @ellingson2020ElectromagneticsVol2[§3.1, Eq. (3.16)].
For E18, $macron(R) approx #rounded(reflected-fraction-te)$.

== Experiment E18 computational domain and refinement
The unit of length is the tangential period $L_x$. The box is periodic in all directions, with
$
  L_z = 48 L_x, quad z_f = 28 L_x, quad z_c = 23 L_x.
$
Vacuum occupies $0 <= z < z_f$ and the dielectric $z_f <= z < L_z$.
The cells are cubes with $h = L_x\/N_x$.
The fields do not vary along $y$, so we use $N_y = 2$.
Then $N_z = 48 N_x$, and the face lies on the cell face between the cells $k = 28 N_x - 1$ and $k = 28 N_x$ on every level.
As in E16 the periodic box has a second face at $z = 0$, which the fields must not reach.

#align(center)[
  #table(
    columns: 4,
    align: center,
    inset: 6pt,
    stroke: 0.5pt,
    table.header([$N_x$], [$N_z$], [cells per wavelength\ in vacuum], [cells per wavelength\ in the dielectric]),
    ..(24, 48, 96, 192).map(cells => (
      [#cells], [#(48 * cells)],
      [#rounded(cells / calc.sqrt(2), digits: 1)], [#rounded(cells / calc.sqrt(2) / n-dielectric, digits: 1)],
    )).flatten(),
  )
]
The wavelength is that of the central plane wave, $L_x\/sqrt(2)$ in vacuum.
The other plane waves have between #rounded(calc.sqrt(2) * kx / measured-waves.at(10).omega, digits: 2)
and #rounded(calc.sqrt(2) * kx / measured-waves.at(0).omega, digits: 2) times this wavelength.
We use eight MPI ranks along $z$ and four OpenMP threads per rank.

=== Duration
We set $h_t = h\/2$ on every level. The target duration is
$
  T = 28 thin L_x\/c_"vac",
$
which is $56 N_x$ steps.
The slowest plane wave of the measurement has the angle $60 degree$ and moves along $z$ with $c_"vac"\/2$.
It travels $14 L_x$ during the run, which is the distance of $5 L_x$ to the face and back, plus the length of the packet.
The fastest plane waves of the packet move along $z$ with less than $0.9 thin c_"vac"$.
After their reflection they stay more than $2 L_x$ away from the second face.

== Experiment E18 initialization
We sample @E18-analytic-potentials at $t = 0$ and at $t = -h_t$ at the staggered position of $A_y$,
which is $(x_i, y_j + h\/2, z_k)$, and we set all other components and $phi$ to zero in both slots.
The discrete gauge holds in exact arithmetic, because neither $A_y$ nor its permittivity varies along $y$.

== What E18 exercises in the update
The only component $A_y$ is tangential to the face and does not vary along its own direction.
Both discrete divergences vanish, and the permittivity of $A_y$ is the cell value.
The update is the discrete wave equation of one component in the $x z$ plane with the cell permittivity.
An update that only scales the wave equation of each component by the cell permittivity gives the same result.
E18 is therefore the easy polarization.
It establishes the packet, the oblique travel on the mesh and the measurement, which E19 uses for the TM case.

== Experiment E18 measured quantities

=== E18.meas01: Reflection magnitude for each normal wavenumber
Let $sans(M)_1$ be the set of the cells with $z_k < z_f$ (side 1),
and let $(x_bold(i), z_bold(i))$ be the position of the sample of $E_y$ that cell $bold(i)$ stores, which is $(x_i, z_k)$.
For a normal wavenumber $q$ define the two sums
$
  I(q) &= sum_(bold(i) in sans(M)_1) E_(y,bold(i)+1/2 bold(e)_y)^("num",-1/2) exp(-i(k_x x_bold(i) + q z_bold(i))), \
  P(q) &= sum_(bold(i) in sans(M)_1) E_(y,bold(i)+1/2 bold(e)_y)^("num",n_T-1/2) exp(-i(k_x x_bold(i) - q z_bold(i))),
$
where $n_T$ is the final step. The first sum uses the initial observation and picks the plane wave that travels toward the face.
The second sum uses the final observation and picks the plane wave that travels away from the face. The measured reflection magnitude is
$
  r^"num" (q) = abs(P(q))/abs(I(q)).
$ <E18-reflection>
We evaluate it at twelve wavenumbers:
$q_j = q_0 + j Q\/10$ for $j = -5, dots, 5$, and the Brewster wavenumber $q_"B" = k_x n_1\/n_2$ of E19.
Their angles of incidence are
#measured-waves.map(wave => rounded(angle-of(wave).deg(), digits: 1) + "°").join(", ", last: " and ").
The reference is $abs(Gamma(q_j))$.

*Purpose:* measure the reflection coefficient of the face as a function of the angle.
The numerical dispersion does not enter this ratio (deduction).
The reflection keeps the tangential wavenumber and the frequency.
On the mesh the reflected partner of the plane wave $(k_x, q)$ is $(k_x, -q)$, as in the continuum,
because the update in vacuum does not distinguish the two directions along $z$.
Between the two observations each spatial Fourier component only changes its phase.
The ratio is valid because side 1 contains only the incident packet at the initial observation
and only the reflected packet at the final observation.

=== E18.meas02: Field energy on each side
E16.meas03 with the cubic cells of E18, sampled every $N_x$ steps and at the final step.
The continuum values after the separation are $rho_1 = macron(R)$ and $rho_2 = 1 - macron(R)$ (@E18-energy-fractions).

*Purpose:* a check of the energy split that does not use Fourier sums.
It is less sharp than E18.meas01, because it averages over all angles.

=== E18.meas03: Field next to the second face
Let $sans(M)_"second"$ be the two layers of cells next to the second face, $k = 0$ and $k = N_z - 1$. Define
$
  b = max_n (max_(bold(i) in sans(M)_"second",c) abs(E_(c,bold(i)+1/2 bold(e)_c)^("num",n-1/2)))
    / (max_(bold(i) in sans(M),c) abs(E_(c,bold(i)+1/2 bold(e)_c)^("num",-1/2))),
$
where $n$ runs over the steps of E18.meas02.

*Purpose:* the continuum solution describes one face. This value shows how little of the field reaches the second face.

== Experiment E18 saved output
Each run writes the following files, and `run.json` with its configuration.

#table(
  columns: (1fr, 4fr),
  inset: 6pt,
  align: left + top,
  stroke: 0.5pt,
  [`reflection.csv`],
  [One row per wavenumber $q_j$: $q_j$, its angle, $abs(I(q_j))$, $abs(P(q_j))$, $r^"num" (q_j)$ and $abs(Gamma(q_j))$ (E18.meas01).],
  [`measurements.csv`],
  [The energies `U_side1` and `U_side2` with step and time (E18.meas02), and `field_at_second_face` (E18.meas03).],
  [`modes.csv`],
  [$q_m$, $omega_m$, $q'_m$, $a_m$ and $Gamma_m$ of every plane wave. The run record contains $macron(R)$.],
  [`profiles.csv`],
  [The numerical and the analytic $E_y$ along one line parallel to $z$ at the initial and at the final observation.],
  [`fields/`],
  [The full arrays `E` and `B_R2` at $0$, $T\/4$, $T\/2$, $3T\/4$ and $T$, for pictures of the $x z$ plane.],
)

== Experiment E18 analysis
We plot $r^"num" (q_j)$ against the angle, on the curve $abs(Gamma)$ from `modes.csv`.
For each $q_j$ and each pair of successive levels we report the observed order of $abs(r^"num" (q_j) - abs(Gamma(q_j)))$,
$
  p_"obs" = log_2 (e(h)\/e(h\/2))
$
(Salari and Knupp @salari2000CodeVerificationMms[§4.2, Eq. (25)]).
We also report $rho_1$ and $rho_2$ at the final step against $macron(R)$ and $1 - macron(R)$.
