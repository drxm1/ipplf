#import "../../docs/manual/preamble.typ": *
#import "guidedModeValues.typ": *
#import "@preview/cetz:0.3.4": canvas, draw

#let rounded(value, digits: 3) = str(calc.round(value, digits: digits))
#let te-modes = guided-modes("TE")
#let tm-modes = guided-modes("TM")
#let all-modes = te-modes + tm-modes
#let angles = all-modes.map(mode => mode.angle.deg())
#let levels = (24, 48, 96, 192) // N_x
#let duration = 30 // T in units of L_x / c_vac
#let safety-factor = 0.99
#let courant = safety-factor / calc.sqrt(3) // h_t / h
#let stride(n) = calc.quo(n, 8) // steps between two observations
// The smallest multiple of the stride whose R2 fields have a time of at least T.
#let last-step(n) = stride(n) * calc.ceil((duration * n / courant + 0.5) / stride(n))
#let long-steps = 1000000
#let long-stride = 64
// Fewest cubic cells per wavelength for an evanescent field with the phase velocity beta c.
#let gamma-beam = 118.417
#let beta-beam = calc.sqrt(1 - 1 / calc.pow(gamma-beam, 2))
#let cells-for-evanescent-field = 2 * calc.pi * gamma-beam * calc.sqrt((1 - calc.pow(courant, 2) * calc.pow(beta-beam, 4)) / 12)

= Experiment E22 detailed specs <experiment-e22>

E22 tests total internal reflection and the evanescent field in the vacuum.
A dielectric slab in the periodic box carries guided modes.
A guided mode is a wave that travels along the slab and is totally reflected at both interfaces.
In the vacuum its field decays with the distance from the slab, which is called an evanescent field.
Each guided mode is an exact solution in the periodic box, with one frequency and one profile.
We load one mode and measure its frequency, the decay constant of its evanescent field and the part of the field that is not the mode.
The six modes have angles of incidence from #rounded(calc.min(..angles), digits: 1)° to #rounded(calc.max(..angles), digits: 1)° inside the slab.
The critical angle is #rounded(critical-angle.deg(), digits: 2)° @ellingson2020ElectromagneticsVol2[§5.11, Eq. (5.260)].

In #link(<experiment-e18>)[E18] and #link(<experiment-e19>)[E19] the wave came from the vacuum, where no total reflection occurs.
E22 has the same combination of fields as one Fourier component of the field of a charge that moves along a dielectric:
a propagating field in the dielectric and an evanescent field in the vacuum.

== Experiment E22 analytic continuum reference
Lengths are in units of the period $L_x$ of the box along $x$, and $c_"vac" = 1$.
The slab has the permittivity $epsilon_r = #epsilon-slab$ and fills $abs(z) <= d\/2$ with $d = L_x$.
The permeability is one, and there is no charge and no current.
The box has the length $L_z = #calc.round(box-length) L_x$ along $z$, so the vacuum between the slab and its periodic image has the width $W = L_z - d$.
No quantity depends on $y$. Every mode has the phase
$
  psi = k_x x - omega t, quad k_x = 2 pi\/L_x,
$
and inside the slab the wave vector $(k_x, 0, plus.minus q)$, so its angle of incidence on the interfaces is $theta = arctan(k_x\/q)$.

=== Profile and mode condition
Let $n = sqrt(epsilon_r)$. For $k_x\/n < omega < k_x$ define
$
  q = sqrt(epsilon_r omega^2 - k_x^2), quad kappa = sqrt(k_x^2 - omega^2).
$ <E22-wavenumbers>
The profile of a mode is even or odd in $z$. For $abs(z) <= L_z\/2$ it is
$
  "even:" quad & f(z) = cases(
    cos(q z) & "if" abs(z) <= d\/2,
    C cosh(kappa (L_z\/2 - abs(z))) & "if" abs(z) > d\/2\,
  ) \
  "odd:" quad & f(z) = cases(
    sin(q z) & "if" abs(z) <= d\/2,
    "sgn"(z) thin C sinh(kappa (L_z\/2 - abs(z))) & "if" abs(z) > d\/2\,
  )
$ <E22-profile>
with $C = cos(q d\/2)\/cosh(kappa W\/2)$ for an even and $C = sin(q d\/2)\/sinh(kappa W\/2)$ for an odd mode.
Then $f$ is continuous, and $f$ and $f'$ are continuous across the periodic seam at $z = plus.minus L_z\/2$.
The frequency follows from the mode condition
$
  sigma thin q tan(q d\/2) = kappa tanh(kappa W\/2) quad "(even)", quad
  -sigma thin q cot(q d\/2) = kappa coth(kappa W\/2) quad "(odd)",
$ <E22-mode-condition>
with $sigma = 1$ for the TE case and $sigma = 1\/epsilon_r$ for the TM case.
The hyperbolic functions contain the periodic images of the slab exactly.
Near the slab, where the image is negligible, the vacuum field is $prop e^(-kappa (abs(z) - d\/2))$ with the decay length $1\/kappa$.

=== Potentials and fields
In both cases $phi = 0$.

*TE case* (the electric field is parallel to the interfaces):
$
  vb(A) = vb(hat(y)) thin f(z)/omega sin psi, quad
  E_y = f cos psi, quad B_x = -f'/omega sin psi, quad B_z = k_x/omega f cos psi.
$ <E22-te>
*TM case* (the magnetic field is parallel to the interfaces), with $f$ as the profile of $B_y$ and $epsilon(z)$ the permittivity at $z$:
$
  A_x & = -f'/(epsilon omega^2) cos psi, quad A_z = -(k_x f)/(epsilon omega^2) sin psi, \
  B_y & = f cos psi, quad E_x = f'/(epsilon omega) sin psi, quad E_z = -(k_x f)/(epsilon omega) cos psi.
$ <E22-tm>
The fields are $vb(E) = -partial_t vb(A)$ and $vb(B) = curl vb(A)$ @chew2014generalizedgauge[§2, Eqs. (5)–(6)].
In the TM case $E_z$ and $A_z$ are normal to the interfaces and jump there by the factor $epsilon_r$.
On an interface the reference value of $E_z$ and of $A_z$ is the mean of the two limits.
This is the continuous product $epsilon E_z$ divided by the harmonic mean of the two permittivities (@eq-O1VS-epsd).

=== The six modes
For this slab @E22-mode-condition has three solutions in each case.
We name a mode by its case and by the number $m$ of zeros of its profile inside the slab.

#figure(
  table(
    columns: 8,
    inset: 5pt,
    align: (left, center, right, right, right, right, right, right),
    stroke: 0.4pt,
    table.header(
      [Mode], [Parity], [$omega\/k_x$], [$theta$], [$theta - theta_"c"$], [$q L_x$], [$kappa L_x$], [$1\/(kappa L_x)$],
    ),
    ..all-modes
      .map(mode => (
        [#mode.name],
        [#mode.parity],
        [#rounded(mode.omega / k-x, digits: 5)],
        [#rounded(mode.angle.deg(), digits: 2)°],
        [#rounded(mode.angle.deg() - critical-angle.deg(), digits: 2)°],
        [#rounded(mode.q)],
        [#rounded(mode.kappa)],
        [#rounded(1 / mode.kappa)],
      ))
      .flatten(),
  ),
  caption: [The guided modes of E22, computed from @E22-mode-condition. $theta_"c"$ is the critical angle and $1\/kappa$ the decay length in the vacuum.],
) <E22-modes>

// Profiles of the three TE modes near the slab, to scale along z. The profiles have the
// amplitude 1 inside the slab.
#let profile-sketch = canvas(length: 1mm, {
  import draw: *
  let (z-scale, lane, amplitude) = (22, 15, 6)
  let z-range = 3.2
  rect((-slab-thickness / 2 * z-scale, -lane / 2), (slab-thickness / 2 * z-scale, 2.5 * lane), fill: luma(88%), stroke: none)
  for (index, mode) in te-modes.enumerate() {
    let base = (2 - index) * lane
    line((-z-range * z-scale, base), (z-range * z-scale, base), stroke: 0.3pt + gray)
    line(
      ..range(321).map(i => {
        let z = -z-range + i * z-range / 160
        (z * z-scale, base + amplitude * mode-profile(mode, z))
      }),
      stroke: 0.7pt + black,
    )
    content((-z-range * z-scale - 2, base), anchor: "east", text(size: 7.5pt)[#mode.name, #rounded(mode.angle.deg(), digits: 1)°])
  }
  for z in (-3, -2, -1, 0, 1, 2, 3) {
    content((z * z-scale, -lane / 2 - 2.5), text(size: 7.5pt)[$#z$])
  }
  content((z-range * z-scale + 2, -lane / 2 - 2.5), anchor: "west", text(size: 7.5pt)[$z\/L_x$])
  content((0, 2.5 * lane + 2.5), text(size: 7.5pt)[slab, $epsilon_r = #epsilon-slab$])
  content((2 * z-scale, 2.5 * lane + 2.5), text(size: 7.5pt)[vacuum])
  content((-2 * z-scale, 2.5 * lane + 2.5), text(size: 7.5pt)[vacuum])
})
#figure(
  profile-sketch,
  caption: [The profiles $f(z)$ of the three TE modes near the slab, to scale.
    The box continues to $z = plus.minus #calc.round(box-length / 2) L_x$. The mode nearest to the critical angle has the longest evanescent field.],
) <E22-profiles>

=== Proof: The modes solve the equations of the periodic box
*Gauge.* In the TE case $div(epsilon vb(A)) = epsilon partial_y A_y = 0$.
In the TM case $epsilon A_x = -(f'\/omega^2) cos psi$ and $epsilon A_z = -(k_x f\/omega^2) sin psi$,
so $div(epsilon vb(A)) = (k_x f'\/omega^2) sin psi - (k_x f'\/omega^2) sin psi = 0$, also across an interface.
The generalized gauge @chew2014generalizedgauge[Eq. (15)] then holds with $phi = 0$.

*Wave equation.* Without the gauge term the vector equation @chew2014generalizedgauge[Eq. (17)] is
$-curl curl vb(A) - epsilon partial_t^2 vb(A) = vb(0)$ in each medium.
For both cases this is $f'' + (epsilon omega^2 - k_x^2) f = 0$.
By @E22-wavenumbers, $f'' = -q^2 f$ in the slab and $f'' = kappa^2 f$ in the vacuum, which @E22-profile satisfies.

*Conditions at the interfaces.* The tangential component of $vb(A)$, the tangential component of $vb(B)$ and $epsilon A_z$ are continuous
@chew2014generalizedgauge[§3, Eqs. (23), (24), (26)].
In the TE case these are the continuity of $f$ and of $f'$.
In the TM case they are the continuity of $f$ and of $f'\/epsilon$.
With @E22-profile the first holds by the choice of $C$, and the second is @E22-mode-condition.

*Periodicity.* The profile is even or odd about the center of the vacuum gap as well,
so $f$ and $f'$ are continuous across $z = plus.minus L_z\/2$. The dependence on $x$ has the period $L_x$.

== Experiment E22 solver configuration, domain and refinement
#table(
  columns: (1fr, 3.2fr),
  inset: 6pt,
  align: left + top,
  stroke: 0.5pt,
  [Solver], [O1VS with $alpha = 1$, Standard stencil, staggered placement, R2, periodic in all three directions],
  [Domain], [$L_x times 2h times #calc.round(box-length) L_x$ with cubic cells, $h = L_x\/N_x$, so $N_x times 2 times #calc.round(box-length) N_x$ cells],
  [Slab],
  [The $N_x$ cell layers with $15 N_x\/2 <= k < 17 N_x\/2$, where $k = 0$ is the lowest layer. Both interfaces lie on cell faces.],
  [Levels], [$N_x = #levels.map(str).join(", ")$],
  [Time step],
  [$h_t = #safety-factor thin h\/sqrt(3)$, selected with `ChdrSolver::setSafetyFactor(0.99)`, as in #link(<experiment-e16s>)[E16-S]],
  [Observations], [At step 0 and at every positive multiple of $N_x\/8$ up to the last step],
  [Last step], [$n_T$, the smallest multiple of $N_x\/8$ whose R2 fields have a time of at least $T = #duration L_x\/c_"vac"$:
    #levels.map(n => [#last-step(n)]).join(", ") steps],
  [Cases], [One run for each mode of @E22-modes on each level],
  [Long runs], [All six modes on $N_x = #levels.first()$ with $10^6$ steps, observed at step 0 and at every #long-stride\th step.
    The mean of E22.meas02 uses every multiple of $N_x\/8$ also in these runs.],
)
The time step is not a binary fraction. The program selects every observation by its step index and computes each time from that index.
The observation times of two levels differ by half a time step, because the R2 fields of step $n$ have the time $(n - 1\/2) h_t$.

The slab has between #rounded(levels.first() / calc.sqrt(epsilon-slab), digits: 1) and #rounded(levels.last() / calc.sqrt(epsilon-slab), digits: 0) cells
for each wavelength of a plane wave with the frequency $k_x$.
On the coarsest level the shortest decay length, that of TE0, is #rounded(levels.first() / te-modes.first().kappa, digits: 1) cells.

== Experiment E22 initialization
As in #link(<experiment-e18>)[E18]: we sample the potentials of @E22-te or @E22-tm at their native positions,
$vb(A)$ at $t = 0$ and at $t = -h_t$. The scalar potential is zero in both slots.

The loaded state is the continuum mode, which is not exactly a mode of the mesh.
It contains the mode of the mesh, other modes of the mesh with the same $k_x$,
and a weak wave that travels along $-x$, because the two slots have the continuum frequency.
In the periodic box no part of it can leave.

== What E22 exercises in the update
- In the TE case only $A_y$ is nonzero. Its update is the wave equation of one component with the permittivity of the cell, as in E18.
- In the TM case $A_x$ and $A_z$ are nonzero. One layer of $A_z$ samples lies on each interface and uses the harmonic mean of the two permittivities (@eq-O1VS-epsd), as in E19.
- New in E22: the angle of incidence is beyond the critical angle, the vacuum carries an evanescent field,
  and the field stays at the two interfaces during the whole run. In E18 and E19 a packet met the interface once.
- The scalar potential stays zero, as in E19.

== Expected behaviour
The statements of this section are deductions. Their numbers come from a calculation of the modes of the update equations
(@eq-O1VS-ThAd) on the four meshes, outside the solver.

*Frequency.* All six modes exist on all four meshes.
The frequency $omega_h$ of a mode on the mesh is lower than $omega$,
by about $2 times 10^(-3) omega$ on the coarsest mesh, and the difference has the order two.

*Decay constant.* In the vacuum every field component obeys the wave equation of the Standard stencil.
For a field $prop e^(i (k_x x - omega_h t)) e^(-kappa_h z)$ its dispersion relation @fallahi2020mithra20fullwavesimulation[§3.1, Eq. (3.17)]
with the wavenumber $i kappa_h$ along $z$ is
$
  sinh^2 (kappa_h h)/2 = sin^2 (k_x h)/2 - (h/h_t)^2 sin^2 (omega_h h_t)/2.
$ <E22-stencil-relation>
To the order $h^2$, with $S = h_t\/h$ and $e_omega = omega_h\/omega - 1$,
$
  (kappa_h^2 - kappa^2)/kappa^2 approx -2 (omega/kappa)^2 e_omega - h^2/(12 kappa^2) (k_x^4 - S^2 omega^4 + kappa^4).
$
The two terms have opposite signs. Near the critical angle each of them is large, and they nearly cancel.
The decay constant of a mode on the mesh differs from $kappa$ by less than one percent on the coarsest mesh, also near the critical angle.

*Residual.* A difference across an interface is of first order, because $f''$ jumps there:
$
  (f(z_f + h\/2) - f(z_f - h\/2))/h - f'(z_f) = h/8 (f''(z_f^+) - f''(z_f^-)) + cal(O)(h^2).
$
The two cell layers at the interfaces are a fraction of the order $h$ of the box.
At step zero the residual of E22.meas03, which is a root of an energy, therefore has the order $3\/2$.
During the run the other modes of the mesh add a part of the order two, which is the larger part on these meshes.

*A field with a prescribed frequency.* For the component with $k_y = 0$ of the field of a charge with the velocity $beta c_"vac"$ along $x$,
$omega = beta k_x$ is prescribed and $kappa = k_x\/gamma$ with $gamma = (1 - beta^2)^(-1\/2)$.
To the order $h^2$, @E22-stencil-relation then has a positive right side only for more than
$2 pi gamma sqrt((1 - S^2 beta^4)\/12)$ cubic cells per wavelength $2 pi\/k_x$, which is about #calc.round(cells-for-evanescent-field) for $gamma = #gamma-beam$.
With fewer cells the mesh has no evanescent field for this component.
E22 does not test a prescribed frequency. It tests @E22-stencil-relation with the frequency of a guided mode.

== Experiment E22 measured quantities and their purpose
After $n$ completed steps both fields have the time $t = (n - 1\/2) h_t$.
Each component is sampled at its own position (@sec-spatial-placement), and the reference at the same position and time.
The complex reference fields $cal(E)$ and $cal(B)$ are @E22-te and @E22-tm with $cos psi$ replaced by $e^(i psi)$ and $sin psi$ by $-i e^(i psi)$,
so that the real fields are their real parts.
$sum_(bold(i),c)$ covers all cells and the three components, and each sample is counted once.
The energy of two fields is
$
  U(vb(E), vb(B)) = h^3/2 sum_(bold(i),c)
    (epsilon_(c,bold(i)+bold(delta)_c^E) abs(E_(c,bold(i)+bold(delta)_c^E))^2 + abs(B_(c,bold(i)+bold(delta)_c^B))^2),
$
with the permittivity that the update uses for each component, as in E19.meas02.

=== E22.meas01: Mode coefficients
As E15.meas01 and E15.meas02 (`E15.typ`, "Measured quantities"). For $F in {E, B}$ and each observation,
$
  c_F = (2 sum_(bold(i),c) F_(c,bold(i)+bold(delta)_c^F)^"num" thin overline(cal(F)_(c,bold(i)+bold(delta)_c^F)))
        / (sum_(bold(i),c) abs(cal(F)_(c,bold(i)+bold(delta)_c^F))^2), quad
  rho_F = abs(c_F), quad delta_F = arg(c_F).
$
The phases are unwrapped in the order of the observations.
A mode with the frequency $omega_h$ gives $delta_F = (omega - omega_h) t + "const"$.
With the observations $(t_j, delta_(E,j))$ and their means $macron(t)$ and $macron(delta)$,
$
  lambda = (sum_j (t_j - macron(t)) (delta_(E,j) - macron(delta))) / (sum_j (t_j - macron(t))^2), quad
  e_omega = (omega_h - omega)/omega = -lambda/omega.
$

*Purpose:* $e_omega$ is the error of the frequency of the guided mode.
It contains the dispersion error of the wave in the slab and in the vacuum and the error of the phase of the total reflection. E22 does not separate them.
$rho_F$ is the amplitude of the projection of the numerical field on the continuum mode.
It falls if the field loses energy, and also if its shape changes or if other modes beat with it.

=== E22.meas02: Decay constant of the evanescent field
Let $F = E_y$ in the TE case and $F = B_y$ in the TM case, and let $(x_i, z_k)$ be the position of its sample in the cell $(i, j, k)$.
At each observation and for each cell layer $k$,
$
  a_k = 2/(N_x N_y) sum_(i,j) F_(i j k)^"num" e^(-i k_x x_i), quad b_k = (a_k e^(i omega t))/c_F,
$
with $c_F$ of E22.meas01 for the same field. The program saves the mean $p_k$ of $b_k$ over step 0 and every multiple of $N_x\/8$ up to the last step, with equal weights.
For the continuum mode $b_k = f(z_k)$ at every time.
For a mode of the mesh $b_k$ is its profile, and the other parts of the field oscillate in $b_k$ and are reduced in the mean.
The mean needs this dense sequence of steps. With every #long-stride\th step it is wrong, because other modes then appear at the frequency zero.

Let $s_k = z_k - d\/2$ be the distance of a layer from the upper interface,
and let $k_1$ and $k_2$ be the first layers with $s_k >= L_x\/8$ and with $s_k >= L_x\/8 + 1\/kappa$.
The decay constant $kappa_h$ on the mesh is the solution of
$
  (op("Re") p_(k_2))/(op("Re") p_(k_1)) = (cosh(kappa_h (W\/2 - s_(k_2))))/(cosh(kappa_h (W\/2 - s_(k_1))))
$
for an even mode, and of the same equation with $sinh$ for an odd mode. Its error is $e_kappa = kappa_h\/kappa - 1$.
We compare it with $e_kappa^"st" = kappa_h^"st"\/kappa - 1$, where $kappa_h^"st"$ is the solution of @E22-stencil-relation with $omega_h$ of E22.meas01.

*Purpose:* $kappa_h$ is the decay constant of the evanescent field on the mesh.
If it agrees with $kappa_h^"st"$, the evanescent field of the mode is that of the vacuum update at the frequency of the mode.
An error of the interfaces can still change this frequency, and with it the decay constant.
Near the critical angle the start excites modes that beat with the guided mode with periods of 100 to 200 $L_x\/c_"vac"$.
The mean over $T$ reduces them only in part, so the comparison is most precise in the long runs.

=== E22.meas03: Residual outside the mode and total energy
As E15.meas03. With the fitted fields subtracted before the sum,
$
  r = sqrt(U(vb(E)^"num" - op("Re")(c_E cal(E)), thin vb(B)^"num" - op("Re")(c_B cal(B))) / U^"ref"), quad
  U^"ref" = U(op("Re") cal(E), op("Re") cal(B)).
$
The change of the total energy is $U(vb(E)^"num", vb(B)^"num")\/U_0 - 1$, with $U_0$ the value at step zero.

*Purpose:* $r$ is the part of the field that remains after the two fits, relative to the energy of the continuum mode.
The two fits are independent, so $r$ does not test the ratio of the electric and the magnetic amplitude.
It contains the difference between the profile on the mesh and @E22-profile, and the other modes of the mesh that the start excites.
We expect the order $3\/2$ at step zero, a decrease under refinement during the run, and no growth over the long runs.
The total energy shows a loss or a growth of the whole field. It cannot show a loss at one reflection.

== Experiment E22 saved output
The inputs are `E22-<mode>-n<N_x>.yaml`, for example `E22-TM2-n96.yaml`, and `E22-<mode>-long-n24.yaml` for the long runs.
Each run writes the following files.

#table(
  columns: (1fr, 4fr),
  inset: 6pt,
  align: left + top,
  stroke: 0.5pt,
  [`observations.csv`],
  [
    One row per observation with the columns `step`, `time`, `c_E_re`, `c_E_im`, `c_B_re`, `c_B_im` (E22.meas01),
    `U`, `U_reference` and `r` (E22.meas03). The time is that of the R2 fields.
  ],
  [`profile.csv`],
  [
    One row per cell layer $k$ with the columns `k`, `z`, `mean_re`, `mean_im` and `reference`:
    the height $z_k$ of the sample above the center of the slab, $p_k$ of E22.meas02 and $f(z_k)$.
    It is written once, at the end of the run.
  ],
  [`run.json`],
  [
    The configuration of the run.
    Under `effective.analytic_reference` it holds `mode`, `parity`, `omega`, `q`, `kappa`, `angle_degrees`, `slab_thickness` and `box_length`.
    Under `effective.observations` it holds `every`, the number of steps between two observations, and `last_step`.
    Under `effective.solver` it holds `time_step`.
  ],
)

== Experiment E22 analysis
For each mode we plot $abs(e_omega)$, $abs(e_kappa)$ with $abs(e_kappa^"st")$, and the largest $r$ against $N_x$.
For each mode we tabulate, on the finest level, $e_omega$, $r$ at step zero, the largest $r$,
the largest $abs(rho_F - 1)$ over both fields and the largest $abs(U\/U_0 - 1)$ over the observations,
and the observed orders of $e_omega$ and of $r$ at step zero between all successive levels.
An order uses the magnitudes of two errors. We state it if one of them is zero or if their signs differ.
We plot $abs(p_k)$ against $z$ on a logarithmic axis, above $abs(f(z_k))$.
For each long run we report $e_kappa$ and $e_kappa^"st"$, the largest sampled $abs(rho_F\/rho_F^0 - 1)$ over both fields,
with $rho_F^0$ the value at step zero, the largest sampled $r$, the ratio of the largest $r$ of the last tenth of the run to that of the first tenth,
and the largest sampled $abs(U\/U_0 - 1)$.

== What E22 does not test
E22 has no incidence from the dielectric below the critical angle, no absorbing boundary and no source.
The vacuum gap is wide. A narrow gap, through which two slabs exchange energy, is a separate experiment.
The six modes have $k_y = 0$ and one $k_x$.
The modes are at least #rounded(calc.min(..angles) - critical-angle.deg(), digits: 1)° from the critical angle,
and the component with $k_y = 0$ of the field of a charge with $gamma = #gamma-beam$ is $0.002 degree$ from it.
