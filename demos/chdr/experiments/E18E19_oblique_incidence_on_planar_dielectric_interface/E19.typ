#import "../../docs/manual/preamble.typ": *
#import "obliquePacketValues.typ": *

= Experiment E19 detailed specs <experiment-e19>

E19 repeats #link(<experiment-e18>)[E18] with the other polarization:
the electric field lies in the plane of incidence (the TM case), and the magnetic field has only the component $B_y$.
The electric field now has a component normal to the face.
The reflection coefficient has a zero at the Brewster angle, which lies inside the range of angles of the packet.
This is the first experiment in which the O1VS update differs from an update that only scales each component by the cell permittivity.

Everything that is not stated otherwise below is the same as in E18:
the plane waves and their phases (@E18-modes, @E18-phases), the amplitudes and the parameters of the packet (@E18-amplitudes),
the domain, the refinement levels, the duration and the saved output.

== Experiment E19 analytic continuum solution
Let $k_(j,m) = n_j omega_m$ be the wavenumber of plane wave $m$ in medium $j$. We define three unit vectors in the $x z$ plane,
$
  vu(e)_("i",m) = 1/k_(1,m) (q_m, 0, -k_x), quad
  vu(e)_("r",m) = 1/k_(1,m) (-q_m, 0, -k_x), quad
  vu(e)_("t",m) = 1/k_(2,m) (q'_m, 0, -k_x).
$ <E19-polarization>
Each is perpendicular to the wave vector of its plane wave,
which is $(k_x, 0, q_m)$, $(k_x, 0, -q_m)$ and $(k_x, 0, q'_m)$.
The potentials are
$
  vb(A)(x,z,t) = sum_(m=-M)^M a_m cases(
    vu(e)_("i",m) cos psi_("i",m) + R_m vu(e)_("r",m) cos psi_("r",m) & "if" z < z_f,
    S_m vu(e)_("t",m) cos psi_("t",m) & "if" z >= z_f
  ), quad phi = 0,
$ <E19-analytic-potentials>
$
  R_m = R(q_m), quad S_m = n_1/n_2 (1 + R_m), quad
  R(q) = (epsilon_2 q - epsilon_1 q')/(epsilon_2 q + epsilon_1 q').
$ <E19-coefficients>
The fields are
$
  vb(E) &= -sum_m a_m omega_m cases(
    vu(e)_("i",m) sin psi_("i",m) + R_m vu(e)_("r",m) sin psi_("r",m),
    S_m vu(e)_("t",m) sin psi_("t",m)
  ), \
  B_y &= -sum_m a_m cases(
    k_(1,m) (sin psi_("i",m) + R_m sin psi_("r",m)),
    k_(2,m) S_m sin psi_("t",m)
  ),
$ <E19-analytic-EB>
with the upper line for $z < z_f$ and the lower line for $z >= z_f$.
$E_y$, $B_x$ and $B_z$ are zero.

$R_m$ is the reflection coefficient of the magnetic field $B_y$, and $1 + R_m$ is its transmission coefficient.
With $eta_j = 1\/n_j$ and the cosines of E18, $R_m = -Gamma_"TM"$ of Ellingson @ellingson2020ElectromagneticsVol2[§5.7, Eqs. (5.184), (5.187)].
The sign differs because Ellingson counts the reflected magnetic field in the opposite direction @ellingson2020ElectromagneticsVol2[§5.7, Eq. (5.166)].
$R$ is zero at the Brewster angle, $tan theta_"B" = n_2\/n_1$ @ellingson2020ElectromagneticsVol2[§5.10, Eq. (5.253)],
which is #rounded(angle-of(brewster-wave).deg(), digits: 1)° here and corresponds to $q_"B" = k_x n_1\/n_2$.
At the central angle of $45 degree$, $R = #rounded(r-tm(plane-wave(q-center)))$.

=== Proof: Each plane wave solves the equations of the two media
*Gauge.* In each medium every term of $vb(A)$ is perpendicular to its wave vector, so $div vb(A) = 0$ there.
At the face $epsilon_r A_z$ is continuous (shown below).
Then $div (epsilon_r vb(A)) = 0$ everywhere, and the generalized gauge @chew2014generalizedgauge[Eq. (15)] holds with $phi = 0$.

*Wave equation.* Without the gauge term the vector equation @chew2014generalizedgauge[Eq. (17)] is
$-curl curl vb(A) - epsilon_j partial_t^2 vb(A) = vb(0)$ in medium $j$.
For a plane wave that is perpendicular to its wave vector, $curl curl vb(A) = k_(j,m)^2 vb(A)$,
and $k_(j,m)^2 = epsilon_j omega_m^2$.

*Conditions at the face.* $A_x$, $B_y = (curl vb(A))_y$ and $epsilon_r A_z$ are continuous
@chew2014generalizedgauge[§3, Eqs. (23), (24), (26)]. With equal phases at the face these are
$
  q_m/k_(1,m) (1 - R_m) = q'_m/k_(2,m) S_m, quad
  k_(1,m) (1 + R_m) = k_(2,m) S_m, quad
  epsilon_1/k_(1,m) (1 + R_m) = epsilon_2/k_(2,m) S_m.
$
The second and the third condition are the same because $k_(j,m)^2 = epsilon_j omega_m^2$.
The first and the second give @E19-coefficients.

=== Energy fractions
As in E18,
$
  macron(R)_"TM" = (sum_m a_m^2 omega_m^2 R_m^2)/(sum_m a_m^2 omega_m^2),
$ <E19-energy-fractions>
and the continuum values after the separation are $rho_1 = macron(R)_"TM"$ and $rho_2 = 1 - macron(R)_"TM"$.
For E19, $macron(R)_"TM" approx #rounded(reflected-fraction-tm)$, against $macron(R) approx #rounded(reflected-fraction-te)$ in E18.

== What E19 exercises in the update
The vector potential has the components $A_x$ and $A_z$.
- $A_z$ is normal to the face.
  Its samples lie at $(x_i, y_j, z_k + h\/2)$, and one layer of them lies on the face.
  There the update uses the harmonic mean of $epsilon_1$ and $epsilon_2$ (`HarmonicFace`).
- The continuum divergence $div (epsilon_r vb(A))$ is zero, but its discrete form is not zero identically.
  The gauge term of the update acts on this discretization error.
- The scalar potential is not exercised.
  Its update does not depend on $vb(A)$, so it stays zero.

An update that only scales the wave equation of each component by the cell permittivity treats $A_x$ and $A_z$ like $A_y$ in E18.
It then reflects each plane wave with $abs(Gamma(q))$ of E18 instead of $abs(R(q))$ (deduction from @E18-coefficients and @E19-coefficients).
At the Brewster angle this is $abs(Gamma) = #rounded(calc.abs(gamma-te(brewster-wave)))$ instead of zero,
and at $45 degree$ it is #rounded(calc.abs(gamma-te(plane-wave(q-center)))) instead of #rounded(calc.abs(r-tm(plane-wave(q-center)))).
E19.meas01 therefore distinguishes the two updates at every angle and on every level.

== Experiment E19 initialization
We sample @E19-analytic-potentials at $t = 0$ and at $t = -h_t$:
$A_x$ at $(x_i + h\/2, y_j, z_k)$ and $A_z$ at $(x_i, y_j, z_k + h\/2)$.
$A_y$ and $phi$ are zero in both slots.
$A_z$ is discontinuous at the face, and a sample of $A_z$ on the face uses the expression for $z < z_f$.
At these two times the field at the face is below $10^(-4)$ of the peak.
The samples of the continuum potentials do not satisfy the discrete gauge exactly.
The residual is a discretization error of the initial state, as in E03.

== Experiment E19 measured quantities

=== E19.meas01: Reflection magnitude for each normal wavenumber
E18.meas01 with the samples of $B_y$ in place of those of $E_y$.
Cell $bold(i)$ stores the R2 sample $B_(y,bold(i)+bold(delta)_y^B)$ at $(x_bold(i), z_bold(i)) = (x_i + h\/2, z_k + h\/2)$:
$
  I(q) &= sum_(bold(i) in sans(M)_1) B_(y,bold(i)+bold(delta)_y^B)^("num",-1/2) exp(-i(k_x x_bold(i) + q z_bold(i))), \
  P(q) &= sum_(bold(i) in sans(M)_1) B_(y,bold(i)+bold(delta)_y^B)^("num",n_T-1/2) exp(-i(k_x x_bold(i) - q z_bold(i))), quad
  r^"num" (q) = abs(P(q))/abs(I(q)).
$
The twelve wavenumbers $q_j$ are those of E18. The reference is $abs(R(q_j))$, which is zero at $q_"B"$.

*Purpose:* as E18.meas01, for the TM coefficient.
We use $B_y$ because it is the only magnetic component and it is continuous at the face.
The samples of $B_y$ with $z_k + h\/2 = z_f$ lie on the face and belong to side 1.

=== E19.meas02: Field energy on each side with the permittivity of each component
Let $epsilon_(c,bold(i)+1/2 bold(e)_c)$ be the permittivity that the update uses for component $c$ at its own position,
which is the harmonic mean of the two cells that share that face (`HarmonicFace`).
E19.meas02 is E18.meas02 with this permittivity in place of the cell permittivity:
$
  U_S^(n-1/2) = h^3/2 sum_(bold(i) in sans(M)_S,c)
    (epsilon_(c,bold(i)+1/2 bold(e)_c) abs(E_(c,bold(i)+1/2 bold(e)_c)^("num",n-1/2))^2
     + abs(B_(c,bold(i)+bold(delta)_c^B)^("num",n-1/2))^2).
$
The two definitions differ only for the samples of $E_z$ on the face, where $E_z$ is discontinuous.
The continuum values after the separation are $rho_1 = macron(R)_"TM"$ and $rho_2 = 1 - macron(R)_"TM"$.

E18.meas03 is used unchanged.

== Experiment E19 saved output and analysis
As in E18.
`reflection.csv` holds $abs(R(q_j))$ as the reference, `modes.csv` lists $R_m$ instead of $Gamma_m$,
and `profiles.csv` holds $B_y$, $E_x$ and $E_z$.
We plot $r^"num" (q_j)$ against the angle on the curve $abs(R)$, together with the curve $abs(Gamma)$ of E18.
We report the observed order of $abs(r^"num" (q_j) - abs(R(q_j)))$ for each $q_j$, including $q_"B"$,
and $rho_1$ and $rho_2$ at the final step.
