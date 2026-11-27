#import "../../docs/manual/preamble.typ": *

= Experiment E20 detailed specs <experiment-e20>

E20 repeats #link(<experiment-e16>)[E16] with a pulse of the potentials that carries no field.
Both potentials are derivatives of one function $Lambda$, so that $vb(E) = vb(B) = vb(0)$.
The scalar potential and the component of $vb(A)$ normal to the face are nonzero.
At the face this pulse is reflected and transmitted with its own coefficients.
The main result is that the reconstructed $vb(E)$ stays zero while the pulse crosses the face.
The second result is the accuracy of the scalar potential behind the face:
its split between the two sides and its error on each side.

Everything that is not stated otherwise below is the same as in E16:
the media, the phases (@E16-phases), the profile $a$ (@E16-profile), the solver configuration,
the domain, the position of the face, the refinement levels, the time step, the duration and the observation times.

== Experiment E20 analytic continuum reference

The potentials are
$
  vb(A) = grad Lambda, quad phi = -partial_t Lambda, quad
  Lambda(z,t) = cases(
    a(s_"i") + Gamma_phi a(s_"r") & "if" z < z_f,
    tau_phi a(s_"t") & "if" z >= z_f
  ),
$ <E20-gauge-function>
with
$
  Gamma_phi = (n_1^3 - n_2^3)/(n_1^3 + n_2^3), quad
  tau_phi = 1 + Gamma_phi = (2 n_1^3)/(n_1^3 + n_2^3).
$ <E20-coefficients>
Written out with @E16-phases, they are
$
  phi(z,t) = v_1 cases(
    a'(s_"i") + Gamma_phi a'(s_"r") & "if" z < z_f,
    tau_phi a'(s_"t") & "if" z >= z_f
  ), quad
  vb(A)(z,t) = vb(hat(z)) cases(
    a'(s_"i") - Gamma_phi a'(s_"r") & "if" z < z_f,
    (n_2\/n_1) tau_phi a'(s_"t") & "if" z >= z_f
  ).
$ <E20-analytic-potentials>
The fields are $vb(E) = -partial_t grad Lambda + grad partial_t Lambda = vb(0)$ and $vb(B) = curl grad Lambda = vb(0)$
@chew2014generalizedgauge[§2, Eqs. (5)–(6)].

#let n-dielectric = calc.sqrt(2.13)
#let n-cubed = calc.pow(n-dielectric, 3)
#let gamma-phi = (1 - n-cubed) / (1 + n-cubed)
#let rounded(value, digits: 4) = str(calc.round(value, digits: digits))

For E20, $Gamma_phi approx #rounded(gamma-phi)$ and $tau_phi approx #rounded(1 + gamma-phi)$.
The reflection is stronger than that of the electric field in E16 (@E16-coefficients).

=== Proof: The potentials solve the equations of the two media
*Gauge and scalar equation.*
With @E20-gauge-function the generalized gauge $div (epsilon_r vb(A)) = -chi partial_t phi$ @chew2014generalizedgauge[Eq. (15)] is
$
  partial_z (epsilon_r partial_z Lambda) = chi partial_t^2 Lambda,
$ <E20-gauge-equation>
and the scalar equation @eq-Chew-scalar without charge is its time derivative.
In medium $j$, where $chi = epsilon_j^2$, this is $partial_z^2 Lambda = epsilon_j partial_t^2 Lambda$.
The three terms of $Lambda$ satisfy it, as the three terms of $A_x$ do in E16.

*Vector equation.*
In @eq-Chew-vector the curl term is zero.
By @E20-gauge-equation the gauge term is
$epsilon_r grad (chi^(-1) div (epsilon_r grad Lambda)) = epsilon_r grad partial_t^2 Lambda = epsilon_r partial_t^2 vb(A)$.

*Conditions at the face.*
$phi$ and $epsilon_r partial_z phi$ are continuous @chew2014generalizedgauge[§3, Eqs. (27), (29)],
and so is the normal product $epsilon_r A_z$ @chew2014generalizedgauge[§3, Eq. (26)].
All three hold if $Lambda$ and $epsilon_r partial_z Lambda$ are continuous.
With equal phases at the face and $epsilon_j = n_j^2$, these two conditions are
$
  1 + Gamma_phi = tau_phi, quad n_1^3 (1 - Gamma_phi) = n_2^3 tau_phi,
$
and @E20-coefficients is their solution.
$A_z$ itself jumps at the face, $A_z (z_f^+) = (epsilon_1\/epsilon_2) A_z (z_f^-)$.

=== The split of the scalar potential
After the two pulses have separated, @E20-analytic-potentials and $dif s_"t" = (n_2\/n_1) dif z$ give
$
  integral_(z < z_f) chi phi^2 dif z = Gamma_phi^2 Q_"i", quad
  integral_(z >= z_f) chi phi^2 dif z = (n_2/n_1)^3 tau_phi^2 Q_"i" = (1 - Gamma_phi^2) Q_"i",
$ <E20-fractions>
where $Q_"i"$ is the same integral of the incident pulse.
For E20, $Gamma_phi^2 approx #rounded(gamma-phi * gamma-phi)$ and $1 - Gamma_phi^2 approx #rounded(1 - gamma-phi * gamma-phi)$.

== Experiment E20 initialization
We build the histories from samples of $Lambda$ at the cell centers,
$Lambda_k^n := Lambda(z_k, n h_t)$ for $n in {0, -1}$.
At these two times only the incident term is nonzero, as in E16.
We set $A_x = A_y = 0$ and
$
  A_(z,k+1/2)^n = (Lambda_(k+1)^n - Lambda_k^n)/h_z quad "for" n in {0, -1},
$
$
  phi_k^(-1/2) = -(Lambda_k^0 - Lambda_k^(-1))/h_t, quad
  phi_k^(1/2) = phi_k^(-1/2) - h_t/(epsilon_1 h_z^2) (Lambda_(k+1)^0 - 2 Lambda_k^0 + Lambda_(k-1)^0).
$ <E20-initialization>
The samples do not depend on the transverse indices.
The last expression is one step of the scalar update in the first medium, applied to $Lambda$ (see the next section).
With these histories the reconstructed electric field is zero in exact arithmetic.
Samples of @E20-analytic-potentials themselves would leave an $E_z$ of the order $h_z^2$ from the start.

== What E20 exercises in the update
#link(<experiment-e19>)[E19] exercised the component of $vb(A)$ normal to the face
and the harmonic mean of two different permittivities (@eq-O1VS-epsd), with $phi = 0$.
Its gauge term acted only on a discretization error.
E20 adds the scalar potential at a face, which #link(<experiment-e03>)[E03] exercised only in a uniform medium.
The gauge term is now the only term of the vector update, because the curl term is zero.
The samples of $A_z$ at $z_k + h_z\/2 = z_f$ lie on the face and use the harmonic mean.
E20 does not test the curl term or a tangential component, and it contains no physical field.

=== Why the electric field stays zero
Let $G f := sum_d D_d^- (epsilon_d D_d^+ f)$ be the operator of the scalar update, with $epsilon_d$ from @eq-O1VS-epsd.
Assume that at step $n$
$
  A_d^n = D_d^+ Lambda^n, quad A_d^(n-1) = D_d^+ Lambda^(n-1), quad
  phi^(n-1/2) = -(Lambda^n - Lambda^(n-1))/h_t, quad
  phi^(n+1/2) = -(Lambda^(n+1) - Lambda^n)/h_t,
$ <E20-discrete-gauge-form>
where
$
  Lambda^(n+1) = 2 Lambda^n - Lambda^(n-1) + h_t^2/chi G Lambda^n.
$ <E20-discrete-wave>
The difference operators commute.
For $A_d = D_d^+ Lambda$ the Standard Laplacian of $A_d$ therefore equals $D_d^+ gamma$ (@eq-O1VS-gamma), and $g = G Lambda$ (@eq-O1VS-g).
With @eq-O1VS-ThAd the vector update is then
$
  A_d^(n+1) = 2 A_d^n - A_d^(n-1) + h_t^2 D_d^+ ((G Lambda^n)/chi) = D_d^+ Lambda^(n+1).
$
The scalar update is linear, and $phi^(n+1/2)$ is the difference of two successive $Lambda$.
With @E20-discrete-wave it gives $phi^(n+3/2) = -(Lambda^(n+2) - Lambda^(n+1))\/h_t$,
where $Lambda^(n+2)$ follows from @E20-discrete-wave one step later.
Thus @E20-discrete-gauge-form holds at every step if it holds at $n = 0$, and @E20-initialization ensures this.
The R2 reconstruction (@sec-time-placement) then gives
$
  E_d^(n-1/2) = -(A_d^n - A_d^(n-1))/h_t - D_d^+ phi^(n-1/2) = 0
$
at every step and on every face, also on the material face.
This is a deduction from the update equations of O1VS.
It uses that the same $epsilon_d$ stands in the scalar flux, in $g$ and in front of $partial_t^2 A_d$.

== Experiment E20 measured quantities and their purpose
The sets $sans(M)_1$ and $sans(M)_2$ are those of E16, and $sans(M)$ is their union.

=== E20.meas01: Electric field relative to its gradient term
The reference of $vb(E)$ is zero, so the relative error of E16.meas01 has no denominator.
The R2 reconstruction (@sec-time-placement) forms $E_z$ from two terms,
$
  E_z^(n-1/2) = -(A_z^n - A_z^(n-1))/h_t - D_z^+ phi^(n-1/2).
$
In E20 both terms are nonzero and they cancel.
We take the size of the second term as the denominator.
For $n$ completed steps, define
$
  N_E^(n-1/2) = sum_(bold(i) in sans(M),c) abs(E_(c,bold(i)+1/2 bold(e)_c)^("num",n-1/2))^2, quad
  S_E^(n-1/2) = sum_(bold(i) in sans(M)) abs((D_z^+ phi^("num",n-1/2))_(bold(i)+1/2 bold(e)_z))^2, quad
  e_E^(n-1/2) = sqrt(N_E^(n-1/2) / S_E^(n-1/2)).
$
Thus $e_E$ is the fraction of the gradient term that remains after the cancellation.
We also save $B_"max"^(n-1/2)$, the largest magnitude of all R2 magnetic samples.
It is zero because $A_z$ does not vary along $x$ and $y$.

*Purpose:* this is the main result of E20.
We expect $e_E$ at the level of the rounding error at every observation, also while the pulse crosses the face.
The rounding error of a difference quotient grows as $1\/h_t$, so $e_E$ can double from one level to the next.
An increase near $T\/2$ would show that the vector update and the scalar update do not agree at the face.

=== E20.meas02: Weighted sum of the scalar potential on each side
With $chi_bold(i) = epsilon_bold(i)^2$ define
$
  Q_S^(n+1/2) = sum_(bold(i) in sans(M)_S) chi_bold(i) abs(phi_bold(i)^("num",n+1/2))^2, quad
  rho_S^(n+1/2) = Q_S^(n+1/2) / (Q_1^(1/2) + Q_2^(1/2)).
$
The scalar potential of the current slot has the time $(n + 1\/2) h_t$.
We evaluate $Q_S$ at the steps selected for the R2 fields.
The denominator is the value at the initial observation.

*Purpose:* compare the split of the scalar potential with @E20-fractions.
After the separation the continuum values are $rho_1 = Gamma_phi^2$ and $rho_2 = 1 - Gamma_phi^2$.
We report the relative deviations $delta_1$ and $delta_2$ of @E16-energy-deviations, with $Gamma_phi$ in place of $Gamma$.
This measures the reflection and the transmission of the scalar update with the harmonic face permittivity.

=== E20.meas03: Relative error of the scalar potential on each side
This is E16.meas01 with the samples of $phi$ in place of those of $vb(E)$. For $S in {1,2}$ define
$
  N_(phi,S)^(n+1/2) = sum_(bold(i) in sans(M)_S) abs(phi_bold(i)^("num",n+1/2) - phi_bold(i)^("ref",n+1/2))^2, quad
  D_(phi,S)^(n+1/2) = sum_(bold(i) in sans(M)_S) abs(phi_bold(i)^("ref",n+1/2))^2, quad
  e_(phi,S)^(n+1/2) = sqrt(N_(phi,S)^(n+1/2) / D_(phi,S)^(n+1/2)).
$
The reference $phi_bold(i)^("ref",n+1/2)$ is @E20-analytic-potentials at the cell center and at the time $(n + 1\/2) h_t$.
We evaluate the sums at the three observations, at the steps of E20.meas02.
We do not form $e_(phi,2)$ while $D_(phi,2) = 0$.

*Purpose:* E20.meas02 uses squared magnitudes, which do not show the sign or the phase of a pulse.
$Gamma_phi$ is negative, and a reflected pulse with the wrong sign gives $e_(phi,1) approx 2$.
The scalar potential has no sample on the face, so the reference needs no rule there.

== Experiment E20 saved output
Each run writes the following files.

#table(
  columns: (1fr, 4fr),
  inset: 6pt,
  align: left + top,
  stroke: 0.5pt,
  [`measurements.csv`],
  [
    One row per quantity and observation, with the target time, the selected step index and the actual time.
    The quantities are `E_square_sum` and `gradient_phi_square_sum` ($N_E$ and $S_E$ of E20.meas01), `B_max`,
    `phi_weighted_square_sum_side1` and `phi_weighted_square_sum_side2` ($Q_1$ and $Q_2$ of E20.meas02),
    and the four sums `phi_error_side1`, `phi_reference_side1`, `phi_error_side2` and `phi_reference_side2`
    ($N_(phi,S)$ and $D_(phi,S)$ of E20.meas03).
    The first five are also written every $N_"pulse"\/16$ steps.
  ],
  [`profiles.csv`],
  [
    At the three observations, the numerical and the analytic profiles of $phi$ and $A_z$
    and the numerical profile of $E_z$ along one line parallel to $z$.
    The components are named `phi`, `Az` and `Ez`.
    On the face, where $A_z$ jumps, the analytic sample is the mean of the two limits.
  ],
  [`run.json`],
  [
    The configuration of the run.
    Under `effective.analytic_reference` it holds `reflection_coefficient` and `transmission_coefficient`
    ($Gamma_phi$ and $tau_phi$), and `reflected_fraction` and `transmitted_fraction` (@E20-fractions).
  ],
)

== Experiment E20 analysis
We report $e_E$ on each level at the three observations and its maximum over the history, and we plot the history.
We report $rho_1$ and $rho_2$ at $T$ with $delta_1$, $delta_2$ and their observed orders, as in E16.
We report $e_(phi,1)$ and $e_(phi,2)$ at $T$ with their observed orders.
We plot the profiles of $phi$ and $A_z$ at $T\/2$ near the face, and the profile of $phi$ at $T$.
