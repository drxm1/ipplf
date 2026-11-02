
#import "../../../../docs/manual/preamble.typ": *

== Chew's Gauge

=== Why a regular lorenz gauge is not enough
The problem with a regular lorenz gauge is the following:
When $epsilon(vb(x)) eq.not "const"$, in the gauge $div vb(A) = - epsilon partial_t phi$, Ampere's and Gauss's laws (cf. TODO) stay coupled:
$
  "scalar:"& quad quad quad &epsilon^2 partial_t^2 phi - partial_t vb(A) dot grad epsilon &= epsilon nabla^2 phi + grad epsilon dot grad phi + rho #<eq-local-lorenz-scalar>\
  "vector:"& quad quad quad &epsilon partial_t^2 vb(A)
  &= nabla^2 vb(A) + vb(J) + (partial_t phi) grad epsilon #<eq-local-lorenz-vector>
$

==== Derivation of the remaining coupling in the $div vb(A) = - epsilon partial_t phi$ gauge.

We can see where the coupling comes from by inserting the local Lorenz
condition into the potential equations before making any homogeneous-medium
simplification. In the following calculation, $epsilon = epsilon_r (vb(x))$
is stationary and sufficiently smooth, and $mu = 1$. The calculation is a
deduction from @chew2014generalizedgauge[Eqs. (7)-(11)]; we extend the
homogeneous condition in Eq. (11) to a spatially varying $epsilon$, as
stated above.

*Scalar equation:* We start from the potential form of Gauss's Law in
@eq-MW4-gauss-inOurMatModel-writtenWithPotentials.
We first expand the weighted divergence and then use the local Lorenz
condition stated above:
$
  div (epsilon vb(A)) & = epsilon div vb(A) + vb(A) dot grad epsilon \
                      & = -epsilon^2 partial_t phi + vb(A) dot grad epsilon.
$
Since $partial_t epsilon = 0$, inserting this into
@eq-MW4-gauss-inOurMatModel-writtenWithPotentials
gives
$
  epsilon^2 partial_t^2 phi
  = div (epsilon grad phi) + rho + (partial_t vb(A)) dot grad epsilon.
$ <draft-local-lorenz-scalar>
Expanding $div (epsilon grad phi) = epsilon nabla^2 phi + grad epsilon dot grad phi$
recovers @eq-local-lorenz-scalar. The term $(partial_t vb(A)) dot grad epsilon$
still couples this equation to $vb(A)$.

*Vector equation:* In Ampère's Law
@eq-MW2-ampere-inOurMatModel-writtenWithPotentials
we use $curl curl vb(A) = grad (div vb(A)) - nabla^2 vb(A)$, as in
@chew2014generalizedgauge[Eq. (10)], and then the local Lorenz condition:
$
  epsilon partial_t^2 vb(A)
  &= nabla^2 vb(A) + vb(J) - grad (div vb(A)) - epsilon grad (partial_t phi) \
  &= nabla^2 vb(A) + vb(J) + grad (epsilon partial_t phi) - epsilon grad (partial_t phi).
$ <draft-local-lorenz-vector>
Expanding the last gradient by the product rule gives @eq-local-lorenz-vector.
It leaves a term containing $phi$, just as the scalar equation
retains a term containing $vb(A)$. Both additional terms vanish when
$grad epsilon = 0$. For varying $epsilon$, the local Lorenz condition therefore
leaves a coupled system.

Note that the product expansions above assume a smooth $epsilon$. At a sharp
dielectric interface they are not ordinary pointwise equations. The interface
conditions for the generalized gauge used below are derived from the
unexpanded equations in @chew2014generalizedgauge[Sec. 3, Eqs. (23)-(31)].

=== Chew's generalized lorenz gauge decouples them

We can use Chew's Generalized Lorentz gauge to decouple the scalar and vector eqs, keeping the interface conditions intact. Chew's Gauge with $alpha=1$ is the following:
$
  quad div (epsilon vb(A)) = -chi partial_t phi, quad quad "with" quad chi = epsilon^2.
$
\
Inserting that gauge into Ampere & Gauss (cf. TODO), we obtain
@chew2014generalizedgauge[Eqs. (16)-(17)]:
$
  "vector:" quad epsilon partial_t^2 vb(A) &= - curl curl vb(A) + epsilon grad (epsilon^(-2) div (epsilon vb(A))) + vb(J) #<eq-Chew-vector>\
  "scalar:" quad epsilon^2 partial_t^2 phi &= div (epsilon grad phi) + rho #<eq-Chew-scalar>
$

=== Rewriting the vector and scalar equations in Chew's gauge to more convenient forms
We can use $curl curl vb(A) = grad (div vb(A)) - nabla^2 vb(A)$ to rewrite the vector eq to look like the homogeneous update + $vb(T(A))$:
$
  partial_t^2 vb(A) = 1/epsilon nabla^2 vb(A) + 1/epsilon vb(J) + 1/epsilon &vb(T(A))\
  &vb(T(A)) := epsilon grad (chi^(-1) div (epsilon vb(A))) - grad (div vb(A)) #<eq-Chew-material-correction>
$
$vb(T(A))$ is zero where $epsilon="const"$.

Similarly, we can rewrite the scalar equation (TODO:IncludeDerivation) to the following form:
$ chi partial_t^2 phi = nabla^2 phi + div ((epsilon - 1) grad phi) + rho $
