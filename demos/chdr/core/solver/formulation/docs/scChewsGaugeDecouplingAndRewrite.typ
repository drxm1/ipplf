
#import "../../../../docs/manual/preamble.typ": *

== Chew's Gauge

=== Why a regular lorenz gauge is not enough
The problem with a regular lorenz gauge is the following:
When $epsilon(vb(x)) eq.not "const"$, in the gauge $div vb(A) = - epsilon partial_t phi$, Ampere's and Gauss's laws (cf. TODO) stay coupled:
$
  "scalar:" quad quad epsilon^2 partial_t^2 phi - partial_t vb(A) dot grad epsilon &= epsilon nabla^2 phi + grad epsilon dot grad phi + rho\
  "vector:" "TODO" &= "TODO"
$

=== Chew's generalized lorenz gauge decouples them

We can use *Chew's Generalized Lorentz gauge* ($alpha=1$) to decouple the scalar and vector eqs, keeping the interface conditions intact: $quad div (epsilon vb(A)) = -chi partial_t phi, quad chi = epsilon^2$.
\
Inserting that gauge into Ampere & Gauss (cf. TODO), we obtain:
$
  "vector:" quad epsilon partial_t^2 vb(A) &= - curl curl vb(A) + epsilon grad (epsilon^(-2) div (epsilon vb(A))) + vb(J)\
  "scalar:" quad epsilon^2 partial_t^2 phi &= div (epsilon grad phi) + rho
$

=== Rewriting the vector and scalar equations in Chew's gauge to more convenient forms
We can use $curl curl vb(A) = grad (div vb(A)) - nabla^2 vb(A)$ to rewrite the vector eq to look like the homogeneous update + $vb(T(A))$:
$
  partial_t^2 vb(A) = 1/epsilon nabla^2 vb(A) + 1/epsilon vb(J) + 1/epsilon &vb(T(A))\
  &vb(T(A)) := epsilon grad (chi^(-1) div (epsilon vb(A))) - grad (div vb(A))
$
$vb(T(A))$ is zero where $epsilon="const"$.

Similarly, we can rewrite the scalar equation (TODO:IncludeDerivation) to the following form:
$ chi partial_t^2 phi = nabla^2 phi + div ((epsilon - 1) grad phi) + rho $
