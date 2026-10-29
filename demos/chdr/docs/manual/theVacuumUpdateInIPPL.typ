#import "./preamble.typ": *

#let chapter = [
  #show: equate.with(sub-numbering: false)



  = The vacuum update in IPPL's existing FDTD Solvers
  Let $epsilon eq.triple epsilon_r (vb(x))$. In this section we describe how the StandardFDTDSolver and the NonStandardFDTDSolver that were already in IPPL when this master's thesis began update the four-potential.

  For updating $(phi, vb(A))$, the relevant MW Eqs. are Ampère and Gauss (cf. @eq-MW2-ampere-in-our-mat-model and @eq-MW4-gauss-in-our-mat-model):
  $
             curl vb(B) & = epsilon partial_t vb(E) + vb(J) \
    div (epsilon vb(E)) & = rho.
  $
  We can rewrite these equations as equations for $phi$ and $vb(A)$:
  $
    "Ampere:" quad quad curl curl vb(A) &= -epsilon partial_t^2 vb(A) - epsilon partial_t grad phi + vb(J) #<eq-MW2-ampere-inOurMatModel-writtenWithPotentials>\
    "Gauss:" quad quad - partial_t div (epsilon vb(A)) &= div (epsilon grad phi) + rho. #<eq-MW4-gauss-inOurMatModel-writtenWithPotentials>
  $
  In a homogeneous lorentz gauge ($div vb(A) = -epsilon partial_t phi$), when $epsilon(vb(x)) = "const"$ they decouple:
  // $
  //   dalembert_epsilon bold(phic) = - rho\/epsilon quad "and" quad dalembert_epsilon Ac = -vb(J) quad quad "where" dalembert_epsilon := (nabla^2 - epsilon partial_t^2)
  // $ -->
  $
    partial_t^2 vb(A) = 1/epsilon nabla^2 vb(A) + vb(J)/epsilon quad quad quad
    partial_t^2 phi = 1/epsilon nabla^2 phi + rho/epsilon^2
  $

  IPPL's FDTD Solvers implement the $epsilon=1$ case componentwise:

  $
    cal(A)^(n+1) = 2 cal(A)^n - cal(A)^(n-1) + h_t^2 (L_vb(h)[cal(A)^n] + cal(S)^n).
  $
  where $n$ is a step label, $cal(A) := (phi, vb(A))$, $cal(S) := (rho, vb(J))$ and
  - $(L_vb(h)^"Standard"u) = (D_"xx"+D_"yy"+D_"zz") u, quad quad h_x^2 (D_"xx"u)_(i j k) = u_(i+1,j k) - 2 u_(i j k) + u_(i-1,j k)$
  - $(L_vb(h)^"NS"u)_(i j k) =
    (D_"zz" u)_(i j k) +
    ((1 - 2 cal(a))Q u)_(i j k) +
    (cal(a)Q u)_(i j, k-1) +
    (cal(a)Q u)_(i j, k+1), quad Q = D_"xx" + D_"yy"$


]

#standalone(title: "The vacuum update in IPPL's existing FDTD Solvers", chapter)
