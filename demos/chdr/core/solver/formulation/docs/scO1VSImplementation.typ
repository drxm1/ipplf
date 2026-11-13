

#import "../../../../docs/manual/preamble.typ": *

== The O1-VS implementation <sec-O1VS-implementation>
This section describes the implemented option O1 (Chew's gauge) on the virtually staggered (staggered spatial convention, see @sec-spatial-placement) spatial convention and the R2 time-convention (see @sec-time-placement).

The following calculation concerns the interior stencil on a uniform Cartesian mesh.
`O1VSUpdate::advance` excludes the outermost layers owned by the Mur conditions because their boundary updates are separately handled.

=== Spatial locations and difference operators <sec-O1VS-operators>

We use the spatial indices defined in @sec-spatial-placement.
Here we specify the locations and differences needed for the O1-VS equations.
Let $bold(i)$ index a cell centre at $vb(x)_bold(i)$, and let $h_d$ be
the mesh spacing in direction $d$.
The quantities $phi_bold(i)$, $rho_bold(i)$, $epsilon_bold(i)$ and
$chi_bold(i) = epsilon_bold(i)^2$ are assigned to that centre.
The components $A_(d,bold(i)+1/2 bold(e)_d)$,
$J_(d,bold(i)+1/2 bold(e)_d)$ and the coefficient
$epsilon_(d,bold(i)+1/2 bold(e)_d)$ are assigned to the positive $d$-face,
at $vb(x)_bold(i) + (h_d/2) vb(e)_d$.
These are the locations used by `LevelContract::positionOffset` and
`O1VSOperator`. Their arrays use the integer index $bold(i)$ for this face,
as described in @sec-spatial-placement.

We use the forward and backward differences of Ryu et al.
@ryu2016potentialfdtd[Appendix A, Eqs. (19)-(23), manuscript pp. 5-6].
For a centre scalar $f$ and a face component $v_d$, these are
$
  (D_d^+ f)_(bold(i)+1/2 bold(e)_d)
  &= (f_(bold(i)+bold(e)_d) - f_bold(i))/h_d, \
  (D_d^- v_d)_bold(i)
  &= (v_(d,bold(i)+1/2 bold(e)_d) - v_(d,bold(i)-1/2 bold(e)_d))/h_d.
$
The first expression is located at the positive face, while the second is at the
centre between the two faces. Composing them gives, by direct substitution,
$
  (D_d^- D_d^+ f)_bold(i)
  = (f_(bold(i)+bold(e)_d) - 2 f_bold(i) + f_(bold(i)-bold(e)_d))/h_d^2.
$
Summing over the three directions therefore gives the Standard Laplacian:
$L_vb(h)^"Standard" = sum_d D_d^- D_d^+$.
This is also the composition documented in `staggeredDifference`.
The same second-difference stencil acts on each shifted face lattice for
$A_d$, using neighbouring values of that component.

=== Permittivity at material interfaces

The update below needs a permittivity $epsilon_d$ at each positive
$d$-face, whereas the material values are stored at cell centres. We first
derive this coefficient for a stationary, planar interface normal to the
$d$-direction, exactly halfway between two neighbouring centres. Write $x$
for the coordinate in that direction, $x_"L"$ and $x_"R"$ for the centres,
and $x_"f"$ for the interface. Thus $x_"R" - x_"L" = h_d$ and each
centre is $h_d/2$ from the interface. The two half-cells have constant,
positive permittivities $epsilon_"L"$ and $epsilon_"R"$.

We use Chew's generalized gauge and
assume that no additional free surface charge, free surface current or
dipole layer is imposed at the interface. Bound polarization charge is
already represented by the permittivity. The conditions needed for the
scalar potential are @chew2014generalizedgauge[Sec. 3, Eqs. (27), (29)]
$
                        phi(x_"f"^-) & = phi(x_"f"^+), \
  epsilon_"L" partial_x phi(x_"f"^-) & = epsilon_"R" partial_x phi(x_"f"^+).
$
Here $-$ and $+$ denote the limits from the left and right. Thus the
potential and the flux $epsilon partial_x phi$ are continuous, while the
normal derivative itself can jump.

Let $phi_"L"$ and $phi_"R"$ be
the stored centre values, and introduce one auxiliary interface value
$phi_"f"$. We approximate $phi$ by a linear function on each half-cell.
Following Eymard, Gallouët and Herbin
@eymard2000FiniteVolume[Sec. 7.1, Eqs. (7.3)-(7.5), Example 7.1,
  manuscript p. 22],
we require the two approximations of the interface flux to agree:
$
  F = epsilon_"L" (phi_"f" - phi_"L")/(h_d/2)
  = epsilon_"R" (phi_"R" - phi_"f")/(h_d/2).
$
Our $F$ approximates $+epsilon partial_x phi$, while Eymard et al. use the
opposite sign for the diffusive flux. Solving for the auxiliary value gives
$
  (epsilon_"L" + epsilon_"R") phi_"f"
  = epsilon_"L" phi_"L" + epsilon_"R" phi_"R",
  quad quad
  phi_"f" = (epsilon_"L" phi_"L" + epsilon_"R" phi_"R")
  /(epsilon_"L" + epsilon_"R").
$
Substituting this into either expression for $F$, we obtain
$
  F = (2 epsilon_"L" epsilon_"R")/(epsilon_"L" + epsilon_"R")
  (phi_"R" - phi_"L")/h_d.
$
We therefore recover the difference between the two centre values, now multiplied by the *harmonic mean*:
$
  epsilon_d = cal(H)(epsilon_"L", epsilon_"R")
  := (2 epsilon_"L" epsilon_"R")/(epsilon_"L" + epsilon_"R")
  = 2/(epsilon_"L"^(-1) + epsilon_"R"^(-1)).
$ <eq-O1VS-harmonic-mean>
With $epsilon_"L" = epsilon_bold(i)$ and
$epsilon_"R" = epsilon_(bold(i)+bold(e)_d)$, this is the coefficient computed by `HarmonicFace<d>`:
$
  epsilon_(d,bold(i)+1/2 bold(e)_d)
  = cal(H)(epsilon_bold(i), epsilon_(bold(i) + bold(e)_d)).
$ <eq-O1VS-epsd>
The scalar flux is then $F = epsilon_d D_d^+ phi$.
If the two cells contain the same material, $cal(H)(epsilon, epsilon) = epsilon$.

The vector update also needs a face coefficient. It forms the product
$epsilon_d A_d$ inside the weighted divergence $g$ in @eq-O1VS-g.
This divergence approximates $div (epsilon vb(A))$ in Chew's gauge @chew2014generalizedgauge[Eq. (15)].
Chew's interface condition requires continuity of the normal product $epsilon A_d$ for the continuous fields
@chew2014generalizedgauge[Sec. 3, Eq. (26)].

O1-VS uses the same $epsilon_d$ in this product, in the scalar flux,
and as the coefficient multiplying the vector potential's second time derivative.
The calculation above establishes the harmonic coefficient for the scalar flux only.

*Open point:* the shared use of $epsilon_d$ in the vector update is not yet justified.
This requires showing that the two discrete potential equations remain compatible with the discrete gauge and with charge continuity.
This derivation has not been done.

=== Spatial discretization

We now construct the spatial terms used in the `O1-VS(R2)` update from
@eq-Chew-vector and @eq-Chew-scalar.
We take $mu = alpha = 1$ in Chew's choice of the gauge coefficient
@chew2014generalizedgauge[Eq. (18)].
Here $partial_d$ denotes the derivative with respect to $x_d$.

==== Scalar equation

The continuous flux component $epsilon partial_d phi$
is represented by $epsilon_d D_d^+ phi$, at the positive $d$-face.
Taking its backward divergence gives the centre quantity
$sum_d D_d^- (epsilon_d D_d^+ phi)$, which represents
$div (epsilon grad phi)$ in @eq-Chew-scalar. Thus
$
  chi_bold(i) partial_t^2 phi_bold(i)
  = sum_d (D_d^- (epsilon_d D_d^+ phi))_bold(i) + rho_bold(i).
$
To relate this to the existing Standard update, write $epsilon_d = 1 +
(epsilon_d - 1)$ inside the flux. Linearity and the composition above give
$
  sum_d (D_d^- (epsilon_d D_d^+ phi))_bold(i)
  & = (L_vb(h)^"Standard" phi)_bold(i)
    + sum_d (D_d^- ((epsilon_d - 1) D_d^+ phi))_bold(i).
$ <eq-O1VS-scalar-flux-split>
The last sum is the correction evaluated by
`O1VSOperator::scalarCorrection`.

==== Vector equation

We now discretize the continuous correction defined in @eq-Chew-material-correction,
obtained from @chew2014generalizedgauge[Eq. (17)].
We define the two discrete divergences
$
  gamma_bold(i) & = sum_q (D_q^- A_q)_bold(i) #<eq-O1VS-gamma> \
      g_bold(i) & = sum_q (D_q^- (epsilon_q A_q))_bold(i). #<eq-O1VS-g>
$

#[
  #set par(justify: false)
  #table(
    columns: (30%, 30%, 40%),
    align: (left + horizon, left + horizon, left + horizon),
    inset: (x: 4pt, y: 7pt),
    stroke: none,
    table.hline(stroke: 0.5pt),
    table.header(
      [*Continuous expression*],
      [*Discrete counterpart*],
      [*Discrete out location*],
    ),
    table.hline(stroke: 0.4pt),
    [$div vb(A)$],
    [$gamma$ (@eq-O1VS-gamma)],
    [Cell centre],
    [$div (epsilon vb(A))$],
    [$g$ (@eq-O1VS-g)],
    [Cell centre],
    [$chi^(-1) div (epsilon vb(A))$],
    [$g/chi$],
    [Cell centre],
    [$partial_d (div vb(A))$],
    [$D_d^+ gamma$],
    [Positive $d$-face],
    [$epsilon partial_d (chi^(-1) div (epsilon vb(A)))$],
    [$epsilon_d D_d^+ (g/chi)$],
    [Positive $d$-face],
    table.hline(stroke: 0.5pt),
  )
]

Both $gamma$ and $g/chi$ are centre scalars, so their forward differences return to the positive faces.
Subtracting the fourth row of the table from the fifth gives the discrete correction
$
  (T_h A)_(d,bold(i)+1/2 bold(e)_d)
  = epsilon_(d,bold(i)+1/2 bold(e)_d) (D_d^+ (g / chi))_(bold(i)+1/2 bold(e)_d)
  - (D_d^+ gamma)_(bold(i)+1/2 bold(e)_d).
$ <eq-O1VS-ThAd>
Note that $chi$ is used at the centres
inside the difference, whereas the factor $epsilon_d$ is at the output face.

With explicit spatial indices, the weighted divergence in @eq-O1VS-g is
$
  g_bold(i) = sum_q
  (epsilon_(q,bold(i)+1/2 bold(e)_q) A_(q,bold(i)+1/2 bold(e)_q)
  - epsilon_(q,bold(i)-1/2 bold(e)_q) A_(q,bold(i)-1/2 bold(e)_q))/h_q.
$
The two products belong to the opposite $q$-faces surrounding the centre
$bold(i)$. The result is assigned to that centre.

Thus the vector equation, with space discretized but time still continuous,
is
$
  epsilon_(d,bold(i)+1/2 bold(e)_d) partial_t^2 A_(d,bold(i)+1/2 bold(e)_d)
  = (L_vb(h)^"Standard" A_d)_(bold(i)+1/2 bold(e)_d)
  + (T_h A)_(d,bold(i)+1/2 bold(e)_d)
  + J_(d,bold(i)+1/2 bold(e)_d).
$
`O1VSOperator::curlCorrection` returns the negative of our material correction.
In the code's notation, $C_h^"code" = -T_h$. `O1VSUpdate::advance` subtracts that result.
The manual uses $T$ and $T_h$ with the same sign in the continuous and discrete equations.
`O1VSOperator::inverseTimeDerivativeCoefficient` returns $1/epsilon_d$ for the vector update and $1/chi$ for the scalar update.

=== Time discretization and implemented update
The following update equations are used when `MaterialFormulation::O1VS` is selected together with `Stencil::Standard`, `SpatialDiscretization::Staggered` and `EBReconstructionConvention::R2`.
`validateMaterialUpdate` rejects every other combination with `O1VS`.
They assume a finite $epsilon >= 1$ in every cell and $alpha = 1$, so that $chi = epsilon^2$.
The code accepts $alpha > 1$ only when every direction is periodic, and then uses $chi = alpha epsilon^2$.
The superscripts are the actual time levels defined in @sec-time-placement, and the subscripts are the mesh locations defined in @sec-spatial-placement.
The solve that advances the vector potential from $t_n$ to $t_(n+1)$ advances the scalar potential from $t_(n+1/2)$ to $t_(n+3/2)$.
The vector update is evaluated at the face $bold(i) + 1/2 bold(e)_d$, and the scalar update at the centre $bold(i)$.
The second equality in the scalar update uses @eq-O1VS-scalar-flux-split.
$
  A_(d,bold(i)+1/2 bold(e)_d)^(n+1)
  &= 2 A_(d,bold(i)+1/2 bold(e)_d)^n - A_(d,bold(i)+1/2 bold(e)_d)^(n-1) #<equate:revoke> \
  &quad + h_t^2/epsilon_(d,bold(i)+1/2 bold(e)_d)
  [(L_vb(h)^"Standard" A_d^n)_(bold(i)+1/2 bold(e)_d)
    + (T_h A^n)_(d,bold(i)+1/2 bold(e)_d)
    + J_(d,bold(i)+1/2 bold(e)_d)^n], \
  phi_bold(i)^(n+3/2)
  &= 2 phi_bold(i)^(n+1/2) - phi_bold(i)^(n-1/2)
  + h_t^2/chi_bold(i) [ sum_d (D_d^- (epsilon_d D_d^+ phi^(n+1/2)))_bold(i) + rho_bold(i)^(n+1/2) ] \
  &= 2 phi_bold(i)^(n+1/2) - phi_bold(i)^(n-1/2) #<equate:revoke> \
  &quad + h_t^2/chi_bold(i) [ (L_vb(h)^"Standard" phi^(n+1/2))_bold(i)
    + sum_d (D_d^- ((epsilon_d - 1) D_d^+ phi^(n+1/2)))_bold(i)
    + rho_bold(i)^(n+1/2) ]. #<equate:revoke>
$
