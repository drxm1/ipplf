#import "../../../docs/manual/preamble.typ": *

== Discrete potential energies <sec-discrete-potential-energies>

`potentialEnergy` in `FieldEnergy.h` computes the discrete wave energies $W_A$ and $W_phi$ of the O1-VS vector- and scalar-potential updates. This implementation requires the Standard stencil, staggered spatial placement and the R2 time convention on a uniform 3D mesh. The material coefficients must come from the solver's `MaterialOwner`. These energies are used to check conservation by the numerical update. The expressions below specialize the leapfrog energy of @kreiss2002SecondOrderWave[Sec. 2, conserved form before Eq. (2.10)] to Chew's potential equations @chew2014generalizedgauge[Eqs. (15)-(18)] (deduction). All variables in this subchapter are normalized solver variables.

=== What is summed

We use the spatial indices of @sec-spatial-placement. The integer vector $bold(i)$ denotes a cell centre. Scalar-potential samples and the coefficients $epsilon_bold(i)$ and $chi_bold(i)$ occupy that centre. The vector-potential component $A_(d,bold(i)+1/2 bold(e)_d)$ and its harmonic coefficient $epsilon_(d,bold(i)+1/2 bold(e)_d)$ occupy the positive $d$-face. A magnetic component in direction $d$ occupies the edge parallel to $bold(e)_d$ at $bold(i)+1/2 bold(e)_p+1/2 bold(e)_q$, where $(d,p,q)$ runs through $(x,y,z)$, $(y,z,x)$ and $(z,x,y)$. These are the staggered locations of @ryu2016potentialfdtd[Appendix A, Eqs. (19)-(24)]. The implementation stores all these samples at integer array index $bold(i)$, and their spatial subscripts here identify their distinct mesh locations.

After $n$ completed advances, define the current and previous vector-potential samples by
$
  a_(d,bold(i)+1/2 bold(e)_d) = A_(d,bold(i)+1/2 bold(e)_d)^n,
  quad b_(d,bold(i)+1/2 bold(e)_d) = A_(d,bold(i)+1/2 bold(e)_d)^(n-1).
$
They represent $t_n$ and $t_(n-1)$, with $t_n=n h_t$. In R2 the scalar entries in the same current and previous slots are $phi_bold(i)^(n+1/2)$ and $phi_bold(i)^(n-1/2)$, at $t_(n+1/2)$ and $t_(n-1/2)$ (see @sec-time-placement). In the sums below we abbreviate them by $phi_bold(i)$ and $psi_bold(i)$, respectively. Thus sharing a storage slot does not imply sharing a sample time. This is the convention assigned by `LevelContract`; the relative half-step shift agrees with the potential updates in @ryu2016potentialfdtd[Appendix A, Eqs. (26)-(27)], after shifting the time origin to our vector-potential levels. The cell volume is $Delta V=h_x h_y h_z$. The harmonic face coefficients are defined in @eq-O1VS-epsd, and $chi_bold(i)=alpha epsilon_bold(i)^2$, with generalized-gauge parameter $alpha$ @chew2014generalizedgauge[Eq. (18)].

The weighted divergence is located at the cell centre. Expanding @eq-O1VS-g with the staggered differences of @ryu2016potentialfdtd[Appendix A, Eq. (23)] gives
$
  g_(a,bold(i)) = sum_d
  (epsilon_(d,bold(i)+1/2 bold(e)_d) a_(d,bold(i)+1/2 bold(e)_d)
  - epsilon_(d,bold(i)-1/2 bold(e)_d) a_(d,bold(i)-1/2 bold(e)_d))/h_d.
$
The definition of $g_(b,bold(i))$ replaces $a$ by $b$. For each cyclic triple $(d,p,q)$, the forward curl gives the edge sample @ryu2016potentialfdtd[Appendix A, Eq. (24)]
$
  B_(a,d,bold(i)+1/2 bold(e)_p+1/2 bold(e)_q)
  &= (a_(q,bold(i)+bold(e)_p+1/2 bold(e)_q)-a_(q,bold(i)+1/2 bold(e)_q))/h_p #<equate:revoke> \
  &quad - (a_(p,bold(i)+bold(e)_q+1/2 bold(e)_p)-a_(p,bold(i)+1/2 bold(e)_p))/h_q.
$
Again, replacing $a$ by $b$ defines $B_(b,d,bold(i)+1/2 bold(e)_p+1/2 bold(e)_q)$. These are the raw curls of the two vector-potential states, evaluated separately. They are not the magnetic field averaged in time for R2 field output. Eq. (24) of the cited manuscript prints $Delta z$ in the second term of the $y$-component. That difference is along $x$, so the formula above uses $h_x$, consistent with its Eq. (19) and `staggeredCurl`.

Let $cal(I)_"owned"$ contain the owned cell indices across all MPI ranks, excluding halo copies. Each rank sums its own indices, and `potentialEnergy` combines the rank-local sums. Applying `energyDetail::Reduction` therefore gives the global energy
$
  W_A &= (Delta V)/2 sum_(bold(i) in cal(I)_"owned") sum_d
  epsilon_(d,bold(i)+1/2 bold(e)_d)
  ((a_(d,bold(i)+1/2 bold(e)_d)-b_(d,bold(i)+1/2 bold(e)_d))/h_t)^2 #<equate:revoke> \
  &quad + (Delta V)/2 sum_(bold(i) in cal(I)_"owned") sum_d
  B_(a,d,bold(i)+1/2 bold(e)_p+1/2 bold(e)_q)
  B_(b,d,bold(i)+1/2 bold(e)_p+1/2 bold(e)_q) #<equate:revoke> \
  &quad + (Delta V)/2 sum_(bold(i) in cal(I)_"owned")
  (g_(a,bold(i)) g_(b,bold(i)))/chi_bold(i).
$ <eq-diag-vector-energy>
The first line sums the squared time differences on faces, weighted by the face permittivity. The second sums products of raw magnetic samples on matching edges. The third sums products of weighted divergences at cell centres. Each owned index contributes its associated positive faces and transverse edges once. No interpolation to a common spatial point is performed; these are separate sums over staggered samples.

For the scalar-potential energy, the forward gradient occupies a face @ryu2016potentialfdtd[Appendix A, Eq. (19)]:
$
  (D_d^+ phi)_(bold(i)+1/2 bold(e)_d)
  = (phi_(bold(i)+bold(e)_d)-phi_bold(i))/h_d,
$
and similarly for $psi$. The scalar counterpart of the same leapfrog energy is
$
  W_phi & = (Delta V)/2 sum_(bold(i) in cal(I)_"owned")
          chi_bold(i) ((phi_bold(i)-psi_bold(i))/h_t)^2 #<equate:revoke> \
        & quad + (Delta V)/2 sum_(bold(i) in cal(I)_"owned") sum_d
          epsilon_(d,bold(i)+1/2 bold(e)_d)
          (D_d^+ phi)_(bold(i)+1/2 bold(e)_d)
          (D_d^+ psi)_(bold(i)+1/2 bold(e)_d).
$ <eq-diag-scalar-energy>
Its time-difference term is summed at cell centres, and its gradient products on matching faces. `FieldEnergy.h` labels $W_A$ at $t_(n-1/2)$ and $W_phi$ at $t_n$, the respective midpoints of the paired potential states in R2.

=== What the values establish

For a source-free periodic run with fixed material coefficients, the discrete update conserves $W_A$ and $W_phi$ in exact arithmetic. This follows from the cited leapfrog identity applied to the two potential equations. A change in either value therefore provides a check of conservation by the numerical update. For open boundaries, the implementation computes the same sums but includes no boundary-flux term; these values need not be constant. This conservation check alone does not establish the accuracy of the electromagnetic fields.

The definition of $W_A$ differs from the electromagnetic field energy $U$ in @units-scale-relation-energy. Its first term contains the time difference of $vb(A)$ without the scalar-potential gradient, its magnetic term multiplies fields from adjacent times, and it includes the weighted-divergence term. It must therefore be distinguished from the field-energy estimate written as `field_energy`.

=== Material contributions and output columns

`vacuum` and `dielectric` partition the terms of $W_A$ by material, and `total` is their sum. Cell-centre terms use the cell's permittivity; face terms use the face coefficient. For the existing vacuum/dielectric presets, an edge term is assigned to the dielectric contribution when any adjoining cell is dielectric. Thus the split follows the location of each term, not just the material of the cell whose array index stores it. The `scalar` column contains $W_phi$ over the whole domain. The `absolute` column sums the magnitudes of the elementary terms in $W_A$, with the same volume factor, before cancellation between positive and negative terms. These are the quantities returned by `potentialEnergy` and written to `energy.csv` by `FieldRunOutput`.
