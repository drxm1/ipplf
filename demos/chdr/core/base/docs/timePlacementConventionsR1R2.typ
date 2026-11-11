#import "../../../docs/manual/preamble.typ": *

== Time-placement interpretation of the stored fields <sec-time-placement>

The electric field $vb(E)$ is usually reconstructed from $(phi, vb(A))$ via $vb(E) = - partial_t vb(A) - grad phi$.
Similarly, the magnetic field $vb(B)$ via $curl vb(A)$.

However we need to be careful about at which times the stored scalar and vector potentials are interpreted to be, and how exactly we compute $partial_t vb(A)$. This becomes important once $phi eq.not 0$, because in that case we need to ensure that both parts of the computed field $E$, $partial_t vb(A)$ and $grad phi$ are expressions at the same time.

In the whole manual we write the time of a discrete quantity explicitly.
The integer $n$ counts the completed solves, and $t_m = m h_t$ for an integer or half-integer $m$.
A superscript on a discrete quantity is its actual time level $m$.
For example $phi^(n+1/2)$ is the scalar potential at $t_(n+1/2)$, and $vb(E)^(n-1/2)$ is the electric field at $t_(n-1/2)$.
A superscript never names a storage slot.
We name the two stored potential histories by their code names instead: `A_n` is the current slot and `A_nm1` is the previous slot.
Which time levels these two slots hold after $n$ solves depends on the time convention, and the table below lists them.

Consider the current such implementation `evaluate_EB()`. The discrete gradient operator in $grad phi^n$ acts only spatially and therefore, if the scalar potential in the current slot is $phi^n$, at the time $t_n$ (which is my interpretation of how the time is assigned to that stored field in the existing FDTD solver's), the output $grad phi^n$ also is at $t_n$. In `evaluate_EB()`, $partial_t vb(A)$ is calculated via the backward difference $(A^n - A^(n-1))/h_t$ (see `FDTDSolverBase::evaluate_EB`), so its output theoretically lives at time $t_(n-1/2)$.
Since this is time is different from $grad phi^n$ at $t_n$ we have a problem as soon as $phi eq.not 0$. This makes the order of $E$ in the time step $h_t$ only 1 and not 2.
To see this, we expand $vb(A)$ in a Taylor series around $t_(n-1/2)$ and around $t_n$:
$
  (vb(A)^n - vb(A)^(n-1))/h_t
  &= partial_t vb(A)(t_(n-1/2)) + h_t^2/24 partial_t^3 vb(A)(t_(n-1/2)) + cal(O)(h_t^4) \
  &= partial_t vb(A)(t_n) - h_t/2 partial_t^2 vb(A)(t_n) + cal(O)(h_t^2).
$
The difference is of second order only at $t_(n-1/2)$.
If we read the computed $vb(E)$ at $t_n$, where $grad phi^n$ lives, the term $partial_t vb(A)$ has the error $-(h_t\/2) partial_t^2 vb(A)$.
If we read it at $t_(n-1/2)$ instead, the same expansion for $phi$ gives $grad phi^n = grad phi(t_(n-1/2)) + (h_t\/2) partial_t grad phi + cal(O)(h_t^2)$.
So with either reading, an error of first order in $h_t$ remains whenever $grad phi$ changes in time.

Therefore it can make sense to use MITHRA's half-step convention, where the scalar potential stored after $n$ solves is $phi^(n+1/2)$, the value at time $t_(n+1/2)$, while the stored vector potential is still $vb(A)^n$, the value at time $t_n$. Then both output terms in the expression for calculating $vb(E)$ align in their time-location, at $t_(n-1/2)$ if one takes the gradient of the previous scalar slot, $grad phi^(n-1/2)$.

M. Winkler's MSc thesis confirms that after $n$ solves the stored scalar and vector potentials are at the same time $t_n$, so that the current slot holds $phi^n$ and $vb(A)^n$.
There a superscript $n$ denotes the time $n h_t$, and the four-potential is one quantity that carries this superscript @winkler2024ChargeConservingFDTD[p. 16].
Its update uses the same two levels $n$ and $n-1$ for all four components, with the source term $(rho^n, vb(J)^n)$ at $n h_t$ @winkler2024ChargeConservingFDTD[Eqs. (3.2), (3.4), p. 17].
The reconstruction of $vb(E)$ and $vb(B)$ in the table is read from the code of `FDTDSolverBase::evaluate_EB`.
The following table therefore lists *R1* as the time-convention (in demos/chdr I sometimes also write `ebReconstructionConvention` instead of time-convention) that the currently-in-ippl FDTD solver assumes, and *R2* is MITHRA's convention where $phi$ and $vb(A)$ are shifted by half a time-step. The different implementation options described in @sec-material-treatment-options have different requirements on which time-level convention is used (see @sec-O1VS-implementation).

#[
  #set text(size: 9pt)
  #table(
    columns: (auto, 1fr, 1fr),
    align: (left + horizon, center + horizon, center + horizon),
    inset: (x: 5pt, y: 4.5pt),
    stroke: none,
    table.hline(stroke: 0.6pt),
    [*after $n$ solves*],
    [*R1* (in IPPL's FDTD Solvers currently)],
    [*R2* (half-step, MITHRA)],
    table.hline(stroke: 0.4pt),
    [current slot `A_n`], [$phi^n$, $vb(A)^n$], [$phi^(n+1\/2)$, $vb(A)^n$],
    [previous slot `A_nm1`], [$phi^(n-1)$, $vb(A)^(n-1)$], [$phi^(n-1\/2)$, $vb(A)^(n-1)$],
    [source for the next solve], [$rho^n$, $vb(J)^n$], [$rho^(n+1\/2)$, $vb(J)^n$],
    table.hline(stroke: 0.4pt + luma(85%)),
    [$vb(E)$],
    [$-(vb(A)^n - vb(A)^(n-1))\/h_t - grad phi^n$\ #text(
        fill: red,
      )[at $t_(n-1\/2)$ and $t_n$]],
    [$vb(E)^(n-1\/2) = -(vb(A)^n - vb(A)^(n-1))\/h_t - grad phi^(n-1\/2)$],
    [$vb(B)$],
    [$vb(B)^n = curl vb(A)^n$],
    [$vb(B)^(n-1\/2) = 1/2 curl (vb(A)^n + vb(A)^(n-1))$],
    table.hline(stroke: 0.4pt + luma(85%)),
    [order of $vb(E)$ ($phi eq.not 0$)],
    [1],
    [2],
    table.hline(stroke: 0.6pt),
  )

  Implementing the time convention *R2* therefore requires changing the time at which the stored scalar potential and charge density are interpreted to be (in `demos/chdr/ChDRSolver` etc). The current implementation in `demos/chdr` allows for the time placement `R2` combined with the spacial placement `staggered` (with its own operators, see @sec-O1VS-operators).

  The Mur boundary conditions need no change for *R2*.
  `AbsorbingBC.h` applies the second-order condition of Mur @mur1981AbsorbingBoundaryFdtd[Sec. IV, Eq. (16)] to each of the four potential components separately.
  For one component it computes the next slot at the outermost owned cell from the next, the current and the previous slot of that same component, at this cell, at its inward neighbour and at their tangential neighbours.
  This recurrence contains no source, no other component and no absolute time.
  Its coefficients depend only on $h_t$, the mesh spacings and the wave speed.
  Shifting the time labels of $phi$ by $h_t\/2$ therefore leaves it unchanged, and under *R2* the same condition holds for $phi$, centred at the half level $t_(n+1\/2)$ of its current slot instead of $t_n$.
  This is a deduction from the constant coefficients of the recurrence.

  The spatial placement enters in the same way.
  The recurrence is centred half a cell inside the outermost owned sample of each component.
  In the staggered placement the component $A_d$ normal to the face lies at $x_i + h_d\/2$ (see @sec-spatial-placement).
  Its Mur plane is therefore shifted by $h_d\/2$ in the positive $d$-direction relative to the plane of $phi$ and of the tangential components, at the lower and at the upper face.

  What *R2* does require at a Mur face are ghost values in both stored slots, because $vb(E)$ takes $grad phi$ from the previous slot and $vb(B)$ averages the curl over both slots.
  After each Mur update the solver exchanges the halos and then fills the exterior ghosts with the quadratic extrapolation $g = 3 u_0 - 3 u_1 + u_2$ (`OpenFaceClosure.h`), where $u_0$ is the outermost owned sample and $u_1$, $u_2$ are the next two samples inward.
  `prepareHistories()` does the same for loaded histories.

  All four components use one Mur speed per face, which is why O1-VS with Mur faces requires $alpha = 1$ (see @sec-O1VS-implementation).
  Applying the condition to each potential component does not guarantee that the gauge condition and Gauss's law hold at the boundary.
  @experiment-e04 measures these defects.

]
