#import "../../docs/manual/preamble.typ": *

= Experiment E09 detailed specs <experiment-e09>

E09 repeats #link(<experiment-e05>)[E05] in vacuum, $epsilon_r = 1$, with three solver configurations on the same meshes.
It compares `ChdrSolver` with the Standard solver of IPPL.
Everything that is not stated otherwise below is the same as in E05:
the Gaussian profile, the box, the initial center, the travel distance $40w$, the ladder $N_"pulse" = 32,dots,1024$,
the time step $h_t = h_z\/2$, the targets, the initialization and the saved output.
E10 to E14 use the same three configurations.

== Experiment E09 solver configurations <E09-configurations>
All three configurations run the same class, `ChdrSolver<StandardFDTDSolver>`, in the same executable.
None of them runs a bare `StandardFDTDSolver`.
A configuration sets three options of `ChdrSolver`:
the update that `step()` executes,
the placement of the stored fields (@sec-spatial-placement),
and the reconstruction of $vb(E)$ and $vb(B)$ with its time labels (@sec-time-placement).
All three use the inherited slot rotation `timeShift()` of `StandardFDTDSolver`.
#table(
  columns: (0.9fr, 1.8fr, 1fr, 1.4fr),
  inset: 5pt,
  align: left + top,
  stroke: 0.5pt,
  table.header([*Name*], [*Update in `step()`*], [*Placement*], [*Reconstruction*]),
  [O1VS],
  [`O1VSUpdate`: effective source of the O1VS formulation ($alpha = 1$), then `StandardFDTDSolver::step()` with $h_t$],
  [staggered],
  [R2, `reconstructStaggered`],

  [Homogeneous],
  [`HomogeneousUpdate`: `StandardFDTDSolver::step()` with $tau = h_t\/sqrt(epsilon_r)$],
  [staggered],
  [R2, `reconstructStaggered`],

  [Standard],
  [`HomogeneousUpdate`, as above],
  [collocated],
  [R1, `evaluate_EB()` of `StandardFDTDSolver`],
)
O1VS is the configuration of E01–E08.
O1VS and Homogeneous differ only in the update.
For $epsilon_r = 1$ the effective source of O1VS is zero and $tau = h_t$, so both execute IPPL's vacuum step (`O1VSUpdate.h`, `HomogeneousUpdate.h`).
Homogeneous and Standard execute the same update, but they store and reconstruct $vb(E)$ and $vb(B)$ differently.

In the Standard configuration every operation of a time step is a method of `StandardFDTDSolver`.
`ChdrSolver` adds only the loading of the initial histories and the time labels of `LevelContract`.
Before the studies we check this: a few steps of the Standard configuration must give the same $vb(A)$, $vb(E)$ and $vb(B)$ as a bare `StandardFDTDSolver` with the same loaded histories, bit for bit.

The staggered R2 reconstruction (`reconstructStaggered`) is
$
  E_d^(n-1/2) = -(A_d^n - A_d^(n-1))\/h_t - D_(+,d) phi^(n-1/2), quad
  vb(B)^(n-1/2) = 1/2 (C_+ vb(A)^n + C_+ vb(A)^(n-1)),
$
with the forward difference $D_(+,d)$ between the two cell centers next to the face of $E_d$ and the forward-difference curl $C_+$.
The collocated R1 reconstruction (`evaluate_EB()` in `FDTDSolverBase.hpp`) is
$
  vb(E)^n = -(vb(A)^n - vb(A)^(n-1))\/h_t - D_c phi^n, quad vb(B)^n = D_c times vb(A)^n,
$
with central differences $D_c$ over $2h$.

== Experiment E09 sample locations and times <E09-samples>
The table gives the location and the time of every sample after $n$ completed steps.
$vb(x)_bold(i)$ is the center of cell $bold(i)$, as in @sec-spatial-placement.
#table(
  columns: (0.7fr, 1.6fr, 1.6fr),
  inset: 5pt,
  align: left + top,
  stroke: 0.5pt,
  table.header([*Quantity*], [*O1VS and Homogeneous (staggered, R2)*], [*Standard (collocated, R1)*]),
  [$phi$], [$vb(x)_bold(i)$, $t_(n+1/2)$], [$vb(x)_bold(i)$, $t_n$],
  [$A_d$], [$vb(x)_bold(i) + (h_d\/2) vb(e)_d$, $t_n$], [$vb(x)_bold(i)$, $t_n$],
  [$E_d$], [$vb(x)_bold(i) + (h_d\/2) vb(e)_d$, $t_(n-1/2)$], [$vb(x)_bold(i)$, $t_(n-1/2)$ for the difference of $vb(A)$, $t_n$ for $D_c phi^n$],
  [$B_d$], [$vb(x)_bold(i) + sum_(q != d) (h_q\/2) vb(e)_q$, $t_(n-1/2)$], [$vb(x)_bold(i)$, $t_n$],
)
These are the locations of `LevelContract::positionOffset` and the times of `LevelContract::labelTime`.
The initial histories are loaded at these locations and times (`writeSamples` in `BoxSetup.h`).
Every error norm, every profile and every frame compares a sample with the continuum value at the location and the time of that sample.
No interpolation enters the error norms.
In E09 and E10 $phi = 0$, so the electric field of R1 belongs to $t_(n-1/2)$ alone.

== Experiment E09 medium and mesh
With $epsilon_r = 1$ the speed is $v = c_"vac" = 1$, so $T = 40w\/v = 40w$.
We use $N_x = N_y = 4$ transverse cells ($h_x = h_y = w\/4$) instead of 8.
E05 showed that the transverse mesh does not change the results, and it halves the cost.
The time step $h_t = h_z\/2$ is below the vacuum Standard bound (Fallahi @fallahi2020mithra20fullwavesimulation[Eq. (3.20)]),
and it is also the step that `StandardFDTDSolver` selects itself, $min(h)\/2$.

== Experiment E09 measurements
Each configuration measures E01.meas01–04 and saves the profiles and frames of E05, at the sample locations and times of @E09-samples.
For Standard these are the cell centers, with $A_x$ and $B_y$ at $t_n$ and $E_x$ at $t_(n-1/2)$.

=== E09.meas01: differences between the configurations
From the saved profiles, compute for each level, target and pair of configurations the largest difference of the numerical $A_x$ and $E_x$,
relative to the largest continuum value of the same component at that target.
The pairs are O1VS–Homogeneous and Homogeneous–Standard.
$A_x$ and $E_x$ have the same $z$ and the same time in all three configurations.
They differ only in $x$: $x_i + h_x\/2$ in the staggered placement and $x_i$ in the collocated one.
The plane pulse does not depend on $x$ and $y$, so this offset does not matter.
The profile line holds all information, since the fields depend only on $z$ and $t$.
$B_y$ is not compared directly: the layouts place it at different positions and times and compute it with different stencils.

*Purpose:* show whether the three configurations compute the same potential and electric field.
A difference at rounding level means the same numbers, while a larger difference identifies the configuration that changes the result.
