#import "../../docs/manual/preamble.typ": *

#pagebreak()
= Experiment E09 detailed specs <experiment-e09>

E09 repeats #link(<experiment-e05>)[E05] in vacuum, $epsilon_r = 1$, with four solver configurations on the same setup.
It compares `ChdrSolver` with the `StandardFDTDSolver`.
Everything that is not stated otherwise below is the same as in E05:
the Gaussian profile, the box, the initial center, the travel distance $40w$, the ladder $N_"pulse" = 32,dots,1024$,
the time step $h_t = h_z\/2$, the targets, the initialization and the saved output.
E10 to E14 use the same four configurations.
The fourth one, StandardR2Prime, is described in @E11-standard-r2prime.

== Experiment E09 solver configurations <E09-configurations>
All configurations run the same class, `ChdrSolver<StandardFDTDSolver>`, in the same executable.
None of them runs a bare `StandardFDTDSolver`.
A configuration sets three options of `ChdrSolver`:
the update that `step()` executes,
the spatial placement of the stored fields (@sec-spatial-placement),
and the time levels R1, R2 or R2′, which give the time of each potential slot and each output field (@sec-time-placement).
The reconstruction function of $vb(E)$ and $vb(B)$ follows from the placement and the time levels.

#table(
  columns: (0.9fr, 1.8fr, 1fr, 1.4fr),
  inset: 5pt,
  align: left + top,
  stroke: 0.5pt,
  table.header(
    [*Name*],
    [*Update in `step()`*],
    [*Placement*],
    [*Time levels, reconstruction*],
  ),
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

  [StandardR2Prime],
  [`HomogeneousUpdate`, as above],
  [collocated],
  [R2′, `evaluate_EB()` for $vb(E)$, average of two central curls for $vb(B)$],
)
O1VS is the configuration of E01–E08.
O1VS and Homogeneous differ only in the update.
For $epsilon_r = 1$ the effective source of O1VS is zero and $tau = h_t$, so both execute `StandardFDTDSolver`'s vacuum step (`O1VSUpdate.h`, `HomogeneousUpdate.h`).
Homogeneous and Standard execute the same update, but they differ in placement and time levels, and so in the reconstruction of $vb(E)$ and $vb(B)$.

In the Standard configuration the field update is `StandardFDTDSolver::step()` and the reconstruction is its `evaluate_EB()`.
Around `step()`, `HomogeneousUpdate` only copies the source, with the charge divided by $epsilon_r$, and sets the time step $tau$.
Otherwise `ChdrSolver` adds only the loading of the initial histories and the time labels of `LevelContract`.
Before the studies we check this at $epsilon_r = 1$: a few steps of the Standard configuration must give the same $vb(A)$, $vb(E)$ and $vb(B)$ as a bare `StandardFDTDSolver` with the same loaded histories, bit for bit.

With the staggered placement and R2, the reconstruction (`reconstructStaggered`) is
$
  E_d^(n-1/2) = -(A_d^n - A_d^(n-1))\/h_t - D_(+,d) phi^(n-1/2), quad
  vb(B)^(n-1/2) = 1/2 (C_+ vb(A)^n + C_+ vb(A)^(n-1)),
$
with the forward difference $D_(+,d)$ between the two cell centers next to the face of $E_d$ and the forward-difference curl $C_+$.
With the collocated placement and R1, the reconstruction (`evaluate_EB()` in `FDTDSolverBase.hpp`) is
$
  vb(E)^n = -(vb(A)^n - vb(A)^(n-1))\/h_t - D_c phi^n, quad vb(B)^n = D_c times vb(A)^n,
$
with central differences $D_c$ over $2h$.
The placement selects the stencils and the positions of the components.
The time levels select the slot of $phi$ that $vb(E)$ reads and whether $vb(B)$ is averaged.
In R2 the slot $n$ holds $phi$ at $t_(n+1/2)$, so $vb(E)$ reads the previous slot, and $vb(B)$ is the average of two curls.
In R1 the slot $n$ holds $phi$ at $t_n$, and $vb(B)$ is not averaged.
In R2′ the slot $n$ holds $phi$ at $t_(n-1/2)$, so $vb(E)$ reads the current slot, and $vb(B)$ is the average of two curls.

== Experiment E09 sample locations and times <E09-samples>
The table gives the location and the time of every sample after $n$ completed steps.
$vb(x)_bold(i)$ is the center of cell $bold(i)$, as in @sec-spatial-placement.
#table(
  columns: (0.6fr, 1.5fr, 1.5fr, 1.1fr),
  inset: 5pt,
  align: left + top,
  stroke: 0.5pt,
  table.header(
    [*Quantity*],
    [*O1VS and Homogeneous (staggered, R2)*],
    [*Standard (collocated, R1)*],
    [*StandardR2Prime (collocated, R2′)*],
  ),
  [$phi$], [$vb(x)_bold(i)$, $t_(n+1/2)$], [$vb(x)_bold(i)$, $t_n$], [$vb(x)_bold(i)$, $t_(n-1/2)$],
  [$A_d$],
  [$vb(x)_bold(i) + (h_d\/2) vb(e)_d$, $t_n$],
  [$vb(x)_bold(i)$, $t_n$],
  [$vb(x)_bold(i)$, $t_n$],

  [$E_d$],
  [$vb(x)_bold(i) + (h_d\/2) vb(e)_d$, $t_(n-1/2)$],
  [$vb(x)_bold(i)$, $t_(n-1/2)$ for the difference of $vb(A)$, $t_n$ for $D_c phi^n$],
  [$vb(x)_bold(i)$, $t_(n-1/2)$],

  [$B_d$],
  [$vb(x)_bold(i) + sum_(q != d) (h_q\/2) vb(e)_q$, $t_(n-1/2)$],
  [$vb(x)_bold(i)$, $t_n$],
  [$vb(x)_bold(i)$, $t_(n-1/2)$],
)
These are the locations of `LevelContract::positionOffset` and the times of `LevelContract::labelTime`.
The initial histories are loaded at these locations and times (`writeSamples` in `BoxSetup.h`).
Every error norm and every profile compares a sample with the continuum value at the location and the time of that sample.
The frames store the samples with these locations and times.
No interpolation enters the error norms.
At each target, $vb(A)$ and $phi$ are measured after the step whose time of that quantity is nearest to the target,
and $vb(E)$ and $vb(B)$ after the step chosen with the time of $vb(E)$, as in E01 and E03.
So $vb(B)$ of Standard ($t_n$) is half a step later than $vb(B)$ of O1VS and Homogeneous ($t_(n-1/2)$).
Each is compared with the continuum at its own time.
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
For StandardR2Prime they are the cell centers too, with $B_y$ at $t_(n-1/2)$.

=== E09.meas01: differences between the configurations
From the saved profiles, compute for each level, target and pair of configurations the largest difference of the numerical $A_x$ and $E_x$,
relative to the largest continuum value of the same component at that target.
Every pair of the four configurations is compared.
$A_x$ and $E_x$ have the same $z$ and the same time in all four configurations.
They differ only in $x$: $x_i + h_x\/2$ in the staggered placement and $x_i$ in the collocated one.
The plane pulse does not depend on $x$ and $y$, so this offset does not matter.
The profile line holds all information, since the fields depend only on $z$ and $t$.
$B_y$ is not compared directly: the layouts place it at different positions and times and compute it with different stencils.

*Purpose:* show whether the four configurations compute the same potential and electric field.
A difference at rounding level means the same numbers, while a larger difference identifies the configuration that changes the result.
