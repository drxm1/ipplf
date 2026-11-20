#import "../../docs/manual/preamble.typ": *

= Experiment E11 detailed specs <experiment-e11>

E11 repeats #link(<experiment-e07>)[E07] in vacuum, $epsilon_r = 1$, with the four solver configurations of E09 (@E09-configurations).
Everything that is not stated otherwise below is the same as in E07:
the Gaussian dipole pulse, the box $[-6w, 6w]^3$, the levels $N = 96,dots,324$, the time step $h_t = h\/2$,
the travel distance $2w$, the initialization, the schedule and the analysis.

== Experiment E11 medium and timing
With $v = 1$ the duration is $T_3 = 2w\/v = 2w$.
The time step is below the vacuum Standard bound $h\/sqrt(3)$ (Fallahi @fallahi2020mithra20fullwavesimulation[Eq. (3.20)]).
The bounds of E07 hold also for $v = 1$:
at the native samples and the startup times the inward-traveling term of the generator is below $10^(-9)$ of the peak values of $vb(E)$, $vb(B)$, $vb(A)$ and $phi$,
and at the selected times every other periodic copy has $abs(s) >= 1.96w$.
The collocated samples are cell centers, which have at least the distance $sqrt(3)h\/2$ from the center.

== Experiment E11 measurements
Each configuration loads the histories and measures E03.meas01–05 at the sample locations and times of @E09-samples.
For Standard these are the cell centers, and $t_n$ for $vb(A)$, $phi$ and $vb(B)$.
The electric field of R1 combines the difference of $vb(A)$, centered at $t_(n-1/2) = (n-1/2)h_t$, with $phi^n$ at $n h_t$.
We compare it with the continuum field at $t_(n-1/2)$.
Its error then contains the term
$
  vb(E)_"R1" - vb(E)(t_(n-1/2)) = -(h_t\/2) nabla partial_t phi(t_(n-1/2)) + O(h_t^2 + h^2)
$
(`LevelContract`, a deduction by Taylor expansion), which is first order in $h_t$.
E03.meas06 is measured only for the two staggered configurations, because its operators are those of the staggered layout.

The two staggered configurations save the $z = 0$ section of E03 at the initial state and at the final field target on every level.
As in E03, the section interpolates every component trilinearly from its own sample locations (@E09-samples) to the points $(x_i, y_j, 0)$,
and the error of a sample is formed before the interpolation.

=== E11.meas01: difference between O1VS and Homogeneous
From the sections, compute the largest difference of $vb(E)$ and of $vb(B)$ between O1VS and Homogeneous at the final target,
relative to the largest $abs(vb(E))$ and $abs(vb(B))$ of the section.
Both configurations use the same positions and the same interpolation to the section, so the interpolation does not enter the difference.
The collocated configurations save no section: their positions and their interpolation differ from those of the staggered ones,
and for Standard also the time of $vb(B)$.
They are compared through the error norms against the continuum.

=== E11.meas02: samples along the central grid lines
With `write_native_lines`, every configuration saves $vb(A)$, $phi$, $vb(E)$ and $vb(B)$ in the cells of the three grid lines
through the cell with the index $N\/2 - 1$ in every direction, for the initial state and for the final target.
Each quantity is taken in the state that its own clock selects, and $vb(B)$ in the state of $vb(E)$, as in E03.meas01–05.
Every row holds one component at its own sample location (@E09-samples), its numerical value and the continuum value at that location and at the comparison time of the quantity.
Nothing is interpolated.
For Standard the continuum value of $vb(E)$ is taken at $t_(n-1/2)$, as in its error norm.
The loaded potentials of the initial state equal their continuum values.

*Purpose:* the error norms show whether `ChdrSolver` is as accurate as the Standard solver for a pulse with both potentials;
the direct difference shows whether O1VS and the Homogeneous update compute the same numbers;
the samples along the grid lines show the difference of every configuration to the continuum point by point.

== Experiment E11 fourth configuration StandardR2Prime <E11-standard-r2prime>
The time levels R1 of Standard are those of the authors of the solver.
Winkler puts the whole four-potential and the source at the time $n h_t$ (Winkler @winkler2024ChargeConservingFDTD[Sec. 3.1, p. 16, and Eqs. (3.1)–(3.4), p. 17]).
`StandardFDTDSolver::step()` updates the four components separately, so the update does not fix the time of $phi$ relative to $vb(A)$.
Only `evaluate_EB()` combines them, and its $vb(E)$ has a single time only if the current slot holds $phi$ half a step before $vb(A)$.
These are the time levels R2′ (`R2Prime` in `LevelContract`, @sec-time-placement).
StandardR2Prime is the Standard configuration with these time levels.
It loads $phi$ at $t_(-1/2)$ and $t_(-3/2)$ instead of $t_0$ and $t_(-1)$, and $vb(A)$ as Standard.
So at a given time its $phi$ has received one update more than that of O1VS, whose slots start at $t_(1/2)$ and $t_(-1/2)$.
Its $vb(E)$ is the output of `evaluate_EB()`, now at $t_(n-1/2)$ alone.
Its $vb(B)$ is the average of the central curls of both slots at $t_(n-1/2)$ (`ChdrSolver::reconstructEB`), not the $vb(B)^n$ of `evaluate_EB()`.
StandardR2Prime measures E03.meas01–05 as Standard.
It runs in all of E09 to E14.
Its first runs used an executable that differs from that of the other configurations only by this configuration.

*Purpose:* separate the first-order term of R1 from the collocated placement.
If that term is the only cause of the slower convergence of the error of $vb(E)$ of Standard, the error of StandardR2Prime converges with second order.
