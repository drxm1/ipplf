#import "../../docs/manual/preamble.typ": *

= Experiment E12 detailed specs <experiment-e12>

E12 repeats #link(<experiment-e08>)[E08] in vacuum, $epsilon_r = 1$, with the four solver configurations of E11 (@E09-configurations, @E11-standard-r2prime).
Everything that is not stated otherwise below is the same as in E08:
the Gaussian dipole pulse, the all-face Mur box $[-6w, 6w]^3$, the levels, the time step $h_t = h\/2$,
the duration $16w\/v$, the event table in units of $eta = v t\/w$ and the late targets.

== Experiment E12 boundary and timing
With $v = 1$ the Mur speed is $1$ on all faces, edges and corners, and the duration is $T_4 = 16w$.
The runs advance through $ceil(T_4\/h_t + 1\/2)$ steps: 257, 385, 577 and 865 on the four meshes.

== Experiment E12 measurements
Each configuration measures E04.meas01 at the sample locations and times of @E09-samples.
E12 compares the energy left in the box.
The energy-weighted field error of E04.meas01 is saved as in E04, but E12 does not compare it.
The energy sums the square of every component at its own sample location, without interpolation, as in E04.
For Standard (R1), $vb(E)$ belongs to $t_(n-1/2)$ as in E11, while $vb(B)^n$ belongs to $t_n$.
To align both in time, the energy uses $(vb(B)^n + vb(B)^(n-1))\/2$ at $t_(n-1/2)$, the average of the R2 time levels.
This is the $vb(B)$ of the collocated placement with R2 (`solver::reconstructMagnetic` with averaging):
the central curl of `evaluate_EB()` applied to both potential slots, at the cell centers.
StandardR2Prime has $vb(E)$ and this averaged $vb(B)$ at $t_(n-1/2)$ already, and its energy uses them.
The monitor flux E04.meas02, the constraint defects E04.meas03 and the periodic control are measured only for the two staggered configurations:
their interpolation and operators are those of the staggered layout.

=== E12.meas01: differences between the configurations
Compare the energy traces of the four configurations at every step, relative to the trace of O1VS.
O1VS and Homogeneous use the same layout and samples, so a difference above rounding comes from the update.

*Purpose:* show whether `ChdrSolver` and the Standard solver leave the same energy in the box after Mur,
including the incidence angle effect at $eta = 10$ (E08).

=== E12.meas02: parts of the initial energy
At the initial state, before the first update, save three sums of the reconstructed fields:
$
  U_E = h^3/2 sum_(bold(i),c) epsilon abs(E_(c,bold(i)))^2, quad
  U_(overline(B)) = h^3/2 sum_(bold(i),c) abs(overline(B)_(c,bold(i)))^2, quad
  U_(B^2) = h^3/2 sum_(bold(i),c) 1/2 (abs(B_(c,bold(i))^0)^2 + abs(B_(c,bold(i))^(-1))^2).
$
Here $vb(B)^0$ and $vb(B)^(-1)$ are the curls of the two potential slots with the curl of the placement,
and $overline(vb(B)) = (vb(B)^0 + vb(B)^(-1))\/2$ is the averaged field of the energy.
Every component is squared at its own sample location (@E09-samples).
$U_E + U_(overline(B))$ is the numerical energy of E04.meas01 at the initial state.
Normalize $U_E$ by the same sum over the continuum $vb(E)$ at $t_(-1/2)$,
and $U_(overline(B))$ and $U_(B^2)$ by the same sum over the continuum $vb(B)$ at $t_(-1/2)$.
Since $abs((a+b)\/2)^2 = (abs(a)^2 + abs(b)^2)\/2 - abs((a-b)\/2)^2$,
the difference $U_(B^2) - U_(overline(B))$ is the energy that the average of the two curls removes.
The measurement is optional (`write_energy_parts`) and runs on the coarsest mesh, $N = 96$.

*Purpose:* show which part of the initial energy deficit comes from the electric reconstruction,
from the curl of the placement and from the average of the two curls.
