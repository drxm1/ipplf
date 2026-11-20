#import "../../docs/manual/preamble.typ": *

#pagebreak()
= Experiment E10 detailed specs <experiment-e10>

E10 repeats #link(<experiment-e06>)[E06] in vacuum, $epsilon_r = 1$, with the four solver configurations of E09 (@E09-configurations).
Everything that is not stated otherwise below is the same as in E06:
the Gaussian profile without periodic wrapping, the $6w$ box with Mur faces at both ends of $z$, the ladder, the duration $7.5w\/v$,
the targets, the measurements E02.meas01–02 and the saved output.
As in E09 we use $4 times 4$ transverse cells and $h_t = h_z\/2$.

== Experiment E10 boundary and timing
With $v = 1$ the Mur speed is $v_"Mur" = 1$ on both faces.
All four configurations use the same IPPL boundary code (`FDTDBoundary.h`, `AbsorbingBC.h`, `OpenFaceClosure.h`) on each potential component at its own position.
For $A_x$ these positions have the same $z$ in both layouts.
The event times of E06 hold in units of $w\/v$, which are now units of $w$.

== Experiment E10 measurements
Each configuration measures E01.meas01–04 before contact and E02.meas01–02, at the sample locations and times of @E09-samples.

=== E10.meas01: differences between the configurations
As E09.meas01, from the saved profiles at every target, before contact and after the pulse has left.
In addition, compare the remainders $d_E$ and $d_B$ of the four configurations.

*Purpose:* show whether the configurations also agree when Mur acts on the pulse.
