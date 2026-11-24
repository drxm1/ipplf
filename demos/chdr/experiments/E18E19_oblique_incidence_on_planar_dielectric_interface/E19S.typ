#import "../../docs/manual/preamble.typ": *
#import "obliquePacketValues.typ": *

= Experiment E19-S detailed specs <experiment-e19s>

E19-S repeats #link(<experiment-e19>)[E19] with the three differences of #link(<experiment-e18s>)[E18-S]:
the time step $h_t = 0.99 thin h_t^"max"$ with the last step $n_T$ of E18-S,
the reflection magnitude of the continuum field through the same sums,
and the Richardson extrapolation from the two finest levels.
Everything that is not stated otherwise below is the same as in E19.

== Experiment E19-S measured quantities
E19.meas01–02 and E18.meas03 are used unchanged.

*E19-S.meas04* is E18-S.meas04 (@E18S-window-reference) with the field and the reference of E19.meas01:
the continuum field $B_y^"ref"$ of @E19-analytic-EB at the R2 samples of $B_y$, which lie at $(x_i + h\/2, z_k + h\/2)$.
The samples on the face, $z_k + h\/2 = z_f$, belong to side 1 and use the expression for $z < z_f$.
Both expressions agree there, because $B_y$ is continuous.
The split @E18S-split uses $abs(R(q))$ in place of $abs(Gamma(q))$.

*E19-S.meas05* is E18-S.meas05 (@E18S-richardson) for $r^"num"$ of E19.meas01 and for $rho_S$ of E19.meas02.
$rho_1^"R"$ is compared with $macron(R)_"TM"$ (@E19-energy-fractions).

At the Brewster wavenumber $q_"B"$ the reference $abs(R)$ is zero,
and $r^"win" (q_"B")$ is the part that the measurement alone produces there.
The extrapolation of the magnitude at $q_"B"$ needs the sign condition of E18-S.meas05.

== Expected behaviour
As in E18-S.
In E19 the observed order of $abs(r^"num" - abs(R))$ is close to two at all twelve wavenumbers on all levels,
so the measurement part is expected to be smaller than in E18 compared with the error of the update.

== Experiment E19-S saved output and analysis
As in E18-S, with the inputs `inputs/E19-S-n<N>.yaml` and the case `E19-S`.
`reflection.csv` holds $abs(R(q_j))$ as the reference, as in E19.
