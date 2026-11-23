// Standalone: typst compile E15.typ
// 2026-11-23: added the independently initialized million-step run to check
// accumulation at fixed resolution. The two projections and their residual
// supply all reported observables.
// 2026-11-23: report the continuum amplitude error and the drift from the first
// sample separately to distinguish reconstruction accuracy from long-run stability.
#import "@preview/physica:0.9.8": *
#set document(title: "E15. Long-run propagation in a periodic dielectric")
#set page(paper: "a4", margin: (x: 22mm, y: 18mm), numbering: "1")
#set text(font: "Libertinus Serif", size: 11pt, lang: "en")
#set par(leading: 0.55em)
#set heading(numbering: none)
#show math.equation: set text(font: "New Computer Modern Math")
#show raw: set text(size: .88em)

= E15. Long-run propagation in a periodic dielectric <experiment-e15>

Does the O1VS update preserve the amplitude of a single travelling wave as its phase slips relative
to the continuum wave? We follow the electric and magnetic fields for 1000 continuum periods on
three meshes, then repeat the coarse run for one million steps to measure drift at fixed resolution.

== Analytic continuum reference
The unit cube is periodic in all three directions and filled with a homogeneous, lossless
dielectric of relative permittivity $epsilon_r=2.13$ and unit permeability.
Charge and current vanish. In solver units the wave speed is $v=1/sqrt(epsilon_r)$.
The wave vector, polarization, angular frequency and period are
$ vb(k)=2pi(1,1,1), quad vu(p)=1/sqrt(2) (1,-1,0), quad
  omega=v abs(vb(k)), quad T=2pi/omega. $
With $theta=vb(k) dot vb(x)-omega t$ and $vu(k)=vb(k)/abs(vb(k))$, take
$ vb(A)=(vu(p)+vu(k))cos theta, quad phi=v cos theta. $
The longitudinal vector potential and the scalar potential give equal and opposite contributions to $vb(E)$.
From $vb(E)=-partial_t vb(A)-grad phi$ and $vb(B)=curl vb(A)$
@chew2014generalizedgauge[§1, Eqs. (5)–(6)], the physical fields are the real parts of
$ vb(cal(E))=i omega vu(p) e^(i theta), quad
  vb(cal(B))=i (vb(k) times vu(p)) e^(i theta). $

== Initialization and runs
Use staggered O1VS with $alpha=1$ and R2 reconstruction. For $N$ cells per axis,
set $h=1/N$, $h_t=h/2$ and $t_j=j h_t$. The angular frequency of this discrete mode is
$ a=pi/N, quad b=arcsin(sqrt(3/epsilon_r) (h_t/h) sin a), quad omega_d=2b/h_t. $
This follows by inserting the diagonal Fourier mode into the three-dimensional leapfrog update,
following Schneider @schneider2010UnderstandingFdtd[§7.4].
Initialize the potentials above with $omega_d$ in the phase to excite a single forward discrete wave.
At the positive $d$-face, load $A_d$ at $t_0$ and $t_(-1)$.
At the cell center, load $phi$ at $t_(1/2)$ and $t_(-1/2)$ in the current and previous slots,
respectively (R2, #raw("LevelContract")).

#align(center, table(
  columns: (auto, auto, auto, auto), align: center, stroke: none, inset: (x: 9pt, y: 3pt),
  table.hline(), table.header([Input], [$N$], [Final step], [Observation stride (steps)]),
  table.hline(stroke: .4pt),
  [#raw("n16.yaml")], [16], [26 964], [27],
  [#raw("n32.yaml")], [32], [53 928], [54],
  [#raw("n64.yaml")], [64], [107 855], [108],
  [#raw("n16-million.yaml")], [16], [1 000 000], [27], table.hline(),
))
Use eight MPI ranks along $z$ and four OpenMP threads per rank.
The three refinement runs stop at the R2 sample nearest to $1000 T$.
Start the extended run independently with the same coarse-grid histories.
Record step zero, every $floor(T/h_t+1/2)$ steps and the final step, once each.
Compare the two coarse histories at their shared step indices.

== Expected behaviour <e15-expected>
The equal wave-vector components and cubic cells keep the discrete gradient parallel to $vb(k)$,
so the longitudinal field cancels also on the mesh. The staggered spatial difference gives
the factor $sin(a)/a$, and averaging the two curls for R2 gives an additional $cos b$ in $vb(B)$
(@schneider2010UnderstandingFdtd[§§7.3–7.4]; native placement in @ryu2016potentialfdtd[Appendix A]).
Thus the predicted amplitude ratios and continuum-relative phase lag are
$ kappa_E=sin(a)/a, quad kappa_B=kappa_E cos b, quad delta_"pred"(t)=(omega-omega_d)t. $

#pagebreak()
== Measured quantities <e15-measurements>
After $n$ completed steps both fields have time $t=(n-1/2)h_t$, including the initial
observation at $-h_t/2$. For a cell centered at $vb(x)_bold(i)$,
sample $E_d$ at $vb(x)_bold(i)+h bold(e)_d/2$ and $B_d$ at
$vb(x)_bold(i)+h(bold(e)_p+bold(e)_q)/2$, where $p,q$ are the two axes other than $d$.
Evaluate each continuum component at the same position and time
(#raw("LevelContract"), @ryu2016potentialfdtd[Appendix A]).
Below, $sum_(bold(i),d)$ covers all owned cells and the three components.

*E15.meas01 and E15.meas02 — electric and magnetic mode coefficients.*
For $F=E,B$, write $F_(bold(i)d)$ for the numerical sample and $cal(F)_(bold(i)d)$
for its complex continuum reference. At each observation form
$ c_F = (2 sum_(bold(i),d) F_(bold(i)d) overline(cal(F)_(bold(i)d))) /
        (sum_(bold(i),d) abs(cal(F)_(bold(i)d))^2), quad
  rho_F=abs(c_F), quad delta_F=arg(c_F). $
The overline denotes complex conjugation. The factor of two follows from the equal norms
of the orthogonal sine and cosine samples. Unwrap each phase in observation order.
For each field the signed continuum amplitude error is $rho_F-1$
and the drift from the first sample is $D_(F) (t) = frac(rho_(F) (t), rho_(F) (-h_t/2)) - 1$.

=== Residual outside the wave
*E15.meas03.* Define the sampled field energy and its continuum reference by
$ U(E,B) = h^3/2 sum_(bold(i),d) (epsilon_r E_(bold(i)d)^2 + B_(bold(i)d)^2), quad
  U_"ref"=U(op("Re")cal(E),op("Re")cal(B)). $
The weights are those of the electric and magnetic energy in the dielectric
@ellingson2020ElectromagneticsVol2[§3.1, Eqs. (3.20), (3.22)].
Subtract the fitted fields pointwise before summing the residual energy,
$ r = sqrt(U(E-op("Re")(c_E cal(E)), B-op("Re")(c_B cal(B)))/U_"ref"). $
In exact arithmetic $r=0$. A longitudinal electric field from incomplete gauge cancellation
contributes to $r$ even when the fitted transverse amplitude is constant.
Orthogonality and equal electric and magnetic reference energies give
$U/U_"ref"=(rho_E^2+rho_B^2)/2+r^2$. For the energy history use $U/U_0-1$,
where $U_0$ is the energy inferred at the initial observation.

== Output and comparison
Each run writes #raw("run.json") with its configuration and #raw("results.json") with
#raw("N"), #raw("period"), #raw("kappa_E"), #raw("kappa_B") and #raw("phase_lag_rate") ($omega-omega_d$).
Each row of #raw("observations") records #raw("n"), #raw("t"), #raw("c_E_re"), #raw("c_E_im"),
#raw("c_B_re"), #raw("c_B_im") and #raw("r"). The middle mesh also writes #raw("slice-<step>.json")
with native $E_x$ on the central $x y$ layer at steps 0, 14, 27 and 41.

At 1000 periods, halving $h$ should roughly quarter the final lag, whose leading term is
proportional to $h^2$. Report $max_t abs(rho_F-1)$ for each field.
For the million-step run also give the largest sampled drift $abs(D_F)$ over both fields,
the peak residual and the largest $abs(U/U_0-1)$.

#show bibliography: set text(size: 9pt)
#bibliography("E15-references.bib", style: "ieee", title: none)
