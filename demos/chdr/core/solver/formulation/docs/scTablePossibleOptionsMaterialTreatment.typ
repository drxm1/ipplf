#import "../../../../docs/manual/preamble.typ": *

== Our table of possible options for the material/interface treatment
With no free surface sources: $[E_t] = 0, thick [epsilon E_n] = 0, thick [vb(B)] = 0, thick epsilon_2 E_(n,2) = epsilon_1 E_(n,1)$ (TODO: Write a section where we derive these). We have some options:
#v(5pt)

#table(
  columns: (87pt, 183pt, 1fr),
  align: (left + horizon, center + horizon, left),
  inset: (x: 4pt, y: 5pt),
  stroke: none,

  [#text(
    // fill: PSIblue,
    weight: "bold",
  )[O0 #h(2pt) Scale by $epsilon$\
    #text(size: 8pt, fill: black)[but stay in lorentz gauge]]],
  [
    $
      partial_t^2 vb(A) & = epsilon^(-1) nabla^2 vb(A) + epsilon^(-1) vb(J) \
        partial_t^2 phi & = epsilon^(-1) nabla^2 phi + epsilon^(-2) rho
    $
  ],
  [Negative control. Works for TE but is #text(fill: red)[incorrect] for example for TM at 45° incidence (TODO:SupportThisClaim).],

  table.hline(stroke: 0.4pt + luma(85%)),
  [#text(
      // fill: PSIblue,
      weight: "bold",
    )[O1 #h(2pt) Chew's gauge]
  ],
  [
    $ div (epsilon vb(A)) = -epsilon^2 partial_t phi $
    // #text(size: 7.5pt)[Generalized Lorenz gauge, $alpha = 1$.]
  ],
  [Modify the $(phi, vb(A))$ update to support this more general gauge. Elegant solution.],

  table.hline(stroke: 0.4pt + luma(85%)),
  [#text(
      // fill: PSIblue,
      weight: "bold",
    )[O4 #h(2pt) Polarization]// \
    // Charge + current
  ],
  [
    $ vb(P) = (epsilon - 1) vb(E) $
    $ rho_b = -div vb(P), quad vb(J)_b = partial_t vb(P) $
  ],
  [Add $(rho_b, vb(J)_b)$ computed from $vb(E)$ to the vacuum potential equations. $vb(P)$ then comes from the reconstructed $vb(E)$.],

  table.hline(stroke: 0.4pt + luma(85%)),
  [#text(
      // fill: PSIblue,
      weight: "bold",
    )[O3 #h(2pt) Yee]// \
    // Update $vb(E), vb(B)$
  ],
  [
    $
      epsilon partial_t vb(E) & = curl vb(B) - vb(J) \
              partial_t vb(B) & = -curl vb(E)
    $
  ],
  [Update $vb(E), vb(B)$ directly./* staggered both in space and time, with $epsilon$ at the same locations as $vb(E)$.*/],
)// END TABLE

*O0* extends the homogeneous Lorez equations across the interface. See @chew2014generalizedgauge[Eqs. (12)-(13)], with $mu=1$ (TODO:VerifySource).
That option is successful for source-free TE / normal-incidence cases, while it failes for oblique incidence TM cases. This was shown in a python benchmark (TODO:DocumentThisLater). It fails to represent the discontinuous jump in TODO.

*O4* TODO: See my 3w pres source file, haus1989Polarization and the py bench document paragraph on this option.
