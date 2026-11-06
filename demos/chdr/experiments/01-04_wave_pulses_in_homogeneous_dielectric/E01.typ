#import "@preview/physica:0.9.8": *
#set math.equation(numbering: "(1)")

= Experiment E01 detailled specs

== Experiment E01 analytic continuum reference
As the continuum analytic reference for E01 we choose a transverse potential pulse.\

We assume a nonmagnetic, lossless, nondispersive medium with $epsilon_r = 2.13$.
The material is homogeneous and spans the whole mesh domain. There is no free charge or current.

We first define a family of exact continuum solutions. Let
$
  vb(A)(z,t) = vb(hat(x)) a(s(z,t)), quad s(z,t) := z - z_0 - v t, quad phi = 0, quad v = c_"vac"\/sqrt(epsilon_r);
$ <E01-analytic-potentials>

The continuum fields $vb(E)$ and $vb(B)$ can be reconstructed via Chew §1, Eqs. (5)-(6):
$
  vb(E) = -partial_t vb(A) - grad phi, quad vb(B) = curl vb(A)
$

The chain rule for the potential in @E01-analytic-potentials gives $partial_t A_x = -v a'(s)$ and $partial_z A_x = a'(s)$. Only $(curl vb(A))_y = partial_z A_x$ is nonzero in the curl. By substitution we obtain the continuum $vb(E), vb(B)$:
$
  vb(E)(z,t) = vb(hat(x)) v pdv(a(s),s)(z,t), quad vb(B)(z,t) = vb(hat(y)) pdv(a(s),s)(z,t).
$ <E01-analytic-EB>

@E01-analytic-potentials and @E01-analytic-EB hold for all $t$, including the negative time that is used to initialize the previous potential `A_nm1`. We call $a$ a _potential profile_.

=== Proof: The potential family solves the homogeneous equations
Since $A_x$ is independent of x, $div vb(A) = 0$, so the homogeneous Lorenz gauge $div vb(A) = -epsilon_r partial_t phi$ TODO:chew[Eq. (11)] is satisfied.

We now verify the continuum wave equations for the scalar and vector potentials
in a homogeneous medium, with zero free charge and current. In the Lorenz
gauge these are
#link("https://arxiv.org/pdf/1406.4780v1")[Chew, §2, Eqs. (12)–(13)]:
$
  nabla^2 phi - epsilon_r partial_t^2 phi &= 0, \
  nabla^2 vb(A) - epsilon_r partial_t^2 vb(A) &= vb(0).
$
Since $A_x (z,t) = a(s(z,t))$ is independent of $x$ and $y$, we obtain $nabla^2 A_x = partial_z^2 A_x$.
If $a in cal(C)^2 ("Dom"(a))$, $partial_z^2 A_x - epsilon_r partial_t^2 A_x = a''(s) - a''(s) = 0$. Therefore this family of potentials provides source free Maxwell solutions.

=== Specification of the profile
We choose a smooth pulse that is exactly zero outside of a finite interval. Let
$
  p(s) = cases(
    A_0 exp(1 - 1/(1 - (s\/w)^2)) "if" abs(s) < w,
    0 "if" abs(s) >= w
  )\
  "with" 0 < 2w < L_z
$
Alternatively we could have chosed a Gaussian potential profile $p(s) = A_0 e^(-s^2\/(2 sigma^2))$, but it has the disadvantage that its tail does not vanish. $L_z$ is the domain length.

Differentiation gives
$
  p'(s) = cases(
    - (2s p(s))/(w^2 (1-(s\/w)^2)^2) & "if" abs(s) < w,
    0 & "if" abs(s) >= w
  ).
$

#import "@preview/cetz:0.3.4": canvas, draw

#let profile-values(u) = {
  if calc.abs(u) >= 1 {
    (0, 0)
  } else {
    let q = 1 - u * u
    let potential = calc.exp(1 - 1 / q)
    (potential, -2 * u * potential / (q * q))
  }
}

#let profile-axes(position, ticks, ordinate) = {
  draw.rect((0, 0), (5.2, 3), stroke: 0.35pt + luma(65%))
  draw.line(position(-1.35, 0), position(1.35, 0), stroke: 0.35pt + luma(75%))
  for u in (-1, 0, 1) {
    let x = position(u, 0).at(0)
    draw.line((x, 0), (x, -0.08), stroke: 0.35pt + luma(65%))
    draw.content((x, -0.16), anchor: "north", [$#u$])
  }
  for value in ticks {
    let y = position(0, value).at(1)
    draw.line((0, y), (-0.08, y), stroke: 0.35pt + luma(65%))
    draw.content((-0.16, y), anchor: "east", [$#value$])
  }
  draw.content((0, 3.16), anchor: "south-west", ordinate)
  draw.content((2.6, -0.54), anchor: "north", [$s\/w$])
}

#let profile-panel(derivative: false) = {
  let limits = if derivative { (-2.5, 2.5) } else { (-0.05, 1.1) }
  let ticks = if derivative { (-2, -1, 0, 1, 2) } else { (0, 0.5, 1) }
  let ordinate = if derivative { [$w p'(s)\/A_0$] } else { [$p(s)\/A_0$] }
  let position(u, value) = (5.2 * (u + 1.35) / 2.7,
    3 * (value - limits.at(0)) / (limits.at(1) - limits.at(0)))
  let points = range(501).map(i => {
    let u = -1.35 + 2.7 * i / 500
    position(u, profile-values(u).at(if derivative { 1 } else { 0 }))
  })
  set text(size: 9pt)
  canvas(length: 10mm, {
    profile-axes(position, ticks, ordinate)
    draw.line(..points, stroke: 1pt + blue)
  })
}

#figure(
  grid(columns: (1fr, 1fr), column-gutter: 8pt, align: center,
    profile-panel(), profile-panel(derivative: true)),
  caption: [Profile, in dimensionless form.],
) <fig-e01-profile>

We define the periodic profile
$
  a(s) = sum_(m in ZZ) p(s-m L_z), quad a(s+L_z) = a(s).
$
Since all derivatives of $p$ vanish at $s = plus.minus w$, the joins are smooth.
Since $2w < L_z$, at most one term in the sum is nonzero.
Therefore we can evaluate the profile $a(s(z,t))$ directly by measuring the displacement from the nearest periodic copy:
$
  r := s - L_z floor(s/L_z + 1/2), quad -L_z/2 <= r < L_z/2,
$
$
  a(s) = p(r), quad a'(s) = p'(r).
$
We return 0 in code for $a, a'$ when $abs(r) >= w$. The jump in $r$ at the periodic join causes no $a$-discontinuity, because since $2w < L_z$, $p$ is zero in a neighbourhood there. Pseudocode:

  #import "@preview/lovelace:0.3.1": pseudocode-list
  #block(breakable: false)[
    #pseudocode-list(
      line-numbering: none,
      line-gap: 0.6em,
      indentation: 1em,
      stroke: 0.4pt + blue,
    )[
      + *input:* $z, t, z_0, v, L_z, w, A_0$, with $0 < 2w < L_z$
      + *output:* $(a(s(z,t)), a'(s(z,t)))$
      + $s <- z-z_0-v t$
      + $r <- s-L_z floor(s/L_z+1/2)$
      + *if* $abs(r) >= w$:
        + *return* $(0,0)$
      + $q <- 1-(r\/w)^2$
      + $"potential" <- A_0 exp(1-1/q)$
      + $"derivative" <- -(2r "potential")/(w^2 q^2)$
      + *return* $("potential", "derivative")$
    ]
  ]

== Experiment E01 initialization/sampling method
We choose initialization option A described in `../experiment_definition_draft.md`,
that is we sample the continuum potential history directly:
$
  vb(A)^0 eq.triple_"stagg" vb(A)_"exact" (0), quad vb(A)^(-1) eq.triple_"stagg" vb(A)_"exact" (- h_t)
$

Here $h_t$ is the time step, and the superscripts are time-level indices. Time-level zero is by convention at $t=0$, and level minus one at $t = -h_t$. In this notation we also use half-levels to denote when a discrete expression lives a half step $h_t\/2$ apart due to it being calculated by a difference operator.

A reader might wonder why the triple equality with subscript "stagg", and what this should mean.\
In order to distinguish the discrete samples stored in the virtually staggered `O1-VS` convention, we do the following:
Let $(x_i, y_j, z_k)$ denote a cell centre, and we write the spacial indices explicitly. The indices $(i,j,k)$ denote that a field is stored in the code such that it corresponds to the cell with the corresponding centre. We additionally add $+1/2$ or $-1/2$ to the index of a field component if that field component physically lives a half-cell shifted in the direction that corresponds to the said index. This is necessary because in the spatially staggered implementation some stored field components physically have a position that is offset from the cell-centre due to the finite difference operators being used to calculate certain derivatives.

In the convention used by the `O1-VS` update, $A_x$ is for example located at $(x_i + h_x\/2, y_j, z_k)$.  Similarly the locations of $A_y$ and $A_z$ are shifted in $y$ and $z$. These offsets, in code, are returned by `LevelContract::positionOffset`. Their staggered grid "basis" (TODO: What does this mean?) is Ryu et al., Appendix A, Eqs. (19)-(24) TODO:ryu2016potentialfdtd.

For E01, sampling the only nonzero component therefore means the following, which is the meaning of the above $eq.triple_"stagg"$-notation (we do not use that notation anymore from here on).
$
  (A_x^n)_(i+1/2,j,k)
  = A_(x,"exact")(z_k, n h_t)
  = a(z_k - z_0 - v n h_t), quad n in {0,-1}.
$

In the code, the components are still stored at integer indices $(i,j,k)$, and `writeSamples` in `core/solver/BoxSetup.h` adds the physical offset before the analytic is evaluated.

The magnetic component $B_y$ has a different location:
$(x_i+h_x/2, y_j, z_k+h_z/2)$. Its continuum comparison must use
$z_k+h_z\/2$. Thus the half-cell shift in $z$ must be considered in the magnetic-field
comparison, but not in the initialization of $A_x$.

We must measure the initial field error as part of the experiment,
since the initialized electric and magnetic fiels will only approximate the analytic field.
For example the initialized electric field will be
$
  vb(E)^(-1/2) = - (vb(A)^0 - vb(A)^(-1))/h_t eq.triple^"approx"_"stagg" - partial_t vb(A)
$

== Experiment E01 measured quantities and their purpose

=== E01.meas01: Relative electric-field error in the global discrete $L^2$ norm
Let $bold(e)_c$ be the unit index vector in direction $c in {x,y,z}$, and let $n$ count
completed solver steps. The numerical sample
$E_(c,bold(i)+1/2 bold(e)_c)^("num",n-1/2)$ is located at
$vb(x)_bold(i) + (h_c/2) bold(e)_c$, where $vb(x)_bold(i)$ is the cell centre.
The reference sample $E_(c,bold(i)+1/2 bold(e)_c)^("ref",n-1/2)$ is the
continuum solution evaluated at that same position and at the R2 observation
time $t_(n-1/2) = (n-1/2)h_t$. Let $sans(M)$ be the set of all cell indices on the computational domain. Define
$
  e_(E,thin"L2_global")^(n-1/2) = sqrt(
    (sum_(bold(i) in sans(M),c) abs(E_(c,bold(i)+1/2 bold(e)_c)^("num",n-1/2) - E_(c,bold(i)+1/2 bold(e)_c)^("ref",n-1/2))^2)
    / (sum_(bold(i) in sans(M),c) abs(E_(c,bold(i)+1/2 bold(e)_c)^("ref",n-1/2))^2)
  ).
$

*Purpose:* measure the global discrepancy from the exact electric field for
the refinement study. For E01, the reference has only an $x$ component, but
the numerator also includes any unwanted numerical $y$ and $z$ components.
An error of $0.01$ means that the global L2 error norm is one percent of the reference global L2 norm. Here, _global_ means the whole computational domain $sans(M)$. So this error does not tell us the discrepancies at each location; instead it is just a measure of the total summed error.

The initial observation at $n=0$ therefore uses $t_(-1/2)=-h_t/2$
(see `LevelContract::halfStepLabel`). This norm combines
errors in pulse position, amplitude and shape, and sums all local contributions. Thus it does not distinguish those causes or specify the largest local error.

=== E01.meas02: Relative magnetic-field error in the global discrete $L^2$ norm
The magnetic component $c$ is shifted by half a cell in the other two
directions. Define its index offset
$
  bold(delta)_c^B := 1/2 sum_(d != c) bold(e)_d,
$
where $d in {x,y,z}$. The sample
$B_(c,bold(i)+bold(delta)_c^B)^("num",n-1/2)$ is located at
$vb(x)_bold(i) + 1/2 sum_(d != c) h_d bold(e)_d$.
For example, $bold(delta)_y^B = 1/2 (bold(e)_x+bold(e)_z)$.

We use the `R2` magnetic output from `solver.magneticField()`,
which averages the staggered curls of the current and previous vector potentials
(#link("https://arxiv.org/pdf/2009.13645v1")[MITHRA 2.0, Eqs. (3.50)–(3.52)],
with the spatial curl of Ryu et al., Appendix A, Eq. (24)).
Its time label is $n-1/2$. The reference $B_(c,bold(i)+bold(delta)_c^B)^("ref",n-1/2)$ is the continuum field at the same component position and $t_(n-1/2) = (n-1/2)h_t$.

Using the same global sums and reference normalization as E01.meas01, define
$
  e_(B,thin"L2_global")^(n-1/2) = sqrt(
    (sum_(bold(i) in sans(M),c) abs(B_(c,bold(i)+bold(delta)_c^B)^("num",n-1/2) - B_(c,bold(i)+bold(delta)_c^B)^("ref",n-1/2))^2)
    / (sum_(bold(i) in sans(M),c) abs(B_(c,bold(i)+bold(delta)_c^B)^("ref",n-1/2))^2)
  ).
$
The discrete norm is as in
#link("https://www.osti.gov/servlets/purl/759450")[Salari and Knupp, §4.1,
Eq. (15)]. In code, each owned cell must contribute only once across all MPI ranks, and ghost cells must be excluded. E01's reference has only a $y$ component, while numerical
$B_x$ and $B_z$ also contribute to the error.

*Purpose:* measure the magnetic-field discrepancy, including errors from
potential evolution, the spatial curl and temporal averaging. The magnetic
reconstruction differs from the electric-field reconstruction, so
electric-field agreement alone does not check this output. At $n=0$, this
measurement includes the initial curl and averaging errors. The raw magnetic
output from `computeRawMagneticField()` has integer time index $n$ and is a
separate output.

=== E01.meas03: Relative electric-field error in the global discrete $L^infinity$ norm
Using the samples, staggered positions and half-step time of E01.meas01,
define
$
  e_(E,thin"Linf_global")^(n-1/2) =
    (max_(bold(i) in sans(M),c) abs(E_(c,bold(i)+1/2 bold(e)_c)^("num",n-1/2) - E_(c,bold(i)+1/2 bold(e)_c)^("ref",n-1/2)))
    / (max_(bold(i) in sans(M),c) abs(E_(c,bold(i)+1/2 bold(e)_c)^("ref",n-1/2))).
$
This is the maximum-error norm of
#link("https://www.osti.gov/servlets/purl/759450")[Salari and Knupp, §4.1,
Eq. (18)], normalized by the largest sampled reference magnitude. The maxima may occur at different locations.

*Purpose:* measure the largest componentwise field discrepancy anywhere on
the grid. An error confined to a small part of the domain can contribute
little to E01.meas01 while remaining large at an individual sample.

A value of $0.01$ means that the largest absolute component error is 1% of the global sampled reference peak (and not 1% of the local field components value).

=== E01.meas04: Relative vector-potential error in the global discrete $L^2$ norm
Define
$
  e_(A,thin"L2_global")^n = sqrt(
    (sum_(bold(i) in sans(M),c) abs(A_(c,bold(i)+1/2 bold(e)_c)^("num",n) - A_(c,bold(i)+1/2 bold(e)_c)^("ref",n))^2)
    / (sum_(bold(i) in sans(M),c) abs(A_(c,bold(i)+1/2 bold(e)_c)^("ref",n))^2)
  ).
$

*Purpose:* check the four potential advanced by O1-VS directly. This is gauge-dependent.
It might help interpret the field errors of $vb(E)$ and $vb(B)$.
