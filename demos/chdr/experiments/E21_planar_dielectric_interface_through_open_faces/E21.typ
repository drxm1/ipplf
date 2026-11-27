#import "../../docs/manual/preamble.typ": *

#let e03 = link(<experiment-e03>)[E03]
#let e04 = link(<experiment-e04>)[E04]

// The numbers of this specification, in units of w and w/c_vac.
#let levels = (8, 12, 18) // cells per w
#let small-extent = (10, 8, 8)
#let margin = 6
#let reference-extent = small-extent.map(length => length + 2 * margin)
#let reference-x = (-4 - margin, 6 + margin) // the reference box along x
#let reference-half-width = 4 + margin // and along y and z
#let interface-x = 3.5
#let half-width = 4
#let comparison-time = 6.5
#let late-time = 40
#let safety-factor = 0.99
#let courant = safety-factor / calc.sqrt(3) // h_t / h
// The first step whose R2 fields have a time of at least `time`, on the level with n cells per w.
#let first-step-at(time, n) = calc.ceil(time * n / courant + 0.5)
#let rounded(value, digits: 2) = str(calc.round(value, digits: digits))
#let cell-counts(extent, n) = extent.map(length => str(length * n)).join($times$)
#let million-cells(extent, n) = rounded(extent.fold(1, (product, length) => product * length * n) / 1e6, digits: 1)
#let junction-nearest = calc.sqrt(calc.pow(interface-x, 2) + calc.pow(half-width, 2))
#let junction-farthest = calc.sqrt(calc.pow(interface-x, 2) + 2 * calc.pow(half-width, 2))

// Sketches of the two boxes in the plane z = 0, to scale. One canvas unit is w.
#import "@preview/cetz:0.3.4": canvas, draw
#let mur-red = rgb(220, 0, 90)
#let dielectric-grey = luma(88%)
#let sketch-label(body) = text(size: 7.5pt, body)
// The pulse at t = 0: the support w < r < 3w around the origin.
#let pulse-ring(stroke-width) = {
  draw.circle((0, 0), radius: 3, fill: luma(70%), stroke: stroke-width + black)
  draw.circle((0, 0), radius: 1, fill: white, stroke: stroke-width + black)
}

#let small-box-sketch = canvas(length: 6.5mm, {
  import draw: *
  rect((interface-x, -half-width), (6, half-width), fill: dielectric-grey, stroke: none)
  pulse-ring(0.4pt)
  line((0, 0), (0.36, 0.72), stroke: 0.6pt + black, mark: (end: ">", fill: black))
  content((-0.12, 0.5), anchor: "east", sketch-label[$vb(d)$])
  line((interface-x, -half-width), (interface-x, half-width), stroke: 0.9pt + black)
  rect((-4, -half-width), (6, half-width), stroke: 1.6pt + mur-red)
  for y in (-half-width, half-width) {
    circle((interface-x, y), radius: 0.14, fill: black, stroke: none)
  }
  content((interface-x - 0.25, half-width - 0.45), anchor: "east", sketch-label[junction line])
  content((-3.75, 3.55), anchor: "west", sketch-label[vacuum])
  content((4.75, 3.55), sketch-label[$epsilon_r = 2.13$])
  content((0, -3.5), sketch-label[pulse at $t = 0$])
  content((-3.75, -3.5), anchor: "west", sketch-label[$Omega_"S"$])
  content((1, half-width + 0.5), text(size: 7.5pt, fill: mur-red)[open faces (Mur rows)])
  for (x, label) in ((-4, $-4w$), (0, $0$), (interface-x, $x_f$), (6, $6w$)) {
    content((x, -half-width - 0.5), sketch-label(label))
  }
  for (y, label) in ((-half-width, $-4w$), (half-width, $4w$)) {
    content((-4.3, y), anchor: "east", sketch-label(label))
  }
})

#let reference-box-sketch = canvas(length: 3.2mm, {
  import draw: *
  let (low, high) = reference-x
  let half = reference-half-width
  rect((interface-x, -half), (high, half), fill: dielectric-grey, stroke: none)
  pulse-ring(0.3pt)
  line((interface-x, -half), (interface-x, half), stroke: 0.6pt + black)
  rect((low, -half), (high, half), stroke: 1.2pt + mur-red)
  rect((-4, -half-width), (6, half-width), stroke: (paint: black, thickness: 0.7pt, dash: "dashed"))
  line((6, 0), (high, 0), stroke: 0.5pt + black, mark: (start: ">", end: ">", fill: black))
  content((6 + margin / 2, 1), sketch-label[$m = #margin w$])
  line((-2.5, half-width), (-2.5, half), stroke: 0.5pt + black, mark: (start: ">", end: ">", fill: black))
  content((-1.9, half-width + margin / 2), anchor: "west", sketch-label[$m$])
  content((-3.5, -3), anchor: "west", sketch-label[$Omega_"S"$])
  content((low + 0.7, half - 1.3), anchor: "west", sketch-label[$Omega_"R"$])
  for (x, label) in ((low, $-#calc.abs(low) w$), (interface-x, $x_f$), (high, $#high w$)) {
    content((x, -half - 1.1), sketch-label(label))
  }
  for (y, label) in ((-half, $-#half w$), (half, $#half w$)) {
    content((low - 0.6, y), anchor: "east", sketch-label(label))
  }
})

= Experiment E21 detailed specs <experiment-e21>

E21 measures the error of the open faces of the box when a planar dielectric interface crosses them.
A dielectric half-space leaves the box through four of its open faces.
Along each of these faces the permittivity jumps at one line. We call these four lines the junction lines.
The dipole pulse of #e03 starts in the vacuum part and leaves through all six faces.
The reference is the same run in a larger box, whose open faces cannot influence the compared samples.
The difference of the two runs is then the effect of the open faces alone.

== Experiment E21 cases
The three cases differ in the dielectric and in the speed that the Mur rows use.

#table(
  columns: (auto, 1fr, 2.2fr),
  inset: 6pt,
  align: left + top,
  stroke: 0.5pt,
  table.header([*Case*], [*Dielectric*], [*Speed in the Mur rows*]),
  [E21], [$epsilon_r = 2.13$ for $x >= x_f$], [`MurSpeed::local()`: each row uses $1\/sqrt(epsilon_r)$ of the cell that it writes],
  [E21-M], [as E21], [`MurSpeed::matched()`: every row uses the vacuum speed],
  [E21-V], [none, $epsilon_r = 1$ everywhere], [the vacuum speed],
)

E21 is the main case.
E21-M shows what the speed of each cell changes.
E21-V gives the error of the same open faces without a dielectric.

All six faces of the box have the second-order Mur rows of #e04, with their edge and corner rows.
For $alpha = 1$ all four potentials have the speed $1\/sqrt(epsilon_r)$ in a uniform medium @chew2014generalizedgauge[§2, text after Eq. (22)].
`MurSpeed::local()` assumes that the material of each boundary cell continues outward.
At a junction line the speed changes from one boundary cell to the next, and the rows contain no term for this change.

== Experiment E21 computational domain and refinement
Lengths are in units of the half-width $w$ of the pulse, and times in units of $w\/c_"vac"$.
@fig-e21-setup shows the two boxes.

#figure(
  grid(
    columns: 2,
    column-gutter: 9mm,
    align: horizon,
    small-box-sketch, reference-box-sketch,
  ),
  caption: [The plane $z = 0$ of E21, to scale.
    Left: the small box $Omega_"S"$ with the support of the pulse at $t = 0$ and the direction $vb(d)$ of the dipole, projected into the plane.
    Right: the reference box $Omega_"R"$, which is $Omega_"S"$ (dashed) with the margin $m$ on every side.
    The dielectric is shaded and continues to the open faces.
    The two dots are junction lines, which are perpendicular to this plane.],
) <fig-e21-setup>

#table(
  columns: (1fr, 3.4fr),
  inset: 6pt,
  align: left + top,
  stroke: 0.5pt,
  [Small box], [$Omega_"S" = [-4w, 6w] times [-4w, 4w]^2$, all six faces open],
  [Reference box],
  [$Omega_"R" = [-#calc.abs(reference-x.at(0)) w, #reference-x.at(1) w] times [-#reference-half-width w, #reference-half-width w]^2$,
    which is $Omega_"S"$ with the margin $m = #margin w$ on every side.
    All six faces are open and use the Mur speed of the case.],
  [Interface], [The plane $x = x_f = #interface-x w$. In both boxes the dielectric fills $x >= x_f$ up to the open faces.],
  [Pulse],
  [The dipole pulse of E03 (@e03-dipole-potentials) with $beta = 8$ (@e03-profile-sharpness), $R_0 = 2w$,
    the center $vb(0)$ and the direction $vb(d)$ of E03, with the vacuum speed $v = 1$],
  [Solver], [As in E03: O1VS with $alpha = 1$, Standard stencil, staggered placement, R2. Cubic cells.],
  [Time step],
  [$h_t = #safety-factor thin h_t^"max"$ with the stability bound $h_t^"max" = h\/sqrt(3)$.
    The program selects it with `ChdrSolver::setSafetyFactor(0.99)`, as in #link(<experiment-e16s>)[E16-S].],
  [Compared interval], [$T_"c" = #comparison-time w$, up to the step $n_"c" = ceil(T_"c"\/h_t + 1\/2)$],
  [Late interval], [$T_"late" = #late-time w$, up to the step $n_"late" = ceil(T_"late"\/h_t + 1\/2)$, in the small box only],
)
$n_"c"$ and $n_"late"$ are the first steps whose R2 fields have a time of at least $T_"c"$ and $T_"late"$.

We use the time step that a production run would use.
For the Standard stencil the largest stable step has the smallest phase error in every direction of propagation,
in the vacuum and in the dielectric
(deduction from the dispersion relation @fallahi2020mithra20fullwavesimulation[§3.1, Eq. (3.17)] and its stability bound, Eq. (3.20)).
E03 and E04 used $h_t = h\/2$. No experiment with open faces has used the step of E21 before.

#align(center)[
  #table(
    columns: 6,
    inset: 5pt,
    align: center,
    stroke: 0.4pt,
    table.header([$w\/h$], [$h_t\/w$], [Cells of $Omega_"S"$], [Cells of $Omega_"R"$], [$n_"c"$], [$n_"late"$]),
    ..levels
      .map(n => (
        [#n],
        [#rounded(courant / n, digits: 4)],
        [#cell-counts(small-extent, n)],
        [#cell-counts(reference-extent, n) (#million-cells(reference-extent, n) million)],
        [#first-step-at(comparison-time, n)],
        [#first-step-at(late-time, n)],
      ))
      .flatten(),
  )
]
On every level the interface and all faces of both boxes lie on cell faces.
The program holds both boxes.
For every cell of $Omega_"S"$ it subtracts the field value of the reference box from the field value of the small box at that cell.
With several MPI ranks each box is divided among the processes in its own way,
so the two values of one cell can lie in the memory of two different processes.
The runs therefore use one MPI rank, with OpenMP threads.
Run the two coarse levels first, and measure the memory and the time of a step before the third level.

The pulse occupies $(1 + t)w < r < (3 + t)w$ in vacuum.
The following times are deductions from this support and the geometry, for the direct path in vacuum.
The reflected and the refracted parts arrive later.

#table(
  columns: (3fr, 1fr),
  inset: 5pt,
  stroke: 0.4pt,
  table.header([Event], [$t\/w$]),
  [The front reaches the interface], [#(interface-x - 3)],
  [The front reaches the nearest open faces], [#(half-width - 3)],
  [The front reaches the junction lines], [#rounded(junction-nearest - 3)],
  [The rear has passed every point of the junction lines], [#rounded(junction-farthest - 1)],
  [End of the compared interval], [#comparison-time],
)

== Experiment E21 initialization
As in E03 (@e03-initialization): we sample the analytic potentials of the pulse at the native positions and times,
in both boxes at the same positions.
At these times the support lies inside $r < 3w + h_t\/2$, which is in the vacuum and $0.5w$ away from the interface.
The pulse of a uniform vacuum is therefore the correct start for all three cases.

== The reference box
Every update and the reconstruction use only neighboring samples (@eq-O1VS-ThAd, @sec-time-placement).
An effect therefore travels at most one cell per step in every direction.
At an open face a Mur row reads the two outermost cell layers,
and the extrapolation of the ghost samples (`OpenFaceClosure.h`) reads the three outermost.
The initial histories are exactly zero outside the support of the profile.
Let $n_"b"$ be the first step at which a potential sample in the three outermost cell layers of $Omega_"R"$ is nonzero.
Then no sample of $Omega_"S"$ depends on the open faces of $Omega_"R"$ at a step
$
  n <= n_"b" + m\/h - 2
$ <E21-independence>
(deduction).
With one cell per step, $n_"b" approx #(1 + margin) w\/h$, so that the bound is about $#(1 + 2 * margin) w\/h$ steps
against $n_"c" approx #rounded(comparison-time / courant, digits: 1) w\/h$.
One cell per step is #rounded(1 / courant) times the vacuum speed, which is why the margin is this wide.
The program records $n_"b"$, and the analysis checks @E21-independence for $n = n_"c"$.

Up to $n_"c"$ the reference is then the solution of the same update on an unbounded mesh.
Both boxes have the same update, the same material map and the same initial samples.
Their difference therefore contains no error of the interior update and no error of the interface.
It is not an error against the continuum solution.

== Experiment E21 measured quantities and their purpose
Let $sans(M)$ be the cells of $Omega_"S"$, $sans(M)_1$ those with $x < x_f$ and $sans(M)_2$ those with $x > x_f$.
E21-V uses the same two sets.
For a cell $bold(i) in sans(M)$ let $bold(i)'$ be the cell of $Omega_"R"$ at the same position, and
$
  Delta E_(c,bold(i)+bold(delta)_c^E)^(n-1/2) = E_(c,bold(i)+bold(delta)_c^E)^("S",n-1/2) - E_(c,bold(i)'+bold(delta)_c^E)^("R",n-1/2),
$
with $Delta B$ defined in the same way from the R2 magnetic samples.
The superscripts S and R denote the small box and the reference box.
The permittivity $epsilon_(c,bold(i)+bold(delta)_c^E)$ of a component and the field energy $U_S$ on a side are those of E19.meas02.
$U_0$ is the sum of $U_1$ and $U_2$ of the reference over $sans(M)$ at step zero.

=== E21.meas01: Energy of the difference field on each side
For $S in {1, 2}$ and $0 <= n <= n_"c"$ define
$
  D_S^(n-1/2) = h^3/2 sum_(bold(i) in sans(M)_S,c)
    (epsilon_(c,bold(i)+bold(delta)_c^E) abs(Delta E_(c,bold(i)+bold(delta)_c^E)^(n-1/2))^2
     + abs(Delta B_(c,bold(i)+bold(delta)_c^B)^(n-1/2))^2), quad
  F_(Delta,S)^(n-1/2) = D_S^(n-1/2) / U_0,
$
and $F_Delta = F_(Delta,1) + F_(Delta,2)$.

*Purpose:* this is the main result.
$F_Delta$ is the energy of the field that the open faces of $Omega_"S"$ add, as a fraction of the energy of the pulse.
It is zero for an ideal open face.
It does not vanish under refinement, because the second-order rows reflect oblique waves also in the continuum (#e04, analysis).
The difference between E21 and E21-V is the effect of the dielectric at the open faces, which includes the junction lines.
The difference between E21-M and E21 is the effect of the speed of each cell.
We expect it on side 2.

=== E21.meas02: Map of the largest difference energy density in the plane $z = 0$
Let $k_0$ be the layer of cells whose centers have $z = h\/2$. For a cell $bold(i) = (i, j, k_0)$ define
$
  u_(Delta,bold(i))^(n-1/2) = 1/2 sum_c
    (epsilon_(c,bold(i)+bold(delta)_c^E) abs(Delta E_(c,bold(i)+bold(delta)_c^E)^(n-1/2))^2
     + abs(Delta B_(c,bold(i)+bold(delta)_c^B)^(n-1/2))^2), quad
  mu_(i j) = max_(1 <= n <= n_"c") u_(Delta,(i,j,k_0))^(n-1/2).
$
We also save $u_0^"max"$, the largest value of the same density of the reference field over $sans(M)$ at step zero.

*Purpose:* the map shows where the added field starts and which way it travels.
The plane cuts two junction lines, at $(x_f, plus.minus #half-width w)$.
A junction line that radiates appears as a pattern that spreads from these two points.

=== E21.meas03: Energy that remains in the small box
For $0 <= n <= n_"late"$ define
$
  F_"rem"^(n-1/2) = (U_1^(n-1/2) + U_2^(n-1/2)) / U_0
$
with the energies of the small box.
On the coarsest level, E21 continues to $10^5$ steps.

*Purpose:* after the pulse has left, $F_"rem"$ is the energy that the open faces returned or that did not leave.
An increase at late times would show that the Mur rows are unstable, for example at the junction lines.
A run of this length shows no long-time stability beyond its last step.

== Experiment E21 saved output
The inputs are `E21-n8.yaml` to `E21-V-n18.yaml`. Each run writes the following files.

#table(
  columns: (1fr, 4fr),
  inset: 6pt,
  align: left + top,
  stroke: 0.5pt,
  [`measurements.csv`],
  [
    // 2026-11-27: the column names of both files are stated, because the report scripts read them.
    One row per quantity and step, with the columns `quantity`, `step`, `time` and `value`.
    The time is that of the R2 fields.
    Up to $n_"c"$: `D_side1` and `D_side2` (E21.meas01),
    and `U_reference_side1` and `U_reference_side2`, the energies of the reference over $sans(M)_1$ and $sans(M)_2$.
    Up to $n_"late"$: `U_side1` and `U_side2` of the small box (E21.meas03).
    In the long run these two are written at every 16th step after $n_"late"$.
  ],
  [`difference_map.csv`],
  [One row per cell of the layer $k_0$ of $Omega_"S"$, with the columns `i`, `j`, `x`, `y` and `value`:
    its indices, the position of its center and $mu_(i j)$ (E21.meas02).],
  [`run.json`],
  [
    The configuration of the run.
    Under `effective.reference_box` it holds `margin_cells` ($m\/h$), `comparison_steps` ($n_"c"$) and `first_nonzero_boundary_step` ($n_"b"$).
    Under `effective.normalization` it holds `initial_energy` ($U_0$) and `peak_initial_energy_density` ($u_0^"max"$).
  ],
)

== Experiment E21 analysis
We check @E21-independence on every run.
For each case and level we report the largest $F_Delta$ of the compared interval, and $F_(Delta,1)$ and $F_(Delta,2)$ at $n_"c"$.
We plot $F_Delta$ against time for the three cases, with the events of the table above.
We show the maps $mu_(i j)\/u_0^"max"$ of the three cases on the finest level.
We report $F_"rem"$ at $20w$ and at $40w$, and for the long run its largest value after $40w$.

== What E21 does not test
The interface has no edge and no corner inside the box, and it lies on cell faces.
The accuracy of the interface itself is the subject of E16 to E19.
E21 has one pulse, one box and one position of the interface, so it gives no reflection coefficient of the open faces as a function of the angle.
