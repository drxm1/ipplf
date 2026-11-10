#import "../../docs/manual/preamble.typ": *

#let chew = [Chew @chew2014generalizedgauge[§2, Eqs. (5)–(13)]]
#let ryu = [Ryu et al. @ryu2016potentialfdtd[Appendix A, Eqs. (19)–(27)]]
#let fallahi = [Fallahi @fallahi2020mithra20fullwavesimulation[§3.1]]
#let salari = [Salari and Knupp @salari2000CodeVerificationMms[§§4.1–4.2]]
#let christlieb = [Christlieb et al. @christlieb2024GaugeConserving[§3.1, Eqs. (25)–(27)]]

= Experiment E03 detailed specification\* <experiment-e03>

E03 measures convergence of a spherical pulse over a radial travel distance $2w$,
before continuum contact with the periodic boundaries.
We refine all three cell widths and the time step together.
The source-free pulse is the field of an electric dipole.
Its scalar and vector potentials are both nonzero.
This three-dimensional study complements the longer axial propagation in #link(<experiment-e01>)[E01]
and supplies the propagation baseline for #link(<experiment-e04>)[E04].

All quantities use the manual's normalized solver units:
$c_"vac"=epsilon_0=mu_0=1$.
The lossless, nondispersive medium fills the domain,
with $epsilon_r=2.13$, $mu_r=1$, no free charge and no current.
Write $epsilon=epsilon_r$ and $v=1/sqrt(epsilon)$.

== Analytic continuum reference

=== A smooth outgoing scalar pulse

Let $vb(x)_0=vb(0)$ be the centre,
$r=norm(vb(x)-vb(x)_0)$ and $vb(n)=(vb(x)-vb(x)_0)/r$ for $r>0$.
Use the same compact profile as E01,
with potential scale $A_0=1$, half-width $w=1$, and initial radius $R_0=2w$:
$
  s=r-R_0-v t, quad
  p(s)=cases(
    A_0 exp(1-1/(1-(s/w)^2)) & abs(s)<w,
    0 & abs(s)>=w
  ).
$ <e03-profile>
Define the scalar generator
$
  f(s)=w^2 p(s), quad psi(r,t)=f(r-R_0-v t)/r.
$ <e03-generator>
The factor $w^2$ makes $psi$ have the dimension of vector potential times length.

We define $psi$ and all its derivatives to be zero outside $abs(s)<w$,
including at $r=0$.
This is smooth for $t>-(R_0-w)/v$,
because an open neighbourhood of the centre then lies outside the support.
Every time used in this experiment, including the negative initialization times, satisfies this restriction.
There is no point source at the centre during the simulated interval.

The construction is a deduction from the homogeneous wave equations in #chew.
For a radial function,
$
  nabla^2 psi=partial_r^2 psi+(2/r)partial_r psi=frac(f'' (s),r),
  quad partial_t^2 psi=v^2 frac(f'' (s),r).
$
Thus $nabla^2 psi-epsilon partial_t^2 psi=0$.
The zero continuation across the support boundaries also solves this equation,
since all derivatives of the compact profile vanish there.

For an implementation without numerical differentiation,
set $u=s/w$ and $q=1-u^2$ inside the support.
Direct differentiation of @e03-profile gives
$
  p'(s)=-2u p(s)/(w q^2), quad
  p''(s)=p(s)/w^2 (4u^2/q^4-2/q^2-8u^2/q^3).
$ <e03-profile-derivatives>
Multiply these expressions by $w^2$ to obtain $f'$ and $f''$.
Test the support condition before evaluating $1/r$, $vb(n)$ or any inverse power of $q$.

*Profile sharpness for the refined study.*
The first study used the profile above and stayed far outside the asymptotic range on all four meshes:
its flanks are too steep for at most 54 cells per pulse width.
The refined study uses the family
$
  p_beta (s)=A_0 exp(beta(1-1/(1-(s\/w)^2))), quad abs(s)<w,
$ <e03-profile-sharpness>
with $p_beta=0$ outside the support and the input parameter `profile_sharpness` $=beta>0$.
The value $beta=1$ is the profile above, and the refined study uses $beta=8$.
Since $beta(1-1\/q)=-beta u^2\/q$, the pulse is close to the Gaussian $A_0 exp(-beta s^2\/w^2)$ near its centre.
Support and smoothness do not change,
so every statement of this section holds for all $beta>0$.
The derivatives become
$
  p_beta '(s)=-2 beta u thin p_beta\/(w q^2), quad
  p_beta ''(s)=p_beta\/w^2 thin (4 beta^2 u^2\/q^4-2 beta\/q^2-8 beta u^2\/q^3).
$
The geometry, the meshes, the time step and all measurements stay as specified.

#block(breakable: false)[
=== The dipole pulse

Choose one fixed unit vector
$
  vb(d)=frac(lr((1,2,3)),sqrt(14)).
$
]
This dipole direction exercises all Cartesian components.
The pulse has spherical support and a dipole angular pattern, not a spherically symmetric field amplitude.
Introduce the following scalar derivatives of @e03-generator:
$
  a=psi_r=f'/r-f/r^2, quad
  b=psi_(r r)=f''/r-2f'/r^2+2f/r^3, quad
  c=f''/r-f'/r^2.
$ <e03-radial-derivatives>
Here and below $f,f',f''$ are evaluated at $s$.
These expressions are used only inside the support.

Define the potentials
$
  vb(A)=vb(d) partial_t psi\/v=-vb(d)f'\/r,
  quad phi=-v vb(d) dot grad psi=-v a(vb(n) dot vb(d)).
$ <e03-dipole-potentials>
Using $vb(E)=-partial_t vb(A)-grad phi$ and $vb(B)=curl vb(A)$ from #chew,
direct differentiation and $psi_(t t)=v^2 nabla^2 psi$ give
$
  vb(E)=v[(b-a\/r)vb(n)(vb(n) dot vb(d))-(b+a\/r)vb(d)],
  quad vb(B)=-c(vb(n) times vb(d)).
$ <e03-dipole-fields>
This is the pulse of an electric dipole along $vb(d)$, with the Hertz potential $vb(d) psi$.
The magnetic field has no radial component, $vb(n) dot vb(B)=0$,
while the electric field generally has one.
In the multipole terminology this is a TM\_r field, where the subscript refers to the radial direction and not to an interface.
@fig-e03-dipole-pulse shows the pulse at $t=0$.

#import "@preview/cetz:0.3.4": canvas, draw

// Fields of the dipole pulse in a plane that contains the dipole axis.
// The point is rho e + zeta d, with a unit vector e normal to d. Lengths are in units of w,
// and `centre` is the radius R_0 + v t of the pulse centre. The result is
// (E dot e, E dot d, B dot (e times d)), with E in units of v A_0/w and B in units of A_0/w.
#let dipole-pulse-fields(rho, zeta, centre: 2.0, sharpness: 8.0) = {
  let r = calc.sqrt(rho * rho + zeta * zeta)
  let u = r - centre
  if calc.abs(u) >= 1 { return (0.0, 0.0, 0.0) }
  let q = 1 - u * u
  let f = calc.exp(sharpness * (1 - 1 / q))
  let slope = -2 * sharpness * u / (q * q)
  let f1 = slope * f
  let f2 = (slope * slope - 2 * sharpness / (q * q) - 8 * sharpness * u * u / (q * q * q)) * f
  let (a, b, c) = ((f1 - f / r) / r, (f2 - 2 * f1 / r + 2 * f / (r * r)) / r, (f2 - f1 / r) / r)
  let (sine, cosine) = (rho / r, zeta / r)
  ((b - a / r) * sine * cosine, (b - a / r) * cosine * cosine - (b + a / r), -c * sine)
}

// Colours of the two maps: white to blue for a magnitude in [0, 1],
// orange through white to blue for a signed value in [-1, 1].
#let magnitude-colour(t) = (255 - 227 * t, 255 - 163 * t, 255 - 84 * t).map(x => int(calc.round(x)))
#let signed-colour(t) = if t >= 0 {
  (255 - 213 * t, 255 - 135 * t, 255 - 41 * t).map(x => int(calc.round(x)))
} else {
  (255 + 20 * t, 255 + 151 * t, 255 + 203 * t).map(x => int(calc.round(x)))
}

// Samples (|E|, B normal to the plane) on a square of half-width 3.4 w around the centre.
#let pulse-half-extent = 3.4
#let pulse-unit = 6.4mm
#let pulse-pixels = 120
#let pulse-samples = range(pulse-pixels * pulse-pixels).map(k => {
  let step = 2 * pulse-half-extent / pulse-pixels
  let rho = -pulse-half-extent + (calc.rem(k, pulse-pixels) + 0.5) * step
  let zeta = pulse-half-extent - (calc.quo(k, pulse-pixels) + 0.5) * step
  let (e-rho, e-zeta, b-out) = dipole-pulse-fields(rho, zeta)
  (calc.sqrt(e-rho * e-rho + e-zeta * e-zeta), b-out)
})

// One map as a raster image, scaled to its largest absolute value.
#let pulse-map(component, colour) = {
  let values = pulse-samples.map(sample => sample.at(component))
  let largest = values.fold(0.0, (m, x) => calc.max(m, calc.abs(x)))
  image(
    bytes(values.map(x => colour(x / largest)).flatten()),
    format: (encoding: "rgb8", width: pulse-pixels, height: pulse-pixels),
    width: 2 * pulse-half-extent * pulse-unit,
  )
}

// Frame, support circles and dipole axis that the two cut panels share.
#let cut-frame(map, title) = {
  let h = pulse-half-extent
  draw.content((0, 0), map)
  draw.rect((-h, -h), (h, h), stroke: 0.35pt + luma(65%))
  for radius in (1, 3) {
    draw.circle((0, 0), radius: radius, stroke: (paint: luma(55%), thickness: 0.35pt, dash: "dashed"))
  }
  draw.line((0, -h), (0, h), stroke: 0.5pt, mark: (end: ">", fill: black, scale: 0.6))
  draw.content((0.18, h - 0.12), anchor: "north-west", [$vb(d)$])
  draw.content((-h, h + 0.16), anchor: "south-west", title)
  for tick in (-3, -2, -1, 0, 1, 2, 3) {
    draw.line((tick, -h), (tick, -h - 0.08), stroke: 0.35pt + luma(65%))
    draw.content((tick, -h - 0.16), anchor: "north", [$#tick$])
  }
  draw.content((0, -h - 0.62), anchor: "north", [distance from the axis $\/ w$])
}

// Arrows of E on three circles: the centre of the pulse and its two side bands.
#let field-arrows() = {
  let largest = pulse-samples.fold(0.0, (m, sample) => calc.max(m, sample.at(0)))
  for radius in (1.57, 2.0, 2.43) {
    for index in range(16) {
      let angle = (index + 0.5) * 22.5deg
      let (rho, zeta) = (radius * calc.sin(angle), radius * calc.cos(angle))
      let (e-rho, e-zeta, _) = dipole-pulse-fields(rho, zeta)
      let (dx, dy) = (0.62 * e-rho / largest, 0.62 * e-zeta / largest)
      if calc.sqrt(dx * dx + dy * dy) > 0.12 {
        draw.line((rho - dx / 2, zeta - dy / 2), (rho + dx / 2, zeta + dy / 2),
          stroke: 0.55pt, mark: (end: ">", fill: black, scale: 0.4))
      }
    }
  }
}

// Orthographic view of the point p = (x, y, z) for the sketch of the geometry.
// The camera stands at the azimuth -35 degrees and the elevation 30 degrees.
#let sketch-camera = (azimuth: -35deg, elevation: 30deg)
#let scalar-product(p, q) = p.zip(q).map(((a, b)) => a * b).sum()
#let sketch-view(p) = {
  let (az, el) = (sketch-camera.azimuth, sketch-camera.elevation)
  (scalar-product(p, (-calc.sin(az), calc.cos(az), 0)),
    scalar-product(p, (-calc.sin(el) * calc.cos(az), -calc.sin(el) * calc.sin(az), calc.cos(el))))
}
#let towards-camera(p) = {
  let (az, el) = (sketch-camera.azimuth, sketch-camera.elevation)
  scalar-product(p, (calc.cos(el) * calc.cos(az), calc.cos(el) * calc.sin(az), calc.sin(el)))
}

// The dipole axis d and two unit vectors normal to it: e spans the cut plane with d,
// and g = d times e spans the equatorial plane with e.
#let sketch-d = (1, 2, 3).map(x => x / calc.sqrt(14))
#let sketch-e = (2, -1, 0).map(x => x / calc.sqrt(5))
#let sketch-g = (3, 6, -5).map(x => x / calc.sqrt(70))
#let sketch-ring(u, v, radius) = range(73).map(i => range(3).map(k =>
  radius * (calc.cos(i * 5deg) * u.at(k) + calc.sin(i * 5deg) * v.at(k))))

// The axis d, the sphere r = 2 w with its equator, and the cut plane of the two maps.
#let geometry-sketch(title) = {
  let h = pulse-half-extent
  let blue = rgb("2a78d6")
  draw.rect((-h, -h), (h, h), stroke: 0.35pt + luma(65%))
  draw.content((-h, h + 0.16), anchor: "south-west", title)
  draw.line(..sketch-ring(sketch-e, sketch-d, 3.3).map(sketch-view), close: true,
    fill: luma(90%), stroke: 0.35pt + luma(60%))
  for radius in (1, 3) {
    draw.circle((0, 0), radius: radius, stroke: (paint: luma(55%), thickness: 0.35pt, dash: "dashed"))
  }
  draw.circle((0, 0), radius: 2, fill: blue.transparentize(80%), stroke: 0.5pt + blue)
  let equator = sketch-ring(sketch-e, sketch-g, 2)
  for (start, end) in equator.slice(0, -1).zip(equator.slice(1)) {
    let front = towards-camera(start) + towards-camera(end) >= 0
    draw.line(sketch-view(start), sketch-view(end),
      stroke: if front { 0.8pt + blue } else { 0.4pt + blue.lighten(45%) })
  }
  for (axis, name) in (((1, 0, 0), $x$), ((0, 1, 0), $y$), ((0, 0, 1), $z$)) {
    let tip = sketch-view(axis.map(x => 2.9 * x))
    draw.line((0, 0), tip, stroke: 0.4pt + luma(40%), mark: (end: ">", fill: luma(40%), scale: 0.45))
    draw.content(sketch-view(axis.map(x => 3.2 * x)), text(fill: luma(30%), name))
  }
  draw.line(sketch-view(sketch-d.map(x => -3.2 * x)), sketch-view(sketch-d.map(x => 3.2 * x)),
    stroke: 0.6pt, mark: (end: ">", fill: black, scale: 0.6))
  draw.content(sketch-view(sketch-d.map(x => 3.2 * x)), anchor: "west", padding: 0.12, [$vb(d)$])
}

#figure(
  grid(columns: 3, column-gutter: 8pt, align: center + top,
    ..(
      { geometry-sketch([Axis $vb(d)$ and cut plane]) },
      {
        cut-frame(pulse-map(0, magnitude-colour), [$abs(vb(E))$ and direction of $vb(E)$])
        field-arrows()
      },
      { cut-frame(pulse-map(1, signed-colour), [$vb(B)$, normal to the plane]) },
    ).map(panel => {
      set text(size: 9pt)
      canvas(length: pulse-unit, panel)
    }),
  ),
  caption: [The dipole pulse at $t=0$ for $beta=8$.
    Its support is $w<r<3w$ (dashed circles), and it is rotationally symmetric about $vb(d)$.
    Left: the axis $vb(d)$, the sphere $r=2w$ with its equator, and a plane that contains $vb(d)$ (grey).
    Centre: $abs(vb(E))$ in this plane, from zero (white) to its largest value (blue).
    The arrows show the direction of $vb(E)$, and their length is proportional to $abs(vb(E))$.
    Right: $vb(B)$ is normal to the plane.
    It points out of the page where the map is blue and into the page where it is orange.],
) <fig-e03-dipole-pulse>

*Why one pulse is enough.*
The first version of this experiment also ran the dual pulse
$vb(A)=grad psi times vb(d)$, $phi=0$ (a magnetic dipole, TE\_r),
whose fields are those above with $vb(E)$ and $vb(B)$ exchanged: $vb(E) -> v vb(B)$, $vb(B) -> -vb(E)\/v$.
That pulse exercises the vector potential only.
The dipole pulse exercises the same update on all three vector components
and in addition the scalar potential, its gradient in $vb(E)$ and the Lorenz condition.
The first study gave errors of $vb(E)$, $vb(B)$ and $vb(A)$ that differ by less than one percent between the two pulses.
The experiment therefore runs the dipole pulse only, and the input has no `polarization` parameter.

=== Why this is a source-free Maxwell solution

The following is a deduction from the potential relations, Lorenz gauge and homogeneous wave equations in #chew.
Spatial and time derivatives commute with the constant-coefficient wave operator,
so every component in @e03-dipole-potentials solves that wave equation.
The Lorenz condition holds because
$
  div vb(A)+epsilon partial_t phi
  =(1/v-epsilon v) vb(d) dot grad partial_t psi=0.
$
The pulse therefore satisfies the homogeneous Lorenz gauge and gives source-free Maxwell fields.
For example,
$
  curl(grad psi times vb(d))
  =(grad grad psi)vb(d)-vb(d)nabla^2 psi,
$
which gives the electric field of @e03-dipole-fields, with $psi_(t t)=v^2 nabla^2 psi$.
No far-field or locally plane-wave approximation is used.
The radial terms proportional to $1/r^2$ and $1/r^3$ must be retained.

== Solver configuration, domain and refinement

#table(
  columns: (1fr, 2fr), inset: 5pt, align: left + top, stroke: 0.4pt,
  table.header([*Setting*], [*Value*]),
  [Formulation], [O1VS, $alpha=1$],
  [Discretization], [Standard stencil, staggered placement, R2 reconstruction],
  [Material and sources], [$epsilon=2.13$, $mu=1$, $rho=0$, $vb(J)=vb(0)$],
  [Boundary conditions], [Periodic in all three directions],
  [Domain], [$Omega=[-6w,6w]^3$],
  [Pulse], [$A_0=w=1$, $R_0=2w$, $vb(d)=frac(lr((1,2,3)),sqrt(14))$],
  [Time step], [$h_t=h/2$, chosen before loading histories],
  [Target duration], [$T_3=2w/v approx 2.919w/c_"vac"$],
)

Initially the support is $w<r<3w$.
At $T_3$, it is $3w<r<5w$, leaving a distance $w$ from the nearest domain face.

The periodic continuum reference is the sum of translated copies of the potentials over the lattice $(12w)ZZ^3$.
During E03's observation interval, only the central copy contributes anywhere in the domain,
because its support radius is below $6w$ and the neighbouring copies have not arrived.
Thus the formulas above are also the exact periodic reference on this interval.
Check this support condition at the actual selected times, including the scalar-potential time.
Periodic numerical tails can still reach a boundary earlier than the continuum support.
The measured field errors include any such contribution.

Let $N=N_x=N_y=N_z$ and $h=h_x=h_y=h_z=12w/N$.
We choose equal cell counts because the physical box is cubic and the pulse varies in every direction.

#align(center)[
  #table(
    columns: 5, inset: 5pt, align: center, stroke: 0.4pt,
    table.header([$N$], [$2w/h$], [$h/w$], [$h_t$], [Owned cells]),
    [96], [16], [$1/8$], [$1/16$], [884,736],
    [144], [24], [$1/12$], [$1/24$], [2,985,984],
    [216], [36], [$1/18$], [$1/36$], [10,077,696],
    [324], [54], [$1/27$], [$1/54$], [34,012,224],
  )
]
The refinement ratio is $3/2$ in every spatial direction and in time.
The physical box, pulse and target duration stay fixed.

The time step satisfies the solver's vacuum Standard bound
$
  h_t=h/2<h/sqrt(3),
$
which is stricter than the homogeneous-dielectric bound here
(#fallahi, Eq. (3.20)).
This combined refinement does not separately identify spatial and temporal orders.

Start with $N=96,144$.
Measure errors, peak memory and time per step before selecting the remaining planned levels.
At $N=324$, one scalar double array occupies about 260 MiB before halos,
and the final field target requires 158 advances, or $5.37 times 10^9$ owned cell advances.
These counts follow from the table and schedule below, not from a runtime measurement.
Retain the geometry and duration when refining.
If the pilot or available resources do not support the remaining levels, report the achieved range and unresolved trend.

== Initialization and sample positions

Sample the exact continuum potential histories directly at their native positions and times.

The R2 time convention in `core/base/LevelContract.h` distinguishes scalar and vector histories.
At completed step $n$, the current vector slot represents $n h_t$,
but the current scalar slot represents $(n+1/2)h_t$.
The initial samples therefore are
$
  vb(A)^0=vb(A)_"ref" (0), quad vb(A)^(-1)=vb(A)_"ref" (-h_t),
$
$
  phi^0=phi_"ref" (h_t/2), quad phi^(-1)=phi_"ref" (-h_t/2).
$ <e03-initialization>
The superscripts here name storage levels, not a common time for all four components.
Loading both scalar and vector slots at integer times would be an error.
The spatial staggering follows #ryu.
The storage times are defined by `LevelContract` and `StaggeredReconstruction.h`,
using the centred temporal formulation of #christlieb after relabelling the history indices.
Obtain every offset and time from the constructed `LevelContract`.

With cell centre $vb(x)_bold(i)$ and Cartesian unit vector $vb(e)_c$,
the O1VS component positions and R2 output times are
#table(
  columns: (0.65fr, 2fr, 1fr), inset: 5pt, stroke: 0.4pt,
  table.header([Quantity], [Position], [Time at step $n$]),
  [$A_c$], [$vb(x)_bold(i)+(h/2)vb(e)_c$], [$n h_t$],
  [$phi$], [$vb(x)_bold(i)$], [$(n+1/2)h_t$],
  [$E_c$], [$vb(x)_bold(i)+(h/2)vb(e)_c$], [$(n-1/2)h_t$],
  [$B_c$], [$vb(x)_bold(i)+(h/2)sum_(d!=c)vb(e)_d$], [$(n-1/2)h_t$],
)
The magnetic output is the average of the curls of the current and previous vector histories,
as returned by `solver.magneticField()`.
The two-curl average follows #fallahi, Eqs. (3.50)–(3.52), with the spatial curl of #ryu.
Every analytic comparison uses the numerical component's own position and time.

Smooth continuum histories do not generally satisfy the discrete Lorenz constraint exactly.
For these curved fields, separately sampling the continuum derivatives does not reproduce a discrete curl identity.
Measure the reconstructed fields and constraints at state zero as startup measurements.
Second-order startup errors are expected from Taylor expansion of the centred staggered differences in #ryu.
The loaded potentials instead agree with the reference sampler by construction.
Check both loaded histories and their time labels on a small mesh, without assigning an order to zero loading errors.

== Measurements and observation times

=== E03.meas01–04: continuum field and potential errors

Use E01.meas01–04 with the complete reference vectors above:
relative global $L^2$ errors for $vb(E)$, R2 $vb(B)$ and $vb(A)$,
and relative global componentwise $L^infinity$ error for $vb(E)$.
For clarity, for a vector quantity $vb(F)$ define
$
  norm(vb(F))_(2,h)^2=h^3 sum_(bold(i) in cal(M),c) abs(F_(c,bold(i)))^2,
  quad norm(vb(F))_(infinity,h)=max_(bold(i) in cal(M),c) abs(F_(c,bold(i))).
$
Each component in this notation is sampled at its own native location.
The relative errors are
$
  e_(F,2)=norm(vb(F)_"num"-vb(F)_"ref")_(2,h)/norm(vb(F)_"ref")_(2,h),
  quad
  e_(E,infinity)=norm(vb(E)_"num"-vb(E)_"ref")_(infinity,h)/norm(vb(E)_"ref")_(infinity,h).
$
Sum or maximize over all owned samples exactly once across MPI ranks, excluding ghosts.
The discrete norms follow #salari, Eqs. (15), (18), with reference normalization defined here.
They combine phase, amplitude and shape errors.
Use the potential error for evolved states and the field errors for startup and evolution.

=== E03.meas05: scalar-potential error

For evolved states, also save the relative scalar $L^2$ error $e_(phi,2)$,
using the scalar version of the preceding norm at #box[$(n+1/2)h_t$].

=== E03.meas06: Gauss and Lorenz-gauge residuals

The homogeneous source-free constraints are $div(epsilon vb(E))=0$ and
$div vb(A)+epsilon partial_t phi=0$ (#chew).
Using the solver's backward staggered divergence $D_h^-$, record
$
  Q_h^(n-1/2)=D_h^- dot (epsilon vb(E)^(n-1/2)),
  quad
  G_h^n=D_h^- dot vb(A)^n+epsilon (phi^n-phi^(n-1))/h_t.
$ <e03-constraints>
The scalar difference in $G_h^n$ is centred at $n h_t$ by @e03-initialization.
Here $Q_h$ is the electric Gauss-law residual and $G_h$ the Lorenz-gauge residual.
Evaluate them on the numerical state with the operators of #ryu.
Synchronize reconstructed E halos before taking its divergence.
`reconstructStaggered` writes owned fields and leaves their halos zero.

For $n>=1$, the homogeneous scalar bulk update and R2 reconstruction imply
$
  D_h^- dot D_h^+ phi^(n-1)
  =epsilon (phi^n-2phi^(n-1)+phi^(n-2))/h_t^2,
  quad Q_h^(n-1/2)=-epsilon (G_h^n-G_h^(n-1))/h_t.
$
The second identity follows by substituting the first into the divergence of
$vb(E)^(n-1/2)=-(vb(A)^n-vb(A)^(n-1))/h_t-D_h^+ phi^(n-1)$
(deduction from #christlieb and the compatible spatial operators in #ryu).
It need not hold for the directly loaded state $n=0$ or where Mur replaces the scalar bulk update.
Keep both startup and selected evolved residuals, without treating them as independent bulk constraints.

Save each dimensional $L^2$ and maximum norm once, with the fixed scales needed to compute
$
  q_2=w norm(Q_h)_(2,h)/(epsilon norm(vb(E)_"ref" (-h_t/2))_(2,h)),
  quad
  g_2=w norm(G_h)_(2,h)/norm(vb(A)_"ref" (0))_(2,h).
$
Use the corresponding initial reference maximum norms for $q_infinity$ and $g_infinity$.
All denominators are fixed for one run, positive, and recorded.
They do not contain an extra factor of the cell width.
Plot the normalized values without duplicating the raw-norm plots.
Do not subtract startup residual norms from evolved ones.
Separately sampling exact fields and applying the stencil is a different, optional truncation check.
The existing vacuum-only `Monitors` interface does not implement these dielectric diagnostics.

=== Schedule and comparison times

Save the initial state and targets $tau in {T_3/4,T_3/2,3T_3/4,T_3}$.
For each positive target select the nearest native time independently:
$
  n_A=floor(tau/h_t+1/2), quad
  n_(E B)=floor(tau/h_t+1), quad
  n_phi=floor(tau/h_t).
$
Evaluate $G_h$ at $n_A$, $Q_h$ at $n_(E B)$ and the TM\_r $e_(phi,2)$ at $n_phi$.
Startup field and constraint rows use storage state $n=0$,
with their distinct time labels from the table above.
Advance through the largest required index and retain each observation at its own step.
Save target time, actual time, index and mesh with every row.
No analytic reference is evaluated at a rounded target when the numerical sample has a different time.

The nearest-time rule leaves offsets of at most $h_t/2$.
Use these recorded times when interpreting order estimates.
If this offset materially changes the trend, repeat the comparison with documented second-order time interpolation to common targets,
using bracketing native samples and applying the same rule on every mesh.

=== Spatial output

On one documented mesh, save a $z=0$ section of E and R2 B at initialization and the final field target.
This plane shows axial and face-diagonal directions, not a body diagonal.
Declare the plotting grid and trilinearly interpolate each component from its native positions.
For error plots, interpolate the native difference $F_"num"-F_"ref"$.
Subtracting the continuum at the plotting points after interpolating $F_"num"$ would also include interpolation error.
Global norms remain uninterpolated.
Add other sections or volume output only when a discrepancy needs localization.

#block(breakable: false)[
== Analysis and evidence required

For every positive error and consecutive pair of meshes, compute
$
  p_"obs"=ln(e(h)/e(2h/3))/ln(3/2).
$
]
This is the usual observed-order estimate (#salari, §4.2).
Plot startup field/constraint errors and evolved field/potential/constraint errors against $h$.
Tabulate the available successive orders.
An error at roundoff or exactly zero has no useful logarithmic order.

For the smooth solution and stable centred scheme,
second-order combined convergence is the expectation from the staggered differences and leapfrog update (#ryu).
It must be established over the measured range.
A finest-pair order far from two, increasing errors, or a residual that does not decrease calls for inspection of sampling,
initialization, time labels and resolution before making an accuracy claim.
Do not assign an arbitrary percentage tolerance or discard a coarse result to force the expected slope.

Before refinement, independently differentiate the potential pair at generic interior points and compare with its field evaluator.
Check the zero continuation and both loaded histories separately.
One small run should confirm the global reductions across an MPI partition.
Save the input, exact source revision and local modifications, build, ranks, threads,
machine, elapsed time, peak memory and all error rows with the study.
