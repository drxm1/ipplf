#import "@preview/physica:0.9.8": *
#set page(paper: "a4", margin: 22mm, numbering: "1")
#set text(size: 10.5pt)
#set par(justify: true)
#set math.equation(numbering: "(1)")
#show heading: set block(above: 1.3em, below: 0.6em)

#let chew = link("https://arxiv.org/abs/1406.4780")[Chew, §2, Eqs. (5)–(13)]
#let ryu = link("https://wcchew.ece.illinois.edu/chew/qem_papers/ryu_chris052216.pdf")[Ryu et al., Appendix A, Eqs. (19)–(27)]
#let fallahi = link("https://arxiv.org/abs/2009.13645")[Fallahi, §3.1]
#let salari = link("https://www.osti.gov/servlets/purl/759450")[Salari and Knupp, §§4.1–4.2]

= Experiment E03 detailed specification\*

E03 measures convergence of a spherical electromagnetic pulse in a homogeneous dielectric,
before it reaches the periodic domain boundaries.
The fields vary in all three spatial directions.
We refine all three cell widths and the time step together.
This complements the long axial propagation in #link("E01.typ")[E01];
it does not repeat E01's long travel distance.

We use two exact source-free solutions, with electric or magnetic field transverse to the radial direction.
Their potentials and fields are given explicitly below.
The second case has a nonzero scalar potential and therefore also tests its time placement and contribution to the electric field.
These are specifications for implementation and measurement, not simulation results.
The corresponding open-boundary experiment is #link("E04.typ")[E04].

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
It is a normalization choice; it does not change the wave equation.

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
This is an outgoing pulse with spherical support and a dipole angular pattern;
its electromagnetic amplitude is not spherically symmetric.

For an implementation without numerical differentiation,
set $u=s/w$ and $q=1-u^2$ inside the support.
Direct differentiation of @e03-profile gives
$
  p'(s)=-2u p(s)/(w q^2), quad
  p''(s)=p(s)/w^2 (4u^2/q^4-2/q^2-8u^2/q^3).
$ <e03-profile-derivatives>
Multiply these expressions by $w^2$ to obtain $f'$ and $f''$.
Test the support condition before evaluating $1/r$, $vb(n)$ or any inverse power of $q$.
Near $abs(u)=1$, evaluate the exponential times the rational factors in a numerically safe order,
or in logarithmic form; do not form an overflowing factor and multiply it by an underflowed zero.

#block(breakable: false)[
=== Two electromagnetic cases

Choose one fixed unit vector
$
  vb(d)=frac(lr((1,2,3)),sqrt(14)).
$
]
This avoids aligning the dipole with a mesh axis or body diagonal.
It exercises all Cartesian components without adding an orientation sweep.
Introduce the following scalar derivatives of @e03-generator:
$
  a=psi_r=f'/r-f/r^2, quad
  b=psi_(r r)=f''/r-2f'/r^2+2f/r^3, quad
  c=f''/r-f'/r^2.
$ <e03-radial-derivatives>
Here and below $f,f',f''$ are evaluated at $s$.
These expressions are used only inside the support.

*Case TE\_r.* Define
$
  vb(A)_"TE"=grad psi times vb(d)=a(vb(n) times vb(d)), quad phi_"TE"=0.
$ <e03-te-potentials>
Using $vb(E)=-partial_t vb(A)-grad phi$ and $vb(B)=curl vb(A)$ from #chew,
direct differentiation gives
$
  vb(E)_"TE"=v c(vb(n) times vb(d)),
$
$
  vb(B)_"TE"=(b-a/r)vb(n)(vb(n) dot vb(d))-(b+a/r)vb(d).
$ <e03-te-fields>
The subscript TE\_r means $vb(n) dot vb(E)_"TE"=0$.
The magnetic field generally has a radial component.
This is not a TE designation relative to any particular Cartesian boundary face.

*Case TM\_r.* Define the second pair of potentials by
$
  vb(A)_"TM"=vb(d) partial_t psi/v=-vb(d)f'/r,
  quad phi_"TM"=-v vb(d) dot grad psi=-v a(vb(n) dot vb(d)).
$ <e03-tm-potentials>
The corresponding fields are
$
  vb(E)_"TM"=v vb(B)_"TE", quad
  vb(B)_"TM"=-vb(E)_"TE"/v.
$ <e03-tm-fields>
Thus $vb(n) dot vb(B)_"TM"=0$,
while the electric field generally has a radial component.
Both cases must be run separately, with the same profile and grid.

=== Why these are source-free Maxwell solutions

The following is a deduction from the potential relations, Lorenz gauge and homogeneous wave equations in #chew.
Spatial and time derivatives commute with the constant-coefficient wave operator,
so every component in @e03-te-potentials and @e03-tm-potentials solves that wave equation.
The TE\_r vector potential is a curl, hence its divergence is zero.
For TM\_r,
$
  div vb(A)_"TM"+epsilon partial_t phi_"TM"
  =(1/v-epsilon v) vb(d) dot grad partial_t psi=0.
$
Both cases therefore satisfy the homogeneous Lorenz gauge and give source-free Maxwell fields.
For example,
$
  curl(grad psi times vb(d))
  =(grad grad psi)vb(d)-vb(d)nabla^2 psi,
$
which gives @e03-te-fields and, with $psi_(t t)=v^2 nabla^2 psi$, @e03-tm-fields.
No far-field or locally plane-wave approximation is used.
The radial terms proportional to $1/r^2$ and $1/r^3$ must be retained.

The dual relations in @e03-tm-fields give equal continuum electromagnetic energies for the two cases,
since $epsilon v^2=1$ (use the energy definition and Ellingson source in #link("E04.typ")[E04]).
They provide a useful independent algebraic check of the evaluator;
equal numerical errors are not required, because the two initializations use different potential histories.

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
The pulse travels one full support width.
The purpose is to test three-dimensional propagation, curvature, reconstruction and constraints over this fixed distance.
E01 remains the longer propagation study.

The periodic continuum reference is the sum of translated copies of the potentials over the lattice $(12w)ZZ^3$.
During E03's observation interval, only the central copy contributes anywhere in the domain,
because its support radius is below $6w$ and the neighbouring copies have not arrived.
Thus the formulas above are also the exact periodic reference on this interval.
Check this support condition at the actual selected times, including the scalar-potential time.
Periodic numerical tails can still reach a boundary earlier than the continuum support;
the measured field errors include any such contribution.

Let $N=N_x=N_y=N_z$ and $h=h_x=h_y=h_z=12w/N$.
We choose equal cell counts because the physical box is cubic and the pulse varies in every direction.
For a rectangular box, equal cell widths would instead require cell counts proportional to its side lengths.
Keeping one transverse direction fixed, as in E01, is unsuitable for this field.

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
Four levels give three successive order estimates.
The coarsest level is a baseline, not an assumed asymptotic resolution.

The time step satisfies the solver's vacuum Standard bound
$
  h_t=h/2<h/sqrt(3),
$
which is stricter than the homogeneous-dielectric bound here
(#fallahi, Eq. (3.20)).
This is a combined spatial and temporal convergence study;
it does not separately identify either order.

The finest level is substantial: one scalar double array occupies about 260 MiB before halos.
Measure peak memory and time per step on the first two levels,
then account for the actual live arrays and output buffers before running 324 cubed.
For fixed physical duration, the nominal cell-update count grows by $(3/2)^4$ per level.
This is an operation-count estimate, not a measured runtime.
If the last two orders have not stabilized, consider $N=486$ only after estimating its cost.
An incomplete or unresolved ladder must be reported as such;
do not reduce one axis or shorten only the finest run to obtain a favourable trend.

== Initialization and sample positions

We retain initialization option A of #link("experiment_definition_draft.md")[the experiment draft]:
sample the exact continuum potential histories directly.
Do not replace the continuum dispersion relation or project the sampled histories onto a discrete constraint.
Such a projection would define a different initialization experiment.

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
For TE\_r the scalar histories are zero.
For TM\_r, loading both scalar and vector slots at integer times would be an error.
The spatial staggering follows #ryu;
the R2 temporal interpretation uses `LevelContract` and #fallahi, Eqs. (3.50)–(3.55).
The experiment uses the actual `LevelContract` to obtain every offset and time.

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
Raw magnetic output has integer time and is a separate quantity
(#fallahi, Eqs. (3.50)–(3.52), with the spatial curl of #ryu).
Every analytic comparison uses the numerical component's own position and time.

Smooth continuum histories do not generally satisfy the discrete Lorenz constraint exactly.
For these curved fields, separately sampling the continuum derivatives does not reproduce a discrete curl identity.
The initial field and constraint errors must therefore be measured.
E01's exact cancellation for a transversely uniform pulse does not carry over.
Second-order startup errors are expected from Taylor expansion of the centred staggered differences in #ryu;
their magnitude and observed order remain measurements.

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
These norms follow #salari, Eqs. (15), (18).
They combine phase, amplitude and shape errors; none alone identifies the cause.

=== E03.meas05: scalar-potential error

For TM\_r, also save the relative scalar $L^2$ error $e_(phi,2)$,
using the scalar version of the preceding norm at #box[$(n+1/2)h_t$].
For TE\_r, save the absolute maximum $max abs(phi_"num")$ instead;
its exact reference is zero, so a relative scalar error is undefined.
This checks the extra history that distinguishes the two cases.

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
These are physical constraint residuals, not residuals of a rearranged update identity.
The difference operators are those of #ryu, evaluated on the actual histories.

Save the dimensional $L^2$ and maximum norms and the following dimensionless values:
$
  q_2=w norm(Q_h)_(2,h)/(epsilon norm(vb(E)_"ref" (-h_t/2))_(2,h)),
  quad
  g_2=w norm(G_h)_(2,h)/norm(vb(A)_"ref" (0))_(2,h).
$
Use the corresponding initial reference maximum norms for $q_infinity$ and $g_infinity$.
All denominators are fixed for one run, positive for both cases, and recorded.
They do not contain an extra factor of the cell width.
These normalizations are definitions for this experiment.
Measure the initial residual separately; do not subtract it from subsequent residual norms.
The existing vacuum-only `Monitors` interface does not implement these dielectric diagnostics.

=== Schedule and comparison times

Save the initial state and targets $tau in {T_3/4,T_3/2,3T_3/4,T_3}$.
For each positive target select the nearest native time independently:
$
  n_A=floor(tau/h_t+1/2), quad
  n_(E B)=floor(tau/h_t+1), quad
  n_phi=floor(tau/h_t).
$
Evaluate $G_h$ at $n_A$, $Q_h$ at $n_(E B)$ and $e_(phi,2)$ at $n_phi$.
At initialization use storage state $n=0$ for every quantity,
with its distinct time label from the table above.
Advance through the largest required index and retain each observation at its own step.
Save target time, actual time, index, mesh and case with every row.
No analytic reference is evaluated at a rounded target when the numerical sample has a different time.

The nearest-time rule leaves offsets of at most $h_t/2$.
Use these recorded times when interpreting order estimates.
If this offset materially changes the trend, repeat the comparison with documented second-order time interpolation to common targets,
using bracketing native samples and applying the same rule on every mesh.

=== Spatial output

At the initial state and the four target checkpoints,
save central-plane sections of every electric and magnetic component and their reference errors.
Use $x=0$, $y=0$ and $z=0$ sections to see face-axis and diagonal propagation.
For plots only, interpolate each native component to the declared plotting points with trilinear interpolation,
and evaluate the reference there at the actual field time.
State this interpolation on the plots; global error norms use native samples without interpolation.

Save full-domain E, R2 B and optional raw B snapshots at these checkpoints when storage permits.
The last two meshes need not retain all full-volume snapshots:
keep the norms and plane sections at every level,
and record which volume frames were omitted.
Record each array's position and time metadata through `LevelContract`.
Do not interpret fields with different native offsets as colocated vectors in a visualization without interpolation.

#block(breakable: false)[
== Analysis and evidence required

For every positive error and consecutive pair of meshes, compute
$
  p_"obs"=ln(e(h)/e(2h/3))/ln(3/2).
$
]
This is the usual observed-order estimate (#salari, §4.2).
Plot the four errors, the TM\_r scalar error and both constraint residuals against $h$,
showing initialization and final-time values separately.
Tabulate the three orders for each case.
An error at roundoff or exactly zero has no useful logarithmic order; report its magnitude instead.

For the smooth solution and stable centred scheme,
second-order combined convergence is the expectation from the staggered differences and leapfrog update (#ryu).
It must be established over the measured range.
A finest-pair order far from two, increasing errors, or a residual that does not decrease calls for inspection of sampling,
initialization, time labels and resolution before making an accuracy claim.
Do not assign an arbitrary percentage tolerance or discard a coarse result to force the expected slope.

The experiment is useful when it provides a resolved error trend for genuinely three-dimensional fields,
including the nonzero scalar-potential case,
and a measured startup/propagation baseline for E04.
It does not establish interface accuracy, absorption, long-time stability or independent spatial and temporal orders.

Before the full ladder, check the analytic evaluator against independent differentiation at several interior points,
the zero continuation, the dual relations and the loaded scalar times.
One small run should confirm the global reductions across an MPI partition.
These are focused implementation checks, not a request for a new large regression suite.
Save the input, exact source revision and local modifications, build, ranks, threads,
machine, elapsed time, peak memory and all error rows with the study.

== Sources

The equations above are either cited continuum/discrete relations or explicit deductions and measurement definitions.
The cited sources are held in the project literature library under these keys:

- `chew2014generalizedgauge`: W. C. Chew, _Vector Potential Electromagnetic Theory with Generalized Gauge for Inhomogeneous Anisotropic Media_ (2014). #chew.
- `ryu2016potentialfdtd`: C. J. Ryu et al., _Finite-Difference Time-Domain Simulation of the Maxwell–Schrödinger System_ (2016). #ryu.
- `fallahi2020mithra20fullwavesimulation`: A. Fallahi, _MITHRA 2.0: A Full-Wave Simulation Tool for Free Electron Lasers_ (2020). #fallahi.
- `salari2000CodeVerificationMms`: K. Salari and P. Knupp, _Code Verification by the Method of Manufactured Solutions_, SAND2000-1444 (2000). #salari.
