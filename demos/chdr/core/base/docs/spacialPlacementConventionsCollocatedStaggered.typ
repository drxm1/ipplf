
#import "../../../docs/manual/preamble.typ": *

== Spacial placement convention of the stored fields <sec-spatial-placement>

*TODO*: This section right here describes the collocated vs staggered interpretation.
It is vastly interconnected with the implementation choices and operators of the collocated and staggered readings.
The staggered reading can have a few benefits which I should document here. Some useful references / maybe its origins are connected to Yee1966MaxwellIsotropicMedia--Published and Ryu's paper.

*TODO*: This will be connected to the options within chew, see 3wPres slides about Options within O1. I will describe here the "virtually" staggered storage and everything that is connected to handling Chew's gauge within that convention belongs to `core/solver/formulations` docs.

*TODO*: What also belongs here is implications on loading initial histories and saving and interpreting outputs.

=== The staggered placement
For staggered fields, spatial subscripts in the manual identify mesh
locations unless explicitly identified as array indices.
An integer vector $bold(i)$ denotes a cell centre at
$vb(x)_bold(i)$. Its positive $d$-face has index
$bold(i) + 1/2 bold(e)_d$ and lies at
$vb(x)_bold(i) + (h_d/2) vb(e)_d$.
Thus a centre scalar is written $phi_bold(i)$, while a vector component
on that face is written $A_(d,bold(i)+1/2 bold(e)_d)$.
This is the staggered placement of Ryu et al.
@ryu2016potentialfdtd[Appendix A, Eqs. (19)-(23)], with the integer
locations assigned to our cell centres.

The implementation still uses integer array indices for every component.
For example, the vector component at array index $(i,j,k)$ represents
$A_(d,bold(i)+1/2 bold(e)_d)$, where $bold(i) = (i,j,k)$.
The face coefficient stored at that same array index is
$epsilon_(d,bold(i)+1/2 bold(e)_d)$.
Expressions without spatial subscripts denote the corresponding fields
and operators. Superscripts specify the actual time level, as defined in
@sec-time-placement.

A sample is therefore written with its location and its time level,
for example $A_(d,bold(i)+1/2 bold(e)_d)^n$.
For a vector quantity $F$ we write $bold(delta)_d^F$ for the index offset of
its component $d$, so that this component is $F_(d,bold(i)+bold(delta)_d^F)$:
$
  bold(delta)_d^A = bold(delta)_d^J = bold(delta)_d^E = 1/2 bold(e)_d,
  quad
  bold(delta)_d^B = 1/2 sum_(q != d) bold(e)_q.
$
The following table lists the samples of the staggered placement in the R2
time convention after $n$ completed solves.
These locations and time levels are the ones returned by
`LevelContract::positionOffset` and `LevelContract::halfStepLabel`.
#table(
  columns: (auto, auto, 1fr, auto),
  inset: 5pt,
  stroke: 0.4pt,
  table.header([*Quantity*], [*Sample*], [*Location*], [*Time*]),
  [Scalar potential, charge density],
  [$phi_bold(i)^(n+1/2)$, $rho_bold(i)^(n+1/2)$],
  [$vb(x)_bold(i)$],
  [$t_(n+1/2)$],
  [Vector potential, current density],
  [$A_(d,bold(i)+1/2 bold(e)_d)^n$, $J_(d,bold(i)+1/2 bold(e)_d)^n$],
  [$vb(x)_bold(i) + (h_d/2) vb(e)_d$],
  [$t_n$],
  [Electric field],
  [$E_(d,bold(i)+1/2 bold(e)_d)^(n-1/2)$],
  [$vb(x)_bold(i) + (h_d/2) vb(e)_d$],
  [$t_(n-1/2)$],
  [Magnetic field, averaged over both slots],
  [$B_(d,bold(i)+bold(delta)_d^B)^(n-1/2)$],
  [$vb(x)_bold(i) + 1/2 sum_(q != d) h_q vb(e)_q$],
  [$t_(n-1/2)$],
  [Raw magnetic field of the current slot],
  [$B_(d,bold(i)+bold(delta)_d^B)^("raw",n)$],
  [$vb(x)_bold(i) + 1/2 sum_(q != d) h_q vb(e)_q$],
  [$t_n$],
)
