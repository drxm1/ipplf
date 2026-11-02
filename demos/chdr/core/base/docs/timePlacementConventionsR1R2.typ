#import "../../../docs/manual/preamble.typ": *

== Time-placement interpretation of the stored fields

The electric field $vb(E)$ is usually reconstructed from $(phi, vb(A))$ via $vb(E) = - partial_t vb(A) - grad phi$.
Similarly, the magnetic field $vb(B)$ via $curl vb(A)$.

However we need to be careful about at which times the quantities $phi^n$ and $vb(A)^n$ (after $n$ solve's) are interpreted to be, and how exactly we compute $partial_t A^n$. This becomes important once $phi eq.not 0$, because in that case we need to ensure that both parts of the computed field $E$, $partial_t vb(A)$ and $grad phi$ are expressions at the same time.

Consider the current such implementation `evaluate_EB()`. The discrete gradient operator in $grad phi^n$ acts only spatially and therefore, if $phi^n$ is at the time $t_n$ (which is my interpretation of how the time is assigned to that stored field in the existing FDTD solver's), the output $grad phi^n$ also is at $t_n$. In `evaluate_EB()`, $partial_t vb(A)$ is calculated via a TODO:ThisHasAName (TODO:VerifyThatsActuallyTheCase) $(A^n - A^(n-1))/h_t$, so its output theoretically lives at time $t_(n+1/2)$.
Since this is time is different from $grad phi^n$ at $t_n$ we have a problem as soon as $phi eq.not 0$. This makes the order of $E$ only 1 and not 2 (TODO:OrderWrtWhat,TODO:SupportClaim).

Therefore it can make sense to use MITHRA's half-step convention, where $phi^n$ (after $n$ solves) is interpreted to be the value at time $t_(n+1/2)$ while $A^n$ is still interpreted to be the value at time $t_n$. Then both output terms in the expression for calculating $vb(E)$ align in their time-location, at $t_(n-1/2)$ if one takes $grad phi^(n-1)$.

TODO:CheckMWinklMScPaperOnTheImplementedTimeConventions. I have to check in MWinkls MSc thesis that the stored $phi^n$ and $vb(A)^n$ really are stored at the same time $t_n$. If that is actually the case then the following table lists *R1* as the time-convention (in demos/chdr I sometimes also write `ebReconstructionConvention` instead of time-convention) that the currently-in-ippl FDTD solver assumes, and *R2* is MITHRA's convention where $phi$ and $vb(A)$ are shifted by half a time-step. The different implementation options described in TODO:SectionDoesntExistYet have different requirements on which time-level convention is used (see TODO).

#[
  #set text(size: 9pt)
  #table(
    columns: (auto, 1fr, 1fr),
    align: (left + horizon, center + horizon, center + horizon),
    inset: (x: 5pt, y: 4.5pt),
    stroke: none,
    table.hline(stroke: 0.6pt),
    [*after $n$ updates*],
    [*R1* (in IPPL's FDTD Solvers currently (TODO:CheckMWinkl))],
    [*R2* (half-step, MITHRA)],
    table.hline(stroke: 0.4pt),
    [$vb(A)^n$, $vb(J)^n$], [$t_n$], [$t_n$],
    [$phi^n$, $rho^n$], [$t_n$], [$t_(n+1\/2)$],
    table.hline(stroke: 0.4pt + luma(85%)),
    [$vb(E)$],
    [$-(vb(A)^n - vb(A)^(n-1))\/h_t - grad phi^n$\ #text(
        fill: red,
      )[at $t_(n-1\/2)$ and $t_n$]],
    [$-(vb(A)^n - vb(A)^(n-1))\/h_t - grad phi^(n-1)$\ at $t_(n-1\/2)$],
    [$vb(B)$],
    [$curl vb(A)^n$ at $t_n$],
    [$1/2 curl (vb(A)^n + vb(A)^(n-1))$ at $t_(n-1\/2)$],
    table.hline(stroke: 0.4pt + luma(85%)),
    [order of $vb(E)$ ($phi eq.not 0$)],
    [1],
    [2],
    table.hline(stroke: 0.6pt),
  )

  Implementing the time convention *R2* therefore requires changing the interpretation at which time $phi^n, rho^n$ are interpreted to be (in `demos/chdr/ChDRSolver` etc). The current implementation in `demos/chdr` allows for the time placement `R2` combined with the spacial placement `staggered` (with its own operators, see TODO).

  TODO:(Describe how this interacts with the existing Mur BCs implementation (it is compatible, describe how)).

]
