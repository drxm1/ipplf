# Potential ghosts beyond absorbing faces

After loading `A_n` and `A_nm1`, call `prepareHistories()` before advancing
the solver or evaluating fields. It exchanges both histories and completes
their physical ghosts without changing owned values or applying Mur.
Solvers require one ghost layer, three global cells on open axes,
and two owned cells per rank and axis. Periodic templates retain their
field boundary conditions, including on nonperiodic layouts.

`OpenFaceClosure` and `OpenFaceClosureMPI` run the same checks on one and
two ranks. Run them with `ctest -R '^OpenFaceClosure' --output-on-failure`.
The independent forty-coefficient polynomial checks the inherited electric
field and current-level magnetic curl on the outermost layer. Each error
is at most `1e-13` of its positive analytic maximum over the whole box.
Zero physical ghosts must give an error of at least `0.1` in both fields.
The polynomial fixture splits along z and compares every owned field component
with a separate serial field, normalized by its nonzero maximum. It also records
binary64 equality, including signed zero. The explicit x-split thin partition
checks the full tensor extension at both time levels.
The periodic test checks both histories against analytically wrapped samples.

The citation keys in `OuterLayerOracles.h` denote:

- `fornberg1988FiniteDifferenceFormulas`: B. Fornberg (1988), *Generation of
  Finite Difference Formulas on Arbitrarily Spaced Grids*, Table 3, p. 704,
  [doi:10.1090/S0025-5718-1988-0935077-0](https://doi.org/10.1090/S0025-5718-1988-0935077-0).
- `fallahi2020mithra20fullwavesimulation`: A. Fallahi (2020), *MITHRA 2.0:
  A Full-Wave Simulation Tool for Free Electron Lasers*, section 3.1.1,
  Eqs. (3.6)–(3.8), [arXiv:2009.13645](https://arxiv.org/abs/2009.13645).
- `ryu2016potentialfdtd`: C. J. Ryu, A. Y. Liu, W. E. I. Sha and W. C. Chew
  (2016), *Finite-Difference Time-Domain Simulation of the Maxwell–Schrödinger
  System*, Appendix A, Eqs. (19)–(24),
  [doi:10.1109/JMMCT.2016.2605378](https://doi.org/10.1109/JMMCT.2016.2605378).
