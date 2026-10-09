# Periodic ChDR solver audit

This module adds a periodic, normalized-vacuum shell over IPPL's Standard and
NonStandard potential solvers. The inherited `step()` remains the update;
`ChdrSolver<Base>` supplies reconstruction of fields and time labels and the optional
time-step setter.

Five residual monitors and Fourier probe are separate objects that might be useful for measurements along the line.
All work is under `demos/chdr/solver`; the existing vacuum demo and `src/MaxwellSolvers` are reference implementations.

The module provides **both spatial readings**; applications select their
sampling convention explicitly. The supported common-time contract is
R2; R2Prime is its relabelled equivalence experiment, and R1 is the inherited
regression/control. Material, particles and open boundaries belong to later
units.

The constructor requires an explicit reconstruction. Select `Reconstruction::R2`
for the common-time contract. Both collocated and staggered routes are retained
as selectable implementations for later material and sourced applications; the
choice for a particular result will follow its benchmarks. The present vacuum
module supports Standard and NonStandard in the collocated reading, and Standard
in the staggered reading. Record the selected `Reading` and stencil with every result.

## Units and labels

Use normalized units \f$c=\epsilon_0=\mu_0=1\f$, with beam direction z and the
future interface normal x, with the dielectric on \f$x<0\f$. Electron charge is
negative in FEL-style source input. For conversion back to SI, \f$c\f$, \f$\epsilon_0\f$
and \f$\mu_0\f$ in the formulas below denote their SI vacuum values. Given a
length scale \f$L_0\f$ and potential scale \f$\Phi_0\f$, the conversions are
\f$\mathbf x_{\rm SI}=L_0\mathbf x\f$,
\f$t_{\rm SI}=L_0t/c\f$,
\f$\phi_{\rm SI}=\Phi_0\phi\f$,
\f$\mathbf A_{\rm SI}=\Phi_0\mathbf A/c\f$,
\f$\mathbf E_{\rm SI}=\Phi_0\mathbf E/L_0\f$,
\f$\mathbf B_{\rm SI}=\Phi_0\mathbf B/(cL_0)\f$,
\f$\rho_{\rm SI}=\epsilon_0\Phi_0\rho/L_0^2\f$,
\f$\mathbf J_{\rm SI}=\epsilon_0c\Phi_0\mathbf J/L_0^2\f$.
These are dimensional substitutions into the potential equations of
[Fallahi, section 3.1.1, Eqs. (3.6)–(3.9)][fallahi]. No SI input adapter is
implemented here.
## Library references

The citation keys used in source comments identify these works:

- `fallahi2020mithra20fullwavesimulation`: Arya Fallahi (2020), [MITHRA 2.0: A Full-Wave Simulation Tool for Free Electron Lasers][fallahi].
- `christlieb2024GaugeConserving`: Andrew J. Christlieb, William A. Sands and Stephen R. White (2025), [A Particle-in-Cell Method for Plasmas with a Generalized Momentum Formulation, Part III: A Family of Gauge Conserving Methods][gauge].
- `chew2014generalizedgauge`: Weng Cho Chew (2014), [Vector Potential Electromagnetic Theory with Generalized Gauge for Inhomogeneous Anisotropic Media][chew].
- `ryu2016potentialfdtd`: Christopher J. Ryu, Aiyin Y. Liu, Wei E. I. Sha and Weng Cho Chew (2016), [Finite-Difference Time-Domain Simulation of the Maxwell–Schrödinger System][ryu].
- `trefethen1996FiniteDifferenceSpectral`: Lloyd N. Trefethen (1996), [Finite Difference and Spectral Methods for Ordinary and Partial Differential Equations][trefethen].
- `clemens2001FiniteIntegration`: Markus Clemens and Thomas Weiland (2001), [Discrete Electromagnetism with the Finite Integration Technique][fit].
- `jeannerod2013ImprovedErrorBounds`: Claude-Pierre Jeannerod and Siegfried M. Rump (2013), [Improved Error Bounds for Inner Products in Floating-Point Arithmetic][roundoff].
- `oskooi2010meep`: Ardavan F. Oskooi et al. (2010), [Meep: A Flexible Free-Software Package for Electromagnetic Simulations by the FDTD Method][meep].
- `heinzel2002SpectrumSpectralDensity`: Gerhard Heinzel, Albrecht Rüdiger and Roland Schilling (2002), [Spectrum and Spectral Density Estimation by the Discrete Fourier Transform (DFT), Including a Comprehensive List of Window Functions and Some New Flat-Top Windows][window].

[fallahi]: https://arxiv.org/abs/2009.13645
[gauge]: https://arxiv.org/abs/2410.18414
[chew]: https://arxiv.org/abs/1406.4780
[ryu]: https://doi.org/10.1109/JMMCT.2016.2605378
[trefethen]: https://people.maths.ox.ac.uk/trefethen/pdetext.html
[fit]: https://doi.org/10.2528/PIER00080103
[roundoff]: https://doi.org/10.1137/120894488
[meep]: https://doi.org/10.1016/j.cpc.2009.11.008
[window]: https://hdl.handle.net/11858/00-001M-0000-0013-557A-5
