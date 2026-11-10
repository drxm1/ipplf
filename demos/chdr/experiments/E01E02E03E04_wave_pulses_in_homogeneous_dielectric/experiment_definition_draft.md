# Experiment definitions
All experiments below consider a wave pulse initialized in a homogeneous dielectric (on the full mesh domain) that has \f$\epsilon_r = 2.13\f$. The experiments are in 3D, not just 2D. For the plane wave cases one might consider reducing the amount of cells in the y-axis (TODO to be determined).

## Experiment 01 (E01)
A plane wave pulse is initialized from an analytic continuum reference, moving in the z-direction.
It is propagated for some amount of time steps (not too few), with the O1-VS solver acting
on a domain that is fully filled with the homogeneous dielectric and has periodic boundary conditions in the x and y axis.
In the z-axis we have also periodic boundary conditions (this is the difference from experiment E02).
The analytic continuum solution of the analytic reference must be known because it serves as comparison.
The wave is propagated until it is close to the z-boundary of the mesh/domain.

Cases: We cover only a single transverse polarization.
Scope: Source-free radiation with zero scalar potential in a homogeneous dielectric.
Study: Convergence study against the continuum solution w.r.t. refinement. What exact measurables/metrics to compare is now described in `./E01.typ`.

All experiment details are described in `./E01.typ`.

### The analytic case of E01 and sampling/initialization
The exact analytic reference case and how we initialize/sample the field storage for E01 is described in `./E01.typ`.

### Measured quantities
Also described in `./E01.typ`.

## Experiment 02 (E02)
E02 uses E01's pulse and solver settings in a shorter domain, with Mur at
the left (-z) and right (+z) faces. The 6w domain retains the existing E01
comparison after 2w of travel, then lets the pulse reach the right face
after another w of travel. E02.typ defines its own short observation schedule.

We reuse E01's saved 2w profiles and error definitions for the propagation comparison.
After exit, measure the residual E/B norms relative to their initial norms
and keep the electric probe trace; no extra reference evaluation or window analysis is required.
The E02 refinement ladder is 32, 64, 128, 256 and 512 cells per pulse.

Additionaly we perform the following study:
1. Let the wave propagate for even longer such that it passes the Mur boundary at the z-domain-face.
2. Compare the residual fields and spatial profiles across the refinement ladder.

See `E02.typ` for draft specs.

## Experiment 03 (E03)
An exact spherical pulse propagates with O1-VS in a homogeneous dielectric,
with periodic boundaries and refinement in all three directions.
The radial TE and TM cases include a nonzero scalar potential in the TM case.
Compare the evolved fields and potentials with the continuum solution before boundary contact,
and measure the startup and evolved Gauss/gauge residuals.

See `E03.typ` for draft specs.

## Experiment 04 (E04)
Use E03's two cases with dielectric-speed Mur on all six faces,
first comparing with E03 before contact, then continuing through exit and return.
Measure energy returning through an interior cube, remaining electromagnetic energy,
and Gauss/gauge residuals in the interior and accessible boundary region.
The all-face dielectric-speed implementation prerequisite is specified in E04.typ.

See `E04.typ` for draft specs.

# Where E01-E04 live
E01 to E04 measurement code will live here in `experiments/01-04_wave_pulses_in_homogeneous_dielectric`.

# Initialization

There are two options on initialization, listed below.
We are choosing option A for these experiments. If it will show to be useful, option B can be implemented on top of it later.

## Initialization option A: direct sampling to a contiuum potential history
Difficulties with that approach:
- Exact initial samples of the potential do not produce exact reconstructed fields on a finite mesh.
  Consequence:
  1. We must measure the initial field error as part of the experiment.
  2. We must additionally check whether this "continuum sampling" satisfies the discrete gauge constraint.

## Initialization option B: construct a history using the discrete-lattice dispersion relation
Difficulties with that approach:
- Must account for its departure from the continuum initial data.
