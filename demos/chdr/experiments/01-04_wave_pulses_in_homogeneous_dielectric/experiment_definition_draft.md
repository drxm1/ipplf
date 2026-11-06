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
Study: Convergence study against the continuum solution w.r.t. refinement. What exact measurables/metrics to compare is to be determined (TODO).

### The analytic case of E01 and sampling/initialization
The exact analytic reference case and how we initialize/sample the field storage for E01 is described in `./analytic_refs/E01.typ`.


## Experiment 02 (E02)
Exactly the same setup as E01, except that the z-boundary now has Mur BCs enabled.

We do the same cases and study as in E01 (to the point where the wave is close to the boudnary) to verify that the results are not affected by the Mur conditions as long as the plane wave pulse doesn't reach the domain-boundary.

Additionaly we perform the following study:
1. Let the wave propagate for even longer such that it passes the Mur boundary at the z-domain-face.
2. Measure the reflection that is due to Mur and probably also the error in Gauss's law. Also in an error convergence study.

## Experiment 03 (E03)
A circular (TODO: and/or spherical) wave pulse is initialized such that it then propagates towards the outer faces of the mesh domain.
Again it is initialized from a known analytic continuum reference, and again we have the exact continuum reference solution vs time.
That wave is propagated with the O1-VS solver.
We have periodic boundary conditions at each domain-boundary face.
The wave is propagated until it is close to the domain boundaries.

Cases: (TODO: if there is also a TE/TM destinction here, cover them otherwise ignore).
Study: Again an error convergence study.

## Experiment 04 (E04)
Exactly the same setup as E03, but with all-mur epsilon=2.13 absorbing domain-boundary faces.
We do the same cases and study as in E01 again to verify that the results are not affected by the Mur conditions as long as the wave does not reach the domain-boundary.

Additionally we perform the following study:
1. Let the wave propagate for even longer such that it passes the Mur boundary at the domain face(s).
2. Measure the reflection that is due to Mur and probably also the error in Gauss's law. Also in an error convergence study.

# Which "app" code this will run on
To be determined. Vacuum-wave is probably not suitable because it contains a lot of legacy comparison boilerplate.
interface-movies might also not be suitable because it has the code for interface movies from earlier.
planar-interface wave initializations might be usable but this case here concerns no interface.
Ideally we write a small standalone app that supports homogeneous dielectric with the four wave cases for the four experiments above,
with configurable domain sizes etc.

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
