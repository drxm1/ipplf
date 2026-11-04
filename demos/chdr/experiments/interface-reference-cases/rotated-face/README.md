# Rotated planar interface

These inputs and analysis address RC-9, defined in the
[public reference-case index](../references.org). Interface tilt is
`interface_tilt_degrees`; incidence relative to that interface is `angle_degrees`.
They are independent. The convention is the x-z plane, with y out of plane.
All coefficients are referred to the requested physical face and use the same
fixed-component TE/TM definitions as the aligned experiments.

The input matrix and successful program exit do not establish accuracy.
The comparisons, rejected bins, spatial fits and isolation controls must be
assessed together. No universal convergence exponent is imposed.

## Geometry and comparison policy

The physical plane has normal `(cos(alpha), 0, sin(alpha))`. The periodic
construction uses a primitive tangent translation of `p` x cells and `q` z
cells. It retains the requested tilt exactly and adjusts the z cell spacing:

\[
 h_z/h_x = p\cos\alpha/(q\lvert\sin\alpha\rvert).
\]

This is the zero-normal-projection condition for that lattice translation,
combined with ordinary periodic phase matching; it is a deduction from the
coordinate rotation and Oskooi–Johnson §4.5, Eq. (4.30), cited in
[references.org](../references.org). The application records the complete
geometry and every retained carrier. The primitive pair is fixed across each
refinement ladder:

| Interface tilt | Primitive p/q | Actual h_z/h_x | Requested shape |
|---|---:|---:|---|
| 17 degrees | 2/3 | 2.1805684123 | (1, 2, 2) |
| 33 degrees | 4/3 | 2.0531532851 | (1, 2, 2) |

`cell_size` is 1, 0.5 and 0.25 for requested resolutions 20, 40 and 80.
The normal-incidence ladders therefore keep the physical central frequency
fixed. Automatic primitive-pair selection is deliberately not used for these
comparisons, because it can change cell shape between meshes.

The normal rows in `aligned-controls/` use the actual adjusted shape, the same material,
polarization, resolution, time-step policy and central frequency. Its aligned
face is the common physical coordinate 500; each measurement transfers phase
to its own requested face. Translation of the origin does not change the
face-referred half-space coefficients; see the phase-transfer derivation in
[references.org](../references.org), following Schneider §5.8.2. A staircase
has no single represented-face displacement, so the tilted record does not
invent one.

The oblique aligned controls retain the actual cell shape, requested dielectric
resolution 20, central continuum incidence 30 degrees and time-step policy.
Their `shape-17`/`shape-33` names identify the corresponding cell shape; their
material interface is grid aligned. These rows omit `face_position` so the C++
resolver selects qualified layer widths. Ordinary aligned periodicity generally
requires a different frequency adjustment from the tilted construction.
Record both actual frequencies, tangential wavenumbers, resolutions and Gaussian
widths; they are retained in `analysis.json` and its aligned-control record.
The comparison table marks unequal central reference parameters explicitly,
and such rows are not overlaid as exact matches in convergence figures. For
the selected nondispersive media the central continuum angle target agrees,
but the finite-grid resolutions and finite pulses are not identical.

Equal requested spectral controls do not make the lattice normal roots equal.
Consequently the actual Gaussian width can differ between an aligned and a
tilted start, even with equal central frequency. `aligned-comparisons.csv`
records both widths and both frequency/wavenumber pairs. It marks whether the
central reference parameters match; it does not claim identical finite pulses.
The dispersion source is Fallahi Eq. (3.17), cited in
[references.org](../references.org).

## Input matrix and controls

`inputs/` contains TE normal-incidence ladders at both tilts, a TM normal
ladder at 33 degrees, a coarse TM normal row at 17 degrees, and coarse TE/TM
rows at 30-degree incidence. The latter qualify
non-normal illumination without requesting the expensive finest oblique mesh.
Every baseline explicitly uses `tail_sigmas: 6`, 41 frequency bins and a
three-sigma measurement band. Vacuum incidence is selected throughout.
The additional 33-degree TM ladder measures its refinement trend directly.
Convergence is evaluated separately for each polarization, using the TE/TM
boundary conditions documented in [references.org](../references.org).

`controls/` separates the following changes at normal incidence:

- `packet-width`: `tail_sigmas: 8` and `band_sigmas: 4`. Changing the tail
  parameter also changes the Gaussian width and layout. This is a combined
  packet-width/truncation check, not a pure tail-error measurement. The band
  factor changes with it to retain the same queried normal-wavenumber range;
  actual frequency equality is still checked from the records.
- `layout`: larger separation, denser timing sampling and a longer gate.
- `period`: a larger pinned period with the baseline period multiplier's parity,
  retaining baseline timing and probe settings. Verify the actual coordinates
  relative to the phase anchor, gates and queried frequencies before comparing.
  This avoids a probe-rounding change that can accompany an unconstrained
  period enlargement. The retained mode grid and incident amplitude still
  change with the domain; this is not an identical finite pulse.
- `probe-distance`: the box centre moves from 12 to 18 x-cell widths from the face.
- `probe-extent`: the Cartesian half-width changes from two to three cells.
- `n80-isolation`: combined layout, distance and extent changes on the finest
  TE meshes and the 33-degree TM mesh, to compare their effect with the
  corresponding discretization error.

The application derives gates with separate directional group velocities.
Those path estimates need measured controls; they are not rigorous bounds
(Trefethen §1, cited in [references.org](../references.org)). The recorded
single-order Floquet exclusion does not eliminate evanescent contamination
near the staircase. Distance/extent changes and fit residuals address that
remaining extraction concern; the periodic-order basis is Oskooi–Johnson §4.5.

At nonzero incidence, changing the periodic box can change its selected tangent
mode and actual central frequency. Do not interpret such a change as an
isolation estimate. Use integer `period_cells_multiplier` and `tangent_mode`
values that preserve the baseline's exact ratio of tangent index to period
multiplier. A common integer multiple is one option; reducing that ratio can
also give smaller compatible integer pairs. Keep the primitive pair, cell
shape and incidence unchanged. This preserves the tangential wavenumber and
central frequency, apart from floating-point roundoff; it follows directly
from the periodic phase condition in Oskooi–Johnson Eq. (4.30). The C++ resolver
checks the directional travel-time estimates for the pinned period. Preserve
the exact resulting control input and verify the recorded frequency and
tangential wavenumber before comparing. `controls.csv` refuses a raw-coefficient isolation
difference when those central reference parameters differ. It compares only
equal recorded frequencies, with no interpolation, and reports the common valid
bin count. Boundary, packet and probe controls remain pending until their
completed measured records have been supplied.

Two representative oblique controls are supplied. TE at 17 degrees uses
`(M, j) = (942, 111)` against baseline `(628, 74)`; TM at 33 degrees uses
`(476, 60)` against `(357, 45)`. Both increase separation, probe distance and
extent, timing sampling and gate padding. Their period ratios are 3/2 and 4/3.
Actual layer widths, source/probe coordinates, gates, modal sampling and fit
residuals remain part of the comparison. Larger doubled-period fallbacks under
`controls/alternatives/` are unexecuted alternatives, not evidence. The other
two oblique baseline combinations have no dedicated enlarged-domain control
in this selected set; do not silently assign them those control results.

## Analysis and retained evidence

Launch explicit inputs with `../scan.py`. Analyse completed tilted runs and
the corresponding completed aligned controls, for example:

```sh
python3 plot.py TILTED_RUN_20 TILTED_RUN_40 TILTED_RUN_80 CONTROL_RUN \
  --aligned-run ALIGNED_RUN_20 --aligned-run ALIGNED_RUN_40 \
  --aligned-run ALIGNED_RUN_80 \
  --reference-binary /path/to/interface-reference --output rotated-analysis
```

Completed `--scan-manifest` arguments are also accepted for the tilted data.
Use exact run paths for a deliberate subset of an unfinished campaign.
The analysis reuses `references.py` to recompute the declared complex fits
from every valid bin's raw point spectra and lifted component coordinates.
The roots are the recorded directional lattice roots; the incident root is
not negated to construct a reflected wave. The C++ reference calculator supplies
continuum coefficients and flux weights only. There is no aligned lattice
interface oracle for the staircase.

Every run retains `comparison.csv`, `roots-and-fits.csv`, `projection-audit.csv`,
band figures, spatial-residual figures and the exact C++ reference invocation.
The angular table distinguishes the continuum reference angle from the incident
lattice phase angle, obtained from the recorded normal/tangent wavevector.
The latter is not a group-velocity direction. Frequency-dependent incidence is
retained; a single central angle does not describe the whole pulse band.

The aggregate outputs are `central.csv`, `geometry.csv`, `orders.csv`,
`controls.csv`, `aligned-comparisons.csv` and `analysis.json`. Complex and
magnitude errors are separate. Reported exponents are two-level observations,
not assertions of an asymptotic order. Convergence figures overlay available
aligned rows only when their actual cell shape, time step, incidence and central
reference parameters match. Fit residuals and raw-to-published audit
agreement are also distinct from physical accuracy. The `1e-10` audit tolerance
checks coefficient serialization/re-projection, while the input fit limit
controls bin usability. Neither is a certified physical error bound.
Invalid bins remain in the tables; the raw spectra must be retained with them.

Keep `run.json`, input/attempt records, `spectra.csv`, `probe-spectra.csv`,
`modes.csv`, generated references and comparison tables for every reported run.
Save field frames only for the representative coarse rows: fine ladders,
aligned controls and isolation controls disable fields. Six Float64 field
components cost approximately 48 bytes per cell per saved frame, before small
VTK metadata overhead; the material array is separate. Estimate storage from
the actual resolved cell count before a long run. Preserve selected coarse
field collections and movies as large binary artifacts, alongside the small
text records required to reconstruct the comparisons.
