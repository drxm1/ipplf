# Planar dielectric interface in a periodic box

Build `planar-interface` with `IPPL_ENABLE_CHDR=ON`. By default lengths use cell units;
the speed of light, vacuum permittivity and vacuum permeability are one. Every
case has an x-normal planar dielectric in a box periodic on every axis. Packet
presets use relative permittivity 2.13; the plane measurement can choose another
permittivity of at least one. The incidence plane is x-z, with y out of plane.

```sh
OMP_NUM_THREADS=1 planar-interface --input half-space.yaml
OMP_NUM_THREADS=1 planar-interface plane
OMP_NUM_THREADS=1 planar-interface half-space --output periodic-half-space
OMP_NUM_THREADS=1 planar-interface slab --output periodic-slab
OMP_NUM_THREADS=1 planar-interface brewster-te --output periodic-brewster-te
OMP_NUM_THREADS=1 planar-interface brewster-tm --output periodic-brewster-tm
python3 ../interface-movies/render_movie.py periodic-half-space
python3 ../interface-movies/render_movie.py periodic-brewster-tm --field B_raw
```

`--input FILE` uses the [common YAML schema](../../core/input/README.md).
[half-space.yaml](half-space.yaml) selects the existing full packet case.
The historical commands remain available and use the same setup and execution
path. Do not combine them with `--input`.

The packet cases accept `parameters: {smoke: true}`, an optional `run.steps`
or `run.duration`, and the observations `fields` and `diagnostics`. Both
observations default to every ten completed steps, including the initialized
state and the final state. Their schedules can be changed independently with
`first`, `every`, `last`, `include_final` and `enabled`. The `diagnostics`
observation writes the existing energy quantities; it adds no new measurement.
Output directories must be new or empty. A successful recorded run includes
`run.json` with the supplied input, resolved schedules and actual solver settings.
Explicitly enabled packet observations require `output.directory`; disabled
observations may be declared without it. Omitting both output and observations
retains the historical packet command's execution without saved files.

For SI output, set `units.system: si` and provide `length_m` and
`charge_density_C_m3` as described in the common schema. The presets retain
their normalized geometry and initial data; those scales specify their physical
interpretation. `run.duration` then uses seconds. Saved coordinates, times and
fields use SI; energy-column conventions are described in the
[interface-movies README](../interface-movies/README.md#output-and-rendering).

## Plane-pulse reflection measurement

The `plane` case uses an exact periodic lattice pulse and two Fourier probes.
TE measures the fixed Cartesian electric component E_y; TM measures H_y. Empty
parameters select the recorded normal-incidence example. YAML parameters resolve
a new row before field allocation, including the mesh, pulse, probe locations,
separate gates and step count. All presets and requested rows use this same path.

```yaml
schema_version: 1
case: plane
units: {system: normalized}
parameters:
  polarization: TM
  incidence: vacuum
  epsilon: 2.13
  angle_degrees: 30
  cells_per_wavelength: 20
  cell_shape_y: 1
  cell_shape_z: 1
output: {directory: plane-30-tm}
observations:
  diagnostics: {every: 10}
```

The two input forms below are exclusive. At fixed cell shape, `angle_degrees`
uses the nearest integer transverse cell count and adjusts the central frequency
to retain the requested continuum incidence angle. The actual cells per
dielectric wavelength can therefore differ from the requested value.
`transverse_period` instead keeps the requested frequency/resolution and derives
the continuum angle; its length must be an integer number of at least two z cells.
Both use one tangential period. This follows from periodic phase matching
[Oskooi and Johnson, section 4.5](https://arxiv.org/abs/1301.5366).

| Parameter | Meaning and default |
|---|---|
| `polarization` | `TE` or `TM`; default `TE`. |
| `incidence` | `vacuum` or `dielectric`; default `vacuum`. |
| `epsilon` | Dielectric relative permittivity, at least one; default 2.13. |
| `angle_degrees` | Continuum phase angle relative to the interface normal, strictly between -90 and 90; default zero when no period is supplied. |
| `transverse_period` | Alternative physical period in normalized code-length units; mutually exclusive with angle. |
| `cells_per_wavelength` | Requested dielectric resolution at the central frequency, greater than two; default 20. |
| `cell_shape_y`, `cell_shape_z` | Cell spacings relative to h_x=1; both default to one and remain fixed during angle resolution. |
| `courant` | Fraction of the Standard vacuum CFL limit, strictly between zero and one; default 0.9. |
| `tail_sigmas` | Gaussian spectral/envelope extent, greater than three; default 8. |
| `cutoff_fraction` | Fraction of the normal-wavenumber interval above the shared propagation cutoff excluded from the useful band; default zero at normal incidence and 0.5 otherwise. Must be in [0,1), and positive at oblique incidence. |
| `probe_cells` | Interface-to-probe offset in x cells, with the electric half-cell offset added for TE; default 3. |
| `timing_samples` | Samples over the useful normal-wavenumber band for travel-time estimates; default 513, minimum 17. |
| `separation_factor` | Factor of at least one enlarging both computed layers; default one. |
| `gate_padding_steps` | Additional samples after the computed gate ends; default zero. Layer sizes are recomputed to exclude returns through the extended gates. |

Only propagating incident/transmitted bands are supported by this measurement.
The material face remains aligned with x. Its position follows the resolved
layer thickness; general shifted/rotated input is not yet implemented.
The layout uses numerical group velocities and checks second-interface returns,
both periodic wraps and the analytic incident's periodic image. It is a sampled
arrival estimate, not a rigorous error bound. Qualify a new configuration with
denser timing samples and enlarged layers/gates before making an accuracy claim.
The dispersion is [Fallahi, section 3.1.3, Eq. (3.17)](https://arxiv.org/abs/2009.13645);
the group velocity is its derivative
[Trefethen, section 1](https://doi.org/10.1137/1024038).

`parameters: {preset: normal}`, `{preset: brewster}` and
`{preset: dielectric}` select the three recorded numerical layouts, including
their explicit gate lengths. A preset cannot be combined with physical parameter
overrides; use a requested row for those. These are ordinary configurations, not
separate evolution paths.

Standard output contains measured central reflection/transmission, effective
angle/resolution and energy drift. Continuum/lattice comparison formulas remain
separate from the measurement. With `output.directory`, `spectra.csv` contains
the three existing frequencies and complex face-referred amplitude ratios;
`run.json` adds the effective plane settings under `effective.plane`.
Those authoritative plane settings use normalized code units even for SI output;
the existing output factors provide the conversion. Spectrum amplitudes are
dimensionless; its frequency column uses inverse code time or rad/s according
to the selected output units. The material solver uses harmonic face
permittivities and the generalized-gauge potential equations
[Chew, 2014, section 2](https://arxiv.org/abs/1406.4780).

The supported observations are `fourier` and `diagnostics`: Fourier probes sample
every completed state and apply their resolved gates internally; diagnostics may
change cadence but must include the initialized and final states. Neither can be
disabled. An explicit duration may extend the run but cannot truncate either
gate. The plane case does not yet produce field frames.

## Packet cases

The packet cases use the x-z plane with y out of plane and keep periodic
boundaries on every axis. Full meshes are 512 by 2 by 512; `--smoke` selects
32 by 2 by 32 and four steps. The packet starts at z=256 for the normal cases
and z=128 for the Brewster cases. Normal cases run 560 steps and Brewster cases
1,000, with output every ten steps.

The material occupies x from zero to 230 in the half-space cases, or from 156
to 230 in the slab. Its second interface lies at the periodic wrap. The
translating Gaussian packet is approximate and intended for qualitative field
pictures. With `--output`, a case writes `wave.pvti.series`, `material.pvti`
and `energy.csv` in the format described in the
[interface-movies README](../interface-movies/README.md), whose `render_movie.py`
renders them.

## Checks

`ctest -R chdr.interface -E movies -LE long` runs the plane-pulse
checks at normal incidence against the independent lattice and Fresnel values,
strict plane input checks and a new oblique row qualified with denser timing
sampling and enlarged layers/gates
(`test_plane_interface.cpp` calls `PlanePulse.h` directly; it does not run the
`planar-interface` executable) and the four-step smoke run of the half-space
case of `planar-interface`, using YAML with separate field and diagnostic
cadences. The Brewster and two-rank oblique dielectric seam checks carry the
`long` label.
