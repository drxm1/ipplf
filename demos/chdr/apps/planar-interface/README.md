# Planar dielectric interface in a periodic box

Build `planar-interface` with `IPPL_ENABLE_CHDR=ON`. By default lengths use cell units;
the speed of light, vacuum permittivity and vacuum permeability are one. Every
case has one planar dielectric of relative permittivity 2.13 with its normal
along x, in a box that is periodic on every axis.

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

The `plane` mode is the exact-lattice measurement configuration: a quasi
one-dimensional box, an exact lattice plane pulse at normal incidence from
vacuum, and two Fourier probes. It prints the complex reflection and
transmission amplitudes, their Fresnel and exact lattice references, and the
energy drift to standard output. The historical `plane` command writes no files.
YAML with `case: plane` and `output.directory` also saves these measured complex
amplitudes in `spectra.csv` and writes `run.json`. Reflection and transmission
amplitudes remain dimensionless; the frequency column uses inverse code time or
rad/s according to the selected units. The material solver uses
harmonic face permittivities and the generalized-gauge potential equations
[Chew, 2014, section 2](https://arxiv.org/abs/1406.4780).

This measurement retains its full preset mesh, duration and Fourier gates.
`smoke: true` and duration changes are rejected. The supported observations are
`fourier` and `diagnostics`: Fourier probes sample every state and apply their
existing gates internally; diagnostics may change their cadence but must include
the initial and final states. Neither observation can be disabled. The plane
case does not produce field frames.

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
checks at normal incidence against the independent lattice and Fresnel values
(`test_plane_interface.cpp` calls `PlanePulse.h` directly; it does not run the
`planar-interface` executable) and the four-step smoke run of the half-space
case of `planar-interface`, using YAML with separate field and diagnostic
cadences. The Brewster and two-rank oblique dielectric seam checks carry the
`long` label.
