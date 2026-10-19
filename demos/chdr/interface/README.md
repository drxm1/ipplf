# Dielectric interface examples and movies

Build `interface-demo` and `interface-movies` with `IPPL_ENABLE_CHDR=ON`.
Lengths use cell units; the speed of light, vacuum permittivity and vacuum
permeability are one. Movie plots use z horizontally, x vertically and y out
of plane. The two y cells are periodic.

## Open-box annulus movies

```sh
OMP_NUM_THREADS=1 interface-movies half-space --output open-half-space
OMP_NUM_THREADS=1 interface-movies brick --output open-brick
OMP_NUM_THREADS=1 interface-movies slab --output open-slab
OMP_NUM_THREADS=1 interface-movies diagonal --output open-diagonal
OMP_NUM_THREADS=1 interface-movies half-space-tm --output open-half-space-tm
OMP_NUM_THREADS=1 interface-movies brick-tm --output open-brick-tm
OMP_NUM_THREADS=1 interface-movies slab-tm --output open-slab-tm
OMP_NUM_THREADS=1 interface-movies diagonal-tm --output open-diagonal-tm
python3 render_movie.py open-half-space
python3 render_movie.py open-half-space-tm --field B_raw
python3 measure_tm_constraints.py open-half-space-tm
```

The output directory must be new or empty. Add `--smoke` for a 64 by 2 by
64 mesh and four steps. `--steps N` overrides the case's duration. Full cases
use a 600 by 2 by 600 mesh with unit spacing and a frame every eight steps.
Two-rank runs can use `mpiexec --bind-to none -n 2 interface-movies ...`.

The relative permittivity is 2.13 inside the following rectangles and one
outside. Each rectangle is uniform along y; all open x/z faces remain in vacuum.

| Case | Material x interval | Material z interval | Steps |
|---|---|---|---|
| half-space | 300 to 564 | 36 to 564 | 1424 |
| brick | 300 to 480 | 36 to 300 | 1120 |
| slab | 300 to 390 | 36 to 564 | 1280 |

The half-space is a finite, vacuum-padded approximation for illustration.
All three starts are centred at x=150, z=300, with annulus radius 72 and
wavelength parameter 36. `AnnulusStart` samples the y component of the
vector potential at two time levels. Its radial profile is an antiderivative
of a Ricker pulse [Schneider, *Understanding the FDTD Method*, section 5.2.3].
The radial amplitude factor is a large-radius approximation to the vacuum
potential equation [Chew, 2014, section 2]; it is not an exact outgoing
discrete solution. The Gaussian tails are cut off outside the prescribed
annular band, and the driver checks that both histories start entirely in
vacuum, clear of the open faces.

These movies demonstrate field evolution. They do not establish interface
or absorber accuracy. The `diagonal` case uses material at x >= z, centred
initial data at (x,z)=(150,450), and 1,672 steps. Its demo-only solver sets
the upper x and lower z Mur face speeds to 1/sqrt(2.13); other face speeds
are one. Existing unit-speed edge and corner updates are retained. This is
an exploratory boundary treatment, not a validated heterogeneous absorber.
The diagonal movie and overview carry this limitation as a visible caption.
The standard O1-VS Mur path requires alpha = 1 and assumes vacuum permittivity
in the outer layers; the solver rejects other alpha values on open axes but
does not inspect the face permittivity. The diagonal demo overrides the boundary
speeds under the limitations stated above.

Appending `-tm` selects TM initial data for any of the four geometries, retaining
that geometry's duration. Unqualified case names retain the TE initial data.
The shared profile parameter 36 gives a TE field peak wavelength of about 36 cells,
but the extra TM spatial derivative shifts the TM field peak to about 29.5 cells
(large-radius prediction: `36/sqrt(1.5)`, about 29.4). The TM pulse is one derivative
higher, so a TE/TM comparison also compares different incident pulse shapes and spectra.
`TMAnnulusStart` forms the in-plane vector potential from the discrete curl of
the radial profile, with zero scalar potential. Both histories start in vacuum;
their weighted divergence cancels initially, following the same commuting-difference
construction as the periodic TM packet [Ryu et al., 2016, Appendix A]. The
initial field has out-of-plane B and in-plane E. This does not establish Gauss
constraint preservation by the later Mur update. `measure_tm_constraints.py`
records the weighted electric divergence at every saved time, separately in
the interior and boundary band, in `constraints.csv`. It excludes the lowest
x/z cell layers because their backward stencils need exterior E samples that
are not stored. Its dimensionless residual uses the initial displacement-field
peak and the smaller x/z spacing as a fixed normalization.
The TM movie carries a visible warning about its nonzero Gauss residual under
Mur; the recorded residual is a limitation of this demonstration, not evidence
of a constraint-preserving open boundary.

## Existing periodic packet examples

```sh
OMP_NUM_THREADS=1 interface-demo plane
OMP_NUM_THREADS=1 interface-demo half-space --output periodic-half-space
OMP_NUM_THREADS=1 interface-demo slab --output periodic-slab
OMP_NUM_THREADS=1 interface-demo brewster-te --output periodic-brewster-te
OMP_NUM_THREADS=1 interface-demo brewster-tm --output periodic-brewster-tm
python3 render_movie.py periodic-half-space
python3 render_movie.py periodic-brewster-tm --field B_raw
```

Packet examples also use the x-z plane with y out of plane, but retain
periodic boundaries on every axis. Full meshes are 512 by 2 by 512;
`--smoke` selects 32 by 2 by 32 and four steps. The packet starts at z=256
for the normal cases and z=128 for Brewster cases. Normal cases run 560
steps and Brewster cases 1,000, with output every ten steps.

The material occupies x from zero to 230 in half-space cases, or from 156
to 230 in the slab. Its second half-space interface lies at the periodic
wrap. The translating Gaussian packet is approximate and intended for
qualitative field pictures.

The separate `plane` mode retains its exact-lattice measurement configuration.
It prints complex reflection/transmission amplitudes, their Fresnel and lattice
references, and energy drift. The material solver uses harmonic face
permittivities and the generalized-gauge potential equations
[Chew, 2014, section 2](https://arxiv.org/abs/1406.4780).

## Output and rendering

ParaView can open `wave.pvti.series` and `material.pvti`. Field metadata records
component offsets, timestep and separate sample times: E is at the half step
and raw B at the integer step. The renderer uses those times and offsets.

The first line of `energy.csv` records the case parameters. The `time` column
is the electric-field time, `(k - 1/2) dt`, shared by the discrete vector and
aligned field energies. The scalar diagnostic is at `k dt`. Its columns contain:

- `vacuum`, `dielectric`, `total`, `scalar`, `absolute`: the existing discrete
  potential-energy diagnostics; `total` is the vector quantity.
- `field_energy`: the owned-cell quadrature of electromagnetic energy using
  reconstructed E and time-aligned B, with harmonic face permittivities.

The discrete vector quantity uses a product of magnetic fields at adjacent
time levels [Kreiss et al., 2002, section 2; Sakkaplangkul and Bokil, 2021,
section 4.1]. The field-energy column uses squared, time-aligned fields
[Ellingson, *Electromagnetics*, volume 2, section 3.1]. These quantities are
distinct; neither column includes the flux through the open faces.

The renderer needs VTK, NumPy, Matplotlib and FFmpeg; GIF output also needs Pillow.
It writes `movie.mp4`,
a representative `movie.png`, four snapshots in `movie-overview.png`, and
`movie-energy.png`. `--output name.gif` selects GIF instead. One colour scale
is computed over all frames. The energy plot marks the first saved frame
whose outermost field samples reach 1% of that global peak.
For periodic packets this marks a signal at the periodic wrap, not an open face.
The moving time cursor is solid black; the boundary-signal marker is dotted grey.

New data declares `plane=zx` in its metadata. Older saved outputs without that
declaration retain the legacy x-y interpretation; `--plane zx` or `--plane xy`
can select it explicitly. The renderer selects `B_raw` when metadata declares
`polarization=TM`, otherwise `E_raw`; `--field` overrides this choice. Use
`--field B_raw` for legacy TM packet outputs without polarization metadata.

Small regressions and plane checks run with
`ctest -R chdr.interface -LE 'long|reserved'`. The reserved configuration remains
disabled. The movie smoke check runs all eight TE/TM cases for four steps.
