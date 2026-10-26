# Planar dielectric interface in a periodic box

Build `planar-interface` with `IPPL_ENABLE_CHDR=ON`. Lengths use cell units;
the speed of light, vacuum permittivity and vacuum permeability are one. Every
case has one planar dielectric of relative permittivity 2.13 with its normal
along x, in a box that is periodic on every axis.

```sh
OMP_NUM_THREADS=1 planar-interface plane
OMP_NUM_THREADS=1 planar-interface half-space --output periodic-half-space
OMP_NUM_THREADS=1 planar-interface slab --output periodic-slab
OMP_NUM_THREADS=1 planar-interface brewster-te --output periodic-brewster-te
OMP_NUM_THREADS=1 planar-interface brewster-tm --output periodic-brewster-tm
python3 ../interface-movies/render_movie.py periodic-half-space
python3 ../interface-movies/render_movie.py periodic-brewster-tm --field B_raw
```

## Plane-pulse reflection measurement

The `plane` mode is the exact-lattice measurement configuration: a quasi
one-dimensional box, an exact lattice plane pulse at normal incidence from
vacuum, and two Fourier probes. It prints the complex reflection and
transmission amplitudes, their Fresnel and exact lattice references, and the
energy drift to standard output; it writes no files. The material solver uses
harmonic face permittivities and the generalized-gauge potential equations
[Chew, 2014, section 2](https://arxiv.org/abs/1406.4780).

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

`ctest -R chdr.interface -E movies -LE 'long|reserved'` runs the plane-pulse
checks at normal incidence against the independent lattice and Fresnel values
(`test_plane_interface.cpp` calls `PlanePulse.h` directly; it does not run the
`planar-interface` executable) and the four-step smoke run of the half-space
case of `planar-interface`. The Brewster check carries the `long` label. The
reserved configuration remains disabled.
