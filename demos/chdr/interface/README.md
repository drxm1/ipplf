# Periodic dielectric interface examples

Build `interface-demo` with `IPPL_ENABLE_CHDR=ON`. All lengths are in cell units;
the speed of light, vacuum permittivity and vacuum permeability are one.

```sh
OMP_NUM_THREADS=1 interface-demo plane
OMP_NUM_THREADS=1 interface-demo half-space --output frames-half-space
OMP_NUM_THREADS=1 interface-demo slab --output frames-slab
OMP_NUM_THREADS=1 interface-demo brewster-te --output frames-brewster-te
OMP_NUM_THREADS=1 interface-demo brewster-tm --output frames-brewster-tm
python3 render_movie.py frames-half-space
python3 render_movie.py frames-brewster-tm --field B_raw
```

Output directories must be new or empty. Add `--smoke` for a 32 by 32 by 4 box
and four steps. Full packets use 512 by 512 by 4 cells and a frame every ten steps.
Normal cases run 560 steps from y=256; Brewster cases run 1,000 steps from y=128.
The permittivity is 2.13. The material occupies x from zero
to 230 in the half-space cases, or from 156 to 230 in the slab. Requested faces
coincide with realized cell boundaries. The opposite half-space face lies at
the periodic wrap. These finite windows precede the returning periodic packet.

The packet examples show fields qualitatively; their translating Gaussian start
is approximate. The `plane` case uses exact lattice modes and prints complex
electric reflection and transmission, their Fresnel and lattice values, and
energy drift. The TM examples use magnetic amplitude conventions. The solver
uses harmonic face permittivities and the generalized-gauge potential equations
[Chew, 2014, section 2](https://arxiv.org/abs/1406.4780).

ParaView can open `wave.pvti.series` and `material.pvti`. Field data records the
component offsets from cell centres. E is at the half step and raw B at the
integer step; their exact times are stored separately. No interpolation occurs.
The first line of `energy.csv` records the case and its run parameters. The table
records the conserved vector-potential quadratic form. Interface faces and their
in-plane edges count as dielectric; node terms follow the cell material.
The invariant's magnetic term couples adjacent time levels
[Kreiss et al., 2002, section 2](https://doi.org/10.1137/S0036142901397435);
it is not the square of the displayed raw magnetic field.

The renderer needs VTK, NumPy, Matplotlib and Pillow. It reads saved frames only,
uses no OpenGL, and writes a GIF, its final frame and an energy figure. Figures
are labelled preliminary. The small CTest regressions and plane checks run with
`ctest -R chdr.interface -LE long`; the reserved two-rank configuration is disabled.
