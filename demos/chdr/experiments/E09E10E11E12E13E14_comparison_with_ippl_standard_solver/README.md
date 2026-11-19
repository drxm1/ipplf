# Comparison with IPPL's Standard solver

E09 to E14 repeat E05 to E08 with four solver configurations on the same meshes:

| Configuration | Formulation | Placement | Reconstruction |
|---|---|---|---|
| `O1VS` | O1VS | staggered | R2 |
| `Homogeneous` | Homogeneous | staggered | R2 |
| `Standard` | Homogeneous | collocated | R1 |
| `StandardR2Prime` | Homogeneous | collocated | R2Prime |

In the `Standard` configuration `ChdrSolver` runs the inherited `step()`, `timeShift()` and
`evaluate_EB()` of IPPL's `StandardFDTDSolver`. At \(\varepsilon_r=1\) it gives the same
potentials and fields as a bare `StandardFDTDSolver` with the same loaded histories, bit for bit.
`StandardR2Prime` runs the same update and the same electric kernel of `evaluate_EB()`, but it loads
\(\phi\) half a step before \(\mathbf A\). Its \(\mathbf E\) then has the single time \(t_{n-1/2}\),
and `ChdrSolver` replaces the \(\mathbf B\) of `evaluate_EB()` by the average of the central curls of
both potential slots, also at \(t_{n-1/2}\).

| Experiment | Repeats | Medium |
|---|---|---|
| E09 | E05 (periodic plane pulse), 4 by 4 transverse cells | vacuum |
| E10 | E06 (plane pulse through Mur faces), 4 by 4 transverse cells | vacuum |
| E11 | E07 (periodic dipole pulse) | vacuum |
| E12 | E08 (dipole pulse through all-face Mur) | vacuum |
| E13 | E07 | \(\varepsilon_r=2.13\) |
| E14 | E08 | \(\varepsilon_r=2.13\) |

The protocols are in [E09.typ](E09.typ) to [E14.typ](E14.typ), the sections of the manual chapter
"Experiment definitions" after E05 to E08. The executable
`experiments-comparison-with-ippl-standard-solver` reuses the run paths, measurements and writers of
`../E01E02E03E04_wave_pulses_in_homogeneous_dielectric/` and the Gaussian profiles of
`../E05E06E07E08_gaussian_wave_pulses_in_homogeneous_dielectric/`. Every input selects its
configuration with `parameters.solver_configuration`. The CSV files keep the measurement names of
E01 to E04. For `Standard`, R1's electric field is compared at \(t_{n-1/2}\). For both collocated
configurations the constraint residuals, the monitor flux and the periodic control are not
measured, because they use the operators of the staggered layout.
The frames of E09 and E10 keep the array name `B_R2` for the reconstructed magnetic field of every
configuration; in `Standard` runs it holds R1's \(\mathbf B^n\), and `B_R2_time` gives its time.
E09 and E10 save frames with the cadences of E05 and E06.
E11 to E14 write no full-volume frames: with \(\phi\neq0\), R1's \(\mathbf E\) has no single time label,
so a `Standard` run with `--fields-every` stops at its first frame.
The \(z=0\) sections of E11 and E13 have one time column, that of \(\mathbf E\); the errors of
\(\mathbf B\) use the time of \(\mathbf B\).

## E09 and E10

From the project root, run one study per configuration on eight ranks with four threads each:

```sh
bash scripts/build-ippl.sh --worktree main build --target experiments-comparison-with-ippl-standard-solver
for study in E09 E10; do for configuration in O1VS Homogeneous Standard; do
  OMP_NUM_THREADS=4 OMP_PLACES=threads OMP_PROC_BIND=spread \
    mpiexec -n 8 --map-by slot:PE=4:HWTCPUS --bind-to hwthread \
    build/ippl-cpu/bin/experiments-comparison-with-ippl-standard-solver \
    --input ippl/demos/chdr/experiments/E09E10E11E12E13E14_comparison_with_ippl_standard_solver/$study-$configuration-study.yaml
done; done
```

The studies write to `build/results/e09-study-O1VS/` and so on.

## E11 to E14

`run_dipole_study.py` runs one experiment and configuration with the arguments of
`run_e03_study.py`, for example

```sh
python3 ippl/demos/chdr/experiments/E09E10E11E12E13E14_comparison_with_ippl_standard_solver/run_dipole_study.py E11 \
  --executable build/ippl-cpu/bin/experiments-comparison-with-ippl-standard-solver \
  --output build/results/e11-study-Standard --ranks 4 --threads 4 --levels 96 144 216 324 \
  --solver-configuration Standard
```

For `O1VS` and `Homogeneous` in E11 and E13 the script always writes the \(z=0\) section at the
initial state and the final field target on every level (E11.meas01, E13.meas01). The collocated
configurations write it only with `--write-plane`. E12 and E14 need no sections.

The directory holds study inputs of E09 and E10 for `O1VS`, `Homogeneous` and `Standard`.
The `StandardR2Prime` studies of E09 and E10 use the `Standard` inputs with
`solver_configuration: StandardR2Prime`.
