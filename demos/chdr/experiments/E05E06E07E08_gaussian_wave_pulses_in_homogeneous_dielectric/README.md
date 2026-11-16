# Gaussian pulses in a homogeneous dielectric

E05 to E08 repeat E01 to E04 with the Gaussian potential profile
\(p(s)=A_0\exp(-s^2/(2\sigma^2))\), \(\sigma=w/4\).
Their protocols are in [E05.typ](E05.typ) to [E08.typ](E08.typ),
the sections of the manual chapter "Experiment definitions" after E01 to E04.
Each file states only what differs from the corresponding experiment of E01 to E04.

The executable `experiments-gaussian-wave-pulses-in-homogeneous-dielectric` reuses the run paths,
measurements and writers of `../E01E02E03E04_wave_pulses_in_homogeneous_dielectric/`.
Their README describes the inputs, outputs and conventions, which apply here unchanged.
This directory adds only:

- `GaussianProfiles.h`: the Gaussian as the profile of the plane-pulse sampler (E05, E06)
  and as the radial profile of the dipole-pulse sampler (E07, E08),
- `E05.h` and `E06.h`: the two plane-pulse protocols, with the ladder 32 to 1024 of E02,
- `run_experiment.cpp`: the entry point, which selects the samplers and the allowed parameters,
- the inputs, the two study files and the launchers of E07 and E08.

The CSV files keep the measurement names of E01 to E04 (for example `E03.meas01` in an E07 run),
because the definitions are the same. `run.json` names the experiment, the executable and the Gaussian profile.

## E05 and E06

From the project root, build and run the two studies on eight ranks with four threads each:

```sh
bash scripts/build-ippl.sh --worktree main build --target experiments-gaussian-wave-pulses-in-homogeneous-dielectric
for study in E05 E06; do
  OMP_NUM_THREADS=4 OMP_PLACES=threads OMP_PROC_BIND=spread \
    mpiexec -n 8 --map-by slot:PE=4:HWTCPUS --bind-to hwthread \
    build/ippl-cpu/bin/experiments-gaussian-wave-pulses-in-homogeneous-dielectric \
    --input ippl/demos/chdr/experiments/E05E06E07E08_gaussian_wave_pulses_in_homogeneous_dielectric/$study-study.yaml
done
```

The studies write to `build/results/e05-study/` and `build/results/e06-study/`.
Unlike E02, the E06 study needs no reference study, because E06.typ drops the comparison with the periodic run.

`inputs/E05-n32-transverse4.yaml` is the check of E05.typ:
the coarsest level with 4 by 4 transverse cells must give the errors and profiles of the 8 by 8 run to rounding.
It is not part of the study.

## E07 and E08

Run the launchers from the project root with the arguments of `run_e03_study.py`, for example

```sh
python3 ippl/demos/chdr/experiments/E05E06E07E08_gaussian_wave_pulses_in_homogeneous_dielectric/run_e07_study.py \
  --executable build/ippl-cpu/bin/experiments-gaussian-wave-pulses-in-homogeneous-dielectric \
  --output build/results/e07-study --ranks 4 --threads 4 --levels 96 144 216 324 --plane-every 1 1 2 2
python3 ippl/demos/chdr/experiments/E05E06E07E08_gaussian_wave_pulses_in_homogeneous_dielectric/run_e08_study.py \
  --executable build/ippl-cpu/bin/experiments-gaussian-wave-pulses-in-homogeneous-dielectric \
  --output build/results/e08-study --ranks 4 --threads 4 --levels 96 144 216 324 --plane-every 2 3 4 7
```

E07 and E08 reject `profile_sharpness`, so do not pass `--profile-sharpness`.

## Video output

E05 and E06 save the full-domain frames of E01 and E02 at the movie cadence on the levels 32, 128 and 512
(`observations.fields.every` in their inputs). Their transverse mesh has only 8 by 8 cells.
E07 and E08 save no full-domain frames. With `--plane-every` (the observation `plane` in the input),
the run writes the \(z=0\) section of E and R2 B to `plane-z0.csv`, the format of E03's `write_plane`.
In E08, the first row and column of a section interpolate with samples outside the Mur box,
so they are not physical. The sections enter no measurement.
