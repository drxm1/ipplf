# E16 and E17. Normal incidence on a planar dielectric interface

[E16.typ](E16.typ) defines the pulse, the mesh refinement and the measurements.
[E17.typ](E17.typ) exchanges the two media.
Both inputs use the same executable.

From the project root, build the experiment and run the continuum reference check.

```sh
cmake --build build/ippl-cpu --target chdr-experiment-e16-e17 test-e16-e17-reference
ctest --test-dir build/ippl-cpu -R '^chdr\.experiments\.e16-e17\.reference$' --output-on-failure
```

Run the ten inputs from an empty output directory.
Set `exe` and `inputs` to the absolute paths of the executable and this directory's `inputs/`.
On the workstation with 32 hardware threads use

```sh
export OMP_NUM_THREADS=4 OMP_PLACES=threads OMP_PROC_BIND=spread
for experiment in E16 E17; do
  for cells in 128 256 512 1024 2048; do
    mpiexec -n 8 --map-by slot:PE=4:HWTCPUS --bind-to hwthread \
      "$exe" --input "$inputs/$experiment-n$cells.yaml" || exit 1
  done
done
```

Each input writes a directory named after the experiment and mesh level.
`measurements.csv` holds the squared field errors, reference norms squared and energies for each side.
`profiles.csv` holds the numerical and continuum samples along one line.
The three frames in `fields/` contain the full electric field, R2 magnetic field and raw magnetic field.
`run.json` records the input, build, MPI decomposition and effective solver configuration.
