# E15. Long-run propagation of a periodic plane wave

[E15.typ](E15.typ) defines the wave, the initialization and the three measurements.
`experiment-e15` uses the O1VS solver with staggered R2 fields.
Build the target in the configured IPPL build directory.

```sh
cmake --build build/ippl-cpu --target experiment-e15
```

Run from a new study directory. Set `exe` to the absolute path of the built
executable and `inputs` to the absolute path of this experiment's `inputs/`.
On the 32-hardware-thread workstation used for E15, the four launches were

```sh
export OMP_NUM_THREADS=4 OMP_PROC_BIND=close OMP_PLACES=threads
for run in n16 n32 n64 n16-million; do
  mpiexec -n 8 --map-by ppr:8:node:PE=4:HWTCPUS --bind-to hwthread \
    "$exe" --input "$inputs/$run.yaml" || break
done
```

Each input writes its own directory with `results.json` (time histories)
and `run.json` (input, solver, build and MPI layout).
The middle mesh also writes four `slice-<step>.json` files with the native
electric field component `Ex`, its coordinates and its R2 sample time.
The million-step run starts independently from the same wave as `n16`.
