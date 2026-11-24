# E18, E19 and their S cases. Oblique incidence on a planar dielectric interface

E18 measures the reflection of the TE packet using the transverse electric field.
E19 uses the transverse magnetic field of the TM packet,
including the wavenumber at the Brewster angle.
The packet, initialization and observation times are defined in [E18.typ](E18.typ) and [E19.typ](E19.typ).
The cases [E18-S](E18S.typ) and [E19-S](E19S.typ) use the solver safety factor 0.99.
They also apply the incident and reflected Fourier sums to the continuum field,
at the same native positions and times as the numerical field.

From the project root in the CPU/MPI development environment,
build the executable and run its reference and measurement checks.

```sh
cmake --build build/ippl-cpu --target chdr-experiment-e18-e19 test-e18-e19-reference test-e18-e19-measurements
ctest --test-dir build/ippl-cpu -R '^chdr\.experiments\.e18-e19\.' --output-on-failure
```

To execute the study on the workstation, pass the executable and a new output directory.
The script runs all eight inputs with eight MPI ranks along z and four OpenMP threads per rank.
It saves the executable checksum, the input and the MPI log with the results.

```sh
bash ippl/demos/chdr/experiments/E18E19_oblique_incidence_on_planar_dielectric_interface/run_study.sh \
    build/ippl-cpu/bin/chdr-experiment-e18-e19 build/results/e18-e19-study
```

Append `--s` to select the eight E18-S and E19-S inputs instead.
These inputs leave field frames disabled.
Their final observation is the reconstructed field sample nearest time 28.
The energy measurements are saved every `cells_per_period` steps and at the final step.

```sh
bash ippl/demos/chdr/experiments/E18E19_oblique_incidence_on_planar_dielectric_interface/run_study.sh \
    build/ippl-cpu/bin/chdr-experiment-e18-e19 build/results/e18-e19-s-study --s
```

Each completed input writes a separate directory.

| File | Contents |
| --- | --- |
| `reflection.csv` | Incident and reflected Fourier magnitudes, their ratio, the incidence angle and the Fresnel magnitude. The Brewster wavenumber is the twelfth row. |
| `measurements.csv` | The two side energies and the running maximum electric field on the two outer layers divided by the initial peak in the whole box. Divide each side energy by the sum of both initial energies to obtain its energy fraction. |
| `modes.csv` | The 327 incident modes, their frequencies, transmitted normal wavenumbers, normalized amplitudes and signed reflection coefficients. |
| `profiles.csv` | Numerical and continuum samples on the line with storage indices i = j = 0, at the initial and final reconstructed times. Each row uses the native position of its field component. |
| `fields/` | E18/E19 save five full-domain frames of `E` and `B_R2` at steps 0, one quarter, one half, three quarters and the end of the run. E18-S/E19-S do not write these frames. |
| `run.json` | Resolution, time step, solver and build identity, MPI partition, initial electric peak and continuum reflected energy fraction. |

For E18-S and E19-S, `reflection.csv` also contains `reference_incident_magnitude`,
`reference_reflected_magnitude` and their ratio `window_reference_magnitude`.
The Richardson extrapolations of E18-S.meas05 and E19-S.meas05 use the archived
`reflection.csv` and `measurements.csv` from the three finest levels.
They are formed during analysis, after all four levels have run.

At the interface E19 initializes the normal vector potential from the vacuum side,
as prescribed in its initialization section.
The continuum electric profile uses the dielectric value of the normal electric field at that point,
as prescribed in its field formula.
The profile and frame data carry the half-step time of the reconstructed fields.
