# Dielectric-side propagating controls

These two TM rows use the existing fixed-Cartesian H_y Fresnel and aligned-lattice
references in [InterfaceOracles.h](../../../tests/apps/planar-interface/InterfaceOracles.h),
evaluated through the C++ `interface-reference` utility. The coefficient, phase
and flux conventions are derived in [references.org](../references.org), using
Ellingson sections 5.6-5.7 and the documented discrete-interface recurrence.
The incident medium is the dielectric and the transmitted medium is vacuum.

| Input | Central incidence | Dielectric permittivity | Dielectric resolution | Cell spacing |
|---|---:|---:|---|---|
| `dielectric-tm-preset-18.yaml` | 18 degrees | 3.4 | Recorded preset value, about 20.0861 cells per wavelength | (1, 0.85, 1.25) |
| `dielectric-tm-30-n20.yaml` | 30 degrees | 2.13 | Requested 20; use the resolved value in `run.json` | (1, 2, 2) |

The first input selects `preset: dielectric` from
[PlanePresets.h](../../../apps/planar-interface/PlanePresets.h). Its physical
parameters, 1828 by 2 by 52 cells, 78 initial modes, and 3789 steps remain those
of the recorded regression. It requests 41 measured frequencies through the
ordinary application path. The second input resolves a new row from its
explicit physical parameters. It uses the same 41-bin band controls.
The differing materials, shapes and resolutions make these independent cases;
they are not two points of a matched angle or convergence scan.

Both central angles are below the corresponding critical angle for transmission
from the dielectric into vacuum; the criterion is Ellingson's section 5.11,
Eq. (5.260), in [Electromagnetics, Vol. 2](https://doi.org/10.21061/electromagnetics-vol-2).
Every measured frequency has its own propagation checks and reference angle,
because the pulse has fixed tangential wavenumber
[Oskooi and Johnson, section 4.5](https://arxiv.org/abs/1301.5366v1).
Invalid bins remain in the tables with their reasons. The analysis selects only
bins with propagating continuum and lattice roots on both sides.

From this directory, run the explicit inputs within an assigned CPU budget:

```sh
python3 ../scan.py inputs/dielectric-tm-preset-18.yaml inputs/dielectric-tm-30-n20.yaml \
  --binary /path/to/immutable/planar-interface --output-root /path/to/build/dielectric-side \
  --threads 1 --jobs 1 --max-cores 1
python3 plot.py --scan-manifest COMPLETED_SCAN.json \
  --reference-binary /path/to/interface-reference --output dielectric-analysis
```

Alternatively, give `plot.py` explicit completed run directories. It reuses the
common per-bin comparison and band figures, with separate output directories
for each run. `central.csv` retains the material, actual and requested resolution,
cell spacings and preset name; `analysis.json` preserves each run/reference
identity. It produces no shared angle curve or observed convergence order.
The reflection and transmission power columns use the documented continuum flux
weight; the solver's conserved diagnostic remains a separate quantity.

Fields are saved every 1500 steps and at the final state. The preset therefore
saves four frames; the new row's count follows its resolved step count. Energy
diagnostics are saved every ten steps. Completing these inputs establishes no
new error tolerance: inspect the reference differences and run-specific
box/gate qualifications before making an accuracy claim.

These are the below-critical controls associated with the public
[RC-3 definition](../references.org#reference-case-labels). Total internal
reflection, evanescent-field validation, the Goos-Hänchen shift and the other
deferred reference families remain outside this experiment. Existing core
oracle formulas support these two comparisons; no new scattering formula is
introduced by the inputs or analysis entry point.
