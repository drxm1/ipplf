# Normal-incidence comparison

The inputs use TE illumination from vacuum onto relative permittivity 2.13.
The three resolutions are 20, 40 and 80 dielectric cells per central wavelength.
Their `cell_size` values are 1, 0.5 and 0.25, with fixed shape `(1, 2, 2)` and
the common physical face at 500. Thus normal-incidence central frequency remains
fixed while spatial and temporal resolution change. Definitions and source
derivations are in [references.org](../references.org), especially “Normal
incidence and mesh refinement” and the homogeneous dispersion relation.

These inputs require the plane application's `cell_size`, `face_position`, band
spectra and field-output support. Each records 41 bins and sparse field frames.
The fine-mesh control increases the layer-separation factor, timing sampling and
gate padding. Its difference from the corresponding baseline estimates the
combined layout/gate sensitivity; it is not a rigorous error bound.

Run `../scan.py` on the explicit files under `inputs/`, and on `controls/` for
the isolation check. It supplies fresh output directories. Use the completed
scan manifests printed by that command:

```sh
python3 plot.py --scan-manifest RUN_SCAN.json --scan-manifest CONTROL_SCAN.json \
  --reference-binary /path/to/build/bin/interface-reference --output normal-analysis
```

`plot.py` retains the complete per-bin comparisons, generated C++ references and
their provenance. It also writes central values, two-level observed orders,
matched-control differences, band figures and three actual saved pulse profiles
per baseline run. The display selection stops before the first saved frame with
a normal-edge signal above 1% of the global sampled peak, using eight cells at
each x edge. Its values and selected source files are retained in
`profile-selection.json`. This visual criterion is not a boundary-error bound;
the coefficient-isolation controls remain necessary. The raw magnetic frame time is kept
distinct from averaged probe time. Figures are saved as PNG and PDF.

The continuum comparison is primary. Complex error and magnitude error are
separate quantities. An observed exponent from three meshes does not establish
an asymptotic order if isolation or extraction effects are comparable. Assess
the control table before making that claim; do not infer it from process exit.
