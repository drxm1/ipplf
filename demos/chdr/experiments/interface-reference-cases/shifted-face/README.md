# Shifted material-interface comparison

The inputs use normal TE illumination at resolutions 20, 40 and 80, with
`cell_size` 1, 0.5 and 0.25 and fixed shape `(1, 2, 2)`. The common aligned interface
is 500. For each mesh, `face_position` is that physical coordinate plus zero,
one quarter, or one half of the x cell width. The phase origin is the requested
material interface. The represented material jump at a mesh cell face must also
be present in the output record. The directory and parameter names retain `face`
for compatibility; it denotes the material interface in those names.

These inputs require the plane application's physical interface-position and cell-size
support. Zero-shift input parameters match the normal ladder; completed matching
normal runs may be reused as the zero-shift data when their provenance is retained.
The half-shift fine-mesh control enlarges separation and timing qualification.

Run the explicit YAML files with `../scan.py`, then analyse completed manifests:

```sh
python3 plot.py --scan-manifest RUN_SCAN.json --scan-manifest CONTROL_SCAN.json \
  --reference-binary /path/to/build/bin/interface-reference \
  --base-face 500 --output shifted-analysis
```

Alternatively supply the exact completed run directories, for example when
reusing normal-ladder zero-shift data. The analysis requests continuum references
only: the unshifted aligned-lattice oracle is not a reference at an arbitrary
requested phase origin. It retains every bin, its validity, both interface positions,
complex and magnitude errors, observed orders and control comparisons.

[references.org](../references.org), “Shifted interface and phase origin”, derives
why interface-registration error can change complex coefficients without changing
their magnitudes. The plots therefore keep those errors separate and do not
impose a universal first- or second-order slope.
