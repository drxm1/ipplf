# Incidence-angle comparison

The inputs cover vacuum-side TE and TM at 0, 15, 30, 45, Brewster, 65 and 75
degrees, with relative permittivity 2.13 and requested dielectric resolution 20.
The Brewster input is 55.58149809343727 degrees, from Ellingson §5.10,
Eq. (5.253); the fixed-magnetic-component convention and derivation are given in
[references.org](../references.org). Cell shape `(1, 2, 2)` is fixed for every
row, including zero incidence. The C++ setup records the frequency adjustment
needed for integer transverse cells and the resulting actual resolution.

Each input asks for 41 band bins and sparse frames. The Brewster and 75-degree
TM controls increase separation, timing sampling and gate padding, with frames
disabled. They test layout/gate sensitivity at important scan locations.

Run `../scan.py` on the explicit inputs and controls, then pass its completed
manifests to the plotting command:

```sh
python3 plot.py --scan-manifest RUN_SCAN.json --scan-manifest CONTROL_SCAN.json \
  --reference-binary /path/to/build/bin/interface-reference --output angle-analysis
```

The output contains per-bin data and band figures, central amplitude and phase
figures, continuum flux estimates, actual resolution labels and matched-control
differences. Phase is omitted when either compared amplitude is below the
recorded phase floor; absolute reflection is retained at Brewster. Connected
central reference values are explicitly labelled “Fresnel samples”. All reference
values and transmission flux weights come from the C++ calculator.

Each frequency bin has its own continuum reference angle, because the pulse has
fixed tangential wavenumber. Invalid bins are retained in tables with reasons and
excluded from primary figures. `R + T` uses the documented continuum flux weights;
it must not be identified with the solver's discrete conserved diagnostic.
