# Larger-period fallbacks

These inputs double the recorded baseline period multiplier and tangent index.
They are unexecuted alternatives to the selected oblique isolation controls,
not additional validation evidence. Use a fallback if the smaller period fails
the application's isolation checks or if its measured variation is appreciable
relative to the continuum discrepancy. Keep the actual execution record if a
fallback is subsequently run.

The selected TE 17-degree control uses `(M, j) = (942, 111)` instead of the
baseline `(628, 74)`; the selected TM 33-degree control uses `(476, 60)` instead
of `(357, 45)`. Both preserve the exact rational value of `j/M`. Their integer
periods remain subject to the C++ setup checks. The corresponding fallback
pairs are `(1256, 148)` and `(714, 90)`.
