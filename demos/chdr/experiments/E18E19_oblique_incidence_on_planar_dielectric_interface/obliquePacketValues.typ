// Values that E18.typ and E19.typ compute from the parameters of the oblique packet.
// Lengths are in units of the tangential period L_x, and c_vac = 1.

#let eps-dielectric = 2.13
#let n-dielectric = calc.sqrt(eps-dielectric)

#let kx = 2 * calc.pi
#let q-center = kx // central angle of incidence 45 degrees
#let q-half-width = 0.85 * kx
#let sharpness = 8
#let train-period = 192 // Lambda = 4 L_z
#let delta-q = 2 * calc.pi / train-period
#let half-count = int(calc.round(q-half-width / delta-q))

// Spectral weight p_beta(u) of E03, used in the wavenumber: zero for |u| >= 1.
#let spectral-weight(u) = if calc.abs(u) >= 1 { 0 } else {
  calc.exp(sharpness * (1 - 1 / (1 - u * u)))
}

// One entry per plane wave: normal wavenumbers in both media, frequency and energy weight.
#let plane-wave(q) = {
  let omega = calc.sqrt(kx * kx + q * q)
  let amplitude = spectral-weight((q - q-center) / q-half-width)
  (
    q: q,
    omega: omega,
    q-second: calc.sqrt(eps-dielectric * omega * omega - kx * kx),
    amplitude: amplitude,
    energy: amplitude * amplitude * omega * omega,
  )
}
#let modes = range(-half-count, half-count + 1).map(m => plane-wave(q-center + m * delta-q))

// Coefficients of one plane wave incident from vacuum.
#let gamma-te(mode) = (mode.q - mode.q-second) / (mode.q + mode.q-second)
#let r-tm(mode) = (eps-dielectric * mode.q - mode.q-second) / (eps-dielectric * mode.q + mode.q-second)
#let angle-of(mode) = calc.atan(kx / mode.q)

// Energy-weighted mean of a quantity over the plane waves of the packet.
#let energy-mean(quantity) = (
  modes.map(mode => mode.energy * quantity(mode)).sum() / modes.map(mode => mode.energy).sum()
)
#let reflected-fraction-te = energy-mean(mode => calc.pow(gamma-te(mode), 2))
#let reflected-fraction-tm = energy-mean(mode => calc.pow(r-tm(mode), 2))

// The normal wavenumbers of the measurement: eleven equally spaced values and the Brewster value.
#let brewster-wave = plane-wave(kx / n-dielectric)
#let measured-waves = range(-5, 6).map(j => plane-wave(q-center + j * q-half-width / 10)) + (brewster-wave,)

#let rounded(value, digits: 4) = str(calc.round(value, digits: digits))
