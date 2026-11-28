// Guided modes of the dielectric slab of E22, computed from the mode conditions of E22.typ.
// Lengths are in units of the period L_x of the box along x, and c_vac = 1.

#let epsilon-slab = 2.13
#let slab-thickness = 1.0 // d
#let box-length = 16.0 // L_z
#let gap = box-length - slab-thickness // W
#let k-x = 2 * calc.pi
#let critical-angle = calc.asin(1 / calc.sqrt(epsilon-slab))

// The left side minus the right side of the mode condition, multiplied by cos(q d/2) for an
// even mode and by sin(q d/2) for an odd mode. This form has no pole.
#let mode-condition(omega, polarization, parity) = {
  let q = calc.sqrt(epsilon-slab * omega * omega - k-x * k-x)
  let kappa = calc.sqrt(k-x * k-x - omega * omega)
  let weight = if polarization == "TE" { q } else { q / epsilon-slab }
  let (sine, cosine) = (calc.sin(q * slab-thickness / 2), calc.cos(q * slab-thickness / 2))
  let hyperbolic = calc.tanh(kappa * gap / 2)
  if parity == "even" { weight * sine - kappa * hyperbolic * cosine } else {
    -weight * cosine - kappa / hyperbolic * sine
  }
}

#let bisect(f, low, high) = {
  let (a, b) = (low, high)
  for _ in range(60) {
    let middle = (a + b) / 2
    if f(a) * f(middle) <= 0 { b = middle } else { a = middle }
  }
  (a + b) / 2
}

// All frequencies of one polarization and parity between k_x/n and k_x.
#let mode-frequencies(polarization, parity) = {
  let low = k-x / calc.sqrt(epsilon-slab) * (1 + 1e-9)
  let high = k-x * (1 - 1e-9)
  let intervals = 400
  let f = omega => mode-condition(omega, polarization, parity)
  let roots = ()
  for i in range(intervals) {
    let a = low + (high - low) * i / intervals
    let b = low + (high - low) * (i + 1) / intervals
    if f(a) * f(b) < 0 { roots.push(bisect(f, a, b)) }
  }
  roots
}

// The modes of one polarization in the order of their frequency. The order m is the number
// of zeros of the profile inside the slab.
#let guided-modes(polarization) = {
  let entries = ("even", "odd")
    .map(parity => mode-frequencies(polarization, parity).map(omega => (omega: omega, parity: parity)))
    .flatten()
    .sorted(key: entry => entry.omega)
  entries
    .enumerate()
    .map(((order, entry)) => {
      let omega = entry.omega
      let q = calc.sqrt(epsilon-slab * omega * omega - k-x * k-x)
      let kappa = calc.sqrt(k-x * k-x - omega * omega)
      (
        name: polarization + str(order),
        polarization: polarization,
        order: order,
        parity: entry.parity,
        omega: omega,
        q: q,
        kappa: kappa,
        angle: calc.atan(k-x / q),
        amplification: calc.pow(omega / kappa, 2),
      )
    })
}

// The profile of a mode at the height z above the center of the slab, for |z| <= L_z/2.
#let mode-profile(mode, z) = {
  let (half, distance) = (slab-thickness / 2, calc.abs(z))
  let sign = if mode.parity == "odd" and z < 0 { -1 } else { 1 }
  if distance <= half {
    if mode.parity == "even" { calc.cos(mode.q * z) } else { calc.sin(mode.q * z) }
  } else if mode.parity == "even" {
    calc.cos(mode.q * half) * calc.cosh(mode.kappa * (box-length / 2 - distance)) / calc.cosh(mode.kappa * gap / 2)
  } else {
    sign * calc.sin(mode.q * half) * calc.sinh(mode.kappa * (box-length / 2 - distance)) / calc.sinh(mode.kappa * gap / 2)
  }
}
