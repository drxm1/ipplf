#import "./preamble.typ": *

#let section = [
  #show: equate.with(sub-numbering: false)

  == SI input and energy <units-si-input-energy>
  The input `units: {system: normalized}` retains normalized values. With `system: si`, we must also supply `length_m` ($L_0$, in metres) and `charge_density_C_m3` ($rho_0$, in coulombs per cubic metre). Both scales and the resulting conversion factors must be finite and positive. `UnitSystem` applies the scale relations of @units at the input and output while the solver still uses normalized variables. The current presets do not expose dimensional geometry or initial-field parameters. The scales assign a physical interpretation to their normalized values. In SI mode, `run.duration` is in seconds and must equal an integer number of the actual numerical time steps, within roundoff. The time step is not adjusted to fit the duration. The implementation uses the CODATA-2022 values of $c$ and $epsilon_0$ @mohr2025Codata2022[Table XXXIII].

  For a linear, nondispersive dielectric with $mu_r = 1$, substituting @units-scale-relation-x, @units-scale-relation-E and @units-scale-relation-B into the field energy gives the following deduction from @fitzpatrick2014ClassicalElectromagnetism[Eqs. (1.109), (4.70)]:
  $
    U & = 1/2 integral_V (epsilon_0 epsilon_r abs(vb(E))^2 + abs(vb(B))^2 / mu_0) thin dd(V) = U_0 tilde(U), \
    U_0 & = epsilon_0 E_0^2 L_0^3 = (rho_0^2 L_0^5) / epsilon_0.
  $ <units-scale-relation-energy>
  Here $tilde(U)$ is the same field-energy integral in normalized variables; the magnetic term has the same prefactor since $B_0 = E_0/c$ and $epsilon_0 mu_0 c^2 = 1$.

  The output contains field-energy estimates and diagnostics of the potential update. The discrete potential energies computed by `FieldEnergy.h` are defined in @sec-discrete-potential-energies. The table lists the quantities written by `FieldRunOutput` to `energy.csv` and by `chdr-vacuum-wave` to standard output. Here, dimensionless means that the number is expressed in normalized solver units. Multiplication by $U_0$ converts a field-energy estimate to joules, and the following table specifies where the programs apply this conversion.

  #table(
    columns: (1.1fr, 2fr, 1fr),
    inset: 5pt,
    stroke: 0.4pt,
    table.header([*Columns*], [*Quantity*], [*Units*]),
    table.cell(colspan: 3)[*`FieldRunOutput`: `energy.csv`*],
    [`vacuum`, `dielectric`, `total`],
    [The vacuum and dielectric contributions to the discrete wave energy $W_A$, and their sum. $W_A$ is defined in @eq-diag-vector-energy.],
    [Dimensionless in both output modes.],
    [`scalar`],
    [The discrete scalar-potential energy $W_phi$ defined in @eq-diag-scalar-energy.],
    [Dimensionless in both output modes.],
    [`absolute`],
    [Sum of the magnitudes of the individual terms contributing to `total`, before positive and negative terms cancel.],
    [Dimensionless in both output modes.],
    [`field_energy`],
    [Estimate of $U$ using electric and magnetic fields at the same sample time. The magnetic field is obtained by averaging its two adjacent time levels.],
    [Dimensionless in both output modes.],
    [`field_energy_J`],
    [`field_energy` multiplied by $U_0$.],
    [Joules. This column exists only in SI mode.],
    table.cell(colspan: 3)[*`chdr-vacuum-wave`: standard output*],
    [`energyInvariant`],
    [Potential diagnostic for the vacuum update.],
    [Dimensionless in both output modes.],
    [`energyE`, `energyB`, `energy`],
    [Electric-field contribution, magnetic-field contribution, and their sum. The first two use their separate raw-field sample times.],
    [Dimensionless in normalized mode; joules in SI mode.],
  )
]
