
#import "./preamble.typ": *

#let chapter = [
  #show: equate.with(sub-numbering: false)

  = Units
  Unless stated otherwise, all equations use normalized variables in which the vacuum values are $c = 1 = epsilon_0 = mu_0$.
  The solver works in these units. We can relate these normalized quantities to SI quantities.

  In the following relations, notilde-quantities are SI quantities and tilde-quantities are dimensionless. The constants $c, epsilon_0, mu_0$ in the conversion factors are their SI-values.

  Choose a length scale $L_0$ and a charge-density scale $rho_0$, then write
  $
    vb(x) & = L_0 vb(tilde(x)), quad quad #<units-scale-relation-x> \
        t & = L_0/c tilde(t), quad quad #<units-scale-relation-t> \
      rho & = rho_0 tilde(rho). #<units-scale-relation-rho>
  $
  The following other relations follow (see @units-scale-relations-derivations):
  $
    vb(J) & = c rho_0 vb(tilde(J)), quad quad #<units-scale-relation-J>\
    phi & = (rho_0 L_0^2)/epsilon_0 tilde(phi), quad quad #<units-scale-relation-phi>\
    vb(A) & = (rho_0 L_0^2)/(epsilon_0 c) vb(tilde(A)), quad quad #<units-scale-relation-A>\
    vb(E) & = (rho_0 L_0)/epsilon_0 vb(tilde(E)), quad quad #<units-scale-relation-E>\
    vb(B) & = (rho_0 L_0)/(epsilon_0 c) vb(tilde(B)). #<units-scale-relation-B>
  $

  == Derivations of our units-scale relations <units-scale-relations-derivations>
  The remaining scales were chosen such that Maxwell's equations and the potential definitions acquire the desired normalized form.

  The derivation for $E_0$ is below, the ones for the other quantities follow by similar arguments. The ones for $D_0$ and $H_0$ follow from $D_0 = epsilon_0 E_0$ and $H_0 = B_0 / mu_0$, using $c^2 epsilon_0 mu_0 = 1$.

  === Derivation of the scale relation for $vb(E)$
  Starting from Gauss's law with the constitutive relation ($vb(D) = epsilon_0 epsilon_r vb(E)$) in SI units
  $
    div (epsilon_0 epsilon_r vb(E)) = rho,
  $ <gauss-law-SI>
  we can substitute $vb(nabla) = L_0^(-1) vb(tilde(nabla)), thick vb(E) = E_0 vb(tilde(E)), thick rho = rho_0 tilde(rho)$ to get
  $
    (epsilon_0 E_0)/(rho_0 L_0) vb(tilde(nabla)) dot (epsilon_r vb(tilde(E))) = tilde(rho),
  $
  and by choosing that the prefactor be equal to one we get $E_0 = (rho_0 L_0)/epsilon_0$, as in @units-scale-relation-E.

  == SI input and energy <units-si-input-energy>
  The input `units: {system: normalized}` retains normalized values. With `system: si`, we must also supply `length_m` ($L_0$, in metres) and `charge_density_C_m3` ($rho_0$, in coulombs per cubic metre). Both scales and the resulting conversion factors must be finite and positive. `UnitSystem` applies the relations above at the input and output while the solver still uses normalized variables. The current presets do not expose dimensional geometry or initial-field parameters. The scales assign a physical interpretation to their normalized values. In SI mode, `run.duration` is in seconds and must equal an integer number of the actual numerical time steps, within roundoff. The time step is not adjusted to fit the duration. The implementation uses the CODATA-2022 values of $c$ and $epsilon_0$ @mohr2025Codata2022[Table XXXIII].

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
