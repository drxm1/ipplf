
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




]

#standalone(title: "Units", chapter)
