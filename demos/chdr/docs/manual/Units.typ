
#import "./preamble.typ": *

#let chapter = [
  #show: equate.with(sub-numbering: false)

  = Units
  Unless stated otherwise, all equations use normalized variables in which the vacuum values are $c = 1 = epsilon_0 = mu_0$.
  The solver works in these units. We can relate these normalized quantities to SI quantities.

  Choose a length scale $L_0$ and a charge-density scale $rho_0$, then write
  $
    vb(x) = L_0 vb(tilde(x)), quad quad
    t = L_0/c tilde(t), quad quad
    rho = rho_0 tilde(rho).
  $
  The following other relations follow from this:
  $
    vb(J) = c rho_0 vb(tilde(J)), quad quad
    phi = (rho_0 L_0^2)/epsilon_0 tilde(phi), quad quad
    vb(A) = (rho_0 L_0^2)/(epsilon_0 c) vb(tilde(A)), quad quad
    vb(E) = (rho_0 L_0)/epsilon_0 vb(tilde(E)), quad quad
    vb(B) = (rho_0 L_0)/(epsilon_0 c) vb(tilde(B)).
  $



]

#standalone(title: "Units", chapter)
