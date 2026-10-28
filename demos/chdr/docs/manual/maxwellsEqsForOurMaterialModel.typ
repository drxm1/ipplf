#import "./preamble.typ": *

#let chapter = [
  #show: equate.with(sub-numbering: false)

  = Maxwell's equations for our material model

  Substitute @constitutive-D and @constitutive-H into Maxwell's equations @chew2014generalizedgauge[Eqs. 1-4].
  In our normalized units we get
  #pad(left: 1em)[
    #show math.equation: set align(left)
    $
      & #text[Faraday’s law:] #h(2em)
        & partial_t vb(B) &= - curl vb(E), #<eq-MW1-faraday-in-our-mat-model> \
      & #text[Ampère’s law:] #h(2em)
        & epsilon_r partial_t vb(E) &= curl vb(B) - vb(J)_f, #<eq-MW2-ampere-in-our-mat-model> \
      & #text[Gauss’s law for magnetism:] #h(2em)
        & div vb(B) &= 0, #<eq-MW3-gauss-magnetism-in-our-mat-model> \
      & #text[Gauss’s law:] #h(2em)
        & div (epsilon_r vb(E)) &= rho_f, #<eq-MW4-gauss-in-our-mat-model>
    $
  ]
  where we could only move $epsilon_r$ outside of $partial_t$ because the material is assumet to be stationary, so $partial_t epsilon_r = 0$. Note however that $epsilon_r$ is allowed to vary in space, which is why we keep $epsilon_r$ inside the divergence in @eq-MW4-gauss-in-our-mat-model.

  The quantities $rho_f$ and $vb(J)_f$ are the _free_ charge and current densities. The polarization response is already represented by @constitutive-D and @constitutive-H /*TODO: See Haus-Melcher 6.2 pp. 7-9: Setting the free sources to zero still allows waves from initial data and a polarized dielectric*/.

  == Charge continuity
  Taking the divergence of Ampère's law (see @eq-MW2-ampere-in-our-mat-model) and inserting Gauss's law (@eq-MW4-gauss-in-our-mat-model) we get the charge continuity equation
  $
    0 = partial_t div (epsilon_r vb(E)) + div vb(J)_f eq.triple partial_t rho_f + div vb(J)_f
  $

]

#standalone(title: "Maxwell's equations for our material model", chapter)
