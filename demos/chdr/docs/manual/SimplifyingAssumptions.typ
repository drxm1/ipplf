#import "./preamble.typ": *

#let chapter = [
  #show: equate.with(sub-numbering: false)



  = Simplifying material assumptions
  We consider the following material model and usually a domain of two piecewise homogeneous materials (dielectric $cal(D)$ and vacuum), where normalized units are used ($c = epsilon_0 = mu_0 = 1$):
  $
    vb(D) &= epsilon_r vb(E) #<constitutive-D> \
    vb(B) &= vb(H) #<constitutive-H> \
    epsilon_r (cancel(omega), cancel(t), vb(x)) &= cases(epsilon_cal(D) "if" vb(x) in cal(D), epsilon_(not cal(D)) "else (usually vacuum)")
  $

  which contains the following assumptions: Linear material response, Isotropic ($epsilon_cal(D)$ a scalar), spatially local response. Nondispersive, lossless, nonconducting and nonmagnetic, stationary material.

  We impose no additional free surface charge or current at the dielectric interface, but bound polarization charge and current can occur as a response to the moving charge's field.
]

#standalone(title: "Simplifying Assumptions", chapter)
