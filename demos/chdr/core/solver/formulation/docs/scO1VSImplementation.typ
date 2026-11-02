

#import "../../../../docs/manual/preamble.typ": *

== The O1-VS implementation
This section describes the implemented option O1 (Chew's gauge) on the virtually staggered (staggered spatial convention, see TODO) spatial convention and the R2 time-convention (see TODO).

=== The virtually staggered spatial storage convention
See TODO:aSubchapterOf`demos/chdr/core/base/levelContract.typ`.

=== Derivation (operators)


=== Implemented update
The following update equations are used when TODO:CodeOfExactSolverConfiguration is selected. They assume TODO.
Here the upper index on $vb(A)$ and on $phi$ denotes the actual time level TODO:(In contradiction to the notation I used in other sections?, I should standardize the notation), e.g. an on $phi$ the upper index $n + 3/2$ denotes that after the step $n$, the new stored $phi$ value is at $t_(n + 3/2)$.
$
  A_d^(n+1) &= 2 A_d^n - A_d^(n-1) + h_t^2/epsilon_d [L_vb(h)^"Standard" A_d^n - (C_h A^n)_d + J_d^n]\
  phi^(n+3/2) &= 2 phi^(n+1/2) - phi^(n-1/2) + h_t^2/chi [ sum_d D_d^- (epsilon_d D_d^+ phi^(n+1/2)) + rho^(n+1/2) ]\
$
where TODO:(Missing a whole lot of derivations that lead to these operators, I should describe that in the Derivation)
$
  (C_h A)_d & = D_d^+ gamma - epsilon_d D_d^+ (g / chi) #<eq-O1VS-ChAd>\
  gamma & = sum_e D_e^- A_e #<eq-O1VS-gamma>\
  g & = sum_e D_e^- (epsilon_e A_e) #<eq-O1VS-g>\
  epsilon_d & = cal(H)(epsilon_bold(i), epsilon_(bold(i) + bold(e)_d)) #<eq-O1VS-epsd>
$

#todo[
  How do these equations @eq-O1VS-ChAd, @eq-O1VS-gamma, @eq-O1VS-g, @eq-O1VS-epsd relate to the following equations (on the collocated layout)?:
  $
    g(bold(i)) &= sum_e 1/h_e ( epsilon_e (bold(i)) A_e (bold(i)) - epsilon_e (bold(i) - bold(e)_e) A_e (bold(i) - bold(e)_e) )\
    "etc" &= "etc"
  $
]
