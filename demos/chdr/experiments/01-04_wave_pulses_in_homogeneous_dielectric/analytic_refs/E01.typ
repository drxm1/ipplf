#import "@preview/physica:0.9.8": *
#set math.equation(numbering: "(1)")

= Experiment E01 analytic continuum reference and initialization

== Experiment E01 analytic continuum reference
As the continuum analytic reference for E01 we choose a transverse potential pulse.\

We assume a nonmagnetic, lossless, nondispersive medium with $epsilon_r = 2.13$.
The material is homogeneous and spans the whole mesh domain. There is no free charge or current.

Let
$
  vb(A)(z,t) = vb(hat(x)) a(s(z,t)), quad s(z,t) := z - z_0 - v t, quad phi = 0, quad v = c_"vac"\/sqrt(epsilon_r).
$ <E01-analytic-potentials>

The continuum fields $vb(E)$ and $vb(B)$ can be reconstructed via Chew §1, Eqs. (5)-(6):
$
  vb(E) = -partial_t vb(A) - grad phi, quad vb(B) = curl vb(A)
$

The chain rule for the potential in @E01-analytic-potentials gives $partial_t A_x = -v a'(s)$ and $partial_z A_x = a'(s)$. Only $(curl vb(A))_y = partial_z A_x$ is nonzero in the curl. By substitution we obtain the continuum $vb(E), vb(B)$:
$
  vb(E)(z,t) = vb(hat(x)) v pdv(a(s),s)(z,t), quad vb(B)(z,t) = vb(hat(y)) pdv(a(s),s)(z,t).
$ <E01-analytic-EB>

@E01-analytic-potentials and @E01-analytic-EB hold for all $t$, including the negative time that is used to initialize the previous potential `A_nm1`. We call $a$ a _potential profile_.



=== Checking that the potential solves the homogeneous equations
Since $A_x$ is independent of x, $div vb(A) = 0$, so the homogeneous Lorenz gauge $div vb(A) = -epsilon_r partial_t phi$ TODO:chew[Eq. (11)] is satisfied.

We now verify the continuum wave equations for the scalar and vector potentials
in a homogeneous medium, with zero free charge and current. In the Lorenz
gauge these are
#link("https://arxiv.org/pdf/1406.4780v1")[Chew, §2, Eqs. (12)–(13)]:
$
  nabla^2 phi - epsilon_r partial_t^2 phi &= 0, \
  nabla^2 vb(A) - epsilon_r partial_t^2 vb(A) &= vb(0).
$
Since $A_x (z,t) = a(s(z,t))$ is independent of $x$ and $y$, we obtain $nabla^2 A_x = partial_z^2 A_x$.
If $a in cal(C)^2 ("Dom"(a))$, $partial_z^2 A_x - epsilon_r partial_t^2 A_x = a''(s) - a''(s) = 0$.

Therefore this potential provides a source free Maxwell solution.


== Experiment E01 initialization/sampling method
We choose initialization option A described in `../experiment_definition_draft.md`,
that is we sample the continuum potential history directly:
$
  vb(A)^0 eq.triple_"stagg" vb(A)_"exact" (0), quad vb(A)^(-1) eq.triple_"stagg" vb(A)_"exact" (- h_t)
$

Here $h_t$ is the time step, and the superscripts are time-level indices. Time-level zero is by convention at $t=0$, and level minus one at $t = -h_t$. In this notation we also use half-levels to denote when a discrete expression lives a half step $h_t\/2$ apart due to it being calculated by a difference operator.

A reader might wonder why the triple equality with subscript "stagg", and what this should mean.\
In order to distinguish the discrete samples stored in the virtually staggered `O1-VS` convention, we do the following:
Let $(x_i, y_j, z_k)$ denote a cell centre, and we write the spacial indices explicitly. The indices $(i,j,k)$ denote that a field is stored in the code such that it corresponds to the cell with the corresponding centre. We additionally add $+1/2$ or $-1/2$ to the index of a field component if that field component physically lives a half-cell shifted in the direction that corresponds to the said index. This is necessary because in the spatially staggered implementation some stored field components physically have a position that is offset from the cell-centre due to the finite difference operators being used to calculate certain derivatives.

In the convention used by the `O1-VS` update, $A_x$ is for example located at $(x_i + h_x\/2, y_j, z_k)$.  Similarly the locations of $A_y$ and $A_z$ are shifted in $y$ and $z$. These offsets, in code, are returned by `LevelContract::positionOffset`. Their staggered grid "basis" (TODO: What does this mean?) is Ryu et al., Appendix A, Eqs. (19)-(24) TODO:ryu2016potentialfdtd.

For E01, sampling the only nonzero component therefore means the following, which is the meaning of the above $eq.triple_"stagg"$-notation (we do not use that notation anymore from here on).
$
  (A_x^n)_(i+1/2,j,k)
  = A_(x,"exact")(z_k, n h_t)
  = a(z_k - z_0 - v n h_t), quad n in {0,-1}.
$

In the code, the components are still stored at integer indices $(i,j,k)$, and `writeSamples` in `core/solver/BoxSetup.h` adds the physical offset before the analytic is evaluated.

The magnetic component $B_y$ has a different location:
$(x_i+h_x/2, y_j, z_k+h_z/2)$. Its continuum comparison must use
$z_k+h_z\/2$. Thus the half-cell shift in $z$ must be considered in the magnetic-field
comparison, but not in the initialization of $A_x$.

We must measure the initial field error as part of the experiment,
since the initialized electric and magnetic fiels will only approximate the analytic field.
For example the initialized electric field will be
$
  vb(E)^(-1/2) = - (vb(A)^0 - vb(A)^(-1))/h_t eq.triple^"approx"_"stagg" - partial_t vb(A)
$

== TODO's
Here I describe TODO's that I wish to specify more clearly in/about the above.
=== TODO 1
We must choose a sufficiently smooth pulse and a consistent periodic extension for EO1.
The travelling-profile construction is explained in Schneider, §6.2.
TODO: Does this change our equations @E01-analytic-potentials, @E01-analytic-EB?

#text(fill: blue)[
The form of @E01-analytic-potentials and @E01-analytic-EB does not change.
Periodicity restricts the choice of $a$: for a domain of length $L_z$, we need
$a(s + L_z) = a(s)$, with matching derivatives at the periodic join.
The travelling-profile construction is explained in
#link("https://eecs.wsu.edu/~schneidj/ufdtd/chap6.pdf")[Schneider, §6.2].
The periodic construction below follows by superposition of translated
solutions of the homogeneous linear wave equation.

For example, starting from a localized potential profile $p$, define
$
  a(s) = sum_(m in ZZ) p(s - m L_z).
$
This is periodic by reindexing the sum. We require enough decay and smoothness
to differentiate the sum twice. A smooth compactly supported profile satisfies
this requirement because only finitely many copies contribute locally; a
Gaussian and its derivatives also converge in this construction.
Use this same periodic $a$ in the potential histories and its derivative in
the analytic fields.

The resulting periodic continuum solution is exact for all times, including
after the pulse crosses a periodic boundary. For a smooth compactly supported
profile whose support is shorter than $L_z$, it agrees with an isolated pulse
where only that copy contributes. A periodized Gaussian includes nonzero tails
from the other copies, so it is not exactly equal to an isolated Gaussian.
Stopping before the main pulse reaches a boundary does not remove those tails.

The particular profile, its width and the domain length remain to be chosen.
]
