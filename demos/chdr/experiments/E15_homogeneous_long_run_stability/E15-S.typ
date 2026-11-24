// Standalone: typst compile E15-S.typ
#set document(title: "E15-S. One million steps at the dielectric Courant limit")
#set page(paper: "a4", margin: (x: 22mm, y: 18mm))
#set text(font: "Libertinus Serif", size: 11pt, lang: "en")
#set par(leading: 0.55em)
#show math.equation: set text(font: "New Computer Modern Math")

= E15-S. One million steps at the dielectric Courant limit

How much phase lag remains after one million steps at the Courant number for which
the E15 wave travels at its continuum speed?
Repeat the $N=16$ run once in the same periodic unit cube with $epsilon_r=2.13$,
the same diagonal mode $bold(k)=2pi(1,1,1)$ and the same staggered O1VS update and R2 reconstruction.
For $h=1/N$ and $v=1/sqrt(epsilon_r)$, set
$ C = (v h_t)/h = 1/sqrt(3), quad h_t = h sqrt(epsilon_r/3). $
For this mode the dispersion relation gives $omega_d=omega$ at this step
(deduction from @fallahi2020mithra20fullwavesimulation[Eqs. (3.17), (3.20)],
as in _Numerical dispersion of O1VS with vacuum and dielectric_, 24 November 2026, §2.1).
The shortest grid mode can grow linearly at this endpoint.
The current O1VS time-step check uses the vacuum bound.
Select the dielectric endpoint specifically for this run before loading the E15 histories
at their new time levels and discrete frequency.

Use eight MPI ranks along $z$ and four threads per rank.
Record $n=0,16,dots,10^6$ at $t=(n-1/2)h_t$
(one sample per continuum period).
Use the continuum reference, native field locations and definitions of *E15.meas01–03*
in `E15.typ`, “Measured quantities”.
Write the complex electric and magnetic mode coefficients and the residual to `results.json`
in the E15 format. Record the case `E15-S` and the effective step in `run.json`.
Archive both files and the input under `runs/<date>/n16-million/` in the E15-S task directory.
Compare the measured phase lag, amplitude errors, amplitude and energy drifts, and residual
with the archived E15 million-step run, using each run’s own final $t/T$.

#show bibliography: set text(size: 9pt)
#bibliography("E15-references.bib", style: "ieee", title: none)
