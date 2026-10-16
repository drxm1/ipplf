# Potential halo-order regression

`HaloOrder` checks the affine-time solution on all-absorbing Standard and
NonStandard solvers. Its serial witness uses `MPI_COMM_SELF`; the world
communicator uses two ranks split along z. `ShippedOrderSolver` restores
the old exchange-before-Mur sequence as the negative control.

The pulse checks use 16x16x48 all-absorbing boxes for both stencils, with
spacing (1,1,1) for Standard and (1,1,0.5) for NonStandard. Both use a
time step of 0.5. A Gaussian
of width 2.5 cells, centered at x=4,z=24 and constant in y, intersects the
side face near the z seam. It compares every owned sample over 256 steps
against a serial witness and reports exact agreement of the binary64 bits.

`OpenBoundaryRegression` checks both stencils for 256 steps. Standard uses
24x24x24 cells with spacing (1,1,1); NonStandard uses 16x16x48 cells with
spacing (1,1,0.5). Both use a time step of 0.5. Periodic cases use the world
communicator; all-absorbing cases use a serial communicator. The test runs
with both one and two world ranks.

The comparison uses deterministic binary64 histories generated in C++
with unsigned integer arithmetic.
The test compares all three potential histories before each time shift,
then E, raw B and the mean of the current and previous B after the shift.
Field comparisons exclude only the outermost layer of absorbing faces.
Each sample must match the old-order control exactly and remain finite.
