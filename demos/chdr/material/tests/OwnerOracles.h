/** @file OwnerOracles.h
 * @brief Expected values for the material owner: the source cell of every allocated cell and
 * the closed forms of the derived coefficient fields.
 * @details Written independently of the owner, so that the halo and derived-field tests
 * compare it with values that share none of its code. Global cell indices run from 0 to
 * N-1 over the domain; the one ghost layer of a rank holds the indices -1 and N on the
 * axes where its block touches the domain boundary. The cell permittivity is sampled at
 * the cell centre, the position of the scalar potential and the charge.
 */
#pragma once
#include <Kokkos_Core.hpp>
#include "Ippl.h"

namespace chdr::material {
    /** @brief Global index of the cell whose value an allocated cell must hold, on one axis.
     * @details On a periodic axis the ghost is the image \f$g \bmod N\f$. On a non-periodic
     * axis a cell outside the domain takes the nearest domain cell,
     * \f$\min(\max(g, 0), N - 1)\f$. A cell inside the domain is its own source; an
     * inter-rank ghost holds the value of the rank that owns it. Deduction: a physical ghost
     * continues the boundary medium outward, as an absorbing boundary assumes for "the
     * homogeneous region that surrounds" the obstacles
     * [mur1981AbsorbingBoundaryFdtd, section II, p. 378]. A zero normal gradient of
     * \f$\varepsilon\f$ is the index clamp; the harmonic face value on the boundary face then
     * equals the boundary cell, so no artificial interface appears. Serves the halo check:
     * every allocated cell equals its source cell bit for bit.
     *
     * sourceIndex() wraps periodic indices and clamps physical-boundary indices.
     */
    KOKKOS_INLINE_FUNCTION int sourceIndex(int global, int cells, bool periodic) {
        if (periodic) {
            const int remainder = global % cells;
            return remainder < 0 ? remainder + cells : remainder;
        }
        return global < 0 ? 0 : global >= cells ? cells - 1 : global;
    }

    /** @brief Source cell of an allocated cell; the axes act independently.
     * @details Faces, edges and corners of the ghost layer combine the per-axis rule, and so
     * do mixed ghosts (periodic or physical on one axis, inter-rank on another), because
     * wrap and clamp act on each coordinate of a product domain separately (Deduction).
     */
    KOKKOS_INLINE_FUNCTION ippl::Vector<int, 3> sourceCell(ippl::Vector<int, 3> global,
                                                           ippl::Vector<int, 3> cells,
                                                           Kokkos::Array<bool, 3> periodic) {
        ippl::Vector<int, 3> source;
        for (unsigned d = 0; d < 3; ++d)
            source[d] = sourceIndex(global[d], cells[d], periodic[d]);
        return source;
    }

    /** @brief Value of the extended permittivity map at any global index.
     * @details \f$\varepsilon_{\mathrm{ext}}(\mathbf g) = \varepsilon(\mathrm{source}(\mathbf
     * g))\f$, the domain map continued by wrap or clamp. The derived closed forms below read this
     * map, so they are defined on every ghost cell and one cell beyond it.
     * @tparam Map Callable `double(ippl::Vector<int, 3>)`, defined on the domain cells.
     */
    template <class Map>
    KOKKOS_INLINE_FUNCTION double extendedValue(const Map& map, ippl::Vector<int, 3> global,
                                                ippl::Vector<int, 3> cells,
                                                Kokkos::Array<bool, 3> periodic) {
        return map(sourceCell(global, cells, periodic));
    }

    /** @brief Centre coefficient \f$Q = \varepsilon^{-2}\f$.
     * @details The generalized gauge weight is \f$\chi = \alpha\varepsilon^2\mu\f$
     * [chew2014generalizedgauge, section 2, Eq. (18), p. 3]. With \f$\mu = 1\f$ the stored
     * field is \f$\chi^{-1}\f$ without \f$\alpha\f$, which is a solver parameter and multiplies
     * exactly in the kernel. Serves the derived-field check at the cell centres.
     */
    KOKKOS_INLINE_FUNCTION double inverseSquare(double eps) {
        return 1 / (eps * eps);
    }

    /** @brief Harmonic mean, exact for equal arguments.
     * @details \f$H(a, b) = 2ab/(a + b)\f$ is the flux-continuity face coefficient when the
     * discontinuity lies midway between two samples
     * [eymard2000FiniteVolume, section 7.1, Eqs. (7.3)-(7.5), Example 7.1, p. 22]. The
     * equal-argument branch returns \f$a\f$ itself: \f$2a\f$ and \f$a + a\f$ are exact in
     * floating point but \f$aa\f$ is not, so \f$\mathrm{fl}(2aa/(a + a))\f$ misses \f$a\f$ by one
     * unit in the last place for mantissas in \f$[\sqrt2, 2)\f$ (Deduction), and a homogeneous
     * region must give face values equal to its cells.
     */
    KOKKOS_INLINE_FUNCTION double harmonicMean(double a, double b) {
        return a == b ? a : 2 * a * b / (a + b);
    }

    /** @brief Face coefficient \f$F_d(\mathbf g) = H(\varepsilon_{\mathrm{ext}}(\mathbf g),
     * \varepsilon_{\mathrm{ext}}(\mathbf g + \mathbf e_d))\f$ on the extended map.
     * @details The value sits on the positive \f$d\f$-face of cell \f$\mathbf g\f$, where the
     * staggered vector potential component \f$A_d(\mathbf g)\f$ is read
     * [ryu2016potentialfdtd, Appendix A, Eqs. (19)-(24)]; the cell-centre samples put every
     * material interface on such a face, midway between two samples. Deduction for the ghost
     * layer of a non-periodic axis: at \f$g_d = -1\f$ the value is
     * \f$H(\varepsilon_0, \varepsilon_0) = \varepsilon_0\f$, not the neighbouring face
     * \f$F_d(0) = H(\varepsilon_0, \varepsilon_1)\f$; at \f$g_d = N\f$ it is
     * \f$\varepsilon_{N-1} = F_d(N - 1)\f$. On a periodic axis \f$F_d(-1) = F_d(N - 1)\f$. So a
     * copy of the neighbouring face is right everywhere except on the lower physical ghost
     * layer of the face's own axis. Serves the derived-field check on every allocated cell.
     */
    template <class Map>
    KOKKOS_INLINE_FUNCTION double faceValue(const Map& map, ippl::Vector<int, 3> global,
                                            unsigned axis, ippl::Vector<int, 3> cells,
                                            Kokkos::Array<bool, 3> periodic) {
        ippl::Vector<int, 3> upper = global;
        upper[axis] += 1;
        return harmonicMean(extendedValue(map, global, cells, periodic),
                            extendedValue(map, upper, cells, periodic));
    }
}  // namespace chdr::material
