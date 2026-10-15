/** @file BoundaryOracles.h
 * @brief Expected values of the open-boundary layout tests: the role of each owned cell, the
 * owned source of each ghost cell, and the weights of the second-order Mur face, edge and
 * corner updates at a general speed.
 * @details Indices are ghost-offset global indices with one ghost layer: along axis \f$d\f$ the
 * owned cells are \f$1 \ldots n_d\f$ and the ghosts are \f$0\f$ and \f$n_d + 1\f$. An axis is
 * periodic or open; both faces of an open axis carry Mur's condition. These functions are
 * written from the definitions below, independently of the solver.
 */
#pragma once
#include <Kokkos_Core.hpp>

#include <array>

namespace open_boundary_oracle {
    /// Role of an owned cell: updated by the interior stencil, or written by a Mur struct.
    enum class CellRole {
        Interior,
        Face,
        Edge,
        Corner
    };

    /** @brief Open faces whose outermost owned layer contains the cell.
     * @details Bit \f$2d\f$ is set when \f$g_d = 1\f$ and bit \f$2d + 1\f$ when \f$g_d = n_d\f$,
     * for open axes only; periodic axes have no Mur layer and set no bit. The bit number equals
     * the face index (\f$2d\f$ lower, \f$2d + 1\f$ upper). The struct that writes a cell follows
     * from the mask: one bit, the face struct of axis \f$d\f$ with the inward sign of that side;
     * two bits, the edge struct whose normal axes are the two axes and whose edge axis is the
     * third; three bits, the corner struct with the side of every axis. Deduction from the face
     * conditions: Mur's condition holds on every open face
     * [mur1981AbsorbingBoundaryFdtd, section III, Eq. (12), p. 379]; where two or three open
     * faces meet, the compatibility conditions of
     * [bamberger1987SecondOrderCorner, section II, Eqs. (2.10)-(2.11), p. 06] replace them.
     * A cell on a periodic face is an interior cell along that axis.
     * @pre \f$n_d \ge 3\f$ on open axes, so that the two layers of an axis are distinct.
     */
    KOKKOS_INLINE_FUNCTION unsigned murFaceMask(Kokkos::Array<int, 3> g, Kokkos::Array<int, 3> n,
                                                Kokkos::Array<bool, 3> periodic) {
        unsigned mask = 0;
        for (unsigned d = 0; d < 3; ++d) {
            if (!periodic[d]) {
                mask |= unsigned(g[d] == 1) << (2 * d);
                mask |= unsigned(g[d] == n[d]) << (2 * d + 1);
            }
        }
        return mask;
    }

    /// Role from the number of set bits of murFaceMask: 0 interior, 1 face, 2 edge, 3 corner.
    KOKKOS_INLINE_FUNCTION CellRole cellRole(unsigned mask) {
        unsigned bits = 0;
        for (unsigned b = 0; b < 6; ++b) {
            bits += (mask >> b) & 1u;
        }
        return bits == 0
                   ? CellRole::Interior
                   : (bits == 1 ? CellRole::Face : (bits == 2 ? CellRole::Edge : CellRole::Corner));
    }

    /** @brief Owned source, along one axis, of the value that the halo exchange puts at index g.
     * @details Owned indices map to themselves (a ghost inside the global box is a rank seam and
     * receives the neighbour's owned cell). On a periodic axis the ghosts wrap,
     * \f$0 \mapsto n\f$ and \f$n + 1 \mapsto 1\f$. On an open axis the ghosts beyond the faces have
     * no source (-1): the exchange does not write them. A ghost cell of a local field is filled
     * when all three axes return a source, and then holds the owned value at that global index.
     * The accumulation is the transpose: the value in such a ghost cell is added to that owned
     * cell, and a ghost with no source adds nowhere. Deduction from the periodic extension of
     * the global index range; the rule does not depend on the rank decomposition.
     */
    KOKKOS_INLINE_FUNCTION int haloSourceIndex(int g, int n, bool periodic) {
        if (g >= 1 && g <= n) {
            return g;
        }
        if (!periodic) {
            return -1;
        }
        return g == 0 ? n : 1;
    }

    /** @brief Weights \f$C_0 \ldots C_4\f$ of the Mur face update at speed v.
     * @details Mur's condition at a face, normal coordinate \f$\nu\f$ into the box,
     * \f[ \bigl(\partial_\nu\partial_t - v^{-1}\partial_t^2
     *     + (v/2)(\partial_1^2 + \partial_2^2)\bigr)\psi = 0 \f]
     * [engquist1977AbsorbingBoundaryWaves, section 1.B, Eq. (1.10), p. 632;
     * mur1981AbsorbingBoundaryFdtd, section III, Eq. (12), p. 379], is centred at the midpoint of
     * layers 0 and 1 and at \f$t^k\f$: the mixed derivative by a difference of the two layers
     * over \f$2\Delta t\f$, the time and side derivatives as means of the two layers
     * [mur1981AbsorbingBoundaryFdtd, section IV, Eq. (16), p. 380]. Solved for
     * \f$u_0^{k+1}\f$ with \f$p = v\Delta t\f$:
     * \f[ u_0^{k+1} = C_0(u_0^{k-1} + u_1^{k+1}) + C_1u_1^{k-1} + C_2(u_0^k + u_1^k)
     *     + C_3\Sigma_1 + C_4\Sigma_2, \f]
     * \f$\Sigma_s\f$ the sum of the four side neighbours along side axis \f$s\f$ on both layers,
     * \f$C_0 = (p - h_\nu)/(p + h_\nu)\f$, \f$C_1 = -1\f$,
     * \f$C_2 = h_\nu(2 - p^2/h_1^2 - p^2/h_2^2)/(p + h_\nu)\f$,
     * \f$C_s = p^2h_\nu/(2h_s^2(p + h_\nu))\f$ for \f$s = 3, 4\f$ (Deduction: multiply the
     * centred equation by \f$2v\Delta t^2h_\nu\f$). Cubic cells give Mur's Eq. (16).
     */
    inline std::array<double, 5> murFaceWeights(double v, double dt, double hNormal, double hSide1,
                                                double hSide2) {
        const double p = v * dt, den = p + hNormal;
        const double r1 = p * p / (hSide1 * hSide1), r2 = p * p / (hSide2 * hSide2);
        return {(p - hNormal) / den, -1.0, hNormal * (2 - r1 - r2) / den, hNormal * r1 / (2 * den),
                hNormal * r2 / (2 * den)};
    }

    /** @brief Weights \f$E_0 \ldots E_4\f$ of the Mur edge update at speed v.
     * @details Adding the face conditions of the normal axes \f$\nu_1, \nu_2\f$ and removing
     * \f$\partial_{\nu_1}^2 + \partial_{\nu_2}^2\f$ with the wave equation gives
     * \f$((\partial_{\nu_1} + \partial_{\nu_2})\partial_t - (3/(2v))\partial_t^2
     * + (v/2)\partial_e^2)\psi = 0\f$ [bamberger1987SecondOrderCorner, section II, Eq. (2.10),
     * p. 06]. Centred in the \f$2\times2\f$ block of columns (00 the edge cell, 11 diagonal) and at
     * \f$t^k\f$, with \f$p = v\Delta t\f$ and \f$D = 2p(1/h_1 + 1/h_2) + 3\f$:
     * \f$u_{00}^{k+1} = E_0(u_{10}^{k+1} + u_{01}^{k-1}) + E_1(u_{10}^{k-1} + u_{01}^{k+1})
     * + E_2(u_{00}^{k-1} + u_{11}^{k+1}) + E_3\sum u^k + E_4\sum_{\pm e}u^k - u_{11}^{k-1}\f$, with
     * \f$E_0 = (2p(1/h_1 - 1/h_2) - 3)/D\f$, \f$E_1 = (2p(1/h_2 - 1/h_1) - 3)/D\f$,
     * \f$E_2 = (2p(1/h_1 + 1/h_2) - 3)/D\f$, \f$E_3 = (6 - 2p^2/h_e^2)/D\f$,
     * \f$E_4 = (p^2/h_e^2)/D\f$ (Deduction: multiply by \f$8v\Delta t^2\f$).
     */
    inline std::array<double, 5> murEdgeWeights(double v, double dt, double hNormal1,
                                                double hNormal2, double hEdge) {
        const double p = v * dt, a1 = 2 * p / hNormal1, a2 = 2 * p / hNormal2;
        const double pe = p * p / (hEdge * hEdge), den = a1 + a2 + 3;
        return {(a1 - a2 - 3) / den, (a2 - a1 - 3) / den, (a1 + a2 - 3) / den, (6 - 2 * pe) / den,
                pe / den};
    }

    /** @brief Corner weights \f$C_0 \ldots C_{16}\f$ at speed v, raw and normalized by \f$C_0\f$.
     * @details The sum of the three face conditions with the wave equation gives
     * \f$((\partial_x + \partial_y + \partial_z)\partial_t - (2/v)\partial_t^2)\psi = 0\f$, the
     * time derivative of [bamberger1987SecondOrderCorner, section II, Eq. (2.11), p. 06] with
     * inward normals. Centred in the \f$2\times2\times2\f$ cube at \f$t^k\f$, the coefficient of
     * \f$u^{k+1}\f$ in the cell with inward signs \f$s_d = \pm1\f$ is
     * \f$(p\sigma - 2)/(8v\Delta t^2)\f$, \f$\sigma = \sum_d s_d/h_d\f$, that of \f$u^{k-1}\f$ is
     * \f$(-p\sigma - 2)/(8v\Delta t^2)\f$ and that of every \f$u^k\f$ is \f$1/(2v\Delta t^2)\f$
     * (Deduction). Order: \f$j = 0 \ldots 7\f$ has positive signs on the axes of
     * plus[j] (bit d for axis d), \f$j + 8\f$ the \f$u^{k-1}\f$ weight of the same cell, 16 the
     * \f$u^k\f$ weight. normalized[j] \f$= C_j/C_0 = (2 - p\sigma_j)/(2 + p\sum_d h_d^{-1})\f$.
     */
    inline void murCornerWeights(double v, double dt, std::array<double, 3> h,
                                 std::array<double, 17>& raw, std::array<double, 17>& normalized) {
        constexpr unsigned plus[8] = {0, 1, 2, 4, 3, 5, 6, 7};
        const double p = v * dt, scale = 8 * v * dt * dt;
        const double total = 2 + p * (1 / h[0] + 1 / h[1] + 1 / h[2]);
        for (unsigned j = 0; j < 8; ++j) {
            double sigma = 0;
            for (unsigned d = 0; d < 3; ++d) {
                sigma += ((plus[j] >> d) & 1u ? 1.0 : -1.0) / h[d];
            }
            raw[j]            = (p * sigma - 2) / scale;
            raw[j + 8]        = (-p * sigma - 2) / scale;
            normalized[j]     = (2 - p * sigma) / total;
            normalized[j + 8] = (2 + p * sigma) / total;
        }
        raw[16]        = 4 / scale;
        normalized[16] = -4 / total;
    }
}  // namespace open_boundary_oracle
