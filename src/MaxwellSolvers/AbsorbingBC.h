#ifndef ABC_H
#define ABC_H
#include <Kokkos_Core.hpp>

#include "Types/Vector.h"

#include "Field/Field.h"

template <unsigned a, unsigned b>
constexpr KOKKOS_INLINE_FUNCTION auto first() {
    return a;
}
template <unsigned a, unsigned b>
constexpr KOKKOS_INLINE_FUNCTION auto second() {
    return b;
}

/* @struct second_order_abc_face
 * @brief Implements a second-order absorbing boundary condition (ABC) for a face of a 3D
 * computational domain.
 *
 * This struct computes and applies the second-order Mur ABC on a specified face of the domain,
 * minimizing reflections for electromagnetic field solvers. The face is defined by its main axis
 * (the normal direction) and two side axes (the tangential directions).
 *
 * The struct precomputes weights based on mesh spacing and time step, and provides an operator()
 * to update the field at the boundary using current, previous, and next time-step values.
 *
 * @tparam _scalar      Scalar type (e.g., float or double).
 * @tparam _main_axis   The main axis (0=x, 1=y, 2=z) normal to the face.
 * @tparam _side_axes   The two side axes (tangential to the face).
 */
template <typename _scalar, unsigned _main_axis, unsigned... _side_axes>
struct second_order_abc_face {
    using scalar = _scalar;  // Define the scalar type for the weights and calculations
    scalar Cweights[5];  // Weights for the second-order ABC, precomputed based on mesh spacing and
                         // time step
    int sign;            // Sign of the main axis (1 for lower boundary, -1 for upper).
    constexpr static unsigned main_axis =
        _main_axis;  // Main axis normal to the face (0=x, 1=y, 2=z)

    /*!
     * @brief Constructor for the second-order ABC face.
     * @param hr Mesh spacing in each dimension.
     * @param dt Time step size.
     * @param _sign Sign of the main axis (1 or -1).
     * @param waveSpeed Positive wave speed in the units of dt and hr; one preserves the default.
     * The weights discretize Mur (1981), section IV, equation (16), extended
     * to unequal tangential spacings [mur1981AbsorbingBoundaryFdtd].
     */
    KOKKOS_FUNCTION second_order_abc_face(ippl::Vector<scalar, 3> hr, scalar dt, int _sign,
                                          scalar waveSpeed = scalar(1))
        : sign(_sign) {
        const scalar c = waveSpeed;

        // Check that the main axis is one of the three axes and the side axes are different from
        // each other and from the main axis
        constexpr unsigned side_axes[2] = {_side_axes...};
        static_assert(
            (main_axis == 0 && first<_side_axes...>() == 1 && second<_side_axes...>() == 2)
            || (main_axis == 1 && first<_side_axes...>() == 0 && second<_side_axes...>() == 2)
            || (main_axis == 2 && first<_side_axes...>() == 0 && second<_side_axes...>() == 1));
        assert(_main_axis != side_axes[0]);
        assert(_main_axis != side_axes[1]);
        assert(side_axes[0] != side_axes[1]);

        // Calculate prefactor for the weights used to calculate A_np1 at the boundary
        scalar d = 1.0 / (2.0 * dt * hr[main_axis]) + 1.0 / (2.0 * c * dt * dt);

        // Calculate the weights for the boundary condition
        Cweights[0] = (1.0 / (2.0 * dt * hr[main_axis]) - 1.0 / (2.0 * c * dt * dt)) / d;
        Cweights[1] = (-1.0 / (2.0 * dt * hr[main_axis]) - 1.0 / (2.0 * c * dt * dt)) / d;
        assert(abs(Cweights[1] + 1) < 1e-6);  // Cweight[1] should always be -1
        Cweights[2] = (1.0 / (c * dt * dt) - c / (2.0 * hr[side_axes[0]] * hr[side_axes[0]])
                       - c / (2.0 * hr[side_axes[1]] * hr[side_axes[1]]))
                      / d;
        Cweights[3] = (c / (4.0 * hr[side_axes[0]] * hr[side_axes[0]])) / d;
        Cweights[4] = (c / (4.0 * hr[side_axes[1]] * hr[side_axes[1]])) / d;
    }

    /*!
     * @brief Applies the second-order ABC to the boundary of the field.
     * @param A_n Current time step field.
     * @param A_nm1 Previous time step field.
     * @param A_np1 Next time step field.
     * @param c Coordinates of the current point in the field.
     * @return Updated value at the boundary.
     */
    template <typename view_type, typename Coords>
    KOKKOS_INLINE_FUNCTION auto operator()(const view_type& A_n, const view_type& A_nm1,
                                           const view_type& A_np1, const Coords& c) const ->
        typename view_type::value_type {
        uint32_t i = c[0];
        uint32_t j = c[1];
        uint32_t k = c[2];
        using ippl::apply;
        constexpr unsigned side_axes[2] = {_side_axes...};

        // Create vectors with 1 in the direction of the specified axes and 0 in the others
        ippl::Vector<uint32_t, 3> side_axis1_onehot =
            ippl::Vector<uint32_t, 3>{side_axes[0] == 0, side_axes[0] == 1, side_axes[0] == 2};
        ippl::Vector<uint32_t, 3> side_axis2_onehot =
            ippl::Vector<uint32_t, 3>{side_axes[1] == 0, side_axes[1] == 1, side_axes[1] == 2};
        // Create a vector in the direction of the main axis, the sign is used to determine the
        // direction of the offset depending on which side of the boundary we are at
        ippl::Vector<uint32_t, 3> mainaxis_off = ippl::Vector<int32_t, 3>{
            (main_axis == 0) * sign, (main_axis == 1) * sign, (main_axis == 2) * sign};

        // Apply the boundary condition
        // The main axis is the one that is normal to the face, so we need to apply the boundary
        // condition in the direction of the main axis, this means we need to offset the coordinates
        // by the sign of the main axis
        return advanceBoundaryS(
            A_nm1(i, j, k), A_n(i, j, k), apply(A_nm1, c + mainaxis_off),
            apply(A_n, c + mainaxis_off), apply(A_np1, c + mainaxis_off),
            apply(A_n, c + side_axis1_onehot + mainaxis_off),
            apply(A_n, c - side_axis1_onehot + mainaxis_off),
            apply(A_n, c + side_axis2_onehot + mainaxis_off),
            apply(A_n, c - side_axis2_onehot + mainaxis_off), apply(A_n, c + side_axis1_onehot),
            apply(A_n, c - side_axis1_onehot), apply(A_n, c + side_axis2_onehot),
            apply(A_n, c - side_axis2_onehot));
    }

    /*!
     * @brief Advances the boundary condition using the precomputed weights.
     */
    template <typename value_type>
    KOKKOS_FUNCTION value_type advanceBoundaryS(const value_type& v1, const value_type& v2,
                                                const value_type& v3, const value_type& v4,
                                                const value_type& v5, const value_type& v6,
                                                const value_type& v7, const value_type& v8,
                                                const value_type& v9, const value_type& v10,
                                                const value_type& v11, const value_type& v12,
                                                const value_type& v13) const noexcept {
        value_type v0 = Cweights[0] * (v1 + v5) + (Cweights[1]) * v3 + (Cweights[2]) * (v2 + v4)
                        + (Cweights[3]) * (v6 + v7 + v10 + v11)
                        + (Cweights[4]) * (v8 + v9 + v12 + v13);
        return v0;
    }
};

/* @struct second_order_abc_edge
 * @brief Implements a second-order absorbing boundary condition (ABC) for an edge of a 3D
 * computational domain.
 *
 * This struct computes and applies the second-order Mur ABC on a specified edge of the domain,
 * minimizing reflections for electromagnetic field solvers. The edge is defined by its axis (the
 * edge direction) and two normal axes (the directions normal to the edge).
 *
 * The struct precomputes weights based on mesh spacing and time step, and provides an operator()
 * to update the field at the boundary using current, previous, and next time-step values.
 *
 * @tparam _scalar        Scalar type (e.g., float or double).
 * @tparam edge_axis      The axis (0=x, 1=y, 2=z) along the edge.
 * @tparam normal_axis1   The first normal axis to the edge.
 * @tparam normal_axis2   The second normal axis to the edge.
 * @tparam na1_zero       Boolean indicating direction for normal_axis1 (true for lower boundary,
 * false for upper).
 * @tparam na2_zero       Boolean indicating direction for normal_axis2 (true for lower boundary,
 * false for upper).
 */
template <typename _scalar, unsigned edge_axis, unsigned normal_axis1, unsigned normal_axis2,
          bool na1_zero, bool na2_zero>
struct second_order_abc_edge {
    using scalar = _scalar;  // Define the scalar type for the weights and calculations
    scalar Eweights[5];  // Weights for the second-order ABC, precomputed based on mesh spacing and
                         // time step

    /*!
     * @brief Constructor for the second-order ABC edge.
     * @param hr Mesh spacing in each dimension.
     * @param dt Time step size.
     * @param waveSpeed Speed shared by the wave equation and the edge condition.
     */
    KOKKOS_FUNCTION second_order_abc_edge(ippl::Vector<scalar, 3> hr, scalar dt,
                                          scalar waveSpeed = scalar(1)) {
        // Check that the edge axis is one of the three axes and the normal axes are different from
        // each other and from the edge axis
        static_assert(normal_axis1 != normal_axis2);
        static_assert(edge_axis != normal_axis2);
        static_assert(edge_axis != normal_axis1);
        // Check that the axes are in the correct order
        static_assert((edge_axis == 2 && normal_axis1 == 0 && normal_axis2 == 1)
                      || (edge_axis == 0 && normal_axis1 == 1 && normal_axis2 == 2)
                      || (edge_axis == 1 && normal_axis1 == 2 && normal_axis2 == 0));

        // Speed of light in the medium.
        // The same speed enters the wave equation and the edge compatibility
        // condition [bamberger1987SecondOrderCorner, section II, equation (2.10)].
        const scalar c0_ = waveSpeed;

        // Calculate prefactor for the weights used to calculate A_np1 at the boundary
        scalar d = (1.0 / hr[normal_axis1] + 1.0 / hr[normal_axis2]) / (4.0 * dt)
                   + 3.0 / (8.0 * c0_ * dt * dt);

        // Calculate the weights for the boundary condition
        Eweights[0] = (-(1.0 / hr[normal_axis2] - 1.0 / hr[normal_axis1]) / (4.0 * dt)
                       - 3.0 / (8.0 * c0_ * dt * dt))
                      / d;
        Eweights[1] = ((1.0 / hr[normal_axis2] - 1.0 / hr[normal_axis1]) / (4.0 * dt)
                       - 3.0 / (8.0 * c0_ * dt * dt))
                      / d;
        Eweights[2] = ((1.0 / hr[normal_axis2] + 1.0 / hr[normal_axis1]) / (4.0 * dt)
                       - 3.0 / (8.0 * c0_ * dt * dt))
                      / d;
        Eweights[3] =
            (3.0 / (4.0 * c0_ * dt * dt) - c0_ / (4.0 * hr[edge_axis] * hr[edge_axis])) / d;
        Eweights[4] = c0_ / (8.0 * hr[edge_axis] * hr[edge_axis]) / d;
    }

    /*!
     * @brief Applies the second-order ABC to the edge of the field.
     * @param A_n Current time step field.
     * @param A_nm1 Previous time step field.
     * @param A_np1 Next time step field.
     * @param c Coordinates of the current point in the field.
     * @return Updated value at the edge.
     */
    template <typename view_type, typename Coords>
    KOKKOS_INLINE_FUNCTION auto operator()(const view_type& A_n, const view_type& A_nm1,
                                           const view_type& A_np1, const Coords& c) const ->
        typename view_type::value_type {
        uint32_t i = c[0];
        uint32_t j = c[1];
        uint32_t k = c[2];
        using ippl::apply;

        // Get the opposite direction of the normal vectors for the two normal axes (direction to
        // the interior of the domain), 1 or -1 in the direction of the normal axis and 0 in the
        // others
        ippl::Vector<int32_t, 3> normal_axis1_onehot =
            ippl::Vector<int32_t, 3>{normal_axis1 == 0, normal_axis1 == 1, normal_axis1 == 2}
            * int32_t(na1_zero ? 1 : -1);
        ippl::Vector<int32_t, 3> normal_axis2_onehot =
            ippl::Vector<int32_t, 3>{normal_axis2 == 0, normal_axis2 == 1, normal_axis2 == 2}
            * int32_t(na2_zero ? 1 : -1);

        // Calculate the coordinates of the neighboring points in the interior of the domain
        ippl::Vector<uint32_t, 3> acc0 = {i, j, k};
        ippl::Vector<uint32_t, 3> acc1 = acc0 + normal_axis1_onehot;
        ippl::Vector<uint32_t, 3> acc2 = acc0 + normal_axis2_onehot;
        ippl::Vector<uint32_t, 3> acc3 = acc0 + normal_axis1_onehot + normal_axis2_onehot;

        // Create a vector in the direction of the edge axis
        ippl::Vector<uint32_t, 3> axisp{edge_axis == 0, edge_axis == 1, edge_axis == 2};

        // Apply the boundary condition
        return advanceEdgeS(A_n(i, j, k), A_nm1(i, j, k), apply(A_np1, acc1), apply(A_n, acc1),
                            apply(A_nm1, acc1), apply(A_np1, acc2), apply(A_n, acc2),
                            apply(A_nm1, acc2), apply(A_np1, acc3), apply(A_n, acc3),
                            apply(A_nm1, acc3), apply(A_n, acc0 - axisp), apply(A_n, acc1 - axisp),
                            apply(A_n, acc2 - axisp), apply(A_n, acc3 - axisp),
                            apply(A_n, acc0 + axisp), apply(A_n, acc1 + axisp),
                            apply(A_n, acc2 + axisp), apply(A_n, acc3 + axisp));
    }

    /*!
     * @brief Advances the edge boundary condition using the precomputed weights.
     */
    template <typename value_type>
    KOKKOS_INLINE_FUNCTION value_type advanceEdgeS(value_type v1, value_type v2, value_type v3,
                                                   value_type v4, value_type v5, value_type v6,
                                                   value_type v7, value_type v8, value_type v9,
                                                   value_type v10, value_type v11, value_type v12,
                                                   value_type v13, value_type v14, value_type v15,
                                                   value_type v16, value_type v17, value_type v18,
                                                   value_type v19) const noexcept {
        value_type v0 = Eweights[0] * (v3 + v8) + Eweights[1] * (v5 + v6) + Eweights[2] * (v2 + v9)
                        + Eweights[3] * (v1 + v4 + v7 + v10)
                        + Eweights[4] * (v12 + v13 + v14 + v15 + v16 + v17 + v18 + v19) - v11;
        return v0;
    }
};

/* @struct second_order_abc_corner
 * @brief Implements a second-order absorbing boundary condition (ABC) for a corner of a 3D
 * computational domain.
 *
 * This struct computes and applies the second-order Mur ABC on a specified corner of the domain,
 * minimizing reflections for electromagnetic field solvers. The corner is defined by its three
 * axes, with each axis having a boolean indicating whether it is at the lower or upper boundary.
 *
 * The struct precomputes weights based on mesh spacing and time step, and provides an operator()
 * to update the field at the boundary using current, previous, and next time-step values.
 *
 * @tparam _scalar  Scalar type (e.g., float or double).
 * @tparam x0       Boolean indicating direction for x-axis (true for lower boundary, false for
 * upper).
 * @tparam y0       Boolean indicating direction for y-axis (true for lower boundary, false for
 * upper).
 * @tparam z0       Boolean indicating direction for z-axis (true for lower boundary, false for
 * upper).
 */
template <typename _scalar, bool x0, bool y0, bool z0>
struct second_order_abc_corner {
    using scalar = _scalar;  // Define the scalar type for the weights and calculations
    scalar Cweights[17];  // Weights for the second-order ABC, precomputed based on mesh spacing and
                          // time step

    /*!
     * @brief Constructor for the second-order ABC corner.
     * @param hr Mesh spacing in each dimension.
     * @param dt Time step size.
     * @param waveSpeed Speed shared by the wave equation and the corner condition.
     */
    KOKKOS_FUNCTION second_order_abc_corner(ippl::Vector<scalar, 3> hr, scalar dt,
                                            scalar waveSpeed = scalar(1)) {
        // Time derivative of the matched-speed corner compatibility condition
        // [bamberger1987SecondOrderCorner, section II, equation (2.11)].
        const scalar c0_ = waveSpeed;
        Cweights[0] =
            (-1.0 / hr[0] - 1.0 / hr[1] - 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        Cweights[1] =
            (1.0 / hr[0] - 1.0 / hr[1] - 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        Cweights[2] =
            (-1.0 / hr[0] + 1.0 / hr[1] - 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        Cweights[3] =
            (-1.0 / hr[0] - 1.0 / hr[1] + 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        Cweights[4] =
            (1.0 / hr[0] + 1.0 / hr[1] - 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        Cweights[5] =
            (1.0 / hr[0] - 1.0 / hr[1] + 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        Cweights[6] =
            (-1.0 / hr[0] + 1.0 / hr[1] + 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        Cweights[7] =
            (1.0 / hr[0] + 1.0 / hr[1] + 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        Cweights[8] =
            -(-1.0 / hr[0] - 1.0 / hr[1] - 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        Cweights[9] =
            -(1.0 / hr[0] - 1.0 / hr[1] - 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        Cweights[10] =
            -(-1.0 / hr[0] + 1.0 / hr[1] - 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        Cweights[11] =
            -(-1.0 / hr[0] - 1.0 / hr[1] + 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        Cweights[12] =
            -(1.0 / hr[0] + 1.0 / hr[1] - 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        Cweights[13] =
            -(1.0 / hr[0] - 1.0 / hr[1] + 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        Cweights[14] =
            -(-1.0 / hr[0] + 1.0 / hr[1] + 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        Cweights[15] =
            -(1.0 / hr[0] + 1.0 / hr[1] + 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        Cweights[16] = 1.0 / (2.0 * c0_ * dt * dt);
    }

    /*!
     * @brief Applies the second-order ABC to the corner of the field.
     * @param A_n Current time step field.
     * @param A_nm1 Previous time step field.
     * @param A_np1 Next time step field.
     * @param c Coordinates of the current point in the field.
     * @return Updated value at the corner.
     */
    template <typename view_type, typename Coords>
    KOKKOS_INLINE_FUNCTION auto operator()(const view_type& A_n, const view_type& A_nm1,
                                           const view_type& A_np1, const Coords& c) const ->
        typename view_type::value_type {
        // Get direction of the vector to the interior of the domain
        constexpr uint32_t xoff = (x0) ? 1 : uint32_t(-1);
        constexpr uint32_t yoff = (y0) ? 1 : uint32_t(-1);
        constexpr uint32_t zoff = (z0) ? 1 : uint32_t(-1);
        using ippl::apply;

        // Calculate the offset vectors for the neighboring points in the interior of the domain
        const ippl::Vector<uint32_t, 3> offsets[8] = {
            ippl::Vector<uint32_t, 3>{0, 0, 0},       ippl::Vector<uint32_t, 3>{xoff, 0, 0},
            ippl::Vector<uint32_t, 3>{0, yoff, 0},    ippl::Vector<uint32_t, 3>{0, 0, zoff},
            ippl::Vector<uint32_t, 3>{xoff, yoff, 0}, ippl::Vector<uint32_t, 3>{xoff, 0, zoff},
            ippl::Vector<uint32_t, 3>{0, yoff, zoff}, ippl::Vector<uint32_t, 3>{xoff, yoff, zoff},
        };

        // Apply the boundary condition
        return advanceCornerS(
            apply(A_n, c), apply(A_nm1, c), apply(A_np1, c + offsets[1]),
            apply(A_n, c + offsets[1]), apply(A_nm1, c + offsets[1]), apply(A_np1, c + offsets[2]),
            apply(A_n, c + offsets[2]), apply(A_nm1, c + offsets[2]), apply(A_np1, c + offsets[3]),
            apply(A_n, c + offsets[3]), apply(A_nm1, c + offsets[3]), apply(A_np1, c + offsets[4]),
            apply(A_n, c + offsets[4]), apply(A_nm1, c + offsets[4]), apply(A_np1, c + offsets[5]),
            apply(A_n, c + offsets[5]), apply(A_nm1, c + offsets[5]), apply(A_np1, c + offsets[6]),
            apply(A_n, c + offsets[6]), apply(A_nm1, c + offsets[6]), apply(A_np1, c + offsets[7]),
            apply(A_n, c + offsets[7]), apply(A_nm1, c + offsets[7]));
    }

    /*!
     * @brief Advances the corner boundary condition using the precomputed weights.
     */
    template <typename value_type>
    KOKKOS_INLINE_FUNCTION value_type
    advanceCornerS(value_type v1, value_type v2, value_type v3, value_type v4, value_type v5,
                   value_type v6, value_type v7, value_type v8, value_type v9, value_type v10,
                   value_type v11, value_type v12, value_type v13, value_type v14, value_type v15,
                   value_type v16, value_type v17, value_type v18, value_type v19, value_type v20,
                   value_type v21, value_type v22, value_type v23) const noexcept {
        return -(v1 * (Cweights[16]) + v2 * (Cweights[8]) + v3 * Cweights[1] + v4 * Cweights[16]
                 + v5 * Cweights[9] + v6 * Cweights[2] + v7 * Cweights[16] + v8 * Cweights[10]
                 + v9 * Cweights[3] + v10 * Cweights[16] + v11 * Cweights[11] + v12 * Cweights[4]
                 + v13 * Cweights[16] + v14 * Cweights[12] + v15 * Cweights[5] + v16 * Cweights[16]
                 + v17 * Cweights[13] + v18 * Cweights[6] + v19 * Cweights[16] + v20 * Cweights[14]
                 + v21 * Cweights[7] + v22 * Cweights[16] + v23 * Cweights[15])
               / Cweights[0];
    }
};

/* @struct second_order_mur_boundary_conditions
 * @brief Applies second-order Mur absorbing boundary conditions (ABC) to all boundaries of a 3D
 * computational domain.
 *
 * This struct coordinates the application of second-order Mur ABCs on faces, edges, and corners of
 * the domain, minimizing reflections for electromagnetic field solvers. It uses the appropriate
 * face, edge, or corner ABC implementation depending on the location of each boundary point.
 *
 * The apply() method detects the type of boundary (face, edge, or corner) for each point and
 * applies the corresponding second-order ABC using precomputed weights based on mesh spacing and
 * time step.
 *
 * @tparam field_type  The field type (must provide getView() and get_mesh().getMeshSpacing()).
 * @tparam dt_type     The type of the time step.
 */
struct second_order_mur_boundary_conditions {
    /// Read-only Mur speed of each local cell, in the memory space of the potential field.
    template <typename field_type, typename dt_type>
    using cell_speed_view =
        Kokkos::View<const dt_type***, typename field_type::view_type::memory_space>;

    /*!
     * @brief Apply second-order Mur ABCs on the selected nonperiodic domain faces.
     * @details An open face here is a global domain face selected for Mur absorption,
     * not a periodic wrap or an internal MPI partition face. Mur updates its outermost
     * owned cells; this routine does not fill the exterior ghost plane.
     * @param FA_n Field at the current time step.
     * @param FA_nm1 Field at the previous time step.
     * @param FA_np1 Field at the next time step.
     * @param dt Time step size.
     * @param true_nr True number of grid points in each dimension.
     * @param lDom Local domain index for parallel processing.
     * @param openFaces Mask selecting Mur faces: bit 2*d is the lower face and
     * bit 2*d+1 the upper face of axis d. Unset faces receive no Mur update.
     * @param speeds Mur wave speeds in kernel units, indexed in the same paired order.
     * Entries for unselected faces are ignored. With Mur faces on only one axis,
     * the two opposing faces may have different validated speeds; they do not intersect.
     * With Mur faces on two or three axes, validateBoundary requires the same
     * speed on each selected face. The edge and corner updates use that speed.
     * @param cellSpeeds Optional Mur speed of each cell in kernel units, with the extents of
     * the local field view including ghosts. If it is allocated, the face, edge or corner row
     * of an owned cell uses the speed of that cell, and speeds is ignored. One speed per cell
     * satisfies the common-speed condition at every intersection. In a medium that varies
     * along the boundary, each row thus uses the local speed. Engquist and Majda construct
     * local absorbing conditions for variable coefficients with the coefficients at the
     * boundary point [Engquist and Majda (1977), section 1.C, Eqs. (1.18)-(1.26)]. Their
     * second approximation adds a term with the normal derivative of the coefficient,
     * which is zero where the medium continues outward. This row omits the derivative
     * terms (deduction; their Eq. (1.18) varies only the tangential coefficient).
     * @pre The caller has checked the face pairs, speeds and storage with validateBoundary.
     * The intersecting-face restriction follows the matched-speed contract documented
     * there [Bamberger, Joly and Roberts (1987), section II, Eqs. (2.10)-(2.11)].
     * The caller has checked that each owned speed in cellSpeeds is admissible.
     */
    template <typename field_type, typename dt_type>
    void apply(field_type& FA_n, field_type& FA_nm1, field_type& FA_np1, dt_type dt,
               ippl::Vector<uint32_t, 3> true_nr, ippl::NDIndex<3> lDom, unsigned openFaces = 63,
               ippl::Vector<dt_type, 6> speeds                        = ippl::Vector<dt_type, 6>(1),
               const cell_speed_view<field_type, dt_type>& cellSpeeds = {}) {
        using scalar = decltype(dt);

        // Get the mesh spacing
        const ippl::Vector<scalar, 3> hr = FA_n.get_mesh().getMeshSpacing();

        // Get views of the fields
        auto A_n   = FA_n.getView();
        auto A_np1 = FA_np1.getView();
        auto A_nm1 = FA_nm1.getView();

        // Get the local number of grid points in each dimension
        ippl::Vector<uint32_t, 3> local_nr{uint32_t(A_n.extent(0)), uint32_t(A_n.extent(1)),
                                           uint32_t(A_n.extent(2))};

        const unsigned activeMask             = toDispatchMask(openFaces);
        const scalar uniformIntersectionSpeed = selectedIntersectionSpeed(openFaces, speeds);
        const bool useCellSpeeds              = cellSpeeds.is_allocated();

        unsigned localFaces = activeMask;
        for (unsigned d = 0; d < 3; ++d) {
            if (lDom[d].first() != 0)
                localFaces &= ~(1u << d);
            if (uint32_t(lDom[d].last()) != true_nr[d] - 3)
                localFaces &= ~(8u << d);
        }

        // Each group contains the owned points with exactly these incident open axes.
        using Policy =
            Kokkos::MDRangePolicy<typename decltype(A_n)::execution_space, Kokkos::Rank<3>>;
        auto applyGroup = [&](unsigned axes) {
            Kokkos::Array<uint32_t, 3> begin{}, end{};
            if (!groupBounds(axes, localFaces, local_nr, begin, end))
                return;
            const auto update = KOKKOS_LAMBDA(uint32_t i, uint32_t j, uint32_t k) {
                ippl::Vector<uint32_t, 3> point{i, j, k};
                const unsigned val = boundaryPoint(axes, local_nr, point);
                i                  = point[0];
                j                  = point[1];
                k                  = point[2];

                // A cell speed replaces the face and intersection speeds at this cell.
                const scalar cellSpeed = useCellSpeeds ? cellSpeeds(i, j, k) : scalar(0);
                const ippl::Vector<scalar, 6> faceSpeeds =
                    useCellSpeeds ? ippl::Vector<scalar, 6>(cellSpeed) : speeds;
                const scalar intersectionSpeed =
                    useCellSpeeds ? cellSpeed : uniformIntersectionSpeed;

                // Exactly one incident Mur face selects a face update, even at a periodic wrap.
                if (Kokkos::popcount(val) == 1) {
                    const bool onLowerX = (val & 1u) != 0;
                    const bool onLowerY = (val & 2u) != 0;
                    const bool onLowerZ = (val & 4u) != 0;
                    const bool onUpperX = (val & 8u) != 0;
                    const bool onUpperY = (val & 16u) != 0;
                    const bool onUpperZ = (val & 32u) != 0;
                    if (onLowerX) {
                        second_order_abc_face<scalar, 0, 1, 2> soa(hr, dt, 1, faceSpeeds[0]);
                        A_np1(i, j, k) = soa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    }
                    if (onLowerY) {
                        second_order_abc_face<scalar, 1, 0, 2> soa(hr, dt, 1, faceSpeeds[2]);
                        A_np1(i, j, k) = soa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    }
                    if (onLowerZ) {
                        second_order_abc_face<scalar, 2, 0, 1> soa(hr, dt, 1, faceSpeeds[4]);
                        A_np1(i, j, k) = soa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    }
                    if (onUpperX) {
                        second_order_abc_face<scalar, 0, 1, 2> soa(hr, dt, -1, faceSpeeds[1]);
                        A_np1(i, j, k) = soa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    }
                    if (onUpperY) {
                        second_order_abc_face<scalar, 1, 0, 2> soa(hr, dt, -1, faceSpeeds[3]);
                        A_np1(i, j, k) = soa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    }
                    if (onUpperZ) {
                        second_order_abc_face<scalar, 2, 0, 1> soa(hr, dt, -1, faceSpeeds[5]);
                        A_np1(i, j, k) = soa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    }
                }

                // An edge update requires two incident Mur faces, not merely two domain faces.
                if (Kokkos::popcount(val) == 2) {
                    const bool onLowerX = (val & 1u) != 0;
                    const bool onLowerY = (val & 2u) != 0;
                    const bool onLowerZ = (val & 4u) != 0;
                    const bool onUpperX = (val & 8u) != 0;
                    const bool onUpperY = (val & 16u) != 0;
                    const bool onUpperZ = (val & 32u) != 0;
                    if (onLowerX && onLowerZ) {
                        second_order_abc_edge<scalar, 1, 2, 0, true, true> soa(hr, dt,
                                                                               intersectionSpeed);
                        A_np1(i, j, k) = soa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    } else if (onLowerX && onLowerY) {
                        second_order_abc_edge<scalar, 2, 0, 1, true, true> soa(hr, dt,
                                                                               intersectionSpeed);
                        A_np1(i, j, k) = soa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    } else if (onLowerY && onLowerZ) {
                        second_order_abc_edge<scalar, 0, 1, 2, true, true> soa(hr, dt,
                                                                               intersectionSpeed);
                        A_np1(i, j, k) = soa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    } else if (onLowerX && onUpperZ) {
                        second_order_abc_edge<scalar, 1, 2, 0, false, true> soa(hr, dt,
                                                                                intersectionSpeed);
                        A_np1(i, j, k) = soa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    } else if (onLowerX && onUpperY) {
                        second_order_abc_edge<scalar, 2, 0, 1, true, false> soa(hr, dt,
                                                                                intersectionSpeed);
                        A_np1(i, j, k) = soa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    } else if (onLowerY && onUpperZ) {
                        second_order_abc_edge<scalar, 0, 1, 2, true, false> soa(hr, dt,
                                                                                intersectionSpeed);
                        A_np1(i, j, k) = soa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    } else if (onUpperX && onLowerZ) {
                        second_order_abc_edge<scalar, 1, 2, 0, true, false> soa(hr, dt,
                                                                                intersectionSpeed);
                        A_np1(i, j, k) = soa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    } else if (onUpperX && onLowerY) {
                        second_order_abc_edge<scalar, 2, 0, 1, false, true> soa(hr, dt,
                                                                                intersectionSpeed);
                        A_np1(i, j, k) = soa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    } else if (onUpperY && onLowerZ) {
                        second_order_abc_edge<scalar, 0, 1, 2, false, true> soa(hr, dt,
                                                                                intersectionSpeed);
                        A_np1(i, j, k) = soa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    } else if (onUpperX && onUpperZ) {
                        second_order_abc_edge<scalar, 1, 2, 0, false, false> soa(hr, dt,
                                                                                 intersectionSpeed);
                        A_np1(i, j, k) = soa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    } else if (onUpperX && onUpperY) {
                        second_order_abc_edge<scalar, 2, 0, 1, false, false> soa(hr, dt,
                                                                                 intersectionSpeed);
                        A_np1(i, j, k) = soa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    } else if (onUpperY && onUpperZ) {
                        second_order_abc_edge<scalar, 0, 1, 2, false, false> soa(hr, dt,
                                                                                 intersectionSpeed);
                        A_np1(i, j, k) = soa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    } else {
                        assert(false);
                    }
                }

                // A corner update requires three incident Mur faces.
                if (Kokkos::popcount(val) == 3) {
                    const bool onLowerX = (val & 1u) != 0;
                    const bool onLowerY = (val & 2u) != 0;
                    const bool onLowerZ = (val & 4u) != 0;
                    const bool onUpperX = (val & 8u) != 0;
                    const bool onUpperY = (val & 16u) != 0;
                    const bool onUpperZ = (val & 32u) != 0;
                    if (onLowerX && onLowerY && onLowerZ) {
                        second_order_abc_corner<scalar, 1, 1, 1> coa(hr, dt, intersectionSpeed);
                        A_np1(i, j, k) = coa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    } else if (onUpperX && onLowerY && onLowerZ) {
                        second_order_abc_corner<scalar, 0, 1, 1> coa(hr, dt, intersectionSpeed);
                        A_np1(i, j, k) = coa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    } else if (onLowerX && onUpperY && onLowerZ) {
                        second_order_abc_corner<scalar, 1, 0, 1> coa(hr, dt, intersectionSpeed);
                        A_np1(i, j, k) = coa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    } else if (onUpperX && onUpperY && onLowerZ) {
                        second_order_abc_corner<scalar, 0, 0, 1> coa(hr, dt, intersectionSpeed);
                        A_np1(i, j, k) = coa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    } else if (onLowerX && onLowerY && onUpperZ) {
                        second_order_abc_corner<scalar, 1, 1, 0> coa(hr, dt, intersectionSpeed);
                        A_np1(i, j, k) = coa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    } else if (onUpperX && onLowerY && onUpperZ) {
                        second_order_abc_corner<scalar, 0, 1, 0> coa(hr, dt, intersectionSpeed);
                        A_np1(i, j, k) = coa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    } else if (onLowerX && onUpperY && onUpperZ) {
                        second_order_abc_corner<scalar, 1, 0, 0> coa(hr, dt, intersectionSpeed);
                        A_np1(i, j, k) = coa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    } else if (onUpperX && onUpperY && onUpperZ) {
                        second_order_abc_corner<scalar, 0, 0, 0> coa(hr, dt, intersectionSpeed);
                        A_np1(i, j, k) = coa(A_n, A_nm1, A_np1, ippl::Vector<uint32_t, 3>{i, j, k});
                    } else {
                        assert(false);
                    }
                }
            };
            Kokkos::parallel_for("Mur boundary", Policy(begin, end), update);
        };

        constexpr unsigned x = 1, y = 2, z = 4;
        // Apply second-order ABC to faces
        applyGroup(x);
        applyGroup(y);
        applyGroup(z);
        Kokkos::fence();

        // Apply second-order ABC to edges
        applyGroup(x | y);
        applyGroup(x | z);
        applyGroup(y | z);
        Kokkos::fence();

        // Apply second-order ABC to corners
        applyGroup(x | y | z);
        Kokkos::fence();
    }

private:
    /** @brief Bounds of the locally owned face, edge or corner group with normal axes in axes.
     * On a normal axis, iterate over the side labels 0 (lower) and 1 (upper) owned by this rank.
     * On a free axis, omit the owned Mur endpoints. This assigns intersections to a single group
     * and retains periodic wrap samples. local_nr includes one ghost cell on each side.
     */
    static bool groupBounds(unsigned axes, unsigned localFaces,
                            const ippl::Vector<uint32_t, 3>& local_nr,
                            Kokkos::Array<uint32_t, 3>& begin, Kokkos::Array<uint32_t, 3>& end) {
        for (unsigned d = 0; d < 3; ++d) {
            const bool lower  = localFaces & (1u << d);
            const bool upper  = localFaces & (8u << d);
            const bool normal = axes & (1u << d);
            begin[d]          = normal ? unsigned(!lower) : 1u + lower;
            end[d]            = normal ? 1u + upper : local_nr[d] - 1u - upper;
            if (begin[d] >= end[d])
                return false;
        }
        return true;
    }

    /** @brief Convert the side labels to owned allocation indices and incident Mur face bits.
     * The resulting mask has the same lower-x,y,z then upper-x,y,z ordering as toDispatchMask().
     */
    KOKKOS_INLINE_FUNCTION static unsigned boundaryPoint(unsigned axes,
                                                         const ippl::Vector<uint32_t, 3>& local_nr,
                                                         ippl::Vector<uint32_t, 3>& point) {
        unsigned val = 0;
        for (unsigned d = 0; d < 3; ++d) {
            if (axes & (1u << d)) {
                const unsigned side = point[d];
                val |= 1u << (d + 3 * side);
                point[d] = side == 0 ? 1 : local_nr[d] - 2;
            }
        }
        return val;
    }

    /** @brief Return a selected speed without reading a periodic-face entry.
     * @details Intersection updates only run with at least two open axes,
     * when validation has established that all selected speeds are equal.
     */
    template <typename Scalar>
    static Scalar selectedIntersectionSpeed(unsigned openFaces,
                                            const ippl::Vector<Scalar, 6>& speeds) {
        for (unsigned face = 0; face < 6; ++face)
            if ((openFaces >> face) & 1u)
                return speeds[face];
        return 1;
    }

    /** @brief Convert paired face bits to the Mur dispatcher's internal ordering.
     * @param openFaces Bits ordered lower x, upper x, lower y, upper y, lower z, upper z.
     * @return The same selected faces ordered lower x, lower y, lower z, upper x,
     * upper y, upper z. Only the six face bits are used; unset bits remain unset.
     */
    static constexpr unsigned toDispatchMask(unsigned openFaces) {
        unsigned activeMask = 0;
        for (unsigned d = 0; d < 3; ++d)
            for (unsigned side = 0; side < 2; ++side)
                activeMask |= ((openFaces >> (2 * d + side)) & 1u) << (d + 3 * side);
        return activeMask;
    }
};
#endif
