/** @file MaterialOwner.h
 * @brief Immutable distributed cell permittivity and lazily derived coefficient fields.
 */
#pragma once
#include "Ippl.h"

#include <algorithm>
#include <cstdint>
#include <stdexcept>

namespace chdr::material {
    using Mesh                 = ippl::UniformCartesian<double, 3>;
    using Layout               = ippl::FieldLayout<3>;
    using ScalarField          = ippl::Field<double, 3, Mesh, Mesh::DefaultCentering>;
    using ScalarFieldConstView = ScalarField::view_type::const_type;
    using CellIndex            = ippl::Vector<int, 3>;
    using Position             = ippl::Vector<double, 3>;
    using Policy               = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;

    namespace detail {
        /** @brief Device-copyable index metadata shared by sampling and verification.
         * @details The centre expression retains its arithmetic order, so a geometry
         * predicate sees exactly the mesh demo's original boundary samples.
         */
        struct Coordinates {
            CellIndex first_m, count_m, lower_m;
            Position origin_m, spacing_m;
            bool periodic_m;
            KOKKOS_INLINE_FUNCTION CellIndex global(int i, int j, int k) const {
                return CellIndex(i + first_m[0] - 1, j + first_m[1] - 1, k + first_m[2] - 1);
            }
            KOKKOS_INLINE_FUNCTION Position centre(CellIndex global) const {
                Position result;
                for (unsigned d = 0; d < 3; ++d)
                    result[d] = origin_m[d] + (global[d] + 0.5) * spacing_m[d];
                return result;
            }
        };

        /// Metadata only: no field allocation or host copy of the distributed map.
        inline Coordinates coordinates(Mesh& mesh, Layout& layout) {
            Coordinates result;
            const auto local = layout.getLocalNDIndex(), domain = layout.getDomain();
            for (unsigned d = 0; d < 3; ++d) {
                result.first_m[d] = local[d].first();
                result.count_m[d] = domain[d].length();
                result.lower_m[d] = domain[d].first();
            }
            result.origin_m   = mesh.getOrigin();
            result.spacing_m  = mesh.getMeshSpacing();
            result.periodic_m = layout.isAllPeriodic_m;
            return result;
        }

        /// Include the one ghost layer, edges and corners, excluding view padding.
        inline Policy allocationPolicy(const ScalarField& field) {
            const auto view = field.getView();
            return Policy({0, 0, 0},
                          {int(view.extent(0)), int(view.extent(1)), int(view.extent(2))});
        }

        /// Refuse insufficient local extents collectively before allocating fields.
        inline Layout& checkedLayout(Layout& layout) {
            const auto local = layout.getLocalNDIndex();
            int minimum =
                std::min({int(local[0].length()), int(local[1].length()), int(local[2].length())});
            MPI_Allreduce(MPI_IN_PLACE, &minimum, 1, MPI_INT, MPI_MIN, layout.comm);
            if (minimum < 2)
                throw std::invalid_argument("Material halos need two owned cells per local axis");
            return layout;
        }

        /** @brief Continue the boundary medium outward after inter-rank exchange.
         * @details Mur assumes a homogeneous region outside the scatterers
         * [mur1981AbsorbingBoundaryFdtd, section II, p. 378]. Deduction: copying
         * the nearest boundary value continues the material without a jump in
         * the outward direction. This fills material ghosts; it does not apply
         * an absorbing boundary condition to the electromagnetic fields.
         *
         * The layout currently distinguishes fully periodic from fully physical
         * global boundaries. Internal rank boundaries still use halo exchange.
         * A local allocation has one ghost layer: along each axis, indices 0
         * and extent - 1 are ghosts; indices 1 through extent - 2 are owned.
         * The global index ranges identify which ranks touch a physical face.
         *
         * Each axis sweep includes the full transverse allocation, including
         * ghosts. For example, at a lower x/y edge with interior k, the x sweep
         * fills (0, 1, k) from (1, 1, k); the y sweep then fills (0, 0, k) from
         * (0, 1, k). The z sweep completes corners in the same way. Intermediate
         * values at unswept physical ghosts are overwritten by later sweeps.
         * At a physical/inter-rank edge the transverse source can instead be a
         * value received from a neighbour, already supplied by fillHalo().
         *
         * Kernels run in the field execution space; no map data reaches the host.
         * Each sweep completes before the next one reads its results.
         * @pre The field has one ghost layer and at least two owned cells per
         * axis, as required by checkedLayout(); fillHalo() has completed.
         */
        inline void extendPhysical(ScalarField& field) {
            const auto& layout = field.getLayout();

            // Periodic layouts need no extension
            if (layout.isAllPeriodic_m)
                return;

            const auto local = layout.getLocalNDIndex(), domain = layout.getDomain();
            const auto view = field.getView();
            for (unsigned axis = 0; axis < 3; ++axis) {
                // Identify ranks touching each global boundary
                const bool low  = local[axis].first() == domain[axis].first();
                const bool high = local[axis].last() == domain[axis].last();
                // Last allocated index, including the upper ghost layer.
                const int end = int(view.extent(axis)) - 1;
                // Include transverse ghosts so edges and corners are filled too.
                Kokkos::parallel_for(
                    "Extend physical material", allocationPolicy(field),
                    KOKKOS_LAMBDA(int i, int j, int k) {
                        CellIndex source(i, j, k);
                        if (low && source[axis] == 0)
                            source[axis] = 1;
                        else if (high && source[axis] == end)
                            source[axis] = end - 1;
                        else
                            // Skip this cell's invocation, not the rest of the axis sweep.
                            return;
                        // Clamp only this axis; retain both transverse indices.
                        view(i, j, k) = view(source[0], source[1], source[2]);
                    });
                // The next axis can read ghosts written by this sweep.
                Kokkos::fence();
            }
        }

        /** @brief Fill inter-rank/periodic ghosts, then extend physical boundaries.
         * @details At a physical/inter-rank edge, extendPhysical() reads a
         * transverse ghost supplied by fillHalo(). Reversing these calls would
         * let the physical sweep copy a value before that exchange supplies it.
         */
        inline void synchronize(ScalarField& field) {
            field.fillHalo();
            extendPhysical(field);
        }
        /// Sample owned cells on the device; the functor receives global coordinates.
        template <class Map>
        void rasterize(ScalarField& field, Map map) {
            const auto coordinates = detail::coordinates(field.get_mesh(), field.getLayout());
            const auto view        = field.getView();
            Kokkos::parallel_for(
                "Rasterize material", field.getFieldRangePolicy(),
                KOKKOS_LAMBDA(int i, int j, int k) {
                    const auto global = coordinates.global(i, j, k);
                    bool marked       = false;
                    view(i, j, k)     = map(global, coordinates.centre(global), marked);
                });
        }

    }  // namespace detail

    /** @brief Own a static distributed cell map on the caller's mesh and layout.
     * @details Construction is collective. The mesh/layout and all retained views
     * must outlive their uses and end before Kokkos finalization
     * [kokkos2026ViewProgrammingGuide, lines 293-319]. Storage accepts any value;
     * the solver formulation specifies which materials its update supports.
     */
    class MaterialOwner {
    public:
        /** @brief Rasterize owned cells, and exchange/extend ghosts.
         * @param mesh Uniform Cartesian geometry; the caller keeps it alive while this
         * owner or any solver using its fields exists.
         * @param layout Distributed ownership and periodicity, with the same lifetime
         * as mesh. Every local axis must contain at least two owned cells.
         * @param map Device-copyable double(CellIndex, Position, bool&) functor. It must
         * return the same bits at the same coordinates and set its marked flag.
         * @details Map storage stays in the field execution space.
         */
        template <class Map>
        MaterialOwner(Mesh& mesh, Layout& layout, const Map& map)
            : epsilon_m(mesh, detail::checkedLayout(layout)) {
            detail::rasterize(epsilon_m, map);
            detail::synchronize(epsilon_m);
        }
        MaterialOwner(const MaterialOwner&)            = delete;
        MaterialOwner& operator=(const MaterialOwner&) = delete;
        MaterialOwner(MaterialOwner&&)                 = delete;
        MaterialOwner& operator=(MaterialOwner&&)      = delete;

        /// Complete immutable allocation, including every face, edge and corner ghost.
        ScalarFieldConstView epsilon() const { return epsilon_m.getView(); }
        Mesh& mesh() const { return epsilon_m.get_mesh(); }
        Layout& layout() const { return epsilon_m.getLayout(); }

        /// Local allocation spans, including halos, in bytes.
        std::uint64_t allocatedBytes() const {
            std::uint64_t spans = epsilon().span();
            return spans * sizeof(double);
        }

    private:
        ScalarField epsilon_m;
    };
}  // namespace chdr::material
