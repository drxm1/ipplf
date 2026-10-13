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

        /// Exchange owned values into all periodic face, edge and corner ghosts.
        inline void synchronize(ScalarField& field) {
            field.fillHalo();
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
        /** @brief Rasterize owned cells and exchange periodic ghosts.
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
