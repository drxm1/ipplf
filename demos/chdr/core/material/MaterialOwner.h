/** @file MaterialOwner.h
 * @brief Immutable distributed cell permittivity and lazily derived coefficient fields.
 */
#pragma once
#include "Ippl.h"

#include <algorithm>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <typeindex>
#include <unordered_map>

namespace chdr::material {
    using Mesh                 = ippl::UniformCartesian<double, 3>;
    using Layout               = ippl::FieldLayout<3>;
    using ScalarField          = ippl::Field<double, 3, Mesh, Mesh::DefaultCentering>;
    using ScalarFieldConstView = ScalarField::view_type::const_type;
    using CellIndex            = ippl::Vector<int, 3>;
    using Position             = ippl::Vector<double, 3>;
    using Policy               = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;

    /** @brief Global statistics of the owned relative-permittivity samples in epsilon_m, excluding
     * ghosts.
     * @details Nonfinite cells are excluded from extrema and clear allFinite.
     * Unsigned checksum sums are commutative, including unsigned integer overflow.
     *
     * MapStatistics holds a summary of the material map.
     *
     * It considers the `relative permittivity` values stored in `epsilon_m`,
     * sampled by a map functor at cell centres.
     */
    struct MapStatistics {
        /// Extrema of finite values, e.g. {1, 2.13} if it contains both vacuum and the 2.13=epsilon
        /// material
        double minimum_m = 0, maximum_m = 0;

        /// Whether every value is finite - no NaN or infinity
        bool allFinite_m = true;

        /// count of cells marked by the map functor
        std::uint64_t marked_m = 0;

        /// fingerprint of marked cell locations
        std::uint64_t markedChecksum_m = 0;

        /// fingerprint incorporating locations and material values
        std::uint64_t mapChecksum_m = 0;

        /// A homogeneous map has finite, bitwise equal extrema.
        bool homogeneous() const {
            return allFinite_m
                   && Kokkos::bit_cast<std::uint64_t>(minimum_m)
                          == Kokkos::bit_cast<std::uint64_t>(maximum_m);
        }
    };

    /** @brief Sample a shape at each cell centre and choose its relative permittivity.
     * @details Geometry supplies device-callable contains(x, y, z). The inside and
     * outside values are stored exactly as supplied; this adapter does not average
     * partly covered cells. Marked records geometric membership even when both
     * material values are equal. This adapter uses position rather than cell index.
     */
    template <class Geometry>
    struct TwoMaterialMap {
        Geometry geometry_m;
        double inside_m, outside_m;
        KOKKOS_INLINE_FUNCTION double operator()(CellIndex, Position centre, bool& marked) const {
            marked = geometry_m.contains(centre[0], centre[1], centre[2]);
            return marked ? inside_m : outside_m;
        }
    };

    namespace detail {
        /// SplitMix64 finalizer, preserving the mesh/ occupied-ID checksum.
        KOKKOS_INLINE_FUNCTION std::uint64_t hash(std::uint64_t id) {
            id += UINT64_C(0x9e3779b97f4a7c15);
            id = (id ^ (id >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
            id = (id ^ (id >> 27)) * UINT64_C(0x94d049bb133111eb);
            return id ^ (id >> 31);
        }

        /** @brief Device-copyable index metadata shared by sampling and verification.
         * @details The centre expression retains its arithmetic order, so a geometry
         * predicate sees exactly the mesh demo's original boundary samples.
         */
        struct Coordinates {
            CellIndex first_m, count_m, lower_m;
            Position origin_m, spacing_m;
            Kokkos::Array<bool, 3> periodic_m;
            KOKKOS_INLINE_FUNCTION CellIndex global(int i, int j, int k) const {
                return CellIndex(i + first_m[0] - 1, j + first_m[1] - 1, k + first_m[2] - 1);
            }
            KOKKOS_INLINE_FUNCTION Position centre(CellIndex global) const {
                Position result;
                for (unsigned d = 0; d < 3; ++d)
                    result[d] = origin_m[d] + (global[d] + 0.5) * spacing_m[d];
                return result;
            }
            /** @brief Find the global owned cell whose value belongs in an allocation cell.
             * @details Owned cells and inter-rank ghosts already have an index inside
             * the global domain, so their index stays unchanged. Only indices beyond
             * a global boundary wrap (periodic) or clamp (physical). Applying the
             * rule independently to all axes also handles edges and corners.
             * @pre The input is at most one ghost layer outside the global domain;
             * adding one domain length then makes the remainder nonnegative.
             */
            KOKKOS_INLINE_FUNCTION CellIndex source(CellIndex global) const {
                for (unsigned d = 0; d < 3; ++d) {
                    const int relative = global[d] - lower_m[d];
                    global[d]          = lower_m[d]
                                + (periodic_m[d] ? (relative + count_m[d]) % count_m[d]
                                                 : Kokkos::clamp(relative, 0, count_m[d] - 1));
                }
                return global;
            }
            KOKKOS_INLINE_FUNCTION std::uint64_t id(CellIndex global) const {
                return (std::uint64_t(global[0] - lower_m[0]) * count_m[1] + global[1] - lower_m[1])
                           * count_m[2]
                       + global[2] - lower_m[2];
            }
        };

        /// Metadata only: no field allocation or host copy of the distributed map.
        inline Coordinates coordinates(Mesh& mesh, Layout& layout) {
            Coordinates result;
            const auto local = layout.getLocalNDIndex(), domain = layout.getDomain();
            for (unsigned d = 0; d < 3; ++d) {
                result.first_m[d]    = local[d].first();
                result.count_m[d]    = domain[d].length();
                result.lower_m[d]    = domain[d].first();
                result.periodic_m[d] = layout.isPeriodic(d);
            }
            result.origin_m  = mesh.getOrigin();
            result.spacing_m = mesh.getMeshSpacing();
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
         * Each axis wraps periodically or extends its physical boundary independently.
         * Internal rank boundaries use halo exchange.
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
                if (layout.isPeriodic(axis))
                    continue;
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
        /** @brief Store cell-centred relative permittivity and summarize owned samples.
         * @details Each rank visits only its owned cells; halos would duplicate samples
         * in the global statistics. The map receives the global cell index and centre.
         * Its value is stored unchanged, including NaN or infinity. Its separate bool
         * output marks cells for the occupied count and occupied-location fingerprint;
         * unmarked cells still contribute to the extrema and full-map fingerprint.
         *
         * The full-map fingerprint combines each global ID with the stored value's
         * exact bits. It can reveal value changes that the occupied mask cannot, but
         * a hash match is not a cell-by-cell equality proof. Unsigned sums use the same
         * wraparound arithmetic on each decomposition, without floating-point sums.
         * Nonfinite values set a separate flag and never enter the extrema reduction.
         * Only the reduced scalars reach the host; MPI makes the result rank identical.
         *
         * rasterize() fills `epsilon_m` and computes its statistics
         */
        template <class Map>
        MapStatistics rasterize(ScalarField& field, Map map) {
            const auto coordinates = detail::coordinates(field.get_mesh(), field.getLayout());
            const auto view        = field.getView();
            std::uint64_t marked = 0, occupiedHash = 0, valueHash = 0, nonfinite = 0;
            double minimum = 0, maximum = 0;
            Kokkos::parallel_reduce(
                "Rasterize material", field.getFieldRangePolicy(),
                KOKKOS_LAMBDA(int i, int j, int k, std::uint64_t& count, std::uint64_t& mask,
                              std::uint64_t& values, std::uint64_t& bad, double& low,
                              double& high) {
                    const auto global  = coordinates.global(i, j, k);
                    bool inside        = false;
                    const double value = map(global, coordinates.centre(global), inside);
                    view(i, j, k)      = value;
                    const auto hash    = detail::hash(coordinates.id(global));
                    if (inside) {
                        ++count;
                        mask += hash;
                    }
                    values += detail::hash(hash ^ Kokkos::bit_cast<std::uint64_t>(value));
                    if (!Kokkos::isfinite(value))
                        ++bad;
                    else {
                        low  = Kokkos::min(low, value);
                        high = Kokkos::max(high, value);
                    }
                },
                marked, occupiedHash, valueHash, nonfinite, Kokkos::Min<double>(minimum),
                Kokkos::Max<double>(maximum));

            // Kokkos initializes Min/Max with their reducer identities, not the
            // zero-valued output variables above. If a rank has no finite sample,
            // its reducer identities remain; the nonfinite flag prevents homogeneity.
            std::uint64_t sums[3]{marked, occupiedHash, valueHash};
            // Negating the maximum lets one MPI_MIN reduce both extrema. A -1
            // flag wins over 0 if any rank encountered a nonfinite sample.
            double extrema[3]{minimum, -maximum, nonfinite ? -1.0 : 0.0};
            MPI_Allreduce(MPI_IN_PLACE, sums, 3, MPI_UINT64_T, MPI_SUM, field.getLayout().comm);
            MPI_Allreduce(MPI_IN_PLACE, extrema, 3, MPI_DOUBLE, MPI_MIN, field.getLayout().comm);
            return {extrema[0], -extrema[1], extrema[2] == 0, sums[0], sums[1], sums[2]};
        }

    }  // namespace detail

    /** @brief Own a static cell map and cache coefficient fields by rule name.
     * @details Construction and the first derive call for a rule are collective.
     * Every rank must use the same map/rule and call order. The caller's mesh and
     * layout outlive this object and any solver using its views. Returned const
     * views retain allocations [kokkos2026ViewProgrammingGuide, lines 293-306];
     * all handles must be destroyed before Kokkos finalization. Map values are
     * unrestricted: a solver formulation, not this storage class, sets validity.
     */
    class MaterialOwner {
    public:
        /** @brief Rasterize owned cells, exchange/extend ghosts, and verify placement.
         * @param mesh Uniform Cartesian geometry; the caller keeps it alive while this
         * owner or any solver using its fields exists.
         * @param layout Distributed ownership and periodicity, with the same lifetime
         * as mesh. Every local axis must contain at least two owned cells.
         * @param map Device-copyable double(CellIndex, Position, bool&) functor. It must
         * return the same bits at the same coordinates and set its marked flag.
         * @details One owned-cell device reduction computes all statistics; only
         * reduced scalars cross to the host. Verification accepts unchanged NaN bits.
         */
        template <class Map>
        MaterialOwner(Mesh& mesh, Layout& layout, const Map& map)
            : epsilon_m(mesh, detail::checkedLayout(layout)) {
            statistics_m = detail::rasterize(epsilon_m, map);
            detail::synchronize(epsilon_m);
            verify(map);
        }
        MaterialOwner(const MaterialOwner&)            = delete;
        MaterialOwner& operator=(const MaterialOwner&) = delete;
        MaterialOwner(MaterialOwner&&)                 = delete;
        MaterialOwner& operator=(MaterialOwner&&)      = delete;

        /// Complete immutable allocation, including every face, edge and corner ghost.
        ScalarFieldConstView epsilon() const { return epsilon_m.getView(); }
        /// Cached global summary of owned epsilon samples; this accessor performs no collective.
        const MapStatistics& statistics() const { return statistics_m; }
        Mesh& mesh() const { return epsilon_m.get_mesh(); }
        Layout& layout() const { return epsilon_m.getLayout(); }

        /** @brief Return the cached common relative permittivity, without communication.
         * @throws std::invalid_argument If the cached extrema fail homogeneous().
         * @details Returning the sampled minimum preserves its bits; no mean is formed.
         */
        double homogeneousValue() const {
            if (!statistics_m.homogeneous())
                throw std::invalid_argument("Material is not homogeneous");
            return statistics_m.minimum_m;
        }

        /** @brief Compute a coefficient field from epsilon once, then reuse its allocation.
         * @details Rule supplies static Name and a device-callable
         * double operator()(ScalarFieldConstView,int,int,int). Indices are local
         * allocation indices; the input already includes valid epsilon ghosts.
         * The rule runs on owned cells before its output ghosts are synchronized.
         * Rule::ForwardAxis is -1 for a cell-centred rule, or 0, 1, 2 for a rule
         * using the cell and its next neighbour along that axis. Forward rules
         * receive the lower-boundary correction described in fixLowerFace().
         *
         * Rules must be stateless: the cache key is Name, with a C++ type check,
         * not the contents of a rule object. Repeating a request returns the same
         * const view without reevaluation. A different type with the same name
         * is rejected instead of silently returning the first type's field.
         *
         * The first request is collective; all ranks must use the same rule and
         * call order. Later requests only consult the local cache. Filling this
         * mutable cache changes neither epsilon nor any previously returned field.
         * Returned views retain their allocation as described by the class contract.
         * @throws std::invalid_argument If Name already belongs to another rule type.
         */
        template <class Rule>
        ScalarFieldConstView derive(const Rule& rule) const {
            if (cache_m.count(Rule::Name))
                return derived<Rule>();
            auto field        = std::make_unique<ScalarField>(mesh(), layout());
            const auto output = field->getView();
            const auto input  = epsilon();
            Kokkos::parallel_for(
                "Derive material coefficient", field->getFieldRangePolicy(),
                KOKKOS_LAMBDA(int i, int j, int k) { output(i, j, k) = rule(input, i, j, k); });
            detail::synchronize(*field);
            fixLowerFace(*field, rule);
            cache_m.emplace(Rule::Name, CacheEntry{typeid(Rule), std::move(field)});
            return derived<Rule>();
        }

        /** @brief Retrieve an already cached coefficient without communication or allocation.
         * @details Unlike derive(), this never builds a missing field. The type
         * guard makes Name a unique association with the rule that created it.
         * @throws std::out_of_range If no field with Rule::Name has been derived.
         * @throws std::invalid_argument If that name belongs to another rule type.
         */
        template <class Rule>
        ScalarFieldConstView derived() const {
            const auto found = cache_m.find(Rule::Name);
            if (found == cache_m.end())
                throw std::out_of_range("Material rule has not been derived");
            if (found->second.type_m != std::type_index(typeid(Rule)))
                throw std::invalid_argument("Material rule name belongs to another type");
            return found->second.field_m->getView();
        }

        /// Local allocation spans, including halos and cached fields, in bytes.
        std::uint64_t allocatedBytes() const {
            std::uint64_t spans = epsilon().span();
            for (const auto& entry : cache_m)
                spans += entry.second.field_m->getView().span();
            return spans * sizeof(double);
        }

    private:
        /// Keep the owning field for halo operations and its rule identity for safe lookup.
        struct CacheEntry {
            std::type_index type_m;
            std::unique_ptr<ScalarField> field_m;
        };
        ScalarField epsilon_m;
        MapStatistics statistics_m;
        mutable std::unordered_map<std::string, CacheEntry> cache_m;

    public:
        // NVCC requires public enclosing functions for extended device lambdas.
        // These construction/cache helpers are not mutation entry points for callers.
        /** @brief Collectively check owned values and every face, edge and corner ghost.
         * @details Resolve each allocation index to its global source cell, then
         * reevaluate the map there using the same centre arithmetic as rasterize().
         * Physical ghosts sample the boundary cell's centre, not a point outside
         * the domain. The map must therefore return the same bits on repeated calls.
         *
         * Compare representations rather than floating-point values: an unchanged
         * NaN passes, while a changed NaN payload fails. This checks placement and
         * copying, not whether the values are suitable for a solver formulation.
         * The marked flag does not enter this check; the field stores only values.
         * One MPI sum makes every rank throw if any rank finds a mismatch.
         * @throws std::runtime_error If any allocated sample has unexpected bits.
         */
        template <class Map>
        void verify(Map map) const {
            const auto coordinates = detail::coordinates(mesh(), layout());
            const auto view        = epsilon();
            std::uint64_t bad      = 0;
            Kokkos::parallel_reduce(
                "Verify complete material", detail::allocationPolicy(epsilon_m),
                KOKKOS_LAMBDA(int i, int j, int k, std::uint64_t& errors) {
                    const auto source     = coordinates.source(coordinates.global(i, j, k));
                    bool marked           = false;
                    const double expected = map(source, coordinates.centre(source), marked);
                    errors += Kokkos::bit_cast<std::uint64_t>(view(i, j, k))
                              != Kokkos::bit_cast<std::uint64_t>(expected);
                },
                bad);
            MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_UINT64_T, MPI_SUM, layout().comm);
            if (bad)
                throw std::runtime_error("Material field/halo verification failed");
        }

        /** @brief Replace the copied lower ghost coefficient with the actual boundary pair.
         * @details A forward coefficient at local index i uses epsilon at i and i+1.
         * At the lower physical boundary, ghost index 0 therefore represents the
         * face between ghost 0 and owned cell 1. synchronize() copied coefficient
         * 1 there, but coefficient 1 used cells 1 and 2, which may differ. Reevaluate
         * the rule at index 0 using the already extended epsilon field.
         * This index correction applies to two-centre face coefficients such as
         * harmonic averaging [eymard2000FiniteVolume, section 7.1, Example 7.1, p. 22].
         *
         * The upper boundary needs no analogous repair: its last owned coefficient
         * already used the boundary cell and its equal-valued outer epsilon ghost.
         * The coefficient extension copies that pair's result farther outward.
         * Transverse ghosts must be included to correct boundary edges and corners
         * too. Only ranks touching this global lower face perform the repair;
         * periodic and inter-rank faces already have their exchanged coefficients.
         * @pre The field was produced by this rule and synchronize() has completed.
         */
        template <class Rule>
        void fixLowerFace(ScalarField& field, Rule rule) const {
            if constexpr (Rule::ForwardAxis >= 0) {
                constexpr unsigned Axis = Rule::ForwardAxis;
                if (layout().isPeriodic(Axis)
                    || layout().getLocalNDIndex()[Axis].first()
                           != layout().getDomain()[Axis].first())
                    return;
                const auto output = field.getView();
                const auto input  = epsilon();
                Kokkos::parallel_for(
                    "Derive lower boundary face", detail::allocationPolicy(field),
                    KOKKOS_LAMBDA(int i, int j, int k) {
                        if (CellIndex(i, j, k)[Axis] == 0)
                            output(i, j, k) = rule(input, i, j, k);
                    });
                Kokkos::fence();
            }
        }
    };
}  // namespace chdr::material
