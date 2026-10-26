/** @file CoreTestSupport.h
 * @brief Small device-reduction and host-snapshot helpers shared by the core audits.
 */
#pragma once
#include <limits>

#include "core/base/SpatialDimension.h"
#include "core/material/MaterialOwner.h"
#include "core/solver/BoxSetup.h"
#include "core/solver/ChdrSolver.h"

static_assert(
    chdr::solver::SpatialDim == 3,
    "CoreTestSupport.h requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver::coreTest {
    /** @brief Construct the requested test path while retaining the caller's field owners.
     * @details epsilon equal to zero is a fixture sentinel for the owner-free vacuum
     * constructor, not a physical material value. Every other value constructs a
     * constant map and passes through the normal material validation; epsilon equal
     * to one therefore exercises a different construction path from the sentinel.
     * HomogeneousUpdate retains only the owner's scalar permittivity and allocates
     * its own source copy, so this temporary owner may die on return. The caller's
     * source, E/B fields, mesh and layout must outlive the returned solver.
     */
    template <class PotentialSolverBase>
    auto makeTestSolver(AJFourFieldsType& source, EMFieldsType& electric, EMFieldsType& magnetic,
                        EBReconstructionConvention reconstruction,
                        SpatialDiscretization reading = SpatialDiscretization::Collocated,
                        double epsilon                = 0) {
        using Solver = ChdrSolver<PotentialSolverBase>;
        if (epsilon == 0)
            return std::make_unique<Solver>(source, electric, magnetic, reconstruction, reading);
        const chdr::material::MaterialOwner material(
            source.get_mesh(), source.getLayout(),
            KOKKOS_LAMBDA(ippl::Vector<int, SpatialDim>, ippl::Vector<double, SpatialDim>,
                          bool& marked) {
                marked = false;
                return epsilon;
            });
        return std::make_unique<Solver>(source, electric, magnetic, reconstruction, reading,
                                        material);
    }
    /** @brief Reduce an already nonnegative cell observation over owned cells only.
     * @details The caller supplies a device-copyable scalar functor. This helper
     * performs no MPI operation, so each test retains its original collective
     * timing. A nonfinite observation becomes infinity and cannot disappear in
     * a maximum comparison. Halo cells are excluded by the field range policy.
     * @tparam FieldType IPPL field used only to select the owned-cell range.
     * @tparam ScalarObservation Device-copyable callable returning one nonnegative
     * scalar from three local array indices, including the halo offset.
     * @param field Field whose owned cells define the reduction domain.
     * @param observation Read-only cell measurement; captured inputs must be ready.
     */
    template <class FieldType, class ScalarObservation>
    double localMaximum(const FieldType& field, ScalarObservation observation) {
        const double infinity = std::numeric_limits<double>::infinity();
        double result         = 0;
        Kokkos::parallel_reduce(
            "core audit maximum", field.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k, double& maximum) {
                const double value = observation(i, j, k);
                maximum = Kokkos::max(maximum, Kokkos::isfinite(value) ? value : infinity);
            },
            Kokkos::Max<double>(result));
        return result;
    }

    /** @brief Preserve a field's current bytes in independent host mirror storage.
     * @details Returning a mirror is safe: its allocation is reference-counted
     * and contains values only. The owning mesh/layout and solver are never
     * copied or moved. Completion of deep_copy makes the snapshot ready to read.
     * @tparam FieldType IPPL field; the stored quantity and component count are generic.
     * @param field Source field whose complete allocation, including halos, is copied.
     */
    template <class FieldType>
    auto copyToHost(const FieldType& field) {
        auto host = field.getHostMirror();
        Kokkos::deep_copy(host, field.getView());
        return host;
    }

    /// @brief Combine an owned-cell maximum after the caller's local work has finished.
    template <class FieldType>
    double globalMaximum(const FieldType& field, double local) {
        MPI_Allreduce(MPI_IN_PLACE, &local, 1, MPI_DOUBLE, MPI_MAX, field.getLayout().comm);
        return local;
    }
}  // namespace chdr::solver::coreTest
