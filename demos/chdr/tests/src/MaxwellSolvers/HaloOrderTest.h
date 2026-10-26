#ifndef IPPL_TEST_HALO_ORDER_H
#define IPPL_TEST_HALO_ORDER_H

#include "Ippl.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "MaxwellSolvers/NonStandardFDTDSolver.h"
#include "MaxwellSolvers/StandardFDTDSolver.h"
#include "ShippedOrderSolver.h"
#include "gtest/gtest.h"

namespace fdtd_test {
    using Mesh      = ippl::UniformCartesian<double, 3>;
    using Potential = ippl::Field<ippl::Vector<double, 4>, 3, Mesh, Mesh::DefaultCentering>;
    using Field     = ippl::Field<ippl::Vector<double, 3>, 3, Mesh, Mesh::DefaultCentering>;
    template <ippl::fdtd_bc BC>
    using Standard = ippl::StandardFDTDSolver<Field, Potential, BC>;
    template <ippl::fdtd_bc BC>
    using NonStandard = ippl::NonStandardFDTDSolver<Field, Potential, BC>;

    /** @brief Own source-free solver storage with decomposition along z only.
     * @details MPI_COMM_SELF constructs an independent serial reference with
     * identical global labels. Every test must keep this storage alive longer
     * than its solver, whose mesh, layout and source pointers refer here.
     */
    struct Fields {
        ippl::NDIndex<3> domain;
        Mesh mesh;
        ippl::FieldLayout<3> layout;
        Potential source;
        Field electric, magnetic;
        Fields(std::array<int, 3> n, ippl::Vector<double, 3> h, bool periodic,
               MPI_Comm comm = MPI_COMM_WORLD)
            : Fields(n, h, std::array<bool, 3>{periodic, periodic, periodic}, comm) {}

        /// Select periodic axis pairs while retaining the same one-axis decomposition.
        Fields(std::array<int, 3> n, ippl::Vector<double, 3> h, std::array<bool, 3> periodic,
               MPI_Comm comm = MPI_COMM_WORLD)
            : domain(ippl::Index(n[0]), ippl::Index(n[1]), ippl::Index(n[2]))
            , mesh(domain, h, ippl::Vector<double, 3>(0))
            , layout(comm, domain, {false, false, true}, periodic)
            , source(mesh, layout)
            , electric(mesh, layout)
            , magnetic(mesh, layout) {
            source = 0;
        }
    };

    /** @brief Visit owned samples in global x,y,z,component order on the host.
     * @details Excludes every ghost. Global labels allow direct comparison with
     * a full serial field regardless of the rank's position in the z split.
     */
    template <class F, class Function>
    void owned(F& field, Function function) {
        const auto a   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), field.getView());
        const auto box = field.getLayout().getLocalNDIndex();
        const int g    = field.getNghost();
        for (std::size_t i = 0; i < box[0].length(); ++i)
            for (std::size_t j = 0; j < box[1].length(); ++j)
                for (std::size_t k = 0; k < box[2].length(); ++k)
                    function(i + box[0].first(), j + box[1].first(), k + box[2].first(),
                             a(i + g, j + g, k + g));
    }

    /** @brief Load all owned components from global labels and exchange them.
     * @details Physical exterior ghosts retain zero. Periodic field boundary
     * conditions are applied explicitly, matching the solver's initialization
     * contract; no potential evolution or Mur update is performed here.
     */
    template <class Function>
    void load(Potential& field, Function function) {
        auto a         = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), field.getView());
        const auto box = field.getLayout().getLocalNDIndex();
        const int g    = field.getNghost();
        for (std::size_t i = 0; i < box[0].length(); ++i)
            for (std::size_t j = 0; j < box[1].length(); ++j)
                for (std::size_t k = 0; k < box[2].length(); ++k)
                    for (unsigned q = 0; q < 4; ++q)
                        a(i + g, j + g, k + g)[q] =
                            function(i + box[0].first(), j + box[1].first(), k + box[2].first(), q);
        Kokkos::deep_copy(field.getView(), a);
        field.fillHalo();
        field.getFieldBC().apply(field);
    }
}  // namespace fdtd_test
#endif
