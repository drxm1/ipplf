/** @file OuterLayerTest.h
 * @brief Host-side sampling for the independent polynomial outer-layer check.
 */
#pragma once
#include "Ippl.h"

#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

#include "MaxwellSolvers/StandardFDTDSolver.h"
#include "OuterLayerOracles.h"
#include "gtest/gtest.h"

namespace outer_layer_test {
    using Mesh       = ippl::UniformCartesian<double, 3>;
    using Potentials = ippl::Field<ippl::Vector<double, 4>, 3, Mesh, Mesh::DefaultCentering>;
    using Fields     = ippl::Field<ippl::Vector<double, 3>, 3, Mesh, Mesh::DefaultCentering>;
    using Standard   = ippl::StandardFDTDSolver<Fields, Potentials, ippl::absorbing>;
    /// Test-only access to the fixed polynomial fixture's time step.
    struct Solver : Standard {
        using Standard::Standard;
        void usePolynomialTimeStep() { this->dt = 0.5029188806044014; }
    };
    inline constexpr std::array<double, 3> spacing = {1., .8, 1.25};

    /// Nonperiodic box with a z split by default; self-communicator copies provide rank witnesses.
    struct Box {
        std::array<int, 3> size;
        ippl::NDIndex<3> domain;
        ippl::FieldLayout<3> layout;
        Mesh mesh;
        Potentials source;
        Fields electric, magnetic;
        explicit Box(std::array<int, 3> n = {8, 6, 10}, std::array<double, 3> h = spacing,
                     unsigned splitAxis = 2, MPI_Comm communicator = MPI_COMM_WORLD)
            : size(n)
            , domain(ippl::Index(n[0]), ippl::Index(n[1]), ippl::Index(n[2]))
            , layout(communicator, domain, {splitAxis == 0, splitAxis == 1, splitAxis == 2}, false)
            , mesh(domain, {h[0], h[1], h[2]}, {0., 0., 0.})
            , source(mesh, layout)
            , electric(mesh, layout)
            , magnetic(mesh, layout) {
            source = 0.;
        }
    };

    /// Visit local indices and their global labels, optionally including the entire ghost layer.
    template <typename Field, typename Function>
    void visit(const Field& field, bool ghosts, Function function) {
        const auto local = field.getLayout().getLocalNDIndex();
        const int start = ghosts ? 0 : 1, extra = ghosts ? 1 : 0;
        for (int i = start; i <= int(local[0].length()) + extra; ++i)
            for (int j = start; j <= int(local[1].length()) + extra; ++j)
                for (int k = start; k <= int(local[2].length()) + extra; ++k)
                    function(i, j, k,
                             std::array<int, 3>{local[0].first() + i - 1, local[1].first() + j - 1,
                                                local[2].first() + k - 1});
    }

    /** @brief Capture owned field components in global x,y,z,component order without arithmetic.
     * @details Exactly one rank owns each cell. Bitwise OR with zero-filled nonowner entries
     * preserves binary64 samples, including signed zeros. The field communicator also makes
     * this work for a rank-zero-only MPI_COMM_SELF witness. Ghosts are excluded.
     */
    template <typename Field>
    std::vector<std::uint64_t> ownedFieldBits(const Field& field) {
        constexpr unsigned components = Field::value_type::dim;
        const auto n                  = field.getLayout().getDomain().length();
        std::vector<std::uint64_t> values(n[0] * n[1] * n[2] * components);
        const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), field.getView());
        visit(field, false, [&](int i, int j, int k, auto index) {
            const auto cell = (index[0] * n[1] + index[1]) * n[2] + index[2];
            for (unsigned q = 0; q < components; ++q)
                values[cell * components + q] = std::bit_cast<std::uint64_t>(host(i, j, k)[q]);
        });
        MPI_Allreduce(MPI_IN_PLACE, values.data(), int(values.size()), MPI_UINT64_T, MPI_BOR,
                      field.getLayout().comm);
        return values;
    }

    /** @brief Require field-array rank agreement against a nonzero serial absolute maximum.
     * @details The assertion compares every owned component within 1e-13 of the serial field's
     * maximum; nonfinite samples fail. Exact bit equality is recorded separately. Only the
     * rank owning the self-communicator reference calls this routine.
     */
    template <typename Field>
    void expectSameField(const std::vector<std::uint64_t>& observed, const Field& reference,
                         const char* channel) {
        const auto expected = ownedFieldBits(reference);
        ASSERT_EQ(observed.size(), expected.size());
        double error = 0, scale = 0;
        for (std::size_t i = 0; i < expected.size(); ++i) {
            const double a          = std::bit_cast<double>(observed[i]),
                         b          = std::bit_cast<double>(expected[i]);
            const double difference = std::abs(a - b);
            error                   = std::max(error, std::isfinite(difference)
                                                          ? difference
                                                          : std::numeric_limits<double>::infinity());
            scale                   = std::max(scale, std::abs(b));
        }
        ASSERT_GT(scale, 0.);
        EXPECT_LE(error / scale, 1.e-13);
        std::cout << "outer-layer " << channel << " rank_error=" << error / scale
                  << " serial_maximum=" << scale << " bitwise=" << (observed == expected) << '\n';
    }

    /// Load only owned samples; zero ghosts provide the shipped reconstruction's negative control.
    inline void load(Potentials& field, double time) {
        field     = 0.;
        auto host = field.getHostMirror();
        Kokkos::deep_copy(host, field.getView());
        const auto coefficients = open_boundary_oracle::polynomialFixture(false);
        visit(field, false, [&](int i, int j, int k, auto index) {
            const auto x = open_boundary_oracle::fieldPosition(index, spacing, false, false, 0);
            for (unsigned q = 0; q < 4; ++q)
                host(i, j, k)[q] =
                    open_boundary_oracle::polynomialPotential(coefficients, q, x, time, -1, false);
        });
        Kokkos::deep_copy(field.getView(), host);
    }

    /// Maximum across ranks; normalizations are analytic global maxima, never cancelling sums.
    inline double maximum(double local) {
        double result;
        MPI_Allreduce(&local, &result, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        return result;
    }

    /** @brief Relative maximum error of the inherited E or B on the outermost owned layer.
     * @details Both fields have time label zero for the linear-time oracle. The inherited B
     * is the raw current-level curl; no separate time-averaged B is exposed by this solver.
     * The positive analytic maximum over the full box defines the normalization.
     */
    inline double fieldError(const Fields& field, const Box& box, bool magnetic) {
        const auto c = open_boundary_oracle::polynomialFixture(false);
        auto host    = field.getHostMirror();
        Kokkos::deep_copy(host, field.getView());
        double error = 0.;
        visit(field, false, [&](int i, int j, int k, auto index) {
            if (index[0] != 0 && index[0] != box.size[0] - 1 && index[1] != 0
                && index[1] != box.size[1] - 1 && index[2] != 0 && index[2] != box.size[2] - 1)
                return;
            const auto x = open_boundary_oracle::fieldPosition(index, spacing, false, magnetic, 0);
            for (unsigned q = 0; q < 3; ++q)
                error = std::max(error, std::abs(host(i, j, k)[q]
                                                 - open_boundary_oracle::polynomialField(
                                                     c, magnetic, q, x, 0.)));
        });
        return maximum(error)
               / open_boundary_oracle::polynomialFieldMaximum(c, box.size, spacing, false, magnetic,
                                                              0.);
    }

    /// Check all owned and ghost samples against the analytic polynomial, including corners.
    inline double historyError(const Potentials& field, double time) {
        const auto c = open_boundary_oracle::polynomialFixture(false);
        auto host    = field.getHostMirror();
        Kokkos::deep_copy(host, field.getView());
        double error = 0., scale = 0.;
        visit(field, true, [&](int i, int j, int k, auto index) {
            const auto x = open_boundary_oracle::fieldPosition(index, spacing, false, false, 0);
            for (unsigned q = 0; q < 4; ++q) {
                const double expected =
                    open_boundary_oracle::polynomialPotential(c, q, x, time, -1, false);
                error = std::max(error, std::abs(host(i, j, k)[q] - expected));
                scale = std::max(scale, std::abs(expected));
            }
        });
        return maximum(error) / maximum(scale);
    }
}  // namespace outer_layer_test
