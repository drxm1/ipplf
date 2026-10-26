#include "Ippl.h"

#include "../MaxwellSolvers/BoundaryOracles.h"
#include "LinearSolvers/Multigrid.h"
#include "TestUtils.h"
#include "gtest/gtest.h"

namespace {
    using Axes  = std::array<bool, 3>;
    using Mesh  = ippl::UniformCartesian<double, 3>;
    using Field = ippl::Field<double, 3, Mesh, Mesh::DefaultCentering>;
    constexpr std::array<int, 3> cells{6, 5, 4};

    /// Every periodicity pattern is checked along each possible decomposition axis.
    Axes periodicAxes(unsigned mask) {
        return {bool(mask & 1), bool(mask & 2), bool(mask & 4)};
    }

    /// A small non-cubic box; every local extent is at least two on two ranks.
    struct Box {
        ippl::NDIndex<3> domain{ippl::Index(cells[0]), ippl::Index(cells[1]),
                                ippl::Index(cells[2])};
        ippl::FieldLayout<3> layout;
        Mesh mesh{domain, ippl::Vector<double, 3>(1), ippl::Vector<double, 3>(0)};
        Field field;

        Box(Axes periodic, unsigned split)
            : layout(MPI_COMM_WORLD, domain, {split == 0, split == 1, split == 2}, periodic)
            , field(mesh, layout) {}
    };

    /** @brief Map a local allocation index to its owned global source.
     * The independent haloSourceIndex oracle returns -1 beyond an open face.
     * Its indices include the lower ghost offset; the linear label is unique
     * for each global owned cell and exactly representable as a double.
     */
    int source(const ippl::NDIndex<3>& local, Axes periodic, int i, int j, int k) {
        const std::array<int, 3> index{i, j, k};
        int label = 0, stride = 1;
        for (unsigned d = 0; d < 3; ++d) {
            int g = open_boundary_oracle::haloSourceIndex(index[d] + local[d].first(), cells[d],
                                                          periodic[d]);
            if (g < 0) {
                return -1;
            }
            label += (g - 1) * stride;
            stride *= cells[d];
        }
        return label;
    }

    /// Fills owned cells with unique integer labels, calls fillHalo(),
    /// then compares every allocated cell with what is expected.
    void checkFill(Box& box, Axes periodic) {
        auto mirror = Kokkos::create_mirror_view(box.field.getView());
        auto local  = box.layout.getLocalNDIndex();
        Kokkos::deep_copy(mirror, -1.0);
        nestedViewLoop(mirror, 1, [&](int i, int j, int k) {
            mirror(i, j, k) = source(local, periodic, i, j, k);
        });
        Kokkos::deep_copy(box.field.getView(), mirror);
        box.field.fillHalo();
        Kokkos::deep_copy(mirror, box.field.getView());
        nestedViewLoop(mirror, 0, [&](int i, int j, int k) {
            EXPECT_EQ(mirror(i, j, k), source(local, periodic, i, j, k))
                << "local index " << i << ',' << j << ',' << k;
        });
    }

    /// Distinct positive integer deposits; owned cells receive no initial deposit.
    double deposit(const ippl::NDIndex<3>& local, int rank, int i, int j, int k) {
        if (i > 0 && i <= int(local[0].length()) && j > 0 && j <= int(local[1].length()) && k > 0
            && k <= int(local[2].length())) {
            return 0;
        }
        return 1 + i + 10 * j + 100 * k + 1000 * rank;
    }

    /** @brief Sum every rank's ghost deposits into the oracle's global destinations.
     * No neighbour list or fill operation enters this reference. In particular,
     * an open-face deposit has no destination even at a wrap/seam corner.
     */
    auto expectedDeposits(Box& box, Axes periodic) {
        std::array<double, cells[0] * cells[1] * cells[2]> sums{};
        for (int rank = 0; rank < box.layout.comm.size(); ++rank) {
            const auto local = box.layout.getLocalNDIndex(rank);
            ippl::detail::nestedLoop<3>(
                [](unsigned) {
                    return 0;
                },
                [&](unsigned d) {
                    return local[d].length() + 2;
                },
                [&](int i, int j, int k) {
                    const int target = source(local, periodic, i, j, k);
                    if (target >= 0) {
                        sums[target] += deposit(local, rank, i, j, k);
                    }
                });
        }
        return sums;
    }

    /// Check exact sums in owned cells, including simultaneous wraps and rank seams.
    void checkAccumulation(Box& box, Axes periodic) {
        const auto expected = expectedDeposits(box, periodic);
        auto local          = box.layout.getLocalNDIndex();
        auto mirror         = Kokkos::create_mirror_view(box.field.getView());
        nestedViewLoop(mirror, 0, [&](int i, int j, int k) {
            mirror(i, j, k) = deposit(local, box.layout.comm.rank(), i, j, k);
        });
        Kokkos::deep_copy(box.field.getView(), mirror);
        box.field.accumulateHalo();
        Kokkos::deep_copy(mirror, box.field.getView());
        nestedViewLoop(mirror, 1, [&](int i, int j, int k) {
            EXPECT_EQ(mirror(i, j, k), expected[source(local, periodic, i, j, k)])
                << "local index " << i << ',' << j << ',' << k;
        });
    }
}  // namespace

TEST(PeriodicAxes, FillEveryGhost) {
    for (unsigned mask = 0; mask < 8; ++mask) {
        for (unsigned split = 0; split < 3; ++split) {
            SCOPED_TRACE(::testing::Message() << "mask " << mask << ", split " << split);
            const auto periodic = periodicAxes(mask);
            Box box(periodic, split);
            EXPECT_EQ(box.layout.isAllPeriodic_m, mask == 7);
            EXPECT_EQ(box.layout.isAnyPeriodic(), mask != 0);
            for (unsigned d = 0; d < 3; ++d) {
                EXPECT_EQ(box.layout.isPeriodic(d), periodic[d]);
            }
            checkFill(box, periodic);
        }
    }
}

TEST(PeriodicAxes, BooleanInterface) {
    for (bool periodic : {false, true}) {
        Box box(periodicAxes(0), 2);
        box.layout =
            ippl::FieldLayout<3>(MPI_COMM_WORLD, box.domain, {false, false, true}, periodic);
        checkFill(box, periodicAxes(periodic ? 7 : 0));
    }
}

TEST(PeriodicAxes, AccumulateEveryGhost) {
    for (unsigned mask = 0; mask < 8; ++mask) {
        for (unsigned split = 0; split < 3; ++split) {
            SCOPED_TRACE(::testing::Message() << "mask " << mask << ", split " << split);
            const auto periodic = periodicAxes(mask);
            Box box(periodic, split);
            checkAccumulation(box, periodic);
        }
    }
}

TEST(PeriodicAxes, MixedMultigridIsRejected) {
    Box box(periodicAxes(6), 2);
    auto op = [](Field& field) {
        return -laplace(field);
    };
    ippl::multigrid_preconditioner<Field, decltype(op)> multigrid(std::move(op));
    EXPECT_THROW(multigrid.init_fields(box.field), IpplException);
}

/// Uniform layouts remain accepted, with enough cells for coarse levels on both rank counts.
TEST(PeriodicAxes, UniformMultigridIsAccepted) {
    const ippl::NDIndex<3> domain(ippl::Index(16), ippl::Index(16), ippl::Index(16));
    for (bool periodic : {false, true}) {
        ippl::FieldLayout<3> layout(MPI_COMM_WORLD, domain, {false, false, true}, periodic);
        Mesh mesh(domain, ippl::Vector<double, 3>(1), ippl::Vector<double, 3>(0));
        Field field(mesh, layout);
        auto op = [](Field& values) {
            return -laplace(values);
        };
        ippl::multigrid_preconditioner<Field, decltype(op)> multigrid(std::move(op));
        EXPECT_NO_THROW(multigrid.init_fields(field));
    }
}

int main(int argc, char* argv[]) {
    ippl::initialize(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    ippl::finalize();
    return result;
}
