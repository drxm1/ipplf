/** @file TestMain.cpp
 * @brief Shared IPPL/GoogleTest lifetime and rank-global exit status.
 * @details Test fixtures release their fields before RUN_ALL_TESTS returns and
 * before IPPL finalizes Kokkos/MPI. All ranks must use the same test-selection,
 * repetition and shuffle options. Unseeded shuffling uses seed 1 on every rank;
 * explicit nonzero seeds are preserved.
 */
#include "Ippl.h"

#include <exception>
#include <gtest/gtest.h>
#include <iostream>

#include "VacuumFailure.h"

int main(int argc, char* argv[]) {
    ippl::initialize(argc, argv);
    int result = 1;
    try {
        ::testing::InitGoogleTest(&argc, argv);
        // Unexpected per-rank exceptions must reach the MPI-abort handler.
        GTEST_FLAG_SET(catch_exceptions, false);
        if (GTEST_FLAG_GET(shuffle) && GTEST_FLAG_GET(random_seed) == 0)
            GTEST_FLAG_SET(random_seed, 1);
        result = RUN_ALL_TESTS();
    } catch (const IpplException& error) {
        result = chdr::vacuum::abortFailure(error.where(), error.what());
    } catch (const std::exception& error) {
        result = chdr::vacuum::abortFailure("RUN_ALL_TESTS", error.what());
    }
    MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_INT, MPI_MAX, *ippl::Comm);
    ippl::finalize();
    return result;
}
