/** @file VacuumFailure.h
 * @brief Shared fatal-error reporting for initialized vacuum MPI jobs.
 * @details Utility/IpplException.h defines a global exception type independent
 * of std::exception. Callers catch both types explicitly and retain IPPL until
 * this handler completes. No recovery or MPI synchronization is attempted.
 */
#pragma once

#include "Ippl.h"

#include <iostream>
#include <string_view>

#include "Utility/IpplException.h"

namespace chdr::vacuum {
    /// @brief Report rank, context and message, then stop the MPI job.
    /// @return Failure status only if MPI_Abort unexpectedly returns.
    inline int abortFailure(std::string_view where, const char* what) {
        std::cerr << "vacuum rank " << ippl::Comm->rank() << " [" << where << "]: " << what
                  << std::endl;
        ippl::Comm->abort(1);
        return 1;
    }
}  // namespace chdr::vacuum
