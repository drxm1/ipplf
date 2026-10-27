/** @file Util.h
 * @brief Small utilities of the ChDR core that belong to no other unit.
 */
#pragma once

#include "Ippl.h"

#include <iostream>
#include <string_view>

#include "Utility/IpplException.h"

namespace chdr::util {
    /** @brief Report rank, context and message, then stop the MPI job.
     * @details Utility/IpplException.h defines a global exception type independent of
     * std::exception; callers catch both types explicitly and retain IPPL until this
     * handler completes. No recovery or MPI synchronization is attempted.
     * @return Failure status only if MPI_Abort unexpectedly returns.
     */
    inline int abortFailure(std::string_view where, const char* what) {
        std::cerr << "rank " << ippl::Comm->rank() << " [" << where << "]: " << what << std::endl;
        ippl::Comm->abort(1);
        return 1;
    }
}  // namespace chdr::util
