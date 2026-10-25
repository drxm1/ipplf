/** @file FieldContainer.cpp
 * @brief Ownership check and storage allocation of the shared field container.
 */
#include "core/base/FieldContainer.h"

#include <algorithm>
#include <stdexcept>

namespace chdr::solver {
    void checkLocalDomain(const ippl::FieldLayout<SpatialDim>& layout) {
        const auto& local    = layout.getLocalNDIndex();
        const int localValid = std::all_of(local.begin(), local.end(), [](const auto& axis) {
            return axis.length() >= 2;
        });
        int globalValid      = 0;
        MPI_Allreduce(&localValid, &globalValid, 1, MPI_INT, MPI_MIN, layout.comm);
        if (!globalValid)
            throw std::runtime_error(
                "Every rank must own at least two cells per local axis, required by IPPL halo "
                "exchange");
    }

    void initializeFieldStorage(ChDRFieldContainer& fields) {
        fields.initializeFields();
        // Scalar assignment zeros every source component in owned and halo cells.
        fields.getJ() = 0.0;
    }
}  // namespace chdr::solver
