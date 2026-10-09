/** @file PeriodicSetup.cpp
 * @brief Construct a periodic mesh/layout owner without allocating its field views.
 */
#include "PeriodicSetup.h"

#include <cmath>
#include <stdexcept>

namespace chdr::solver {

    ChDRFieldContainer makePeriodicFieldContainer(std::array<int, SpatialDim> cells,
                                                  std::array<double, SpatialDim> lengths,
                                                  unsigned splitAxis) {
        // The axis along which MPI decomposition happens must be one of x,y,z.
        if (splitAxis >= SpatialDim) {
            throw std::invalid_argument("Split axis must be 0, 1 or 2");
        }

        ippl::NDIndex<SpatialDim> domain;
        ippl::Vector<double, SpatialDim> spacing, lower(0.0), upper;

        // Define along which axis the MPI decomposition happens
        std::array<bool, SpatialDim> decomposition{false, false, false};
        decomposition[splitAxis] = true;

        for (unsigned d = 0; d < SpatialDim; ++d) {
            // Validate each global cell count and the length of the whole box
            if (cells[d] < 2 || !std::isfinite(lengths[d]) || lengths[d] <= 0.0) {
                throw std::invalid_argument("Periodic box needs positive lengths and >= 2 cells");
            }
            // Index domain
            domain[d] = ippl::Index(cells[d]);
            upper[d]  = lengths[d];
            // Grid spacing
            spacing[d] = lengths[d] / cells[d];
        }
        // Leave the field allocation to the caller:
        // allocated fields retain pointers to this container's mesh and layout.
        return ChDRFieldContainer(spacing, lower, upper, decomposition, domain, lower, true);
    }
}  // namespace chdr::solver
