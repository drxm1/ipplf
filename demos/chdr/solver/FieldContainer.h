/** @file FieldContainer.h
 * @brief Field container types, the ownership check and the storage allocation shared by
 * every ChDR field problem.
 * @details The types are those of the FEL demo: FELFieldContainer owns one mesh, its
 * layout and the E, B and J storage; VField_t and SourceField_t are its three- and
 * four-component cell fields. This header is the only core header that names the FEL
 * headers, so a change of the container touches one place.
 */
#pragma once

#include "Ippl.h"

#include "FELFieldContainer.hpp"
#include "SpatialDimension.h"

namespace chdr::solver {
    /// Caller-owned mesh, layout, source and output storage container.
    using ChDRFieldContainer = FELFieldContainer<double, SpatialDim>;

    /// Three-component electric or magnetic field storage type on an IPPL mesh.
    using EMFieldsType = VField_t<double, SpatialDim>;

    /// Four-component scalar/vector potential or charge/current storage type.
    /// Component zero stores \f$\phi\f$ or \f$\rho\f$; components 1 through SpatialDim
    /// store the Cartesian components of \f$\mathbf A\f$ or \f$\mathbf J\f$.
    using AJFourFieldsType = SourceField_t<double, SpatialDim>;

    /** @brief Require at least two owned cells per axis on every rank.
     * @param layout Initialized, consistent 3D layout on an MPI intracommunicator.
     * @pre All ranks in layout.comm call this in matching collective order.
     * @throws std::runtime_error If any rank has fewer than two cells on an axis.
     * @details Checks extents for the default one-cell halo, not arbitrary stencils.
     * Any grid size/split axes; see HaloCells::exchangeBoundaries in src/Field/HaloCells.hpp.
     */
    void checkLocalDomain(const ippl::FieldLayout<SpatialDim>& layout);

    /** @brief Allocate E/B/J storage and zero the charge/current source.
     * @param fields Container owning the mesh, layout and resulting field storage.
     * @pre Fresh setup after checkLocalDomain; keep fields at this address afterward.
     * @details FELFieldContainer::initializeFields calls Field::initialize(mesh, layout)
     * without a halo argument, using its nghost=1 default (src/Field/Field.h).
     * Potential histories belong to the solver; evaluate_EB initializes E and B.
     */
    void initializeFieldStorage(ChDRFieldContainer& fields);
}  // namespace chdr::solver
