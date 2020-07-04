// Copyright (C) 2020 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// NOTE: This file was automatically generated using the IRL
// tool `polyhedron_creator`

#ifndef SRC_IRL_GVM_STELLATED_ICOSAHEDRON_H_
#define SRC_IRL_GVM_STELLATED_ICOSAHEDRON_H_

#include "src/geometry/general/geometry_type_traits.h"
#include "src/geometry/general/stored_vertex_access.h"
#include "src/geometry/half_edge_structures/half_edge_polyhedron.h"
#include "src/geometry/polyhedrons/general_polyhedron.h"
#include "src/geometry/polyhedrons/tet.h"
#include "src/parameters/defined_types.h"

namespace IRL {

template <class Derived, class VertexType>
class StellatedIcosahedronSpecialization
    : public BasePolyhedron<Derived, VertexType, ProxyTet<Derived>> {
public:
  HalfEdgePolyhedron<VertexType> generateHalfEdgeVersion(void) const;

  template <class HalfEdgePolyhedronType>
  void setHalfEdgeVersion(HalfEdgePolyhedronType *a_half_edge_version) const;

  static constexpr UnsignedIndex_t getNumberOfSimplicesInDecomposition(void);

  static constexpr std::array<UnsignedIndex_t, 4>
  getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet);

  ProxyTet<Derived>
  getSimplexFromDecomposition(const UnsignedIndex_t a_tet) const;
};

template <class VertexType>
class StoredStellatedIcosahedron
    : public StoredVertexAccess<StoredStellatedIcosahedron<VertexType>,
                                VertexType, 32>,
      public StellatedIcosahedronSpecialization<
          StoredStellatedIcosahedron<VertexType>, VertexType> {
  friend StoredVertexAccess<StoredStellatedIcosahedron<VertexType>, VertexType,
                            32>;

public:
  using StoredVertexAccess<StoredStellatedIcosahedron<VertexType>, VertexType,
                           32>::StoredVertexAccess;

  StoredStellatedIcosahedron(void) = default;
};

// Predefined types
using StellatedIcosahedron = StoredStellatedIcosahedron<Pt>;

template <class VertexType>
struct is_polyhedron<StoredStellatedIcosahedron<VertexType>> : std::true_type {
};

} // namespace IRL

#include "stellated_icosahedron.tpp"
#endif // SRC_IRL_GVM_STELLATED_ICOSAHEDRON_H_
