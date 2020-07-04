// Copyright (C) 2020 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// NOTE: This file was automatically generated using the IRL
// tool `polyhedron_creator`

#ifndef SRC_IRL_GVM_STELLATED_ICOSAHEDRON_TPP_
#define SRC_IRL_GVM_STELLATED_ICOSAHEDRON_TPP_

#include <cassert>

#include "src/geometry/general/moment_calculation_through_simplices.h"
#include "src/geometry/half_edge_structures/half_edge.h"

namespace IRL {

namespace stellated_icosahedron_triangulation {
static constexpr UnsignedIndex_t datum_index = 0;
static constexpr std::array<std::array<UnsignedIndex_t, 3>, 50>
    face_triangle_decomposition{
        {{8, 4, 12},  {5, 10, 13}, {2, 4, 14},  {4, 9, 14},  {9, 2, 14},
         {2, 11, 15}, {11, 5, 15}, {5, 2, 15},  {1, 6, 16},  {6, 8, 16},
         {8, 1, 16},  {1, 10, 17}, {10, 7, 17}, {7, 1, 17},  {3, 9, 18},
         {9, 6, 18},  {6, 3, 18},  {3, 7, 19},  {7, 11, 19}, {11, 3, 19},
         {10, 8, 20}, {1, 8, 21},  {8, 10, 21}, {10, 1, 21}, {2, 9, 22},
         {9, 11, 22}, {11, 2, 22}, {3, 11, 23}, {11, 9, 23}, {9, 3, 23},
         {4, 2, 24},  {2, 5, 25},  {6, 1, 26},  {1, 3, 26},  {3, 6, 26},
         {7, 3, 27},  {3, 1, 27},  {1, 7, 27},  {8, 6, 28},  {6, 4, 28},
         {4, 8, 28},  {9, 4, 29},  {4, 6, 29},  {6, 9, 29},  {10, 5, 30},
         {5, 7, 30},  {7, 10, 30}, {11, 7, 31}, {7, 5, 31},  {5, 11, 31}}};
} // namespace stellated_icosahedron_triangulation

template <class Derived, class VertexType>
HalfEdgePolyhedron<VertexType> StellatedIcosahedronSpecialization<
    Derived, VertexType>::generateHalfEdgeVersion(void) const {
  HalfEdgePolyhedron<VertexType> half_edge_version;
  this->setHalfEdgeVersion(&half_edge_version);
  return half_edge_version;
}

template <class Derived, class VertexType>
template <class HalfEdgePolyhedronType>
void StellatedIcosahedronSpecialization<Derived, VertexType>::
    setHalfEdgeVersion(HalfEdgePolyhedronType *a_half_edge_version) const {
  using HalfEdgeType = typename HalfEdgePolyhedronType::half_edge_type;

  a_half_edge_version->resize(180, 32, 60);

  for (UnsignedIndex_t v = 0; v < 32; ++v) {
    a_half_edge_version->getVertex(v).setLocation((*this)[v]);
  }
  static constexpr std::array<UnsignedIndex_t, 180> ending_vertex_mapping{
      {8,  12, 0,  4,  12, 8,  0,  12, 4,  5,  13, 0,  10, 13, 5,  0,  13, 10,
       4,  14, 2,  9,  14, 4,  2,  14, 9,  11, 15, 2,  5,  15, 11, 2,  15, 5,
       6,  16, 1,  8,  16, 6,  1,  16, 8,  10, 17, 1,  7,  17, 10, 1,  17, 7,
       9,  18, 3,  6,  18, 9,  3,  18, 6,  7,  19, 3,  11, 19, 7,  3,  19, 11,
       10, 20, 0,  8,  20, 10, 0,  20, 8,  8,  21, 1,  10, 21, 8,  1,  21, 10,
       9,  22, 2,  11, 22, 9,  2,  22, 11, 11, 23, 3,  9,  23, 11, 3,  23, 9,
       2,  24, 4,  0,  24, 2,  4,  24, 0,  0,  25, 5,  2,  25, 0,  5,  25, 2,
       1,  26, 6,  3,  26, 1,  6,  26, 3,  3,  27, 7,  1,  27, 3,  7,  27, 1,
       6,  28, 8,  4,  28, 6,  8,  28, 4,  4,  29, 9,  6,  29, 4,  9,  29, 6,
       5,  30, 10, 7,  30, 5,  10, 30, 7,  7,  31, 11, 5,  31, 7,  11, 31, 5}};
  static constexpr std::array<UnsignedIndex_t, 180> previous_half_edge_mapping{
      {2,   0,   1,   5,   3,   4,   8,   6,   7,   11,  9,   10,  14,  12,
       13,  17,  15,  16,  20,  18,  19,  23,  21,  22,  26,  24,  25,  29,
       27,  28,  32,  30,  31,  35,  33,  34,  38,  36,  37,  41,  39,  40,
       44,  42,  43,  47,  45,  46,  50,  48,  49,  53,  51,  52,  56,  54,
       55,  59,  57,  58,  62,  60,  61,  65,  63,  64,  68,  66,  67,  71,
       69,  70,  74,  72,  73,  77,  75,  76,  80,  78,  79,  83,  81,  82,
       86,  84,  85,  89,  87,  88,  92,  90,  91,  95,  93,  94,  98,  96,
       97,  101, 99,  100, 104, 102, 103, 107, 105, 106, 110, 108, 109, 113,
       111, 112, 116, 114, 115, 119, 117, 118, 122, 120, 121, 125, 123, 124,
       128, 126, 127, 131, 129, 130, 134, 132, 133, 137, 135, 136, 140, 138,
       139, 143, 141, 142, 146, 144, 145, 149, 147, 148, 152, 150, 151, 155,
       153, 154, 158, 156, 157, 161, 159, 160, 164, 162, 163, 167, 165, 166,
       170, 168, 169, 173, 171, 172, 176, 174, 175, 179, 177, 178}};
  static constexpr std::array<UnsignedIndex_t, 180> next_half_edge_mapping{
      {1,   2,   0,   4,   5,   3,   7,   8,   6,   10,  11,  9,   13,  14,
       12,  16,  17,  15,  19,  20,  18,  22,  23,  21,  25,  26,  24,  28,
       29,  27,  31,  32,  30,  34,  35,  33,  37,  38,  36,  40,  41,  39,
       43,  44,  42,  46,  47,  45,  49,  50,  48,  52,  53,  51,  55,  56,
       54,  58,  59,  57,  61,  62,  60,  64,  65,  63,  67,  68,  66,  70,
       71,  69,  73,  74,  72,  76,  77,  75,  79,  80,  78,  82,  83,  81,
       85,  86,  84,  88,  89,  87,  91,  92,  90,  94,  95,  93,  97,  98,
       96,  100, 101, 99,  103, 104, 102, 106, 107, 105, 109, 110, 108, 112,
       113, 111, 115, 116, 114, 118, 119, 117, 121, 122, 120, 124, 125, 123,
       127, 128, 126, 130, 131, 129, 133, 134, 132, 136, 137, 135, 139, 140,
       138, 142, 143, 141, 145, 146, 144, 148, 149, 147, 151, 152, 150, 154,
       155, 153, 157, 158, 156, 160, 161, 159, 163, 164, 162, 166, 167, 165,
       169, 170, 168, 172, 173, 171, 175, 176, 174, 178, 179, 177}};
  static constexpr std::array<UnsignedIndex_t, 180> face_mapping{
      {0,  0,  0,  1,  1,  1,  2,  2,  2,  3,  3,  3,  4,  4,  4,  5,  5,  5,
       6,  6,  6,  7,  7,  7,  8,  8,  8,  9,  9,  9,  10, 10, 10, 11, 11, 11,
       12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, 16, 16, 16, 17, 17, 17,
       18, 18, 18, 19, 19, 19, 20, 20, 20, 21, 21, 21, 22, 22, 22, 23, 23, 23,
       24, 24, 24, 25, 25, 25, 26, 26, 26, 27, 27, 27, 28, 28, 28, 29, 29, 29,
       30, 30, 30, 31, 31, 31, 32, 32, 32, 33, 33, 33, 34, 34, 34, 35, 35, 35,
       36, 36, 36, 37, 37, 37, 38, 38, 38, 39, 39, 39, 40, 40, 40, 41, 41, 41,
       42, 42, 42, 43, 43, 43, 44, 44, 44, 45, 45, 45, 46, 46, 46, 47, 47, 47,
       48, 48, 48, 49, 49, 49, 50, 50, 50, 51, 51, 51, 52, 52, 52, 53, 53, 53,
       54, 54, 54, 55, 55, 55, 56, 56, 56, 57, 57, 57, 58, 58, 58, 59, 59, 59}};
  static constexpr std::array<UnsignedIndex_t, 180> opposite_half_edge_mapping{
      {78,  5,   7,   150, 8,   1,   114, 2,   4,   117, 14,  16,  162, 17,
       10,  72,  11,  13,  108, 23,  25,  153, 26,  19,  90,  20,  22,  96,
       32,  34,  177, 35,  28,  123, 29,  31,  126, 41,  43,  144, 44,  37,
       81,  38,  40,  87,  50,  52,  168, 53,  46,  141, 47,  49,  105, 59,
       61,  159, 62,  55,  132, 56,  58,  135, 68,  70,  171, 71,  64,  99,
       65,  67,  15,  77,  79,  84,  80,  73,  0,   74,  76,  42,  86,  88,
       75,  89,  82,  45,  83,  85,  24,  95,  97,  102, 98,  91,  27,  92,
       94,  69,  104, 106, 93,  107, 100, 54,  101, 103, 18,  113, 115, 120,
       116, 109, 6,   110, 112, 9,   122, 124, 111, 125, 118, 33,  119, 121,
       36,  131, 133, 138, 134, 127, 60,  128, 130, 63,  140, 142, 129, 143,
       136, 51,  137, 139, 39,  149, 151, 156, 152, 145, 3,   146, 148, 21,
       158, 160, 147, 161, 154, 57,  155, 157, 12,  167, 169, 174, 170, 163,
       48,  164, 166, 66,  176, 178, 165, 179, 172, 30,  173, 175}};
  for (UnsignedIndex_t n = 0;
       n < static_cast<UnsignedIndex_t>(ending_vertex_mapping.size()); ++n) {
    HalfEdgeType &current_half_edge = a_half_edge_version->getHalfEdge(n);
    current_half_edge = HalfEdgeType(
        &a_half_edge_version->getVertex(ending_vertex_mapping[n]),
        &a_half_edge_version->getHalfEdge(previous_half_edge_mapping[n]),
        &a_half_edge_version->getHalfEdge(next_half_edge_mapping[n]),
        &a_half_edge_version->getFace(face_mapping[n]));
    current_half_edge.setOppositeHalfEdge(
        &a_half_edge_version->getHalfEdge(opposite_half_edge_mapping[n]));
    current_half_edge.getFace()->setStartingHalfEdge(&current_half_edge);
    current_half_edge.getVertex()->setHalfEdge(&current_half_edge);
  }
}

template <class Derived, class VertexType>
constexpr UnsignedIndex_t StellatedIcosahedronSpecialization<
    Derived, VertexType>::getNumberOfSimplicesInDecomposition(void) {
  return static_cast<UnsignedIndex_t>(
      stellated_icosahedron_triangulation::face_triangle_decomposition.size());
}

template <class Derived, class VertexType>
constexpr std::array<UnsignedIndex_t, 4>
StellatedIcosahedronSpecialization<Derived, VertexType>::
    getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet) {
  assert(a_tet < StellatedIcosahedronSpecialization::
                     getNumberOfSimplicesInDecomposition());
  return {
      stellated_icosahedron_triangulation::face_triangle_decomposition[a_tet]
                                                                      [0],
      stellated_icosahedron_triangulation::face_triangle_decomposition[a_tet]
                                                                      [1],
      stellated_icosahedron_triangulation::face_triangle_decomposition[a_tet]
                                                                      [2],
      stellated_icosahedron_triangulation::datum_index};
}

template <class Derived, class VertexType>
ProxyTet<Derived> StellatedIcosahedronSpecialization<Derived, VertexType>::
    getSimplexFromDecomposition(const UnsignedIndex_t a_tet) const {
  assert(a_tet < this->getNumberOfSimplicesInDecomposition());
  return {static_cast<const Derived &>(*this),
          this->getSimplexIndicesFromDecomposition(a_tet)};
}

} // namespace IRL
#endif // SRC_IRL_GVM_STELLATED_ICOSAHEDRON_TPP_
