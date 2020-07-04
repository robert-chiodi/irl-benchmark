// Copyright (C) 2020 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// NOTE: This file was automatically generated using the IRL
// tool `polyhedron_creator`

#ifndef SRC_IRL_GVM_STELLATED_DODECAHEDRON_TPP_
#define SRC_IRL_GVM_STELLATED_DODECAHEDRON_TPP_

#include <cassert>

#include "src/geometry/general/moment_calculation_through_simplices.h"
#include "src/geometry/half_edge_structures/half_edge.h"

namespace IRL {

namespace stellated_dodecahedron_triangulation {
static constexpr UnsignedIndex_t datum_index = 0;
static constexpr std::array<std::array<UnsignedIndex_t, 3>, 54>
    face_triangle_decomposition{
        {{8, 9, 20},   {9, 4, 20},   {4, 16, 20},  {16, 17, 21}, {17, 2, 21},
         {2, 12, 21},  {12, 2, 22},  {2, 10, 22},  {10, 3, 22},  {3, 13, 22},
         {13, 12, 22}, {9, 5, 23},   {5, 15, 23},  {15, 14, 23}, {14, 4, 23},
         {4, 9, 23},   {3, 19, 24},  {19, 18, 24}, {18, 1, 24},  {1, 13, 24},
         {13, 3, 24},  {7, 11, 25},  {11, 6, 25},  {6, 14, 25},  {14, 15, 25},
         {15, 7, 25},  {12, 13, 26}, {13, 1, 26},  {1, 8, 26},   {8, 1, 27},
         {1, 18, 27},  {18, 5, 27},  {5, 9, 27},   {9, 8, 27},   {16, 4, 28},
         {4, 14, 28},  {14, 6, 28},  {6, 17, 28},  {17, 16, 28}, {6, 11, 29},
         {11, 10, 29}, {10, 2, 29},  {2, 17, 29},  {17, 6, 29},  {7, 15, 30},
         {15, 5, 30},  {5, 18, 30},  {18, 19, 30}, {19, 7, 30},  {7, 19, 31},
         {19, 3, 31},  {3, 10, 31},  {10, 11, 31}, {11, 7, 31}}};
} // namespace stellated_dodecahedron_triangulation

template <class Derived, class VertexType>
HalfEdgePolyhedron<VertexType> StellatedDodecahedronSpecialization<
    Derived, VertexType>::generateHalfEdgeVersion(void) const {
  HalfEdgePolyhedron<VertexType> half_edge_version;
  this->setHalfEdgeVersion(&half_edge_version);
  return half_edge_version;
}

template <class Derived, class VertexType>
template <class HalfEdgePolyhedronType>
void StellatedDodecahedronSpecialization<Derived, VertexType>::
    setHalfEdgeVersion(HalfEdgePolyhedronType *a_half_edge_version) const {
  using HalfEdgeType = typename HalfEdgePolyhedronType::half_edge_type;

  a_half_edge_version->resize(180, 32, 60);

  for (UnsignedIndex_t v = 0; v < 32; ++v) {
    a_half_edge_version->getVertex(v).setLocation((*this)[v]);
  }

  static constexpr std::array<UnsignedIndex_t, 180> ending_vertex_mapping{
      {8,  20, 0,  9,  20, 8,  4,  20, 9,  16, 20, 4,  0,  20, 16, 16, 21, 0,
       17, 21, 16, 2,  21, 17, 12, 21, 2,  0,  21, 12, 2,  22, 12, 10, 22, 2,
       3,  22, 10, 13, 22, 3,  12, 22, 13, 5,  23, 9,  15, 23, 5,  14, 23, 15,
       4,  23, 14, 9,  23, 4,  19, 24, 3,  18, 24, 19, 1,  24, 18, 13, 24, 1,
       3,  24, 13, 11, 25, 7,  6,  25, 11, 14, 25, 6,  15, 25, 14, 7,  25, 15,
       12, 26, 0,  13, 26, 12, 1,  26, 13, 8,  26, 1,  0,  26, 8,  1,  27, 8,
       18, 27, 1,  5,  27, 18, 9,  27, 5,  8,  27, 9,  4,  28, 16, 14, 28, 4,
       6,  28, 14, 17, 28, 6,  16, 28, 17, 11, 29, 6,  10, 29, 11, 2,  29, 10,
       17, 29, 2,  6,  29, 17, 15, 30, 7,  5,  30, 15, 18, 30, 5,  19, 30, 18,
       7,  30, 19, 19, 31, 7,  3,  31, 19, 10, 31, 3,  11, 31, 10, 7,  31, 11}};
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
      {102, 5,   13,  117, 8,   1,   57,  11,  4,   120, 14,  7,   15,  2,
       10,  12,  20,  28,  132, 23,  16,  144, 26,  19,  30,  29,  22,  90,
       17,  25,  24,  35,  43,  141, 38,  31,  171, 41,  34,  72,  44,  37,
       93,  32,  40,  114, 50,  58,  153, 53,  46,  84,  56,  49,  123, 59,
       52,  6,   47,  55,  168, 65,  73,  159, 68,  61,  108, 71,  64,  96,
       74,  67,  39,  62,  70,  177, 80,  88,  135, 83,  76,  126, 86,  79,
       51,  89,  82,  150, 77,  85,  27,  95,  103, 42,  98,  91,  69,  101,
       94,  105, 104, 97,  0,   92,  100, 99,  110, 118, 66,  113, 106, 156,
       116, 109, 45,  119, 112, 3,   107, 115, 9,   125, 133, 54,  128, 121,
       81,  131, 124, 147, 134, 127, 18,  122, 130, 78,  140, 148, 174, 143,
       136, 33,  146, 139, 21,  149, 142, 129, 137, 145, 87,  155, 163, 48,
       158, 151, 111, 161, 154, 63,  164, 157, 165, 152, 160, 162, 170, 178,
       60,  173, 166, 36,  176, 169, 138, 179, 172, 75,  167, 175}};
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
constexpr UnsignedIndex_t StellatedDodecahedronSpecialization<
    Derived, VertexType>::getNumberOfSimplicesInDecomposition(void) {
  return static_cast<UnsignedIndex_t>(
      stellated_dodecahedron_triangulation::face_triangle_decomposition.size());
}

template <class Derived, class VertexType>
constexpr std::array<UnsignedIndex_t, 4>
StellatedDodecahedronSpecialization<Derived, VertexType>::
    getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet) {
  assert(a_tet < StellatedDodecahedronSpecialization::
                     getNumberOfSimplicesInDecomposition());
  return {
      stellated_dodecahedron_triangulation::face_triangle_decomposition[a_tet]
                                                                       [0],
      stellated_dodecahedron_triangulation::face_triangle_decomposition[a_tet]
                                                                       [1],
      stellated_dodecahedron_triangulation::face_triangle_decomposition[a_tet]
                                                                       [2],
      stellated_dodecahedron_triangulation::datum_index};
}

template <class Derived, class VertexType>
ProxyTet<Derived> StellatedDodecahedronSpecialization<Derived, VertexType>::
    getSimplexFromDecomposition(const UnsignedIndex_t a_tet) const {
  assert(a_tet < this->getNumberOfSimplicesInDecomposition());
  return {static_cast<const Derived &>(*this),
          this->getSimplexIndicesFromDecomposition(a_tet)};
}

} // namespace IRL
#endif // SRC_IRL_GVM_STELLATED_DODECAHEDRON_TPP_
