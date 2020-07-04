// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/irl/irl_intersection.h"

#include <omp.h>

// Timing IRL directory
#include "src/irl/stellated_dodecahedron.h"
#include "src/irl/stellated_icosahedron.h"

// IRL source directory
#include "src/generic_cutting/generic_cutting.h"
#include "src/generic_cutting/half_edge_cutting/half_edge_cutting.tpp"
#include "src/geometry/general/plane.h"
#include "src/geometry/polyhedrons/dodecahedron.h"
#include "src/geometry/polyhedrons/octahedron.h"
#include "src/geometry/polyhedrons/rectangular_cuboid.h"
#include "src/geometry/polyhedrons/symmetric_decompositions/symmetric_hexahedron.h"
#include "src/geometry/polyhedrons/symmetric_decompositions/symmetric_triangular_prism.h"
#include "src/geometry/polyhedrons/triangular_prism.h"

void irl_prismByPlanes(const double *a_prism_pts, const int a_number_of_planes,
                       const double *a_planes, double *a_volume,
                       double *a_times) {
  const auto prism = IRL::TriangularPrism::fromRawDoublePointer(6, a_prism_pts);

  IRL::PlanarLocalizer localizer;
  localizer.setNumberOfPlanes(a_number_of_planes);
  for (IRL::UnsignedIndex_t n = 0;
       n < static_cast<IRL::UnsignedIndex_t>(a_number_of_planes); ++n) {
    localizer[n] =
        IRL::Plane(IRL::Normal(a_planes[n * 4 + 0], a_planes[n * 4 + 1],
                               a_planes[n * 4 + 2]),
                   a_planes[n * 4 + 3]);
  }
  double start = omp_get_wtime();
  auto &half_edge = IRL::setHalfEdgeStructure(prism);
  auto segmented = half_edge.generateSegmentedPolyhedron();
  double end = omp_get_wtime();
  a_times[0] = end - start;

  start = omp_get_wtime();
  for (const auto &plane : localizer) {
    IRL::truncateHalfEdgePolytope(&segmented, &half_edge, plane);
  }
  end = omp_get_wtime();
  a_times[1] = end - start;

  start = omp_get_wtime();
  (*a_volume) = segmented.calculateVolume();
  end = omp_get_wtime();
  a_times[2] = end - start;

  start = omp_get_wtime();
  IRL::updatePolytopeStorage<IRL::TriangularPrism>();
  end = omp_get_wtime();
  a_times[0] += end - start;
}

void irl_unitCubeByPlanes(const double *a_cube_pts,
                          const int a_number_of_planes, const double *a_planes,
                          double *a_volume, double *a_times) {
  const auto cube = IRL::RectangularCuboid::fromBoundingPts(
      IRL::Pt(a_cube_pts[0], a_cube_pts[1], a_cube_pts[2]),
      IRL::Pt(a_cube_pts[3], a_cube_pts[4], a_cube_pts[5]));

  IRL::PlanarLocalizer localizer;
  localizer.setNumberOfPlanes(a_number_of_planes);
  for (IRL::UnsignedIndex_t n = 0;
       n < static_cast<IRL::UnsignedIndex_t>(a_number_of_planes); ++n) {
    localizer[n] =
        IRL::Plane(IRL::Normal(a_planes[n * 4 + 0], a_planes[n * 4 + 1],
                               a_planes[n * 4 + 2]),
                   a_planes[n * 4 + 3]);
  }
  double start = omp_get_wtime();
  auto &half_edge = IRL::setHalfEdgeStructure(cube);
  auto segmented = half_edge.generateSegmentedPolyhedron();
  double end = omp_get_wtime();
  a_times[0] = end - start;

  start = omp_get_wtime();
  for (const auto &plane : localizer) {
    IRL::truncateHalfEdgePolytope(&segmented, &half_edge, plane);
  }
  end = omp_get_wtime();
  a_times[1] = end - start;

  start = omp_get_wtime();
  (*a_volume) = segmented.calculateVolume();
  end = omp_get_wtime();
  a_times[2] = end - start;

  start = omp_get_wtime();
  IRL::updatePolytopeStorage<IRL::RectangularCuboid>();
  end = omp_get_wtime();
  a_times[0] += end - start;
}

void irl_triPrismByPlanes(const double *a_tri_prism_pts,
                          const int a_number_of_planes, const double *a_planes,
                          double *a_volume, double *a_times) {
  // IRL Octahedron is a Triangulated Triangular Prism
  const auto octahedron =
      IRL::Octahedron::fromRawDoublePointer(6, a_tri_prism_pts);

  IRL::PlanarLocalizer localizer;
  localizer.setNumberOfPlanes(a_number_of_planes);
  for (IRL::UnsignedIndex_t n = 0;
       n < static_cast<IRL::UnsignedIndex_t>(a_number_of_planes); ++n) {
    localizer[n] =
        IRL::Plane(IRL::Normal(a_planes[n * 4 + 0], a_planes[n * 4 + 1],
                               a_planes[n * 4 + 2]),
                   a_planes[n * 4 + 3]);
  }
  double start = omp_get_wtime();
  auto &half_edge = IRL::setHalfEdgeStructure(octahedron);
  auto segmented = half_edge.generateSegmentedPolyhedron();
  double end = omp_get_wtime();
  a_times[0] = end - start;

  start = omp_get_wtime();
  for (const auto &plane : localizer) {
    IRL::truncateHalfEdgePolytope(&segmented, &half_edge, plane);
  }
  end = omp_get_wtime();
  a_times[1] = end - start;

  start = omp_get_wtime();
  (*a_volume) = segmented.calculateVolume();
  end = omp_get_wtime();
  a_times[2] = end - start;

  start = omp_get_wtime();
  IRL::updatePolytopeStorage<IRL::Octahedron>();
  end = omp_get_wtime();
  a_times[0] += end - start;
}

void irl_triHexByPlanes(const double *a_tri_hex_pts,
                        const int a_number_of_planes, const double *a_planes,
                        double *a_volume, double *a_times) {
  // IRL Dodecahedron is a Triangulated Hexahedron
  const auto tri_hex =
      IRL::Dodecahedron::fromRawDoublePointer(8, a_tri_hex_pts);

  IRL::PlanarLocalizer localizer;
  localizer.setNumberOfPlanes(a_number_of_planes);
  for (IRL::UnsignedIndex_t n = 0;
       n < static_cast<IRL::UnsignedIndex_t>(a_number_of_planes); ++n) {
    localizer[n] =
        IRL::Plane(IRL::Normal(a_planes[n * 4 + 0], a_planes[n * 4 + 1],
                               a_planes[n * 4 + 2]),
                   a_planes[n * 4 + 3]);
  }
  double start = omp_get_wtime();
  auto &half_edge = IRL::setHalfEdgeStructure(tri_hex);
  auto segmented = half_edge.generateSegmentedPolyhedron();
  double end = omp_get_wtime();
  a_times[0] = end - start;

  start = omp_get_wtime();
  for (const auto &plane : localizer) {
    IRL::truncateHalfEdgePolytope(&segmented, &half_edge, plane);
  }
  end = omp_get_wtime();
  a_times[1] = end - start;

  start = omp_get_wtime();
  (*a_volume) = segmented.calculateVolume();
  end = omp_get_wtime();
  a_times[2] = end - start;

  start = omp_get_wtime();
  IRL::updatePolytopeStorage<IRL::Dodecahedron>();
  end = omp_get_wtime();
  a_times[0] += end - start;
}

void irl_symPrismByPlanes(const double *a_sym_prism_pts,
                          const int a_number_of_planes, const double *a_planes,
                          double *a_volume, double *a_times) {
  const auto sym_prism =
      IRL::SymmetricTriangularPrism::fromRawDoublePointer(11, a_sym_prism_pts);

  IRL::PlanarLocalizer localizer;
  localizer.setNumberOfPlanes(a_number_of_planes);
  for (IRL::UnsignedIndex_t n = 0;
       n < static_cast<IRL::UnsignedIndex_t>(a_number_of_planes); ++n) {
    localizer[n] =
        IRL::Plane(IRL::Normal(a_planes[n * 4 + 0], a_planes[n * 4 + 1],
                               a_planes[n * 4 + 2]),
                   a_planes[n * 4 + 3]);
  }
  double start = omp_get_wtime();
  auto &half_edge = IRL::setHalfEdgeStructure(sym_prism);
  auto segmented = half_edge.generateSegmentedPolyhedron();
  double end = omp_get_wtime();
  a_times[0] = end - start;

  start = omp_get_wtime();
  for (const auto &plane : localizer) {
    IRL::truncateHalfEdgePolytope(&segmented, &half_edge, plane);
  }
  end = omp_get_wtime();
  a_times[1] = end - start;

  start = omp_get_wtime();
  (*a_volume) = segmented.calculateVolume();
  end = omp_get_wtime();
  a_times[2] = end - start;

  start = omp_get_wtime();
  IRL::updatePolytopeStorage<IRL::SymmetricTriangularPrism>();
  end = omp_get_wtime();
  a_times[0] += end - start;
}

void irl_symHexByPlanes(const double *a_sym_hex_pts,
                        const int a_number_of_planes, const double *a_planes,
                        double *a_volume, double *a_times) {
  const auto sym_hex =
      IRL::SymmetricHexahedron::fromRawDoublePointer(14, a_sym_hex_pts);

  IRL::PlanarLocalizer localizer;
  localizer.setNumberOfPlanes(a_number_of_planes);
  for (IRL::UnsignedIndex_t n = 0;
       n < static_cast<IRL::UnsignedIndex_t>(a_number_of_planes); ++n) {
    localizer[n] =
        IRL::Plane(IRL::Normal(a_planes[n * 4 + 0], a_planes[n * 4 + 1],
                               a_planes[n * 4 + 2]),
                   a_planes[n * 4 + 3]);
  }
  double start = omp_get_wtime();
  auto &half_edge = IRL::setHalfEdgeStructure(sym_hex);
  auto segmented = half_edge.generateSegmentedPolyhedron();
  double end = omp_get_wtime();
  a_times[0] = end - start;

  start = omp_get_wtime();
  for (const auto &plane : localizer) {
    IRL::truncateHalfEdgePolytope(&segmented, &half_edge, plane);
  }
  end = omp_get_wtime();
  a_times[1] = end - start;

  start = omp_get_wtime();
  (*a_volume) = segmented.calculateVolume();
  end = omp_get_wtime();
  a_times[2] = end - start;

  start = omp_get_wtime();
  IRL::updatePolytopeStorage<IRL::SymmetricHexahedron>();
  end = omp_get_wtime();
  a_times[0] += end - start;
}

void irl_stelDodecahedronByPlanes(const double *a_stel_dodecahedron_pts,
                                  const int a_number_of_planes,
                                  const double *a_planes, double *a_volume,
                                  double *a_times) {
  const auto stel_dodecahedron =
      IRL::StellatedDodecahedron::fromRawDoublePointer(32,
                                                       a_stel_dodecahedron_pts);

  IRL::PlanarLocalizer localizer;
  localizer.setNumberOfPlanes(a_number_of_planes);
  for (IRL::UnsignedIndex_t n = 0;
       n < static_cast<IRL::UnsignedIndex_t>(a_number_of_planes); ++n) {
    localizer[n] =
        IRL::Plane(IRL::Normal(a_planes[n * 4 + 0], a_planes[n * 4 + 1],
                               a_planes[n * 4 + 2]),
                   a_planes[n * 4 + 3]);
  }
  double start = omp_get_wtime();
  auto &half_edge = IRL::setHalfEdgeStructure(stel_dodecahedron);
  auto segmented = half_edge.generateSegmentedPolyhedron();
  double end = omp_get_wtime();
  a_times[0] = end - start;

  start = omp_get_wtime();
  for (const auto &plane : localizer) {
    IRL::truncateHalfEdgePolytope(&segmented, &half_edge, plane);
  }
  end = omp_get_wtime();
  a_times[1] = end - start;

  start = omp_get_wtime();
  (*a_volume) = segmented.calculateVolume();
  end = omp_get_wtime();
  a_times[2] = end - start;

  start = omp_get_wtime();
  IRL::updatePolytopeStorage<IRL::StellatedDodecahedron>();
  end = omp_get_wtime();
  a_times[0] += end - start;
}

void irl_stelIcosahedronByPlanes(const double *a_stel_icosahedron_pts,
                                 const int a_number_of_planes,
                                 const double *a_planes, double *a_volume,
                                 double *a_times) {
  const auto stel_icosahedron = IRL::StellatedIcosahedron::fromRawDoublePointer(
      32, a_stel_icosahedron_pts);

  IRL::PlanarLocalizer localizer;
  localizer.setNumberOfPlanes(a_number_of_planes);
  for (IRL::UnsignedIndex_t n = 0;
       n < static_cast<IRL::UnsignedIndex_t>(a_number_of_planes); ++n) {
    localizer[n] =
        IRL::Plane(IRL::Normal(a_planes[n * 4 + 0], a_planes[n * 4 + 1],
                               a_planes[n * 4 + 2]),
                   a_planes[n * 4 + 3]);
  }
  double start = omp_get_wtime();
  auto &half_edge = IRL::setHalfEdgeStructure(stel_icosahedron);
  auto segmented = half_edge.generateSegmentedPolyhedron();
  double end = omp_get_wtime();
  a_times[0] = end - start;

  start = omp_get_wtime();
  for (const auto &plane : localizer) {
    IRL::truncateHalfEdgePolytope(&segmented, &half_edge, plane);
  }
  end = omp_get_wtime();
  a_times[1] = end - start;

  start = omp_get_wtime();
  (*a_volume) = segmented.calculateVolume();
  end = omp_get_wtime();
  a_times[2] = end - start;

  start = omp_get_wtime();
  IRL::updatePolytopeStorage<IRL::StellatedIcosahedron>();
  end = omp_get_wtime();
  a_times[0] += end - start;
}
