// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/irl_gvm/irl_gvm_intersection.h"

#include <omp.h>

#include <fstream>

// Timing IRL GVM directory
#include "src/irl_gvm/stellated_dodecahedron.h"
#include "src/irl_gvm/stellated_icosahedron.h"

// IRL source directory
#include "src/generic_cutting/generic_cutting.h"
#include "src/geometry/general/plane.h"
#include "src/geometry/polyhedrons/general_polyhedron.h"
#include "src/geometry/polyhedrons/polyhedron_connectivity.h"
#include "src/geometry/polyhedrons/rectangular_cuboid.h"
#include "src/geometry/polyhedrons/triangular_prism.h"
#include "src/helpers/mymath.h"
#include "src/parameters/defined_types.h"

// If true, will print out each object  to .vtu files
static bool constexpr PRINT_OBJECTS = false;

// This function sets up the cubic mesh. In IRL,
// this means generating PlanarLocalizer objects
// to represent the cell faces (and serve as edges
// in the graph) and connecting them together to form
// LocalizerLink objects. Together, the LocalizerLink objects
// for the graph of the mesh over which we will distribute.
static void setupCubicMesh(const std::array<IRL::Pt, 2> a_bounding_pts,
                           const std::array<int, 3> a_ncells,
                           std::vector<IRL::PlanarLocalizer> *a_localizers,
                           std::vector<IRL::LocalizerLink> *a_links) {
  std::array<double, 3> dx{{(a_bounding_pts[1][0] - a_bounding_pts[0][0]) /
                                static_cast<double>(a_ncells[0]),
                            (a_bounding_pts[1][1] - a_bounding_pts[0][1]) /
                                static_cast<double>(a_ncells[1]),
                            (a_bounding_pts[1][2] - a_bounding_pts[0][2]) /
                                static_cast<double>(a_ncells[2])}};
  a_localizers->resize(a_ncells[0] * a_ncells[1] * a_ncells[2]);
  a_links->resize(a_ncells[0] * a_ncells[1] * a_ncells[2]);
  std::array<IRL::Pt, 2> local_bb;
  for (int k = 0; k < a_ncells[2]; ++k) {
    for (int j = 0; j < a_ncells[1]; ++j) {
      for (int i = 0; i < a_ncells[0]; ++i) {
        const int index = i + j * a_ncells[0] + k * a_ncells[1] * a_ncells[0];
        local_bb[0] =
            a_bounding_pts[0] + IRL::Pt(static_cast<double>(i) * dx[0],
                                        static_cast<double>(j) * dx[1],
                                        static_cast<double>(k) * dx[2]);
        local_bb[1] = local_bb[0] + IRL::Pt(dx[0], dx[1], dx[2]);
        // Set up reconstructions
        (*a_localizers)[index].setNumberOfPlanes(6);
        (*a_localizers)[index][0] =
            IRL::Plane(IRL::Normal(1.0, 0.0, 0.0), local_bb[1][0]);
        (*a_localizers)[index][1] =
            IRL::Plane(IRL::Normal(-1.0, 0.0, 0.0), -local_bb[0][0]);
        (*a_localizers)[index][2] =
            IRL::Plane(IRL::Normal(0.0, 1.0, 0.0), local_bb[1][1]);
        (*a_localizers)[index][3] =
            IRL::Plane(IRL::Normal(0.0, -1.0, 0.0), -local_bb[0][1]);
        (*a_localizers)[index][4] =
            IRL::Plane(IRL::Normal(0.0, 0.0, 1.0), local_bb[1][2]);
        (*a_localizers)[index][5] =
            IRL::Plane(IRL::Normal(0.0, 0.0, -1.0), -local_bb[0][2]);

        // Now set link connections
        (*a_links)[index] = IRL::LocalizerLink(&(*a_localizers)[index]);
        (*a_links)[index].setId(index);
        (*a_links)[index].setEdgeConnectivity(
            0, i != a_ncells[0] - 1 ? &(*a_links)[index + 1] : nullptr);
        (*a_links)[index].setEdgeConnectivity(1, i != 0 ? &(*a_links)[index - 1]
                                                        : nullptr);
        (*a_links)[index].setEdgeConnectivity(
            2,
            j != a_ncells[1] - 1 ? &(*a_links)[index + a_ncells[0]] : nullptr);
        (*a_links)[index].setEdgeConnectivity(
            3, j != 0 ? &(*a_links)[index - a_ncells[0]] : nullptr);
        (*a_links)[index].setEdgeConnectivity(
            4, k != a_ncells[2] - 1
                   ? &(*a_links)[index + a_ncells[0] * a_ncells[1]]
                   : nullptr);
        (*a_links)[index].setEdgeConnectivity(
            5,
            k != 0 ? &(*a_links)[index - a_ncells[0] * a_ncells[1]] : nullptr);
      }
    }
  }

  if (PRINT_OBJECTS) {
    for (int k = 0; k < a_ncells[2]; ++k) {
      for (int j = 0; j < a_ncells[1]; ++j) {
        for (int i = 0; i < a_ncells[0]; ++i) {
          local_bb[0] =
              a_bounding_pts[0] + IRL::Pt(static_cast<double>(i) * dx[0],
                                          static_cast<double>(j) * dx[1],
                                          static_cast<double>(k) * dx[2]);
          local_bb[1] = local_bb[0] + IRL::Pt(dx[0], dx[1], dx[2]);
          const int index = i + j * a_ncells[0] + k * a_ncells[1] * a_ncells[0];
          auto poly =
              IRL::RectangularCuboid::fromBoundingPts(local_bb[0], local_bb[1]);
          auto &half_edge = IRL::setHalfEdgeStructure(poly);
          auto segmented = half_edge.generateSegmentedPolyhedron();
          std::string filename = "cube_mesh_" + std::to_string(index) + ".vtu";
          std::ofstream myfile;
          myfile.open(filename);
          myfile << segmented;
          myfile.close();
        }
      }
    }
  }
}

void irl_gvm_cubeOntoCubicMesh(const double *a_cube_pts, double *a_volume,
                               double *a_times, std::size_t *a_entered_cells) {
  std::array<IRL::Pt, 2> mesh_box{
      {IRL::Pt(-3.0, -3.0, -3.0), IRL::Pt(3.0, 3.0, 3.0)}};
  std::array<int, 3> ncells{{3, 3, 3}};
  static std::vector<IRL::PlanarLocalizer> localizers;
  static std::vector<IRL::LocalizerLink> links;
  static bool mesh_made = false;
  if (!mesh_made) {
    mesh_made = true;
    setupCubicMesh(mesh_box, ncells, &localizers, &links);
  }

  // This will return the volume from each cell, tagged with that cells unique
  // Id set during setupCubicMesh
  const auto cube = IRL::RectangularCuboid::fromBoundingPts(
      IRL::Pt(a_cube_pts[0], a_cube_pts[1], a_cube_pts[2]),
      IRL::Pt(a_cube_pts[3], a_cube_pts[4], a_cube_pts[5]));
  double start = omp_get_wtime();
  const int mid_localizer = ncells[0] / 2 + ncells[1] / 2 * ncells[0] +
                            ncells[2] / 2 * ncells[0] * ncells[1];
  auto tagged_volumes =
      IRL::getVolumeMoments<IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>(
          cube, links[mid_localizer]);
  double end = omp_get_wtime();
  a_times[0] = end - start;

  // Return volume that is the sum of the distributed volume
  *a_volume = 0.0;
  for (const auto &volume : tagged_volumes) {
    *a_volume += volume.volume_moments_m;
  }
  *a_entered_cells = static_cast<std::size_t>(tagged_volumes.size());
}

void irl_gvm_stelIcosahedronOntoCubicMesh(const double *a_stel_icosahedron_pts,
                                          double *a_volume, double *a_times,
                                          std::size_t *a_entered_cells) {
  std::array<IRL::Pt, 2> mesh_box{
      {IRL::Pt(-3.0, -3.0, -3.0), IRL::Pt(3.0, 3.0, 3.0)}};
  std::array<int, 3> ncells{{3, 3, 3}};
  static std::vector<IRL::PlanarLocalizer> localizers;
  static std::vector<IRL::LocalizerLink> links;
  static bool mesh_made = false;
  if (!mesh_made) {
    mesh_made = true;
    setupCubicMesh(mesh_box, ncells, &localizers, &links);
  }

  // This will return the volume from each cell, tagged with that cells unique
  // Id set during setupCubicMesh
  const auto stel_icosahedron = IRL::StellatedIcosahedron::fromRawDoublePointer(
      32, a_stel_icosahedron_pts);

  double start = omp_get_wtime();
  const int mid_localizer = ncells[0] / 2 + ncells[1] / 2 * ncells[0] +
                            ncells[2] / 2 * ncells[0] * ncells[1];
  auto tagged_volumes =
      IRL::getVolumeMoments<IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>(
          stel_icosahedron, links[mid_localizer]);
  double end = omp_get_wtime();
  a_times[0] = end - start;

  // Return volume that is the sum of the distributed volume
  *a_volume = 0.0;
  for (const auto &volume : tagged_volumes) {
    *a_volume += volume.volume_moments_m;
  }
  *a_entered_cells = static_cast<std::size_t>(tagged_volumes.size());
}

// Tetrahedral mesh-type domain with a cuboid decomposed into
// 24 tets using face centroids and the cell center.
// This function sets up the Tet  mesh. In IRL,
// this means generating PlanarLocalizer objects
// to represent the cell faces (and serve as edges
// in the graph) and connecting them together to form
// LocalizerLink objects. Together, the LocalizerLink objects
// for the graph of the mesh over which we will distribute.
static void setupTetMesh(const std::array<IRL::Pt, 2> a_bounding_pts,
                         std::vector<IRL::PlanarLocalizer> *a_localizers,
                         std::vector<IRL::LocalizerLink> *a_links) {
  const std::array<IRL::Pt, 8> cuboid_pts{
      {IRL::Pt(a_bounding_pts[1][0], a_bounding_pts[0][1],
               a_bounding_pts[0][2]),
       IRL::Pt(a_bounding_pts[1][0], a_bounding_pts[1][1],
               a_bounding_pts[0][2]),
       IRL::Pt(a_bounding_pts[1][0], a_bounding_pts[1][1],
               a_bounding_pts[1][2]),
       IRL::Pt(a_bounding_pts[1][0], a_bounding_pts[0][1],
               a_bounding_pts[1][2]),
       IRL::Pt(a_bounding_pts[0][0], a_bounding_pts[0][1],
               a_bounding_pts[0][2]),
       IRL::Pt(a_bounding_pts[0][0], a_bounding_pts[1][1],
               a_bounding_pts[0][2]),
       IRL::Pt(a_bounding_pts[0][0], a_bounding_pts[1][1],
               a_bounding_pts[1][2]),
       IRL::Pt(a_bounding_pts[0][0], a_bounding_pts[0][1],
               a_bounding_pts[1][2])}};

  const std::array<std::array<int, 4>, 6> face_indices{{{0, 1, 2, 3},
                                                        {0, 4, 5, 1},
                                                        {1, 5, 6, 2},
                                                        {2, 6, 7, 3},
                                                        {3, 7, 4, 0},
                                                        {4, 7, 6, 5}}};

  const std::array<IRL::Normal, 6> face_normals{
      {IRL::Normal(1.0, 0.0, 0.0), IRL::Normal(0.0, 0.0, -1.0),
       IRL::Normal(0.0, 1.0, 0.0), IRL::Normal(0.0, 0.0, 1.0),
       IRL::Normal(0.0, -1.0, 0.0), IRL::Normal(-1.0, 0.0, 0.0)}};
  std::array<IRL::Pt, 6> face_centroids;
  for (std::size_t f = 0; f < face_centroids.size(); ++f) {
    face_centroids[f] = IRL::Pt::fromScalarConstant(0.0);
    for (const auto &index : face_indices[f]) {
      face_centroids[f] += cuboid_pts[index];
    }
    face_centroids[f] /= static_cast<double>(face_indices[f].size());
  }

  auto volume_centroid = IRL::Pt::fromScalarConstant(0.0);
  for (const auto vertex : cuboid_pts) {
    volume_centroid += vertex;
  }
  volume_centroid /= static_cast<double>(cuboid_pts.size());

  std::array<int, 24> shared_edge_link{{7,  11, 15, 19, 18, 23, 8,  0,
                                        6,  22, 12, 1,  10, 21, 16, 2,
                                        14, 20, 4,  3,  17, 13, 9,  5}};

  a_localizers->resize(24);
  a_links->resize(24);
  for (std::size_t n = 0; n < a_localizers->size(); ++n) {
    (*a_links)[n] = IRL::LocalizerLink(&(*a_localizers)[n]);
  }
  std::array<IRL::Pt, 4> tet;
  tet[3] = volume_centroid;
  for (std::size_t f = 0; f < face_indices.size(); ++f) {
    tet[2] = face_centroids[f];
    const auto &face = face_indices[f];
    for (std::size_t e = 0; e < face.size(); ++e) {
      tet[0] = cuboid_pts[face[e]];
      tet[1] = cuboid_pts[face[(e + 1) % face.size()]];
      const std::size_t index = e + f * 4;
      (*a_localizers)[index].setNumberOfPlanes(4);

      // Make first plane from on cuboid face
      IRL::Normal normal = face_normals[f];
      (*a_localizers)[index][0] =
          IRL::Plane(normal, normal * face_centroids[f]);

      // Make second plane for face having edge on original cuboid
      normal = IRL::crossProductNormalized(IRL::Pt(tet[1] - tet[3]),
                                           IRL::Pt(tet[0] - tet[3]));
      (*a_localizers)[index][1] = IRL::Plane(normal, normal * tet[3]);

      // Make third plane for tet[0] and face_centroid
      normal = IRL::crossProductNormalized(IRL::Pt(tet[0] - tet[3]),
                                           IRL::Pt(tet[2] - tet[3]));
      (*a_localizers)[index][2] = IRL::Plane(normal, normal * tet[3]);

      // Make fourth plane for tet[1] and face_centroid
      normal = IRL::crossProductNormalized(IRL::Pt(tet[2] - tet[3]),
                                           IRL::Pt(tet[1] - tet[3]));
      (*a_localizers)[index][3] = IRL::Plane(normal, normal * tet[3]);

      // Now link up localizer to rest of the mesh
      (*a_links)[index].setId(static_cast<IRL::UnsignedIndex_t>(index));
      (*a_links)[index].setEdgeConnectivity(0, nullptr);
      (*a_links)[index].setEdgeConnectivity(
          1, &(*a_links)[shared_edge_link[index]]);
      std::size_t neighbor = (e + 1) % face.size() + f * 4;
      (*a_links)[neighbor].setEdgeConnectivity(2, &(*a_links)[index]);
      (*a_links)[index].setEdgeConnectivity(3, &(*a_links)[neighbor]);
    }
  }

  if (PRINT_OBJECTS) {
    tet[3] = volume_centroid;
    for (std::size_t f = 0; f < face_indices.size(); ++f) {
      tet[2] = face_centroids[f];
      const auto &face = face_indices[f];
      for (std::size_t e = 0; e < face.size(); ++e) {
        tet[0] = cuboid_pts[face[e]];
        tet[1] = cuboid_pts[face[(e + 1) % face.size()]];

        const std::size_t index = e + f * 4;
        auto poly = IRL::Tet::fromRawPtPointer(4, tet.data());
        auto &half_edge = IRL::setHalfEdgeStructure(poly);
        auto segmented = half_edge.generateSegmentedPolyhedron();
        std::string filename = "tet_mesh_" + std::to_string(index) + ".vtu";
        std::ofstream myfile;
        myfile.open(filename);
        myfile << segmented;
        myfile.close();
      }
    }
  }
}

void irl_gvm_cubeOntoTetMesh(const double *a_cube_pts, double *a_volume,
                             double *a_times, std::size_t *a_entered_cells) {
  std::array<IRL::Pt, 2> mesh_box{
      {IRL::Pt(-3.0, -3.0, -3.0), IRL::Pt(3.0, 3.0, 3.0)}};
  static std::vector<IRL::PlanarLocalizer> localizers;
  static std::vector<IRL::LocalizerLink> links;
  static bool mesh_made = false;
  if (!mesh_made) {
    mesh_made = true;
    setupTetMesh(mesh_box, &localizers, &links);
  }

  // This will return the volume from each cell, tagged with that cells unique
  // Id set during setupCubicMesh
  const auto cube = IRL::RectangularCuboid::fromBoundingPts(
      IRL::Pt(a_cube_pts[0], a_cube_pts[1], a_cube_pts[2]),
      IRL::Pt(a_cube_pts[3], a_cube_pts[4], a_cube_pts[5]));
  double start = omp_get_wtime();
  auto tagged_volumes =
      IRL::getVolumeMoments<IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>(
          cube, links[0]);
  double end = omp_get_wtime();
  a_times[0] = end - start;

  // Return volume that is the sum of the distributed volume
  *a_volume = 0.0;
  for (const auto &volume : tagged_volumes) {
    *a_volume += volume.volume_moments_m;
  }
  *a_entered_cells = static_cast<std::size_t>(tagged_volumes.size());
}

void irl_gvm_stelIcosahedronOntoTetMesh(const double *a_stel_icosahedron_pts,
                                        double *a_volume, double *a_times,
                                        std::size_t *a_entered_cells) {
  std::array<IRL::Pt, 2> mesh_box{
      {IRL::Pt(-3.0, -3.0, -3.0), IRL::Pt(3.0, 3.0, 3.0)}};
  static std::vector<IRL::PlanarLocalizer> localizers;
  static std::vector<IRL::LocalizerLink> links;
  static bool mesh_made = false;
  if (!mesh_made) {
    mesh_made = true;
    setupTetMesh(mesh_box, &localizers, &links);
  }

  // This will return the volume from each cell, tagged with that cells unique
  // Id set during setupCubicMesh
  const auto stel_icosahedron = IRL::StellatedIcosahedron::fromRawDoublePointer(
      32, a_stel_icosahedron_pts);
  double start = omp_get_wtime();
  auto tagged_volumes =
      IRL::getVolumeMoments<IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>(
          stel_icosahedron, links[0]);
  double end = omp_get_wtime();
  a_times[0] = end - start;

  // Return volume that is the sum of the distributed volume
  *a_volume = 0.0;
  for (const auto &volume : tagged_volumes) {
    *a_volume += volume.volume_moments_m;
  }
  *a_entered_cells = static_cast<std::size_t>(tagged_volumes.size());
}

static IRL::Pt convertToCartesian(const IRL::Pt &a_spherical_pt) {
  return {a_spherical_pt[0] * std::sin(a_spherical_pt[2]) *
              std::cos(a_spherical_pt[1]),
          a_spherical_pt[0] * std::sin(a_spherical_pt[2]) *
              std::sin(a_spherical_pt[1]),
          a_spherical_pt[0] * std::cos(a_spherical_pt[2])};
}

// This function sets up the Spherical Cartesian mesh. In IRL,
// this means generating PlanarLocalizer objects
// to represent the cell faces (and serve as edges
// in the graph) and connecting them together to form
// LocalizerLink objects. Together, the LocalizerLink objects
// for the graph of the mesh over which we will distribute.
static void
setupSphericalCartesianMesh(const double a_radius, const int a_radial_cells,
                         std::vector<IRL::PlanarLocalizer> *a_localizers,
                         std::vector<IRL::LocalizerLink> *a_links) {

  assert(a_radial_cells > 0);

  constexpr static int a_theta_cells = 3;
  constexpr static int a_phi_cells = 3;
  IRL::Pt center(0.0, 0.0, 0.0);

  constexpr static double d_theta =
      2.0 * M_PI / static_cast<double>(a_theta_cells);
  constexpr static double d_phi = 1.0 * M_PI / static_cast<double>(a_phi_cells);
  const double d_radius = a_radius / static_cast<double>(a_radial_cells);

  a_localizers->resize(1 + (a_radial_cells - 1) * a_theta_cells * a_phi_cells);
  a_links->resize(1 + (a_radial_cells - 1) * a_theta_cells * a_phi_cells);
  for (std::size_t n = 0; n < a_localizers->size(); ++n) {
    (*a_links)[n] = IRL::LocalizerLink(&(*a_localizers)[n]);
  }
  (*a_links)[0].setId(0);

  // First reconstruction is a convex sphere like object
  (*a_localizers)[0].setNumberOfPlanes(a_theta_cells * a_phi_cells);
  for (int p = 0; p < a_phi_cells; ++p) {
    for (int t = 0; t < a_theta_cells; ++t) {
      IRL::Normal normal = convertToCartesian(
          IRL::Pt(1.0, (static_cast<double>(t) + 0.5) * d_theta,
                  (static_cast<double>(p) + 0.5) * d_phi));
      (*a_localizers)[0][t + p * a_theta_cells] = IRL::Plane(normal, d_radius);
    }
  }

  // Rest of reconstructions wrap around initial one, propagating outwards
  // towards the sphere radius
  for (int r = 1; r < a_radial_cells; ++r) {
    for (int p = 0; p < a_phi_cells; ++p) {
      for (int t = 0; t < a_theta_cells; ++t) {
        const int index =
            1 + t + p * a_theta_cells + (r - 1) * a_theta_cells * a_phi_cells;
        IRL::Normal normal = convertToCartesian(
            IRL::Pt(1.0, (static_cast<double>(t) + 0.5) * d_theta,
                    (static_cast<double>(p) + 0.5) * d_phi));
        (*a_localizers)[index].setNumberOfPlanes(6);
        // Bottom face pointing towards center
        (*a_localizers)[index][0] =
            IRL::Plane(-normal, -d_radius * static_cast<double>(r));
        // Top face pointing towards outside
        (*a_localizers)[index][1] =
            IRL::Plane(normal, d_radius * static_cast<double>(r + 1));

        // Side face traveling in -theta direction
        double angle = static_cast<double>(t) * d_theta;
        normal = -IRL::Normal(-std::sin(angle), std::cos(angle), 0.0);
        (*a_localizers)[index][2] = IRL::Plane(normal, 0.0);

        // Side face traveling in +theta direction
        angle = static_cast<double>(t) * d_theta;
        normal = IRL::Normal(-std::sin(angle), std::cos(angle), 0.0);
        (*a_localizers)[index][3] = IRL::Plane(normal, 0.0);

        // Side face traveling in -phi direction
        angle = static_cast<double>(p) * d_phi;
        normal = -IRL::Normal(std::cos(angle), 0.0, -std::sin(angle));
        (*a_localizers)[index][4] = IRL::Plane(normal, 0.0);

        // Side face traveling in +phi direction
        angle = static_cast<double>(p) * d_phi;
        normal = IRL::Normal(std::cos(angle), 0.0, -std::sin(angle));
        (*a_localizers)[index][5] = IRL::Plane(normal, 0.0);

        (*a_links)[index].setId(index);
        if (r != 1) {
          (*a_links)[index].setEdgeConnectivity(
              0, &(*a_links)[index - a_theta_cells * a_phi_cells]);
        }
        if (r < a_radial_cells - 1) {
          (*a_links)[index].setEdgeConnectivity(
              1, &(*a_links)[index + a_theta_cells * a_phi_cells]);
        }
        int neighbor = (t + 1) % a_theta_cells + 1 + p * a_theta_cells +
                       (r - 1) * a_theta_cells * a_phi_cells;
        (*a_links)[neighbor].setEdgeConnectivity(2, &(*a_links)[index]);
        (*a_links)[index].setEdgeConnectivity(3, &(*a_links)[neighbor]);
        neighbor = 1 + t + ((p + 1) % a_phi_cells) * a_theta_cells +
                   (r - 1) * a_theta_cells * a_phi_cells;
        (*a_links)[neighbor].setEdgeConnectivity(4, &(*a_links)[index]);
        (*a_links)[index].setEdgeConnectivity(5, &(*a_links)[neighbor]);
      }
    }
  }
  // Have more than one localizer
  if (a_radial_cells > 1) {
    for (int p = 0; p < a_phi_cells; ++p) {
      for (int t = 0; t < a_theta_cells; ++t) {
        const int index = 1 + t + p * a_theta_cells;
        (*a_links)[0].setEdgeConnectivity(index - 1, &(*a_links)[index]);
        (*a_links)[index].setEdgeConnectivity(0, &(*a_links)[0]);
      }
    }
  }

  if (PRINT_OBJECTS) {
    std::array<IRL::Pt, 8> pts;
    std::vector<std::vector<IRL::UnsignedIndex_t>> face_brep(6);
    face_brep[0] = std::vector<IRL::UnsignedIndex_t>({0, 3, 2, 1});
    face_brep[1] = std::vector<IRL::UnsignedIndex_t>({4, 5, 6, 7});
    face_brep[2] = std::vector<IRL::UnsignedIndex_t>({0, 4, 7, 3});
    face_brep[3] = std::vector<IRL::UnsignedIndex_t>({1, 2, 6, 5});
    face_brep[4] = std::vector<IRL::UnsignedIndex_t>({7, 6, 2, 3});
    face_brep[5] = std::vector<IRL::UnsignedIndex_t>({4, 0, 1, 5});
    IRL::PolyhedronConnectivity connectivity(face_brep);
    for (int r = 1; r < a_radial_cells; ++r) {
      for (int p = 0; p < a_phi_cells; ++p) {
        for (int t = 0; t < a_theta_cells; ++t) {
          pts[0] =
              convertToCartesian(IRL::Pt(static_cast<double>(r) * d_radius,
                                         (static_cast<double>(t)) * d_theta,
                                         (static_cast<double>(p)) * d_phi));
          pts[1] =
              convertToCartesian(IRL::Pt(static_cast<double>(r) * d_radius,
                                         (static_cast<double>(t + 1)) * d_theta,
                                         (static_cast<double>(p)) * d_phi));
          pts[2] =
              convertToCartesian(IRL::Pt(static_cast<double>(r) * d_radius,
                                         (static_cast<double>(t + 1)) * d_theta,
                                         (static_cast<double>(p + 1)) * d_phi));
          pts[3] =
              convertToCartesian(IRL::Pt(static_cast<double>(r) * d_radius,
                                         (static_cast<double>(t)) * d_theta,
                                         (static_cast<double>(p + 1)) * d_phi));
          pts[4] =
              convertToCartesian(IRL::Pt(static_cast<double>(r + 1) * d_radius,
                                         (static_cast<double>(t)) * d_theta,
                                         (static_cast<double>(p)) * d_phi));
          pts[5] =
              convertToCartesian(IRL::Pt(static_cast<double>(r + 1) * d_radius,
                                         (static_cast<double>(t + 1)) * d_theta,
                                         (static_cast<double>(p)) * d_phi));
          pts[6] =
              convertToCartesian(IRL::Pt(static_cast<double>(r + 1) * d_radius,
                                         (static_cast<double>(t + 1)) * d_theta,
                                         (static_cast<double>(p + 1)) * d_phi));
          pts[7] =
              convertToCartesian(IRL::Pt(static_cast<double>(r + 1) * d_radius,
                                         (static_cast<double>(t)) * d_theta,
                                         (static_cast<double>(p + 1)) * d_phi));

          const int index =
              1 + t + p * a_theta_cells + (r - 1) * a_theta_cells * a_phi_cells;
          IRL::GeneralPolyhedron poly(pts, &connectivity);
          auto &half_edge = IRL::setHalfEdgeStructure(poly);
          auto segmented = half_edge.generateSegmentedPolyhedron();
          std::string filename =
              "spherical_mesh_" + std::to_string(index) + ".vtu";
          std::ofstream myfile;
          myfile.open(filename);
          myfile << segmented;
          myfile.close();
        }
      }
    }
  }
}

void irl_gvm_cubeOntoSphericalCartesianMesh(const double *a_cube_pts,
                                         double *a_volume, double *a_times,
                                         std::size_t *a_entered_cells) {
  const double radius = 0.5 * std::sqrt(3.0 * 6.0 * 6.0);
  static constexpr bool print_to_files = false;
  static std::vector<IRL::PlanarLocalizer> localizers;
  static std::vector<IRL::LocalizerLink> links;
  static bool mesh_made = false;
  if (!mesh_made) {
    mesh_made = true;
    setupSphericalCartesianMesh(radius, 4, &localizers, &links);
  }

  // This will return the volume from each cell, tagged with that cells unique
  // Id set during setupCubicMesh
  const auto cube = IRL::RectangularCuboid::fromBoundingPts(
      IRL::Pt(a_cube_pts[0], a_cube_pts[1], a_cube_pts[2]),
      IRL::Pt(a_cube_pts[3], a_cube_pts[4], a_cube_pts[5]));

  double start = omp_get_wtime();
  auto tagged_volumes =
      IRL::getVolumeMoments<IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>(
          cube, links[0]);
  double end = omp_get_wtime();
  a_times[0] = end - start;

  // Return volume that is the sum of the distributed volume
  *a_volume = 0.0;
  for (const auto &volume : tagged_volumes) {
    *a_volume += volume.volume_moments_m;
  }
  *a_entered_cells = static_cast<std::size_t>(tagged_volumes.size());
}

void irl_gvm_stelIcosahedronOntoSphericalCartesianMesh(
    const double *a_stel_icosahedron_pts, double *a_volume, double *a_times,
    std::size_t *a_entered_cells) {
  const double radius = 0.5 * std::sqrt(3.0 * 6.0 * 6.0);
  static constexpr bool print_to_files = false;
  static std::vector<IRL::PlanarLocalizer> localizers;
  static std::vector<IRL::LocalizerLink> links;
  static bool mesh_made = false;
  if (!mesh_made) {
    mesh_made = true;
    setupSphericalCartesianMesh(radius, 4, &localizers, &links);
  }

  // This will return the volume from each cell, tagged with that cells unique
  // Id set during setupCubicMesh
  const auto stel_icosahedron = IRL::StellatedIcosahedron::fromRawDoublePointer(
      32, a_stel_icosahedron_pts);

  double start = omp_get_wtime();
  auto tagged_volumes =
      IRL::getVolumeMoments<IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>(
          stel_icosahedron, links[0]);
  double end = omp_get_wtime();
  a_times[0] = end - start;

  // Return volume that is the sum of the distributed volume
  *a_volume = 0.0;
  for (const auto &volume : tagged_volumes) {
    *a_volume += volume.volume_moments_m;
  }
  *a_entered_cells = static_cast<std::size_t>(tagged_volumes.size());
}
