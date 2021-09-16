// Copyright (C) 2020 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/timing_comp/timing_comp.h"

#include <omp.h>

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>

#include "src/timing_comp/distribution_timing.h"
#include "src/timing_comp/files.h"
#include "src/timing_comp/intersection_timing.h"

// IRL Includes
#include "src/geometry/general/normal.h"
#include "src/geometry/general/plane.h"
#include "src/geometry/general/pt.h"

int main(int argc, char** argv) {
  printf("Timing with precision %16.8e\n", omp_get_wtick());

  if (argc != 5) {
    std::cout << "Inproper usage of command-line arguments" << std::endl;
    std::cout << "Four arguments should be supplied. They are:" << std::endl;
    std::cout << "1 -- Results to generate (chosen by integer): " << std::endl;
    std::cout << "     Sample file of interface planes (0)" << std::endl;
    std::cout << "     Plane intersections with polyhedra (1)" << std::endl;
    std::cout << "     Distribution of polyhedron onto meshes (2)" << std::endl;
    std::cout << "2 -- Number of trials per test (integer >1000)" << std::endl;
    std::cout
        << "3 -- Max number of planes to intersect at one time (integer >0)\n"
        << "     Note: This only has an effect if the first input on CLI is 1"
        << std::endl;
    std::cout
        << "4 -- Whether to produce section timings (0), total timings (1), or both (2)\n"
        << "     Note: This only has an effect if the first input on CLI is 1"
	<< std::endl;
    return -1;
  }

  // Set control variables
  const int case_number = std::stoi(std::string(argv[1]));
  const int number_of_trials = std::stoi(std::string(argv[2]));
  const int max_planes = std::stoi(std::string(argv[3]));
  const int timings_to_produce =  std::stoi(std::string(argv[4]));

  if (number_of_trials < 1000) {
    std::cout << "Requires number_of_trials set to >= 1000" << std::endl;
    return -1;
  }
  if (max_planes < 1) {
    std::cout << "Requires max_planes set to >= 1 " << std::endl;
    return -1;
  }

  switch (case_number) {
    // Just export randomly generated planes for display
    // with python script sample_example.py
    case 0: {
      // Write out a small sampling of planes generated
      // the same way they are for the tests.
      const int sample_size = number_of_trials;
      FILE* sample_file = fopen("samples.csv", "w");
      for (int n = 0; n < sample_size; ++n) {
        IRL::Plane plane;
        plane.normal() = randomNormal();
        plane.distance() = randomDouble() * 0.3;
        fprintf(sample_file, "%20.12e,%20.12e,%20.12e,%20.12e\n",
                plane.normal()[0], plane.normal()[1], plane.normal()[2],
                plane.distance());
      }
      fclose(sample_file);
      break;
    }

    // Perform intersections of sets of [1:max_planes] random planes
    // with different polyhedra. Writes out results to
    // irl_timing.txt, r3d_timing.txt, and voftools_timing.txt
    // Results can be displayed in aggregate by running
    // the python script postprocess.py
    case 1: {
      // Intersection of polyhedra with random planes
      auto output_files =
          Files("irl_timing.txt", "r3d_timing.txt", "voftools_timing.txt");

      output_files.writeToFiles(std::to_string(number_of_trials) + " " +
                                std::to_string(max_planes) + "\n\n");

      std::cout << "Intersecting Prism by Planes" << std::endl;
      intersectPrismByPlanes(output_files, number_of_trials, max_planes, timings_to_produce);

      output_files.writeToFiles("\n");

      std::cout << "Intersecting Unit Cube by Planes" << std::endl;
      intersectUnitCubeByPlanes(output_files, number_of_trials, max_planes, timings_to_produce);

      output_files.writeToFiles("\n");

      std::cout << "Intersecting Triangulated Prism by Planes" << std::endl;
      intersectTriPrismByPlanes(output_files, number_of_trials, max_planes, timings_to_produce);

      output_files.writeToFiles("\n");

      std::cout << "Intersecting Triangulated Hexahedron by Planes"
                << std::endl;
      intersectTriHexByPlanes(output_files, number_of_trials, max_planes, timings_to_produce);

      output_files.writeToFiles("\n");

      std::cout << "Intersecting Symmetric Prism by Planes" << std::endl;
      intersectSymPrismByPlanes(output_files, number_of_trials, max_planes, timings_to_produce);

      output_files.writeToFiles("\n");

      std::cout << "Intersecting Symmetric Hexahedron by Planes" << std::endl;
      intersectSymHexByPlanes(output_files, number_of_trials, max_planes, timings_to_produce);

      output_files.writeToFiles("\n");

      std::cout << "Intersecting Stellated Dodecahedron by Planes" << std::endl;
      intersectStelDodecahedronByPlanes(output_files, number_of_trials,
                                        max_planes, timings_to_produce);

      output_files.writeToFiles("\n");

      std::cout << "Intersecting Stellated Icosahedron by Planes" << std::endl;
      intersectStelIcosahedronByPlanes(output_files, number_of_trials,
                                       max_planes, timings_to_produce);

      break;
    }

    // Perform distribution of polyhedron over representative
    // meshes using IRL. Results are exported to
    // distribute_timing.txt.
    case 2: {
      FILE* distribute_file = fopen("distribute_timing.txt", "w");
      std::cout << "Distribute Cube onto Cubic Mesh" << std::endl;
      distributeCubeOntoCubicMesh(distribute_file, number_of_trials);

      fprintf(distribute_file, "\n");

      std::cout << "Distribute Stellated Icosahedron onto Cubic Mesh"
                << std::endl;
      distributeStelIcosahedronOntoCubicMesh(distribute_file, number_of_trials);

      fprintf(distribute_file, "\n");

      std::cout << "Distribute Cube onto Tet Mesh" << std::endl;
      distributeCubeOntoTetMesh(distribute_file, number_of_trials);

      fprintf(distribute_file, "\n");

      std::cout << "Distribute Stellated Icosahedron onto Tet Mesh"
                << std::endl;
      distributeStelIcosahedronOntoTetMesh(distribute_file, number_of_trials);

      fprintf(distribute_file, "\n");

      std::cout << "Distribute Cube onto Spherical Cartesian Mesh" << std::endl;
      distributeCubeOntoSphericalCartesianMesh(distribute_file,
                                               number_of_trials);

      fprintf(distribute_file, "\n");

      std::cout
          << "Distribute Stellated Icosahedron onto Spherical Cartesian Mesh"
          << std::endl;
      distributeStelIcosahedronOntoSphericalCartesianMesh(distribute_file,
                                                          number_of_trials);

      fclose(distribute_file);
      break;
    }

    default: {
      std::cout << "Unknown case switch (first CLI argument) of " << case_number
                << std::endl;
      std::cout << "Run timing_comp with no arguments to display help screen."
                << std::endl;
      return -1;
    }
  }
  return 0;
}

double randomDouble(void) {
  static std::random_device rd;
  static std::mt19937 gen(rd());
  static std::uniform_real_distribution<double> dist(-1.0, 1.0);
  return dist(gen);
}

IRL::Normal randomNormal(void) {
  static std::random_device rd;
  static std::mt19937 gen(rd());
  static std::uniform_real_distribution<double> dist(0.0, 1.0);
  const double theta = 2.0 * M_PI * dist(gen);
  const double phi = std::acos(1.0 - 2.0 * dist(gen));
  return {std::sin(phi) * std::cos(theta), std::sin(phi) * std::sin(theta),
          std::cos(phi)};
}

void setRandomPlanes(std::vector<double>* a_planes,
                     const int a_number_of_planes, const IRL::Pt& a_centroid) {
  IRL::Plane plane;
  for (int rp = 0; rp < a_number_of_planes; ++rp) {
    plane.normal() = randomNormal();
    plane.distance() = plane.normal() * a_centroid + randomDouble() * 0.3;
    (*a_planes)[rp * 4 + 0] = plane.normal()[0];
    (*a_planes)[rp * 4 + 1] = plane.normal()[1];
    (*a_planes)[rp * 4 + 2] = plane.normal()[2];
    (*a_planes)[rp * 4 + 3] = plane.distance();
  }
}

bool sameVolumesFound(double a_total_volume, double a_found_volume) {
  static constexpr double tolerance = 1.0e-14;

  if (std::fabs(1.0 - a_found_volume / a_total_volume) > tolerance) {
    printf("Found volume different from total volume!\n");
    printf("Found: %20.15e\n", a_found_volume);
    printf("Total: %20.15e\n", a_total_volume);
    fflush(stdout);
    return false;
  }
  return true;
}

bool sameVolumesFound(double a_irl_volume, double a_r3d_volume,
                      double a_voftools_volume, double a_scale) {
  a_irl_volume /= a_scale;
  a_r3d_volume /= a_scale;
  a_voftools_volume /= a_scale;
  static constexpr double tolerance = 1.0e-14;

  if (std::fabs(a_irl_volume - a_r3d_volume) > tolerance ||
      std::fabs(a_irl_volume - a_voftools_volume) > tolerance ||
      std::fabs(a_r3d_volume - a_voftools_volume) > tolerance) {
    printf("Different volumes returned!\n");
    printf("IRL (scaled): %20.15e\n", a_irl_volume);
    printf("R3D (scaled): %20.15e\n", a_r3d_volume);
    printf("VOFTools (scaled): %20.15e\n", a_voftools_volume);
    fflush(stdout);
    return false;
  }
  return true;
}

void writeTimes(FILE* a_file, const int a_number_of_planes,
                const Times<4>& a_times) {
  fprintf(a_file, "%19.13e %19.13e %19.13e %19.13e %19.13e\n",
          static_cast<double>(a_number_of_planes), a_times[0], a_times[1],
          a_times[2], a_times[3]);
}

std::array<double, 96> getStelDodecahedronPts(void) {
  // A Stellated Icosahedron. Object is non-convex.
  // 32 points, ordered Pt 1 X/Y/Z, Pt 2 X/Y/Z, ...
  // Matches that used in VOFTools
  std::array<IRL::Pt, 32> pts;
  static constexpr double A = 1.0 / std::sqrt(3.0);
  static constexpr double B = std::sqrt((3.0 - std::sqrt(5.0)) / 6.0);
  static constexpr double C = std::sqrt((3.0 + std::sqrt(5.0)) / 6.0);
  pts[0] = IRL::Pt(A, A, A);
  pts[1] = IRL::Pt(A, A, -A);
  pts[2] = IRL::Pt(A, -A, A);
  pts[3] = IRL::Pt(A, -A, -A);
  pts[4] = IRL::Pt(-A, A, A);
  pts[5] = IRL::Pt(-A, A, -A);
  pts[6] = IRL::Pt(-A, -A, A);
  pts[7] = IRL::Pt(-A, -A, -A);
  pts[8] = IRL::Pt(B, C, 0.0);
  pts[9] = IRL::Pt(-B, C, 0.0);
  pts[10] = IRL::Pt(B, -C, 0.0);
  pts[11] = IRL::Pt(-B, -C, 0.0);
  pts[12] = IRL::Pt(C, 0.0, B);
  pts[13] = IRL::Pt(C, 0.0, -B);
  pts[14] = IRL::Pt(-C, 0.0, B);
  pts[15] = IRL::Pt(-C, 0.0, -B);
  pts[16] = IRL::Pt(0.0, B, C);
  pts[17] = IRL::Pt(0.0, -B, C);
  pts[18] = IRL::Pt(0.0, B, -C);
  pts[19] = IRL::Pt(0.0, -B, -C);

  // Base faces, which will be used to place Stellation point on each face
  std::array<std::array<std::size_t, 5>, 12> base_faces{{{0, 8, 9, 4, 16},
                                                         {0, 16, 17, 2, 12},
                                                         {12, 2, 10, 3, 13},
                                                         {9, 5, 15, 14, 4},
                                                         {3, 19, 18, 1, 13},
                                                         {7, 11, 6, 14, 15},
                                                         {0, 12, 13, 1, 8},
                                                         {8, 1, 18, 5, 9},
                                                         {16, 4, 14, 6, 17},
                                                         {6, 11, 10, 2, 17},
                                                         {7, 15, 5, 18, 19},
                                                         {7, 19, 3, 10, 11}}};

  const double extrusion_length = 0.5;
  IRL::UnsignedIndex_t new_vertex = 20;
  for (const auto& face : base_faces) {
    IRL::Normal face_normal =
        IRL::crossProductNormalized(IRL::Pt(pts[face[1]] - pts[face[0]]),
                                    IRL::Pt(pts[face[2]] - pts[face[0]]));
    IRL::Pt average_pt(0.0, 0.0, 0.0);
    for (const auto& vertex : face) {
      average_pt += pts[vertex];
    }
    average_pt /= static_cast<double>(face.size());
    pts[new_vertex] = average_pt + extrusion_length * face_normal;
    ++new_vertex;
  }

  std::array<double, 96> stel_dodecahedron_pts;
  for (std::size_t n = 0; n < pts.size(); ++n) {
    for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
      stel_dodecahedron_pts[n * 3 + d] = pts[n][d];
    }
  }
  return stel_dodecahedron_pts;
}

std::array<double, 96> getStelIcosahedronPts(void) {
  // A Stellated Icosahedron. Object is non-convex.
  // Note: Matches VOFtools NCICOSAMESH object
  // 32 points, ordered Pt 1 X/Y/Z, Pt 2 X/Y/Z, ...
  std::array<IRL::Pt, 32> pts;
  double T = (1.0 + std::sqrt(5.0)) / 2.0;
  const double A = std::sqrt(1.0 + T * T);
  T /= A;
  const double Ainv = 1.0 / A;
  pts[0] = IRL::Pt(T, Ainv, 0.0);
  pts[1] = IRL::Pt(-T, Ainv, 0.0);
  pts[2] = IRL::Pt(T, -Ainv, 0.0);
  pts[3] = IRL::Pt(-T, -Ainv, 0.0);
  pts[4] = IRL::Pt(Ainv, 0.0, T);
  pts[5] = IRL::Pt(Ainv, 0.0, -T);
  pts[6] = IRL::Pt(-Ainv, 0.0, T);
  pts[7] = IRL::Pt(-Ainv, 0.0, -T);
  pts[8] = IRL::Pt(0.0, T, Ainv);
  pts[9] = IRL::Pt(0.0, -T, Ainv);
  pts[10] = IRL::Pt(0.0, T, -Ainv);
  pts[11] = IRL::Pt(0.0, -T, -Ainv);

  // Base faces, which will be used to place Stellation point on each face
  std::array<std::array<std::size_t, 3>, 20> base_faces{
      {{0, 8, 4},  {0, 5, 10}, {2, 4, 9},  {2, 11, 5}, {1, 6, 8},
       {1, 10, 7}, {3, 9, 6},  {3, 7, 11}, {0, 10, 8}, {1, 8, 10},
       {2, 9, 11}, {3, 11, 9}, {4, 2, 0},  {5, 0, 2},  {6, 1, 3},
       {7, 3, 1},  {8, 6, 4},  {9, 4, 6},  {10, 5, 7}, {11, 7, 5}}};

  const double extrusion_length = 1.0;
  IRL::UnsignedIndex_t new_vertex = 12;
  for (const auto& face : base_faces) {
    IRL::Normal face_normal =
        IRL::crossProductNormalized(IRL::Pt(pts[face[1]] - pts[face[0]]),
                                    IRL::Pt(pts[face[2]] - pts[face[0]]));
    IRL::Pt average_pt(0.0, 0.0, 0.0);
    for (const auto& vertex : face) {
      average_pt += pts[vertex];
    }
    average_pt /= static_cast<double>(face.size());
    pts[new_vertex] = average_pt + extrusion_length * face_normal;
    ++new_vertex;
  }

  std::array<double, 96> stel_icosahedron_pts;
  for (std::size_t n = 0; n < pts.size(); ++n) {
    for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
      stel_icosahedron_pts[n * 3 + d] = pts[n][d];
    }
  }
  return stel_icosahedron_pts;
}
