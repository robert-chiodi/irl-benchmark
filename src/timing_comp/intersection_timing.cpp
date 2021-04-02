// Copyright (C) 2020 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/timing_comp/intersection_timing.h"

#include <omp.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <vector>

#include "src/timing_comp/files.h"
#include "src/timing_comp/times.h"
#include "src/timing_comp/timing_comp.h"

// IRL Includes
#include "src/geometry/general/plane.h"
#include "src/geometry/general/pt.h"
#include "src/geometry/polyhedrons/dodecahedron.h"
#include "src/geometry/polyhedrons/octahedron.h"
#include "src/geometry/polyhedrons/rectangular_cuboid.h"
#include "src/geometry/polyhedrons/symmetric_decompositions/symmetric_hexahedron.h"
#include "src/geometry/polyhedrons/symmetric_decompositions/symmetric_triangular_prism.h"
#include "src/geometry/polyhedrons/triangular_prism.h"

// IRL Timing includes
#include "src/irl/irl_intersection.h"

// IRL GVM Timing includes
#include "src/irl_gvm/irl_gvm_intersection.h"
#include "src/irl_gvm/stellated_dodecahedron.h"
#include "src/irl_gvm/stellated_icosahedron.h"

// R3D Timing includes
extern "C" {
#include "src/r3d/r3d_intersection.h"
}

// VOFTools Timing includes
extern "C" {
#include "src/voftools/voftools_intersection.h"
}

void intersectPrismByPlanes(const Files& a_output_files,
                            const int a_number_of_trials,
                            const int a_max_planes) {
  // A Triangular Prism
  // 6 points, ordered Pt 1 X/Y/Z, Pt 2 X/Y/Z, ...
  std::array<double, 18> prism_pts{{1.0, 0.0, -1.0, 1.0, 1.0, 0.0, 1.0, 0.0,
                                    1.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0,
                                    0.0, 1.0}};

  // Get centroid to translate random planes to
  const auto prism =
      IRL::TriangularPrism::fromRawDoublePointer(6, prism_pts.data());
  const IRL::Pt centroid = prism.calculateCentroid();

  // Volume of object to scale by when comparing results for accuracy
  const double scale = prism.calculateVolume();

  // Will pass plane as Normx, Normy, Normz, Dist,
  // planes stacked contiguously, starting from 0
  std::vector<double> planes(a_max_planes * 4);
  IRL::Plane plane;
  for (int p = 1; p <= a_max_planes; ++p) {
    Times<4> irl_times;
    Times<4> r3d_times;
    Times<4> voftools_times;
    for (int n = 0; n < a_number_of_trials; ++n) {
      setRandomPlanes(&planes, p, centroid);
      Times<4> irl_trial_time;
      Times<4> r3d_trial_time;
      Times<4> voftools_trial_time;
      double irl_volume, r3d_volume, voftools_volume;

      irl_prismByPlanes(prism_pts.data(), p, planes.data(), &irl_volume,
                        irl_trial_time.data());

      r3d_prismByPlanes(prism_pts.data(), p, planes.data(), &r3d_volume,
                        r3d_trial_time.data());

      c_voftools_prismByPlanes(prism_pts.data(), p, planes.data(),
                               &voftools_volume, voftools_trial_time.data());

      if (!sameVolumesFound(irl_volume, r3d_volume, voftools_volume, scale)) {
        std::cout << "Planes are: \n";
        for (int rp = 0; rp < p; ++rp) {
          std::cout << "Normal : (" << planes[rp * 4 + 0] << " "
                    << planes[rp * 4 + 1] << " " << planes[rp * 4 + 2]
                    << ")\n  Distance : " << planes[rp * 4 + 3] << '\n'
                    << std::endl;
        }
        std::exit(-1);
      }

      irl_gvm_prismByPlanes(prism_pts.data(), p, planes.data(), &irl_volume,
                            irl_trial_time.data() + 3);

      r3d_prismByPlanes_total(prism_pts.data(), p, planes.data(), &r3d_volume,
                              r3d_trial_time.data() + 3);

      c_voftools_prismByPlanes_total(prism_pts.data(), p, planes.data(),
                                     &voftools_volume,
                                     voftools_trial_time.data() + 3);

      if (!sameVolumesFound(irl_volume, r3d_volume, voftools_volume, scale)) {
        std::cout << "Planes are: \n";
        for (int rp = 0; rp < p; ++rp) {
          std::cout << "Normal : (" << planes[rp * 4 + 0] << " "
                    << planes[rp * 4 + 1] << " " << planes[rp * 4 + 2]
                    << ")\n  Distance : " << planes[rp * 4 + 3] << '\n'
                    << std::endl;
        }
        std::exit(-1);
      }
      irl_times += irl_trial_time;
      r3d_times += r3d_trial_time;
      voftools_times += voftools_trial_time;
    }
    // Write out time in seconds
    writeTimes(a_output_files.irl, p, irl_times);
    writeTimes(a_output_files.r3d, p, r3d_times);
    writeTimes(a_output_files.voftools, p, voftools_times);
  }
}

void intersectUnitCubeByPlanes(const Files& a_output_files,
                               const int a_number_of_trials,
                               const int a_max_planes) {
  // Pass cube as lower and upper bounding box points
  std::array<double, 6> cube_pts{{-0.5, -0.5, -0.5, 0.5, 0.5, 0.5}};

  // Get centroid to translate random planes to
  const auto cube = IRL::RectangularCuboid::fromBoundingPts(
      IRL::Pt(cube_pts[0], cube_pts[1], cube_pts[2]),
      IRL::Pt(cube_pts[3], cube_pts[4], cube_pts[5]));
  const IRL::Pt centroid = cube.calculateCentroid();

  // Volume of object to scale by when comparing results for accuracy
  const double scale = cube.calculateVolume();

  // Will pass plane as Normx, Normy, Normz, Dist,
  // planes stacked contiguously, starting from 0
  std::vector<double> planes(a_max_planes * 4);
  IRL::Plane plane;
  for (int p = 1; p <= a_max_planes; ++p) {
    Times<4> irl_times;
    Times<4> r3d_times;
    Times<4> voftools_times;
    for (int n = 0; n < a_number_of_trials; ++n) {
      setRandomPlanes(&planes, p, centroid);
      Times<4> irl_trial_time;
      Times<4> r3d_trial_time;
      Times<4> voftools_trial_time;
      double irl_volume, r3d_volume, voftools_volume;

      irl_unitCubeByPlanes(cube_pts.data(), p, planes.data(), &irl_volume,
                           irl_trial_time.data());

      r3d_unitCubeByPlanes(cube_pts.data(), p, planes.data(), &r3d_volume,
                           r3d_trial_time.data());

      c_voftools_unitCubeByPlanes(cube_pts.data(), p, planes.data(),
                                  &voftools_volume, voftools_trial_time.data());

      if (!sameVolumesFound(irl_volume, r3d_volume, voftools_volume, scale)) {
        std::cout << "Planes are: \n";
        for (int rp = 0; rp < p; ++rp) {
          std::cout << "Normal : (" << planes[rp * 4 + 0] << " "
                    << planes[rp * 4 + 1] << " " << planes[rp * 4 + 2]
                    << ")\n  Distance : " << planes[rp * 4 + 3] << '\n'
                    << std::endl;
        }
        std::exit(-1);
      }

      irl_gvm_unitCubeByPlanes(cube_pts.data(), p, planes.data(), &irl_volume,
                               irl_trial_time.data() + 3);

      r3d_unitCubeByPlanes_total(cube_pts.data(), p, planes.data(), &r3d_volume,
                                 r3d_trial_time.data() + 3);

      c_voftools_unitCubeByPlanes_total(cube_pts.data(), p, planes.data(),
                                        &voftools_volume,
                                        voftools_trial_time.data() + 3);

      if (!sameVolumesFound(irl_volume, r3d_volume, voftools_volume, scale)) {
        std::cout << "Planes are: \n";
        for (int rp = 0; rp < p; ++rp) {
          std::cout << "Normal : (" << planes[rp * 4 + 0] << " "
                    << planes[rp * 4 + 1] << " " << planes[rp * 4 + 2]
                    << ")\n  Distance : " << planes[rp * 4 + 3] << '\n'
                    << std::endl;
        }
        std::exit(-1);
      }

      irl_times += irl_trial_time;
      r3d_times += r3d_trial_time;
      voftools_times += voftools_trial_time;
    }
    // Write out time in seconds
    writeTimes(a_output_files.irl, p, irl_times);
    writeTimes(a_output_files.r3d, p, r3d_times);
    writeTimes(a_output_files.voftools, p, voftools_times);
  }
}

void intersectTriPrismByPlanes(const Files& a_output_files,
                               const int a_number_of_trials,
                               const int a_max_planes) {
  // A Triangular Prism  with each quad-face triangulated across a diagonal.
  // Points perturbed to make this case non-convex.
  // 6 points, ordered Pt 1 X/Y/Z, Pt 2 X/Y/Z, ...
  std::array<double, 18> tri_prism_pts{{
      1.0,
      0.0,
      -1.0,
      1.0,
      1.0,
      0.0,
      1.0,
      -0.5,
      1.0,
      0.0,
      -0.5,
      -1.0,
      0.0,
      1.0,
      0.0,
      0.0,
      0.0,
      1.0,
  }};

  // Get centroid to translate random planes to
  const auto tri_prism =
      IRL::Octahedron::fromRawDoublePointer(6, tri_prism_pts.data());
  const IRL::Pt centroid = tri_prism.calculateCentroid();

  // Volume of object to scale by when comparing results for accuracy
  const double scale = tri_prism.calculateVolume();

  // Will pass plane as Normx, Normy, Normz, Dist,
  // planes stacked contiguously, starting from 0
  std::vector<double> planes(a_max_planes * 4);
  IRL::Plane plane;
  for (int p = 1; p <= a_max_planes; ++p) {
    Times<4> irl_times;
    Times<4> r3d_times;
    Times<4> voftools_times;
    for (int n = 0; n < a_number_of_trials; ++n) {
      setRandomPlanes(&planes, p, centroid);
      Times<4> irl_trial_time;
      Times<4> r3d_trial_time;
      Times<4> voftools_trial_time;
      double irl_volume, r3d_volume, voftools_volume;

      irl_triPrismByPlanes(tri_prism_pts.data(), p, planes.data(), &irl_volume,
                           irl_trial_time.data());

      r3d_triPrismByPlanes(tri_prism_pts.data(), p, planes.data(), &r3d_volume,
                           r3d_trial_time.data());

      c_voftools_triPrismByPlanes(tri_prism_pts.data(), p, planes.data(),
                                  &voftools_volume, voftools_trial_time.data());

      if (!sameVolumesFound(irl_volume, r3d_volume, voftools_volume, scale)) {
        std::cout << "Planes are: \n";
        for (int rp = 0; rp < p; ++rp) {
          std::cout << "Normal : (" << planes[rp * 4 + 0] << " "
                    << planes[rp * 4 + 1] << " " << planes[rp * 4 + 2]
                    << ")\n  Distance : " << planes[rp * 4 + 3] << '\n'
                    << std::endl;
        }
        std::exit(-1);
      }

      irl_gvm_triPrismByPlanes(tri_prism_pts.data(), p, planes.data(),
                               &irl_volume, irl_trial_time.data() + 3);

      r3d_triPrismByPlanes_total(tri_prism_pts.data(), p, planes.data(),
                                 &r3d_volume, r3d_trial_time.data() + 3);

      c_voftools_triPrismByPlanes_total(tri_prism_pts.data(), p, planes.data(),
                                        &voftools_volume,
                                        voftools_trial_time.data() + 3);

      if (!sameVolumesFound(irl_volume, r3d_volume, voftools_volume, scale)) {
        std::cout << "Planes are: \n";
        for (int rp = 0; rp < p; ++rp) {
          std::cout << "Normal : (" << planes[rp * 4 + 0] << " "
                    << planes[rp * 4 + 1] << " " << planes[rp * 4 + 2]
                    << ")\n  Distance : " << planes[rp * 4 + 3] << '\n'
                    << std::endl;
        }
        std::exit(-1);
      }
      irl_times += irl_trial_time;
      r3d_times += r3d_trial_time;
      voftools_times += voftools_trial_time;
    }
    // Write out time in seconds
    writeTimes(a_output_files.irl, p, irl_times);
    writeTimes(a_output_files.r3d, p, r3d_times);
    writeTimes(a_output_files.voftools, p, voftools_times);
  }
}

void intersectTriHexByPlanes(const Files& a_output_files,
                             const int a_number_of_trials,
                             const int a_max_planes) {
  // A Hexahedron with each face triangulated across a diagonal.
  // Four points are moved in Y to make this case non-convex.
  // 8 points, ordered Pt 1 X/Y/Z, Pt 2 X/Y/Z, ...
  std::array<double, 24> tri_hex_pts{{
      0.5,  -0.5,  -0.5, 0.5,  0.5,  -0.5, 0.5,  0.25, 0.5, 0.5,  -0.25, 0.5,
      -0.5, -0.25, -0.5, -0.5, 0.25, -0.5, -0.5, 0.5,  0.5, -0.5, -0.5,  0.5,
  }};

  // Get centroid to translate random planes to
  const auto tri_hex =
      IRL::Dodecahedron::fromRawDoublePointer(8, tri_hex_pts.data());
  const IRL::Pt centroid = tri_hex.calculateCentroid();

  // Volume of object to scale by when comparing results for accuracy
  const double scale = tri_hex.calculateVolume();

  // Will pass plane as Normx, Normy, Normz, Dist,
  // planes stacked contiguously, starting from 0
  std::vector<double> planes(a_max_planes * 4);
  IRL::Plane plane;
  for (int p = 1; p <= a_max_planes; ++p) {
    Times<4> irl_times;
    Times<4> r3d_times;
    Times<4> voftools_times;
    for (int n = 0; n < a_number_of_trials; ++n) {
      setRandomPlanes(&planes, p, centroid);
      Times<4> irl_trial_time;
      Times<4> r3d_trial_time;
      Times<4> voftools_trial_time;
      double irl_volume, r3d_volume, voftools_volume;

      irl_triHexByPlanes(tri_hex_pts.data(), p, planes.data(), &irl_volume,
                         irl_trial_time.data());

      r3d_triHexByPlanes(tri_hex_pts.data(), p, planes.data(), &r3d_volume,
                         r3d_trial_time.data());

      c_voftools_triHexByPlanes(tri_hex_pts.data(), p, planes.data(),
                                &voftools_volume, voftools_trial_time.data());

      if (!sameVolumesFound(irl_volume, r3d_volume, voftools_volume, scale)) {
        std::cout << "Planes are: \n";
        for (int rp = 0; rp < p; ++rp) {
          std::cout << "Normal : (" << planes[rp * 4 + 0] << " "
                    << planes[rp * 4 + 1] << " " << planes[rp * 4 + 2]
                    << ")\n  Distance : " << planes[rp * 4 + 3] << '\n'
                    << std::endl;
        }
        std::exit(-1);
      }

      irl_gvm_triHexByPlanes(tri_hex_pts.data(), p, planes.data(), &irl_volume,
                             irl_trial_time.data() + 3);

      r3d_triHexByPlanes_total(tri_hex_pts.data(), p, planes.data(),
                               &r3d_volume, r3d_trial_time.data() + 3);

      c_voftools_triHexByPlanes_total(tri_hex_pts.data(), p, planes.data(),
                                      &voftools_volume,
                                      voftools_trial_time.data() + 3);

      if (!sameVolumesFound(irl_volume, r3d_volume, voftools_volume, scale)) {
        std::cout << "Planes are: \n";
        for (int rp = 0; rp < p; ++rp) {
          std::cout << "Normal : (" << planes[rp * 4 + 0] << " "
                    << planes[rp * 4 + 1] << " " << planes[rp * 4 + 2]
                    << ")\n  Distance : " << planes[rp * 4 + 3] << '\n'
                    << std::endl;
        }
        std::exit(-1);
      }

      irl_times += irl_trial_time;
      r3d_times += r3d_trial_time;
      voftools_times += voftools_trial_time;
    }
    // Write out time in seconds
    writeTimes(a_output_files.irl, p, irl_times);
    writeTimes(a_output_files.r3d, p, r3d_times);
    writeTimes(a_output_files.voftools, p, voftools_times);
  }
}

void intersectSymPrismByPlanes(const Files& a_output_files,
                               const int a_number_of_trials,
                               const int a_max_planes) {
  // A Prism with each face triangulated to a face-internal point
  // Indentation of face points make this case non-convex.
  // 11 points, ordered Pt 1 X/Y/Z, Pt 2 X/Y/Z, ...
  static constexpr double one_over_sqrt2 = 1.0 / std::sqrt(2.0);
  static constexpr double indent = 0.3;
  std::array<double, 33> sym_prism_pts{{1.0,
                                        0.0,
                                        -1.0,
                                        1.0,
                                        1.0,
                                        0.0,
                                        1.0,
                                        0.0,
                                        1.0,
                                        0.0,
                                        0.0,
                                        -1.0,
                                        0.0,
                                        1.0,
                                        0.0,
                                        0.0,
                                        0.0,
                                        1.0,
                                        1.0 - indent,
                                        1.0 / 3.0,
                                        0.0,
                                        0.5,
                                        0.5 - indent * one_over_sqrt2,
                                        -0.5 + indent * one_over_sqrt2,
                                        0.5,
                                        0.5 - indent * one_over_sqrt2,
                                        0.5 - indent * one_over_sqrt2,
                                        0.5,
                                        indent,
                                        0.0,
                                        indent,
                                        1.0 / 3.0,
                                        0.0}};

  // Get centroid to translate random planes to
  const auto sym_prism = IRL::SymmetricTriangularPrism::fromRawDoublePointer(
      11, sym_prism_pts.data());
  const IRL::Pt centroid = sym_prism.calculateCentroid();

  // Volume of object to scale by when comparing results for accuracy
  const double scale = sym_prism.calculateVolume();

  // Will pass plane as Normx, Normy, Normz, Dist,
  // planes stacked contiguously, starting from 0
  std::vector<double> planes(a_max_planes * 4);
  IRL::Plane plane;
  for (int p = 1; p <= a_max_planes; ++p) {
    Times<4> irl_times;
    Times<4> r3d_times;
    Times<4> voftools_times;
    for (int n = 0; n < a_number_of_trials; ++n) {
      setRandomPlanes(&planes, p, centroid);
      Times<4> irl_trial_time;
      Times<4> r3d_trial_time;
      Times<4> voftools_trial_time;
      double irl_volume, r3d_volume, voftools_volume;

      irl_symPrismByPlanes(sym_prism_pts.data(), p, planes.data(), &irl_volume,
                           irl_trial_time.data());

      r3d_symPrismByPlanes(sym_prism_pts.data(), p, planes.data(), &r3d_volume,
                           r3d_trial_time.data());

      c_voftools_symPrismByPlanes(sym_prism_pts.data(), p, planes.data(),
                                  &voftools_volume, voftools_trial_time.data());

      if (!sameVolumesFound(irl_volume, r3d_volume, voftools_volume, scale)) {
        std::cout << "Planes are: \n";
        for (int rp = 0; rp < p; ++rp) {
          std::cout << "Normal : (" << planes[rp * 4 + 0] << " "
                    << planes[rp * 4 + 1] << " " << planes[rp * 4 + 2]
                    << ")\n  Distance : " << planes[rp * 4 + 3] << '\n'
                    << std::endl;
        }
        std::exit(-1);
      }

      irl_gvm_symPrismByPlanes(sym_prism_pts.data(), p, planes.data(),
                               &irl_volume, irl_trial_time.data() + 3);

      r3d_symPrismByPlanes_total(sym_prism_pts.data(), p, planes.data(),
                                 &r3d_volume, r3d_trial_time.data() + 3);

      c_voftools_symPrismByPlanes_total(sym_prism_pts.data(), p, planes.data(),
                                        &voftools_volume,
                                        voftools_trial_time.data() + 3);

      if (!sameVolumesFound(irl_volume, r3d_volume, voftools_volume, scale)) {
        std::cout << "Planes are: \n";
        for (int rp = 0; rp < p; ++rp) {
          std::cout << "Normal : (" << planes[rp * 4 + 0] << " "
                    << planes[rp * 4 + 1] << " " << planes[rp * 4 + 2]
                    << ")\n  Distance : " << planes[rp * 4 + 3] << '\n'
                    << std::endl;
        }
        std::exit(-1);
      }
      irl_times += irl_trial_time;
      r3d_times += r3d_trial_time;
      voftools_times += voftools_trial_time;
    }
    // Write out time in seconds
    writeTimes(a_output_files.irl, p, irl_times);
    writeTimes(a_output_files.r3d, p, r3d_times);
    writeTimes(a_output_files.voftools, p, voftools_times);
  }
}

void intersectSymHexByPlanes(const Files& a_output_files,
                             const int a_number_of_trials,
                             const int a_max_planes) {
  // A Hexahedron with each face triangulated to a face-internal point
  // Each face is indented by 30% of the width (0.3), making this non-convex.
  // 14 points, ordered Pt 1 X/Y/Z, Pt 2 X/Y/Z, ...
  std::array<double, 42> sym_hex_pts{
      {0.5,  -0.5, -0.5, 0.5,  0.5,  -0.5, 0.5,  0.5,  0.5, 0.5, -0.5,
       0.5,  -0.5, -0.5, -0.5, -0.5, 0.5,  -0.5, -0.5, 0.5, 0.5, -0.5,
       -0.5, 0.5,  0.2,  0.0,  0.0,  0.0,  0.0,  -0.2, 0.0, 0.2, 0.0,
       0.0,  0.0,  0.2,  0.0,  -0.2, 0.0,  -0.2, 0.0,  0.0}};

  // Get centroid to translate random planes to
  const auto sym_hex =
      IRL::SymmetricHexahedron::fromRawDoublePointer(14, sym_hex_pts.data());
  const IRL::Pt centroid = sym_hex.calculateCentroid();

  // Volume of object to scale by when comparing results for accuracy
  const double scale = sym_hex.calculateVolume();

  // Will pass plane as Normx, Normy, Normz, Dist,
  // planes stacked contiguously, starting from 0
  std::vector<double> planes(a_max_planes * 4);
  IRL::Plane plane;
  for (int p = 1; p <= a_max_planes; ++p) {
    Times<4> irl_times;
    Times<4> r3d_times;
    Times<4> voftools_times;
    for (int n = 0; n < a_number_of_trials; ++n) {
      setRandomPlanes(&planes, p, centroid);
      Times<4> irl_trial_time;
      Times<4> r3d_trial_time;
      Times<4> voftools_trial_time;
      double irl_volume, r3d_volume, voftools_volume;

      irl_symHexByPlanes(sym_hex_pts.data(), p, planes.data(), &irl_volume,
                         irl_trial_time.data());

      r3d_symHexByPlanes(sym_hex_pts.data(), p, planes.data(), &r3d_volume,
                         r3d_trial_time.data());

      c_voftools_symHexByPlanes(sym_hex_pts.data(), p, planes.data(),
                                &voftools_volume, voftools_trial_time.data());

      if (!sameVolumesFound(irl_volume, r3d_volume, voftools_volume, scale)) {
        std::cout << "Planes are: \n";
        for (int rp = 0; rp < p; ++rp) {
          std::cout << "Normal : (" << planes[rp * 4 + 0] << " "
                    << planes[rp * 4 + 1] << " " << planes[rp * 4 + 2]
                    << ")\n  Distance : " << planes[rp * 4 + 3] << '\n'
                    << std::endl;
        }
        std::exit(-1);
      }

      irl_gvm_symHexByPlanes(sym_hex_pts.data(), p, planes.data(), &irl_volume,
                             irl_trial_time.data() + 3);

      r3d_symHexByPlanes_total(sym_hex_pts.data(), p, planes.data(),
                               &r3d_volume, r3d_trial_time.data() + 3);

      c_voftools_symHexByPlanes_total(sym_hex_pts.data(), p, planes.data(),
                                      &voftools_volume,
                                      voftools_trial_time.data() + 3);

      if (!sameVolumesFound(irl_volume, r3d_volume, voftools_volume, scale)) {
        std::cout << "Planes are: \n";
        for (int rp = 0; rp < p; ++rp) {
          std::cout << "Normal : (" << planes[rp * 4 + 0] << " "
                    << planes[rp * 4 + 1] << " " << planes[rp * 4 + 2]
                    << ")\n  Distance : " << planes[rp * 4 + 3] << '\n'
                    << std::endl;
        }
        std::exit(-1);
      }

      irl_times += irl_trial_time;
      r3d_times += r3d_trial_time;
      voftools_times += voftools_trial_time;
    }
    // Write out time in seconds
    writeTimes(a_output_files.irl, p, irl_times);
    writeTimes(a_output_files.r3d, p, r3d_times);
    writeTimes(a_output_files.voftools, p, voftools_times);
  }
}

void intersectStelDodecahedronByPlanes(const Files& a_output_files,
                                       const int a_number_of_trials,
                                       const int a_max_planes) {
  // A Stellated Icosahedron. Object is non-convex.
  // 32 points, ordered Pt 1 X/Y/Z, Pt 2 X/Y/Z, ...
  // Matches definition in VOFTools
  const auto stel_dodecahedron_pts = getStelDodecahedronPts();

  // Get centroid to translate random planes to
  const auto stel_dodecahedron =
      IRL::StellatedDodecahedron::fromRawDoublePointer(
          32, stel_dodecahedron_pts.data());
  const IRL::Pt centroid = stel_dodecahedron.calculateCentroid();

  // Volume of object to scale by when comparing results for accuracy
  const double scale = stel_dodecahedron.calculateVolume();

  // Will pass plane as Normx, Normy, Normz, Dist,
  // planes stacked contiguously, starting from 0
  std::vector<double> planes(a_max_planes * 4);
  IRL::Plane plane;
  for (int p = 1; p <= a_max_planes; ++p) {
    Times<4> irl_times;
    Times<4> r3d_times;
    Times<4> voftools_times;
    for (int n = 0; n < a_number_of_trials; ++n) {
      setRandomPlanes(&planes, p, centroid);
      Times<4> irl_trial_time;
      Times<4> r3d_trial_time;
      Times<4> voftools_trial_time;
      double irl_volume, r3d_volume, voftools_volume;

      irl_stelDodecahedronByPlanes(stel_dodecahedron_pts.data(), p,
                                   planes.data(), &irl_volume,
                                   irl_trial_time.data());

      r3d_stelDodecahedronByPlanes(stel_dodecahedron_pts.data(), p,
                                   planes.data(), &r3d_volume,
                                   r3d_trial_time.data());

      c_voftools_stelDodecahedronByPlanes(stel_dodecahedron_pts.data(), p,
                                          planes.data(), &voftools_volume,
                                          voftools_trial_time.data());

      if (!sameVolumesFound(irl_volume, r3d_volume, voftools_volume, scale)) {
        std::cout << "Planes are: \n";
        for (int rp = 0; rp < p; ++rp) {
          std::cout << "Normal : (" << planes[rp * 4 + 0] << " "
                    << planes[rp * 4 + 1] << " " << planes[rp * 4 + 2]
                    << ")\n  Distance : " << planes[rp * 4 + 3] << '\n'
                    << std::endl;
        }
        std::exit(-1);
      }

      irl_gvm_stelDodecahedronByPlanes(stel_dodecahedron_pts.data(), p,
                                       planes.data(), &irl_volume,
                                       irl_trial_time.data() + 3);

      r3d_stelDodecahedronByPlanes_total(stel_dodecahedron_pts.data(), p,
                                         planes.data(), &r3d_volume,
                                         r3d_trial_time.data() + 3);

      c_voftools_stelDodecahedronByPlanes_total(stel_dodecahedron_pts.data(), p,
                                                planes.data(), &voftools_volume,
                                                voftools_trial_time.data() + 3);

      if (!sameVolumesFound(irl_volume, r3d_volume, voftools_volume, scale)) {
        std::cout << "Planes are: \n";
        for (int rp = 0; rp < p; ++rp) {
          std::cout << "Normal : (" << planes[rp * 4 + 0] << " "
                    << planes[rp * 4 + 1] << " " << planes[rp * 4 + 2]
                    << ")\n  Distance : " << planes[rp * 4 + 3] << '\n'
                    << std::endl;
        }
        std::exit(-1);
      }

      irl_times += irl_trial_time;
      r3d_times += r3d_trial_time;
      voftools_times += voftools_trial_time;
    }
    // Write out time in seconds
    writeTimes(a_output_files.irl, p, irl_times);
    writeTimes(a_output_files.r3d, p, r3d_times);
    writeTimes(a_output_files.voftools, p, voftools_times);
  }
}

void intersectStelIcosahedronByPlanes(const Files& a_output_files,
                                      const int a_number_of_trials,
                                      const int a_max_planes) {
  // A Stellated Icosahedron. Object is non-convex.
  // Note: Matches VOFtools NCICOSAMESH object
  // 32 points, ordered Pt 1 X/Y/Z, Pt 2 X/Y/Z, ...
  const auto stel_icosahedron_pts = getStelIcosahedronPts();

  // Get centroid to translate random planes to
  const auto stel_icosahedron = IRL::StellatedIcosahedron::fromRawDoublePointer(
      32, stel_icosahedron_pts.data());
  const IRL::Pt centroid = stel_icosahedron.calculateCentroid();

  // Volume of object to scale by when comparing results for accuracy
  const double scale = stel_icosahedron.calculateVolume();

  // Will pass plane as Normx, Normy, Normz, Dist,
  // planes stacked contiguously, starting from 0
  std::vector<double> planes(a_max_planes * 4);
  for (int p = 1; p <= a_max_planes; ++p) {
    Times<4> irl_times;
    Times<4> r3d_times;
    Times<4> voftools_times;
    for (int n = 0; n < a_number_of_trials; ++n) {
      setRandomPlanes(&planes, p, centroid);
      Times<4> irl_trial_time;
      Times<4> r3d_trial_time;
      Times<4> voftools_trial_time;
      double irl_volume, r3d_volume, voftools_volume;

      irl_stelIcosahedronByPlanes(stel_icosahedron_pts.data(), p, planes.data(),
                                  &irl_volume, irl_trial_time.data());

      r3d_stelIcosahedronByPlanes(stel_icosahedron_pts.data(), p, planes.data(),
                                  &r3d_volume, r3d_trial_time.data());

      c_voftools_stelIcosahedronByPlanes(stel_icosahedron_pts.data(), p,
                                         planes.data(), &voftools_volume,
                                         voftools_trial_time.data());

      if (!sameVolumesFound(irl_volume, r3d_volume, voftools_volume, scale)) {
        std::cout << "Planes are: \n";
        for (int rp = 0; rp < p; ++rp) {
          std::cout << "Normal : (" << planes[rp * 4 + 0] << " "
                    << planes[rp * 4 + 1] << " " << planes[rp * 4 + 2]
                    << ")\n  Distance : " << planes[rp * 4 + 3] << '\n'
                    << std::endl;
        }
        std::exit(-1);
      }

      irl_gvm_stelIcosahedronByPlanes(stel_icosahedron_pts.data(), p,
                                      planes.data(), &irl_volume,
                                      irl_trial_time.data() + 3);

      r3d_stelIcosahedronByPlanes_total(stel_icosahedron_pts.data(), p,
                                        planes.data(), &r3d_volume,
                                        r3d_trial_time.data() + 3);

      c_voftools_stelIcosahedronByPlanes_total(stel_icosahedron_pts.data(), p,
                                               planes.data(), &voftools_volume,
                                               voftools_trial_time.data() + 3);

      if (!sameVolumesFound(irl_volume, r3d_volume, voftools_volume, scale)) {
        std::cout << "Planes are: \n";
        for (int rp = 0; rp < p; ++rp) {
          std::cout << "Normal : (" << planes[rp * 4 + 0] << " "
                    << planes[rp * 4 + 1] << " " << planes[rp * 4 + 2]
                    << ")\n  Distance : " << planes[rp * 4 + 3] << '\n'
                    << std::endl;
        }
        std::exit(-1);
      }

      irl_times += irl_trial_time;
      r3d_times += r3d_trial_time;
      voftools_times += voftools_trial_time;
    }
    // Write out time in seconds
    writeTimes(a_output_files.irl, p, irl_times);
    writeTimes(a_output_files.r3d, p, r3d_times);
    writeTimes(a_output_files.voftools, p, voftools_times);
  }
}
