// Copyright (C) 2020 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/timing_comp/distribution_timing.h"

#include <omp.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>

#include "src/timing_comp/files.h"
#include "src/timing_comp/times.h"
#include "src/timing_comp/timing_comp.h"

// IRL GVM Timing includes
#include "src/irl_gvm/irl_gvm_distribution.h"
#include "src/irl_gvm/stellated_dodecahedron.h"
#include "src/irl_gvm/stellated_icosahedron.h"

void distributeCubeOntoCubicMesh(FILE* a_distribute_file,
                                 const int a_number_of_trials) {
  // Pass cube as lower and upper bounding box points
  const std::array<double, 6> cube{{-1.5, -1.5, -1.5, 1.5, 1.5, 1.5}};
  const auto total_volume = 27.0;

  // Will distribute onto a region covering [-1.5,1.5]^3
  Times<1> irl_gvm_times;
  std::size_t cells_entered = 0;
  double volume_error = 0.0;
  double abs_volume_error = -DBL_MAX;
  for (int n = 0; n < a_number_of_trials; ++n) {
    const double x_shift = randomDouble();
    const double y_shift = randomDouble();
    const double z_shift = randomDouble();
    std::array<double, 6> shifted_cube;
    for (std::size_t v = 0; v < shifted_cube.size(); v += 3) {
      shifted_cube[v] = cube[v] + x_shift;
      shifted_cube[v + 1] = cube[v + 1] + y_shift;
      shifted_cube[v + 2] = cube[v + 2] + z_shift;
    }
    Times<1> irl_gvm_trial_time;
    std::size_t trial_cells_entered;
    double irl_gvm_volume;

    irl_gvm_cubeOntoCubicMesh(shifted_cube.data(), &irl_gvm_volume,
                              irl_gvm_trial_time.data(), &trial_cells_entered);


    const double error = std::fabs(1.0 - irl_gvm_volume / total_volume);
    volume_error += error;
    abs_volume_error = std::max(abs_volume_error, error);    

    if (!sameVolumesFound(total_volume, irl_gvm_volume)) {
      printf("Shift of %20.12e %20.12e %20.12e\n", x_shift, y_shift, z_shift);
      printf(
          "Shifted cube defined by: \n (%20.12e %20.12e %20.12e), (%20.12e "
          "%20.12e %20.12e)\n",
          shifted_cube[0], shifted_cube[1], shifted_cube[2], shifted_cube[3],
          shifted_cube[4], shifted_cube[5]);
      std::exit(-1);
    }
    irl_gvm_times += irl_gvm_trial_time;
    cells_entered += trial_cells_entered;
  }
  // Write out time in seconds
  fprintf(a_distribute_file, "%19.13e %19.13e %19.13e %19.13e\n",
          static_cast<double>(cells_entered) /
              static_cast<double>(a_number_of_trials),
          volume_error / static_cast<double>(a_number_of_trials), abs_volume_error,
          irl_gvm_times[0]);
}

void distributeStelIcosahedronOntoCubicMesh(FILE* a_distribute_file,
                                            const int a_number_of_trials) {
  // A Stellated Icosahedron. Object is non-convex.
  // Note: Matches VOFtools NCICOSAMESH object
  // 32 points, ordered Pt 1 X/Y/Z, Pt 2 X/Y/Z, ...
  const auto stel_icosahedron_pts = getStelIcosahedronPts();

  // Set volume of whole object
  const auto stel_icosahedron = IRL::StellatedIcosahedron::fromRawDoublePointer(
      32, stel_icosahedron_pts.data());
  const double total_volume = stel_icosahedron.calculateVolume();

  // Will distribute onto a region covering [-3,3]^3
  Times<1> irl_gvm_times;
  std::size_t cells_entered = 0;
  double volume_error = 0.0;
  double abs_volume_error = -DBL_MAX;  
  for (int n = 0; n < a_number_of_trials; ++n) {
    const double x_shift = randomDouble();
    const double y_shift = randomDouble();
    const double z_shift = randomDouble();
    std::array<double, 96> shifted_stel_icosahedron_pts;
    for (std::size_t v = 0; v < stel_icosahedron_pts.size(); v += 3) {
      shifted_stel_icosahedron_pts[v] = stel_icosahedron_pts[v] + x_shift;
      shifted_stel_icosahedron_pts[v + 1] =
          stel_icosahedron_pts[v + 1] + y_shift;
      shifted_stel_icosahedron_pts[v + 2] =
          stel_icosahedron_pts[v + 2] + z_shift;
    }

    Times<1> irl_gvm_trial_time;
    std::size_t trial_cells_entered;
    double irl_gvm_volume;

    irl_gvm_stelIcosahedronOntoCubicMesh(
        shifted_stel_icosahedron_pts.data(), &irl_gvm_volume,
        irl_gvm_trial_time.data(), &trial_cells_entered);

    const double error = std::fabs(1.0 - irl_gvm_volume / total_volume);
    volume_error += error;
    abs_volume_error = std::max(abs_volume_error, error);

    if (!sameVolumesFound(total_volume, irl_gvm_volume)) {
      printf("Shift of %20.12e %20.12e %20.12e\n", x_shift, y_shift, z_shift);
      std::exit(-1);
    }
    irl_gvm_times += irl_gvm_trial_time;
    cells_entered += trial_cells_entered;
  }
  // Write out time in seconds
  fprintf(a_distribute_file, "%19.13e %19.13e %19.13e %19.13e\n",
          static_cast<double>(cells_entered) /
              static_cast<double>(a_number_of_trials),
          volume_error / static_cast<double>(a_number_of_trials), abs_volume_error,
          irl_gvm_times[0]);
}

void distributeCubeOntoTetMesh(FILE* a_distribute_file,
                               const int a_number_of_trials) {
  // Pass cube as lower and upper bounding box points
  const std::array<double, 6> cube{{-1.5, -1.5, -1.5, 1.5, 1.5, 1.5}};
  const auto total_volume = 27.0;

  // Will distribute onto a region covering [-1.5,1.5]^3
  Times<1> irl_gvm_times;
  std::size_t cells_entered = 0;
  double volume_error = 0.0;
  double abs_volume_error = -DBL_MAX;  
  for (int n = 0; n < a_number_of_trials; ++n) {
    const double x_shift = randomDouble();
    const double y_shift = randomDouble();
    const double z_shift = randomDouble();
    std::array<double, 6> shifted_cube;
    for (std::size_t v = 0; v < shifted_cube.size(); v += 3) {
      shifted_cube[v] = cube[v] + x_shift;
      shifted_cube[v + 1] = cube[v + 1] + y_shift;
      shifted_cube[v + 2] = cube[v + 2] + z_shift;
    }
    Times<1> irl_gvm_trial_time;
    std::size_t trial_cells_entered;
    double irl_gvm_volume;

    irl_gvm_cubeOntoTetMesh(shifted_cube.data(), &irl_gvm_volume,
                            irl_gvm_trial_time.data(), &trial_cells_entered);

    const double error = std::fabs(1.0 - irl_gvm_volume / total_volume);
    volume_error += error;
    abs_volume_error = std::max(abs_volume_error, error);

    if (!sameVolumesFound(total_volume, irl_gvm_volume)) {
      printf("Shift of %20.12e %20.12e %20.12e\n", x_shift, y_shift, z_shift);
      printf(
          "Shifted cube defined by: \n (%20.12e %20.12e %20.12e), (%20.12e "
          "%20.12e %20.12e)\n",
          shifted_cube[0], shifted_cube[1], shifted_cube[2], shifted_cube[3],
          shifted_cube[4], shifted_cube[5]);
      std::exit(-1);
    }
    irl_gvm_times += irl_gvm_trial_time;
    cells_entered += trial_cells_entered;
  }
  // Write out time in seconds
  fprintf(a_distribute_file, "%19.13e %19.13e %19.13e %19.13e\n",
          static_cast<double>(cells_entered) /
              static_cast<double>(a_number_of_trials),
          volume_error / static_cast<double>(a_number_of_trials), abs_volume_error,
          irl_gvm_times[0]);
}

void distributeStelIcosahedronOntoTetMesh(FILE* a_distribute_file,
                                          const int a_number_of_trials) {
  // A Stellated Icosahedron. Object is non-convex.
  // Note: Matches VOFtools NCICOSAMESH object
  // 32 points, ordered Pt 1 X/Y/Z, Pt 2 X/Y/Z, ...
  const auto stel_icosahedron_pts = getStelIcosahedronPts();

  // Set volume of whole object
  const auto stel_icosahedron = IRL::StellatedIcosahedron::fromRawDoublePointer(
      32, stel_icosahedron_pts.data());
  const double total_volume = stel_icosahedron.calculateVolume();

  // Will distribute onto a region covering [-3,3]^3
  Times<1> irl_gvm_times;
  std::size_t cells_entered = 0;
  double volume_error = 0.0;
  double abs_volume_error = -DBL_MAX;  
  for (int n = 0; n < a_number_of_trials; ++n) {
    const double x_shift = randomDouble();
    const double y_shift = randomDouble();
    const double z_shift = randomDouble();
    std::array<double, 96> shifted_stel_icosahedron_pts;
    for (std::size_t v = 0; v < stel_icosahedron_pts.size(); v += 3) {
      shifted_stel_icosahedron_pts[v] = stel_icosahedron_pts[v] + x_shift;
      shifted_stel_icosahedron_pts[v + 1] =
          stel_icosahedron_pts[v + 1] + y_shift;
      shifted_stel_icosahedron_pts[v + 2] =
          stel_icosahedron_pts[v + 2] + z_shift;
    }
    Times<1> irl_gvm_trial_time;
    std::size_t trial_cells_entered;
    double irl_gvm_volume;

    irl_gvm_stelIcosahedronOntoTetMesh(
        shifted_stel_icosahedron_pts.data(), &irl_gvm_volume,
        irl_gvm_trial_time.data(), &trial_cells_entered);

    const double error = std::fabs(1.0 - irl_gvm_volume / total_volume);
    volume_error += error;
    abs_volume_error = std::max(abs_volume_error, error);

    if (!sameVolumesFound(total_volume, irl_gvm_volume)) {
      printf("Shift of %20.12e %20.12e %20.12e\n", x_shift, y_shift, z_shift);
      std::exit(-1);
    }
    irl_gvm_times += irl_gvm_trial_time;
    cells_entered += trial_cells_entered;
  }
  // Write out time in seconds
  fprintf(a_distribute_file, "%19.13e %19.13e %19.13e %19.13e\n",
          static_cast<double>(cells_entered) /
              static_cast<double>(a_number_of_trials),
          volume_error / static_cast<double>(a_number_of_trials), abs_volume_error,
          irl_gvm_times[0]);
}

void distributeCubeOntoSphericalCartesianMesh(FILE* a_distribute_file,
                                              const int a_number_of_trials) {
  // Pass cube as lower and upper bounding box points
  const std::array<double, 6> cube{{-1.5, -1.5, -1.5, 1.5, 1.5, 1.5}};
  const auto total_volume = 27.0;

  // Will distribute onto a region covering [-1.5,1.5]^3
  Times<1> irl_gvm_times;
  std::size_t cells_entered = 0;
  double volume_error = 0.0;
  double abs_volume_error = -DBL_MAX;  
  for (int n = 0; n < a_number_of_trials; ++n) {
    const double x_shift = randomDouble();
    const double y_shift = randomDouble();
    const double z_shift = randomDouble();
    std::array<double, 6> shifted_cube;
    for (std::size_t v = 0; v < shifted_cube.size(); v += 3) {
      shifted_cube[v] = cube[v] + x_shift;
      shifted_cube[v + 1] = cube[v + 1] + y_shift;
      shifted_cube[v + 2] = cube[v + 2] + z_shift;
    }
    Times<1> irl_gvm_trial_time;
    std::size_t trial_cells_entered;
    double irl_gvm_volume;

    irl_gvm_cubeOntoSphericalCartesianMesh(shifted_cube.data(), &irl_gvm_volume,
                                           irl_gvm_trial_time.data(),
                                           &trial_cells_entered);

    const double error = std::fabs(1.0 - irl_gvm_volume / total_volume);
    volume_error += error;
    abs_volume_error = std::max(abs_volume_error, error);

    if (!sameVolumesFound(total_volume, irl_gvm_volume)) {
      printf("Shift of %20.12e %20.12e %20.12e\n", x_shift, y_shift, z_shift);
      printf(
          "Shifted cube defined by: \n (%20.12e %20.12e %20.12e), (%20.12e "
          "%20.12e %20.12e)\n",
          shifted_cube[0], shifted_cube[1], shifted_cube[2], shifted_cube[3],
          shifted_cube[4], shifted_cube[5]);
      std::exit(-1);
    }
    irl_gvm_times += irl_gvm_trial_time;
    cells_entered += trial_cells_entered;
  }
  // Write out time in seconds
  fprintf(a_distribute_file, "%19.13e %19.13e %19.13e %19.13e\n",
          static_cast<double>(cells_entered) /
              static_cast<double>(a_number_of_trials),
          volume_error / static_cast<double>(a_number_of_trials), abs_volume_error,
          irl_gvm_times[0]);
}

void distributeStelIcosahedronOntoSphericalCartesianMesh(
    FILE* a_distribute_file, const int a_number_of_trials) {
  // A Stellated Icosahedron. Object is non-convex.
  // Note: Matches VOFtools NCICOSAMESH object
  // 32 points, ordered Pt 1 X/Y/Z, Pt 2 X/Y/Z, ...
  const auto stel_icosahedron_pts = getStelIcosahedronPts();

  // Set volume of whole object
  const auto stel_icosahedron = IRL::StellatedIcosahedron::fromRawDoublePointer(
      32, stel_icosahedron_pts.data());
  const double total_volume = stel_icosahedron.calculateVolume();

  // Will distribute onto a region covering [-3.0,3.0]^3
  Times<1> irl_gvm_times;
  std::size_t cells_entered = 0;
  double volume_error = 0.0;
  double abs_volume_error = -DBL_MAX;
  for (int n = 0; n < a_number_of_trials; ++n) {
    const double x_shift = randomDouble();
    const double y_shift = randomDouble();
    const double z_shift = randomDouble();
    std::array<double, 96> shifted_stel_icosahedron_pts;
    for (std::size_t v = 0; v < stel_icosahedron_pts.size(); v += 3) {
      shifted_stel_icosahedron_pts[v] = stel_icosahedron_pts[v] + x_shift;
      shifted_stel_icosahedron_pts[v + 1] =
          stel_icosahedron_pts[v + 1] + y_shift;
      shifted_stel_icosahedron_pts[v + 2] =
          stel_icosahedron_pts[v + 2] + z_shift;
    }
    Times<1> irl_gvm_trial_time;
    std::size_t trial_cells_entered;
    double irl_gvm_volume;

    irl_gvm_stelIcosahedronOntoSphericalCartesianMesh(
        shifted_stel_icosahedron_pts.data(), &irl_gvm_volume,
        irl_gvm_trial_time.data(), &trial_cells_entered);

    const double error = std::fabs(1.0 - irl_gvm_volume / total_volume);
    volume_error += error;
    abs_volume_error = std::max(abs_volume_error, error);

    if (!sameVolumesFound(total_volume, irl_gvm_volume)) {
      printf("Shift of %20.12e %20.12e %20.12e\n", x_shift, y_shift, z_shift);
      std::exit(-1);
    }
    irl_gvm_times += irl_gvm_trial_time;
    cells_entered += trial_cells_entered;
  }
  // Write out time in seconds
  fprintf(a_distribute_file, "%19.13e %19.13e %19.13e %19.13e\n",
          static_cast<double>(cells_entered) /
              static_cast<double>(a_number_of_trials),
          volume_error / static_cast<double>(a_number_of_trials), abs_volume_error,
          irl_gvm_times[0]);
}
