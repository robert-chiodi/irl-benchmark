// Copyright (C) 2020 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_TIMING_COMP_TIMING_COMP_H_
#define SRC_TIMING_COMP_TIMING_COMP_H_

#include <array>
#include <vector>

#include "src/timing_comp/times.h"

// IRL Includes
#include "src/geometry/general/normal.h"
#include "src/geometry/general/plane.h"
#include "src/geometry/general/pt.h"

int main(int argc, char** argv);

std::vector<IRL::Plane> generateSampledPlanes(void);

// Generate a random double in the range [-1,1]
// using the Mersenne Twister algorithm
double randomDouble(void);

// Draw normal from uniform distribution on a sphere.
IRL::Normal randomNormal(void);

// Helper function to create and store sets of planes.
void setRandomPlanes(std::vector<double>* a_planes,
                     const int a_number_of_planes, const IRL::Pt& a_centroid);

// Confirm no volume lost when distributing.
bool sameVolumesFound(double a_total_volume, double a_found_volume);

// Confirm all packages found same volume
bool sameVolumesFound(double a_irl_volume, double a_r3d_volume,
                      double a_voftools_volume, double a_scale);

// Helper function to write times to file in consistent format.
void writeTimes(FILE* a_file, const int a_number_of_planes,
                const Times<4>& a_times);

// Returns points for Stellated Dodecahedron object
// to be constructed by each package.
std::array<double, 96> getStelDodecahedronPts(void);

// Returns points for Stellated Icosahedron object
// to be constructed by each package.
std::array<double, 96> getStelIcosahedronPts(void);

#endif  // SRC_TIMING_COMP_TIMING_COMP_H_
