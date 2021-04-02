// Copyright (C) 2020 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_TIMING_COMP_INTERSECTION_TIMING_H_
#define SRC_TIMING_COMP_INTERSECTION_TIMING_H_

#include "src/timing_comp/files.h"

// Time intersection with Triangular Prism
void intersectPrismByPlanes(const Files& a_output_files,
                            const int a_number_of_trials,
                            const int a_max_planes);

// Time intersection with Unit Cube
void intersectUnitCubeByPlanes(const Files& a_output_files,
                               const int a_number_of_trials,
                               const int a_max_planes);

// Time intersection with Triangulated Triangular Prism
void intersectTriPrismByPlanes(const Files& a_output_files,
                               const int a_number_of_trials,
                               const int a_max_planes);

// Time intersection with Triangulated Hexahdron
void intersectTriHexByPlanes(const Files& a_output_files,
                             const int a_number_of_trials,
                             const int a_max_planes);

// Time intersection with Symmetric Triangular Prism
void intersectSymPrismByPlanes(const Files& a_output_files,
                               const int a_number_of_trials,
                               const int a_max_planes);

// Time intersection with Symmetric Hexahedron
void intersectSymHexByPlanes(const Files& a_output_files,
                             const int a_number_of_trials,
                             const int a_max_planes);

// Time intersection with Stellated Dodecahedron
void intersectStelDodecahedronByPlanes(const Files& a_output_files,
                                       const int a_number_of_trials,
                                       const int a_max_planes);

// Time intersection with Stellated Icosahedron
void intersectStelIcosahedronByPlanes(const Files& a_output_files,
                                      const int a_number_of_trials,
                                      const int a_max_planes);

#endif  // SRC_TIMING_COMP_INTERSECTION_TIMING_H_
