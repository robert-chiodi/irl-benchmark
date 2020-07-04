// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This file (and the function definitions in irl_intersection.cpp)
// covers the plane-polyhedron intersection testing for IRL
// when individually timing the initialization, intersection,
// and volume computation.
// NOTE: This is an atypical way to use IRL. Normally, all are performed
// in the single function getVolumeMoments<MomentsType>(polytope,
// reconstruction). For examples of this, and the code that generated the
// "Total" time entries for IRL, please see
// src/irl_gvm/irl_gvm_intersection.h/.cpp

#ifndef SRC_IRL_IRL_INTERSECTION_H_
#define SRC_IRL_IRL_INTERSECTION_H_

// Perform intersection of plane(s) with a Triangular Prism
void irl_prismByPlanes(const double* a_prism_pts, const int a_number_of_planes,
                       const double* a_planes, double* a_volume,
                       double* a_times);

// Perform intersection of plane(s) with a Unit Cube
void irl_unitCubeByPlanes(const double* a_cube_pts,
                          const int a_number_of_planes, const double* a_planes,
                          double* a_volume, double* a_times);

// Perform intersection of plane(s) with a Triangulated Triangular Prism
void irl_triPrismByPlanes(const double* a_tri_prism_pts,
                          const int a_number_of_planes, const double* a_planes,
                          double* a_volume, double* a_times);

// Perform intersection of plane(s) with a Triangulated Hexahedron
void irl_triHexByPlanes(const double* a_tri_hex_pts,
                        const int a_number_of_planes, const double* a_planes,
                        double* a_volume, double* a_times);

// Perform intersection of plane(s) with a Symmetric Triangular Prism
void irl_symPrismByPlanes(const double* a_sym_prism_pts,
                          const int a_number_of_planes, const double* a_planes,
                          double* a_volume, double* a_times);

// Perform intersection of plane(s) with a Symmetric Hexahedron
void irl_symHexByPlanes(const double* a_sym_hex_pts,
                        const int a_number_of_planes, const double* a_planes,
                        double* a_volume, double* a_times);

// Perform intersection of plane(s) with a Stellated Dodecahedron
void irl_stelDodecahedronByPlanes(const double* a_stel_dodecahedron_pts,
                                  const int a_number_of_planes,
                                  const double* a_planes, double* a_volume,
                                  double* a_times);

// Perform intersection of plane(s) with a Stellated Icosahedron
void irl_stelIcosahedronByPlanes(const double* a_stel_icosahedron_pts,
                                 const int a_number_of_planes,
                                 const double* a_planes, double* a_volume,
                                 double* a_times);

#endif  // SRC_IRL_IRL_INTERSECTION_H_
