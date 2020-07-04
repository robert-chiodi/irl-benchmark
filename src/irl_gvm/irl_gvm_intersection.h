// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This directory uses the more idiomatic way of using IRL to obtain
// volume moments from volume-plane intersections using the
// function getVolumeMoments. This function handles conversion from
// IRL polyhedron object to the half edge, performs the truncation,
// and returns the moments directly. This is more efficient, but
// prevents the timing of the individual parts (initialization, intersection,
// and volume calculations) as done in the src/irl directory. This is why
// the two separate versions are provided.

#ifndef SRC_IRL_GVM_IRL_GVM_INTERSECTION_H_
#define SRC_IRL_GVM_IRL_GVM_INTERSECTION_H_

#include <cstddef>

// Perform intersection of plane(s) with a Triangular Prism
void irl_gvm_prismByPlanes(const double *a_prism_pts,
                           const int a_number_of_planes, const double *a_planes,
                           double *a_volume, double *a_times);

// Perform intersection of plane(s) with a Unit Cube
void irl_gvm_unitCubeByPlanes(const double *a_cube_pts,
                              const int a_number_of_planes,
                              const double *a_planes, double *a_volume,
                              double *a_times);

// Perform intersection of plane(s) with a Triangulated Triangular Prism
void irl_gvm_triPrismByPlanes(const double *a_tri_prism_pts,
                              const int a_number_of_planes,
                              const double *a_planes, double *a_volume,
                              double *a_times);

// Perform intersection of plane(s) with a Triangulated Hexahedron
void irl_gvm_triHexByPlanes(const double *a_tri_hex_pts,
                            const int a_number_of_planes,
                            const double *a_planes, double *a_volume,
                            double *a_times);

// Perform intersection of plane(s) with a Symmetric Triangular Prism
void irl_gvm_symPrismByPlanes(const double *a_sym_prism_pts,
                              const int a_number_of_planes,
                              const double *a_planes, double *a_volume,
                              double *a_times);

// Perform intersection of plane(s) with a Symmetric Hexahedron
void irl_gvm_symHexByPlanes(const double *a_sym_hex_pts,
                            const int a_number_of_planes,
                            const double *a_planes, double *a_volume,
                            double *a_times);

// Perform intersection of plane(s) with a Stellated Dodecahedron
void irl_gvm_stelDodecahedronByPlanes(const double *a_stel_dodecahedron_pts,
                                      const int a_number_of_planes,
                                      const double *a_planes, double *a_volume,
                                      double *a_times);

// Perform intersection of plane(s) with a Stellated Icosahedron
void irl_gvm_stelIcosahedronByPlanes(const double *a_stel_icosahedron_pts,
                                     const int a_number_of_planes,
                                     const double *a_planes, double *a_volume,
                                     double *a_times);

#endif // SRC_IRL_GVM_IRL_GVM_INTERSECTION_H_
