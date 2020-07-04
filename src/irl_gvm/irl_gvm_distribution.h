// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This directory uses the idiomatic way of using IRL to 
// distribute a polyhedron over a mesh to get moments using the
// function getVolumeMoments. This function handles conversion from
// IRL polyhedron object to the half edge, performs the distribution,
// and returns the moments directly, tagged with individual cell IDs.

#ifndef SRC_IRL_GVM_IRL_GVM_DISTRIBUTION_H_
#define SRC_IRL_GVM_IRL_GVM_DISTRIBUTION_H_

#include <cstddef>

// Distribute a Cube onto a Cubic Mesh
void irl_gvm_cubeOntoCubicMesh(const double *a_cube_pts, double *a_volume,
                               double *a_times, std::size_t *a_entered_cells);

// Distribute a Stellated Icosahedron onto a Cubic Mesh
void irl_gvm_stelIcosahedronOntoCubicMesh(const double *a_stel_icosahedron_pts,
                                          double *a_volume, double *a_times,
                                          std::size_t *a_entered_cells);

// Distribute a Cube onto a Tet Mesh
void irl_gvm_cubeOntoTetMesh(const double *a_cube_pts, double *a_volume,
                             double *a_times, std::size_t *a_entered_cells);

// Distribute a Stellated Icosahedron onto a Tet Mesh
void irl_gvm_stelIcosahedronOntoTetMesh(const double *a_cube_pts,
                                        double *a_volume, double *a_times,
                                        std::size_t *a_entered_cells);

// Distribute a Cube onto a Spherical Cartesian Mesh
void irl_gvm_cubeOntoSphericalCartesianMesh(const double *a_cube_pts,
                                         double *a_volume, double *a_times,
                                         std::size_t *a_entered_cells);

// Distribute a Stellated Icosahedron onto a Spherical Cartesian Mesh
void irl_gvm_stelIcosahedronOntoSphericalCartesianMesh(
    const double *a_cube_pts, double *a_volume, double *a_times,
    std::size_t *a_entered_cells);

#endif // SRC_IRL_GVM_IRL_GVM_DISTRIBUTION_H_
