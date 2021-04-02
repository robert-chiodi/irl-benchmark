// Copyright (C) 2020 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_TIMING_COMP_DISTRIBUTION_TIMING_H_
#define SRC_TIMING_COMP_DISTRIBUTION_TIMING_H_

#include "src/timing_comp/files.h"

// Distributes Unit Cube onto Cubic Mesh
void distributeCubeOntoCubicMesh(FILE* a_distribute_file,
                                 const int a_number_of_trials);

// Distributes Stellated Icosahedron onto Cubic Mesh
void distributeStelIcosahedronOntoCubicMesh(FILE* a_distribute_file,
                                            const int a_number_of_trials);

// Distriutes Unit Cube onto Tet Mesh
void distributeCubeOntoTetMesh(FILE* a_distribute_file,
                               const int a_number_of_trials);

// Distributes Stellated Icosahedron onto Tet Mesh
void distributeStelIcosahedronOntoTetMesh(FILE* a_distribute_file,
                                          const int a_number_of_trials);

// Distributes Unit Cube onto Spherical Cartesian Mesh
void distributeCubeOntoSphericalCartesianMesh(FILE* a_distribute_file,
                                              const int a_number_of_trials);

// Distributes Stellated Icosahedron onto Spherical Cartesian Mesh
void distributeStelIcosahedronOntoSphericalCartesianMesh(
    FILE* a_distribute_file, const int a_number_of_trials);

#endif  // SRC_TIMING_COMP_DISTRIBUTION_TIMING_H_
