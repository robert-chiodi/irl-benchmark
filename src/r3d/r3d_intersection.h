// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_R3D_R3D_INTERSECTION_H_
#define SRC_R3D_R3D_INTERSECTION_H_

// R3D Include file
#include "r3d.h"

// Initial prism in R3D format
void r3d_init_prism(r3d_poly *poly,  r3d_rvec3* verts);

// Perform intersection of plane(s) with a Triangular Prism
void r3d_prismByPlanes(const double *a_prism_pts,
		       const int a_number_of_planes,
		       const double *a_planes, double *a_volume,
		       double *a_times);

// Perform intersection of plane(s) with a Unit Cube
void r3d_unitCubeByPlanes(const double *a_cube_pts,
                          const int a_number_of_planes,
			  const double *a_planes, double *a_volume,
			  double *a_times);

// Perform intersection of plane(s) with a Triangulated Triangular Prism
void r3d_triPrismByPlanes(const double *a_tri_prism_pts,
			  const int a_number_of_planes,
			  const double *a_planes, double *a_volume,
			  double *a_times);

// Perform intersection of plane(s) with a Triangulated Hexahedron
void r3d_triHexByPlanes(const double *a_tri_hex_pts,
			const int a_number_of_planes,
			const double *a_planes, double *a_volume,
			double *a_times);

// Perform intersection of plane(s) with a Symmetric Triangular Prism
void r3d_symPrismByPlanes(const double *a_sym_prism_pts,			 
			  const int a_number_of_planes,
			  const double *a_planes, double *a_volume,
			  double *a_times);

// Perform intersection of plane(s) with a Symmetric Hexahedron
void r3d_symHexByPlanes(const double *a_sym_hex_pts,
			const int a_number_of_planes,
			const double *a_planes, double *a_volume,
			double *a_times);

// Perform intersection of plane(s) with a Stellated Dodecahedron
void r3d_stelDodecahedronByPlanes(const double *a_stel_dodecahedron_pts,
				  const int a_number_of_planes,
				  const double *a_planes, double *a_volume,
				  double *a_times);

// Perform intersection of plane(s) with a Stellated Icosahedron
void r3d_stelIcosahedronByPlanes(const double *a_stel_icosahedron_pts,
				 const int a_number_of_planes,
				 const double *a_planes, double *a_volume,
				 double *a_times);

// Same function as equivalent above, but only total timed for less overhead
void r3d_prismByPlanes_total(const double *a_prism_pts,
			     const int a_number_of_planes,
			     const double *a_planes, double *a_volume,
			     double *a_times);

// Same function as equivalent above, but only total timed for less overhead
void r3d_unitCubeByPlanes_total(const double *a_cube_pts,
				const int a_number_of_planes,
				const double *a_planes, double *a_volume,
				double *a_times);

// Same function as equivalent above, but only total timed for less overhead
void r3d_triPrismByPlanes_total(const double *a_tri_prism_pts,
				const int a_number_of_planes,
				const double *a_planes, double *a_volume,
				double *a_times);

// Same function as equivalent above, but only total timed for less overhead
void r3d_triHexByPlanes_total(const double *a_tri_hex_pts,
			      const int a_number_of_planes,
			      const double *a_planes, double *a_volume,
			      double *a_times);

// Same function as equivalent above, but only total timed for less overhead
void r3d_symPrismByPlanes_total(const double *a_sym_prism_pts,
				const int a_number_of_planes,
				const double *a_planes, double *a_volume,
				double *a_times);

// Same function as equivalent above, but only total timed for less overhead
void r3d_symHexByPlanes_total(const double *a_sym_hex_pts,
			      const int a_number_of_planes,
			      const double *a_planes, double *a_volume,
			      double *a_times);

// Same function as equivalent above, but only total timed for less overhead
void r3d_stelDodecahedronByPlanes_total(const double *a_stel_dodecahedron_pts,
					const int a_number_of_planes,
					const double *a_planes, double *a_volume,
					double *a_times);

// Same function as equivalent above, but only total timed for less overhead
void r3d_stelIcosahedronByPlanes_total(const double *a_stel_icosahedron_pts,
				       const int a_number_of_planes,
				       const double *a_planes, double *a_volume,
				       double *a_times);


#endif // SRC_R3D_R3D_INTERSECTION_H_
