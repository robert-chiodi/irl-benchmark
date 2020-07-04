// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This code is simply a C interface wrapper to Fortran routines
// written for VOFTools. The VOFTools testing is performed in its
// native Fortran, but drive by the C++ driver functions in
// src/timing_comp/

// NOTE: To ease the use of VOFTools, a 'polyhedron' custom type
// was created, and wrapper functions written for the VOFTools
// functions. Care was taken to not introduce extra costs
// by doing this in order to still accurately gauge the
// performance of VOFTools.

#ifndef SRC_VOFTOOLS_VOFTOOLS_INTERSECTION_H_
#define SRC_VOFTOOLS_VOFTOOLS_INTERSECTION_H_

// Perform intersection of plane(s) with a Triangular Prism
void c_voftools_prismByPlanes(const double *a_prism_pts,
			    const int a_number_of_planes,
			    const double *a_planes, double *a_volume,
			    double *a_times);

// Perform intersection of plane(s) with a Unit Cube
void c_voftools_unitCubeByPlanes(const double *a_cube_pts,
				    const int a_number_of_planes,
				    const double *a_planes, double *a_volume,
				    double *a_times);

// Perform intersection of plane(s) with a Triangulated Triangular Prism
void c_voftools_triPrismByPlanes(const double *a_tri_hex_pts,
			       const int a_number_of_planes,
			       const double *a_planes, double *a_volume,
			       double *a_times);

// Perform intersection of plane(s) with a Triangulated Hexahedron
void c_voftools_triHexByPlanes(const double *a_tri_hex_pts,
			       const int a_number_of_planes,
			       const double *a_planes, double *a_volume,
			       double *a_times);

// Perform intersection of plane(s) with a Symmetric Triangular Prism
void c_voftools_symPrismByPlanes(const double *a_sym_prism_pts,
			       const int a_number_of_planes,
			       const double *a_planes, double *a_volume,
			       double *a_times);

// Perform intersection of plane(s) with a Symmetric Hexahedron
void c_voftools_symHexByPlanes(const double *a_sym_hex_pts,
			       const int a_number_of_planes,
			       const double *a_planes, double *a_volume,
			       double *a_times);

// Perform intersection of plane(s) with a Stellated Dodecahedron
void c_voftools_stelDodecahedronByPlanes(const double *a_stel_dodecahedron_pts,
					 const int a_number_of_planes,
					 const double *a_planes, double *a_volume,
					 double *a_times);

// Perform intersection of plane(s) with a Stellated Icosahedron
void c_voftools_stelIcosahedronByPlanes(const double *a_stel_icosahedron_pts,
					const int a_number_of_planes,
					const double *a_planes, double *a_volume,
					double *a_times);

// Same function as equivalent above, but only total timed for less overhead
void c_voftools_prismByPlanes_total(const double *a_prism_pts,
				    const int a_number_of_planes,
				    const double *a_planes, double *a_volume,
				    double *a_times);

// Same function as equivalent above, but only total timed for less overhead
void c_voftools_unitCubeByPlanes_total(const double *a_cube_pts,
				       const int a_number_of_planes,
				       const double *a_planes, double *a_volume,
				       double *a_times);

// Same function as equivalent above, but only total timed for less overhead
void c_voftools_triPrismByPlanes_total(const double *a_tri_hex_pts,
				       const int a_number_of_planes,
				       const double *a_planes, double *a_volume,
				       double *a_times);

// Same function as equivalent above, but only total timed for less overhead
void c_voftools_triHexByPlanes_total(const double *a_tri_hex_pts,
				     const int a_number_of_planes,
				     const double *a_planes, double *a_volume,
				     double *a_times);

// Same function as equivalent above, but only total timed for less overhead
void c_voftools_symPrismByPlanes_total(const double *a_sym_prism_pts,
				       const int a_number_of_planes,
				       const double *a_planes, double *a_volume,
				       double *a_times);

// Same function as equivalent above, but only total timed for less overhead
void c_voftools_symHexByPlanes_total(const double *a_sym_hex_pts,
				     const int a_number_of_planes,
				     const double *a_planes, double *a_volume,
				     double *a_times);

// Same function as equivalent above, but only total timed for less overhead
void c_voftools_stelDodecahedronByPlanes_total(const double *a_stel_dodecahedron_pts,
					       const int a_number_of_planes,
					       const double *a_planes, double *a_volume,
					       double *a_times);

// Same function as equivalent above, but only total timed for less overhead
void c_voftools_stelIcosahedronByPlanes_total(const double *a_stel_icosahedron_pts,
					      const int a_number_of_planes,
					      const double *a_planes, double *a_volume,
					      double *a_times);

#endif // SRC_VOFTOOLS_VOFTOOLS_INTERSECTION_H_
