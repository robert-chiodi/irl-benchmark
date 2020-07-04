// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/r3d/r3d_intersection.h"

#include <omp.h>

// R3D Include file
#include "r3d.h"

void r3d_init_prism(r3d_poly *poly,  r3d_rvec3* verts) {
	// direct access to vertex buffer
	r3d_vertex *vertbuffer = poly->verts;
	r3d_int *nverts = &poly->nverts;
	*nverts = 6;
	vertbuffer[0].pnbrs[0] = 1;
	vertbuffer[0].pnbrs[1] = 2;
	vertbuffer[0].pnbrs[2] = 3;
	vertbuffer[1].pnbrs[0] = 2;
	vertbuffer[1].pnbrs[1] = 0;
	vertbuffer[1].pnbrs[2] = 4;
	vertbuffer[2].pnbrs[0] = 0;
	vertbuffer[2].pnbrs[1] = 1;
	vertbuffer[2].pnbrs[2] = 5;
	vertbuffer[3].pnbrs[0] = 4;
	vertbuffer[3].pnbrs[1] = 0;
	vertbuffer[3].pnbrs[2] = 5;
	vertbuffer[4].pnbrs[0] = 1;
	vertbuffer[4].pnbrs[1] = 3;
	vertbuffer[4].pnbrs[2] = 5;
	vertbuffer[5].pnbrs[0] = 3;
	vertbuffer[5].pnbrs[1] = 2;
	vertbuffer[5].pnbrs[2] = 4;
	for (r3d_int v = 0; v < 6; ++v){
	  vertbuffer[v].pos = verts[v];
	}
}


void r3d_prismByPlanes(const double *a_prism_pts,
		       const int a_number_of_planes,
		       const double *a_planes, double *a_volume,
		       double *a_times) {

  // Full BREP for Triangular Prism
  r3d_int nvert = 6;
  r3d_int nface = 5;
  
  r3d_int face_flat[18] = {
			   0, 1, 2,
			   0, 3, 4, 1,
			   0, 2, 5, 3,
			   1, 4, 5, 2,
			   3, 5, 4
  };
  r3d_int verts_per_face[5] = {
			       3,4,4,4,3
  };
  
  r3d_rvec3 verts[nvert];
  for(r3d_int p = 0; p < nvert; ++p){
    for(r3d_int d = 0; d < 3; ++d){
      verts[p].xyz[d] = a_prism_pts[p*3+d];
    }
  }
  
  r3d_int **faces = malloc(nface*sizeof(r3d_int*)); 
  for (r3d_int f = 0; f < nface; ++f){
    faces[f] = malloc(verts_per_face[f]*sizeof(r3d_int));
  }
  r3d_int ind = 0;
  for(r3d_int f = 0; f < nface; ++f){
    for(r3d_int v = 0; v < verts_per_face[f]; ++v){
      faces[f][v] = face_flat[ind];
      ++ind;
    }
  }  

  // Multiply normal by -1.0 due to difference in convention.
  // Makes consistent with IRL convention     
  r3d_plane planes[a_number_of_planes];
  for(r3d_int n = 0; n < a_number_of_planes; ++n){
    planes[n].n.x = -a_planes[n*4+0];
    planes[n].n.y = -a_planes[n*4+1];
    planes[n].n.z = -a_planes[n*4+2];
    planes[n].d = a_planes[n*4+3];    
  }

  double start = omp_get_wtime();  
  r3d_poly poly;
  r3d_init_prism(&poly, verts);
  double end = omp_get_wtime();
  a_times[0] = end - start;  

  start = omp_get_wtime();    
  r3d_clip(&poly, planes, a_number_of_planes);
  end = omp_get_wtime();
  a_times[1] = end-start;

  start = omp_get_wtime();
  r3d_reduce(&poly, a_volume, 0);  
  end = omp_get_wtime();
  a_times[2] = end-start;
}

void r3d_unitCubeByPlanes(const double *a_cube_pts,
			  const int a_number_of_planes,
			  const double *a_planes, double *a_volume,
			  double *a_times) {
  r3d_rvec3 pts[2];
  pts[0].x = a_cube_pts[0];
  pts[0].y = a_cube_pts[1];
  pts[0].z = a_cube_pts[2];
  pts[1].x = a_cube_pts[3];
  pts[1].y = a_cube_pts[4];
  pts[1].z = a_cube_pts[5];  

  // Multiply normal by -1.0 due to difference in convention.
  // Makes consistent with IRL convention       
  r3d_plane planes[a_number_of_planes];
  for(r3d_int n = 0; n < a_number_of_planes; ++n){
    planes[n].n.x = -a_planes[n*4+0];
    planes[n].n.y = -a_planes[n*4+1];
    planes[n].n.z = -a_planes[n*4+2];
    planes[n].d = a_planes[n*4+3];    
  }

  double start = omp_get_wtime();  
  r3d_poly poly;  
  r3d_init_box(&poly, pts);
  double end = omp_get_wtime();
  a_times[0] = end - start;  

  start = omp_get_wtime();    
  r3d_clip(&poly, planes, a_number_of_planes);
  end = omp_get_wtime();
  a_times[1] = end-start;

  start = omp_get_wtime();
  r3d_reduce(&poly, a_volume, 0);  
  end = omp_get_wtime();
  a_times[2] = end-start;
}

void r3d_triPrismByPlanes(const double *a_tri_prism_pts,
			const int a_number_of_planes,
			const double *a_planes, double *a_volume,
			double *a_times) {
  // Full BREP for Triangulated Triangular Prism
  r3d_int nvert = 6;
  r3d_int nface = 8;
  
  r3d_int face_flat[24] = {
			   0, 1, 2,
			   4, 3, 5,
			   4, 5, 2,
			   4, 2, 1,
			   4, 1, 0,
			   4, 0, 3,
			   0, 2, 5,
			   0, 5, 3			 
  };
  r3d_int verts_per_face[8] = {
				3,3,3,3,3,3,
				3,3
  };
  
  r3d_rvec3 verts[nvert];
  for(r3d_int p = 0; p < nvert; ++p){
    for(r3d_int d = 0; d < 3; ++d){
      verts[p].xyz[d] = a_tri_prism_pts[p*3+d];
    }
  }
  
  r3d_int **faces = malloc(nface*sizeof(r3d_int*)); 
  for (r3d_int f = 0; f < nface; ++f){
    faces[f] = malloc(verts_per_face[f]*sizeof(r3d_int));
  }
  r3d_int ind = 0;
  for(r3d_int f = 0; f < nface; ++f){
    for(r3d_int v = 0; v < verts_per_face[f]; ++v){
      faces[f][v] = face_flat[ind];
      ++ind;
    }
  }
  
  
  // Multiply normal by -1.0 due to difference in convention.
  // Makes consistent with IRL convention     
  r3d_plane planes[a_number_of_planes];
  for(r3d_int n = 0; n < a_number_of_planes; ++n){
    planes[n].n.x = -a_planes[n*4+0];
    planes[n].n.y = -a_planes[n*4+1];
    planes[n].n.z = -a_planes[n*4+2];
    planes[n].d = a_planes[n*4+3];    
  }

  double start = omp_get_wtime();  
  r3d_poly poly;
  r3d_init_poly(&poly, verts, nvert, faces, verts_per_face, nface);
  double end = omp_get_wtime();
  a_times[0] = end - start;  

  start = omp_get_wtime();    
  r3d_clip(&poly, planes, a_number_of_planes);
  end = omp_get_wtime();
  a_times[1] = end-start;

  start = omp_get_wtime();
  r3d_reduce(&poly, a_volume, 0);  
  end = omp_get_wtime();
  a_times[2] = end-start;

  for (r3d_int f = 0; f < nface; ++f){
    free(faces[f]);
  }
  free(faces);
}


void r3d_triHexByPlanes(const double *a_tri_hex_pts,
			const int a_number_of_planes,
			const double *a_planes, double *a_volume,
			double *a_times) {
  // Full BREP for Triangulated Hexahedron
  r3d_int nvert = 8;
  r3d_int nface = 12;

  
  r3d_int face_flat[36] = {
			   5, 7, 6,
			   5, 4, 7,
			   3, 0, 1,
			   3, 1, 2,
			   4, 3, 7,
			   4, 0, 3,
			   2, 5, 6,
			   2, 1, 5,
			   0, 5, 1,
			   0, 4, 5,
			   3, 6, 7,
			   3, 2, 6
  };
  r3d_int verts_per_face[12] = {
				3,3,3,3,3,3,
				3,3,3,3,3,3
  };
  
  r3d_rvec3 verts[nvert];
  for(r3d_int p = 0; p < nvert; ++p){
    for(r3d_int d = 0; d < 3; ++d){
      verts[p].xyz[d] = a_tri_hex_pts[p*3+d];
    }
  }
  
  r3d_int **faces = malloc(nface*sizeof(r3d_int*)); 
  for (r3d_int f = 0; f < nface; ++f){
    faces[f] = malloc(verts_per_face[f]*sizeof(r3d_int));
  }
  r3d_int ind = 0;
  for(r3d_int f = 0; f < nface; ++f){
    for(r3d_int v = 0; v < verts_per_face[f]; ++v){
      faces[f][v] = face_flat[ind];
      ++ind;
    }
  }
  
  
  // Multiply normal by -1.0 due to difference in convention.
  // Makes consistent with IRL convention     
  r3d_plane planes[a_number_of_planes];
  for(r3d_int n = 0; n < a_number_of_planes; ++n){
    planes[n].n.x = -a_planes[n*4+0];
    planes[n].n.y = -a_planes[n*4+1];
    planes[n].n.z = -a_planes[n*4+2];
    planes[n].d = a_planes[n*4+3];    
  }

  double start = omp_get_wtime();  
  r3d_poly poly;
  r3d_init_poly(&poly, verts, nvert, faces, verts_per_face, nface);
  double end = omp_get_wtime();
  a_times[0] = end - start;  

  start = omp_get_wtime();    
  r3d_clip(&poly, planes, a_number_of_planes);
  end = omp_get_wtime();
  a_times[1] = end-start;

  start = omp_get_wtime();
  r3d_reduce(&poly, a_volume, 0);  
  end = omp_get_wtime();
  a_times[2] = end-start;

  for (r3d_int f = 0; f < nface; ++f){
    free(faces[f]);
  }
  free(faces);  
}

void r3d_symPrismByPlanes(const double *a_sym_prism_pts,
			  const int a_number_of_planes,
			  const double *a_planes, double *a_volume,
			  double *a_times) {
  // Full BREP for Symmetric PRism
  r3d_int nvert = 11;
  r3d_int nface = 18;

  r3d_int face_flat[54] = {
			   6, 0, 1,
			   6, 1, 2,
			   6, 2, 0,
			   7, 1, 0,
			   7, 4, 1,
			   7, 3, 4,
			   7, 0, 3,
			   8, 2, 1,
			   8, 1, 4,
			   8, 4, 5,
			   8, 5, 2,
			   9, 0, 2,
			   9, 3, 0,
			   9, 5, 3,
			   9, 2, 5,
			   10, 4, 3,
			   10, 3, 5,
			   10, 5, 4
  };
  r3d_int verts_per_face[18] = {
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3							
  };
  
  r3d_rvec3 verts[nvert];
  for(r3d_int p = 0; p < nvert; ++p){
    for(r3d_int d = 0; d < 3; ++d){
      verts[p].xyz[d] = a_sym_prism_pts[p*3+d];
    }
  }
  
  r3d_int **faces = malloc(nface*sizeof(r3d_int*)); 
  for (r3d_int f = 0; f < nface; ++f){
    faces[f] = malloc(verts_per_face[f]*sizeof(r3d_int));
  }
  r3d_int ind = 0;
  for(r3d_int f = 0; f < nface; ++f){
    for(r3d_int v = 0; v < verts_per_face[f]; ++v){
      faces[f][v] = face_flat[ind];
      ++ind;
    }
  }
  
  
  // Multiply normal by -1.0 due to difference in convention.
  // Makes consistent with IRL convention     
  r3d_plane planes[a_number_of_planes];
  for(r3d_int n = 0; n < a_number_of_planes; ++n){
    planes[n].n.x = -a_planes[n*4+0];
    planes[n].n.y = -a_planes[n*4+1];
    planes[n].n.z = -a_planes[n*4+2];
    planes[n].d = a_planes[n*4+3];    
  }

  double start = omp_get_wtime();  
  r3d_poly poly;
  r3d_init_poly(&poly, verts, nvert, faces, verts_per_face, nface);
  double end = omp_get_wtime();
  a_times[0] = end - start;  

  start = omp_get_wtime();    
  r3d_clip(&poly, planes, a_number_of_planes);
  end = omp_get_wtime();
  a_times[1] = end-start;

  start = omp_get_wtime();
  r3d_reduce(&poly, a_volume, 0);  
  end = omp_get_wtime();
  a_times[2] = end-start;

  for (r3d_int f = 0; f < nface; ++f){
    free(faces[f]);
  }
  free(faces);  
}


void r3d_symHexByPlanes(const double *a_sym_hex_pts,
			const int a_number_of_planes,
			const double *a_planes, double *a_volume,
			double *a_times) {
  // Full BREP for Symmetric Hex
  r3d_int nvert = 14;
  r3d_int nface = 24;
  
  r3d_int face_flat[72] = {
			   0, 1, 8,
			   1, 2, 8,
			   2, 3, 8,
			   3, 0, 8,
			   5, 1, 9,
			   1, 0, 9,
			   0, 4, 9,
			   4, 5, 9,
			   1, 5, 10,
			   5, 6, 10,
			   6, 2, 10,
			   2, 1, 10,
			   2, 6, 11,
			   6, 7, 11,
			   7, 3, 11,
			   3, 2, 11,
			   0, 3, 12,
			   3, 7, 12,
			   7, 4, 12,
			   4, 0, 12,
			   5, 4, 13,
			   4, 7, 13,
			   7, 6, 13,
			   6, 5, 13
  };
  r3d_int verts_per_face[24] = {
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3   
  };
  
  r3d_rvec3 verts[nvert];
  for(r3d_int p = 0; p < nvert; ++p){
    for(r3d_int d = 0; d < 3; ++d){
      verts[p].xyz[d] = a_sym_hex_pts[p*3+d];
    }
  }
  
  r3d_int **faces = malloc(nface*sizeof(r3d_int*)); 
  for (r3d_int f = 0; f < nface; ++f){
    faces[f] = malloc(verts_per_face[f]*sizeof(r3d_int));
  }
  r3d_int ind = 0;
  for(r3d_int f = 0; f < nface; ++f){
    for(r3d_int v = 0; v < verts_per_face[f]; ++v){
      faces[f][v] = face_flat[ind];
      ++ind;
    }
  }
  
  
  // Multiply normal by -1.0 due to difference in convention.
  // Makes consistent with IRL convention     
  r3d_plane planes[a_number_of_planes];
  for(r3d_int n = 0; n < a_number_of_planes; ++n){
    planes[n].n.x = -a_planes[n*4+0];
    planes[n].n.y = -a_planes[n*4+1];
    planes[n].n.z = -a_planes[n*4+2];
    planes[n].d = a_planes[n*4+3];    
  }

  double start = omp_get_wtime();  
  r3d_poly poly;
  r3d_init_poly(&poly, verts, nvert, faces, verts_per_face, nface);
  double end = omp_get_wtime();
  a_times[0] = end - start;  

  start = omp_get_wtime();    
  r3d_clip(&poly, planes, a_number_of_planes);
  end = omp_get_wtime();
  a_times[1] = end-start;

  start = omp_get_wtime();
  r3d_reduce(&poly, a_volume, 0);  
  end = omp_get_wtime();
  a_times[2] = end-start;

  for (r3d_int f = 0; f < nface; ++f){
    free(faces[f]);
  }
  free(faces);  
}

void r3d_stelDodecahedronByPlanes(const double *a_stel_dodecahedron_pts,
				  const int a_number_of_planes,
				  const double *a_planes, double *a_volume,
				  double *a_times) {
  // Full BREP for Stellated Dodecahedron
  r3d_int nvert = 32;
  r3d_int nface = 60;
  
  r3d_int face_flat[180] = {
			    0, 8, 20,
			    8, 9, 20,
			    9, 4, 20,
			    4, 16, 20,
			    16, 0, 20,
			    0, 16, 21,
			    16, 17, 21,
			    17, 2, 21,
			    2, 12, 21,
			    12, 0, 21,
			    12, 2, 22,
			    2, 10, 22,
			    10, 3, 22,
			    3, 13, 22,
			    13, 12, 22,
			    9, 5, 23,
			    5, 15, 23,
			    15, 14, 23,
			    14, 4, 23,
			    4, 9, 23,
			    3, 19, 24,
			    19, 18, 24,
			    18, 1, 24,
			    1, 13, 24,
			    13, 3, 24,
			    7, 11, 25,
			    11, 6, 25,
			    6, 14, 25,
			    14, 15, 25,
			    15, 7, 25,
			    0, 12, 26,
			    12, 13, 26,
			    13, 1, 26,
			    1, 8, 26,
			    8, 0, 26,
			    8, 1, 27,
			    1, 18, 27,
			    18, 5, 27,
			    5, 9, 27,
			    9, 8, 27,
			    16, 4, 28,
			    4, 14, 28,
			    14, 6, 28,
			    6, 17, 28,
			    17, 16, 28,
			    6, 11, 29,
			    11, 10, 29,
			    10, 2, 29,
			    2, 17, 29,
			    17, 6, 29,
			    7, 15, 30,
			    15, 5, 30,
			    5, 18, 30,
			    18, 19, 30,
			    19, 7, 30,
			    7, 19, 31,
			    19, 3, 31,
			    3, 10, 31,
			    10, 11, 31,
			    11, 7, 31
  };
  r3d_int verts_per_face[60] = {
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3				
  };
  
  r3d_rvec3 verts[nvert];
  for(r3d_int p = 0; p < nvert; ++p){
    for(r3d_int d = 0; d < 3; ++d){
      verts[p].xyz[d] = a_stel_dodecahedron_pts[p*3+d];
    }
  }
  
  r3d_int **faces = malloc(nface*sizeof(r3d_int*)); 
  for (r3d_int f = 0; f < nface; ++f){
    faces[f] = malloc(verts_per_face[f]*sizeof(r3d_int));
  }
  r3d_int ind = 0;
  for(r3d_int f = 0; f < nface; ++f){
    for(r3d_int v = 0; v < verts_per_face[f]; ++v){
      faces[f][v] = face_flat[ind];
      ++ind;
    }
  }
  
  
  // Multiply normal by -1.0 due to difference in convention.
  // Makes consistent with IRL convention     
  r3d_plane planes[a_number_of_planes];
  for(r3d_int n = 0; n < a_number_of_planes; ++n){
    planes[n].n.x = -a_planes[n*4+0];
    planes[n].n.y = -a_planes[n*4+1];
    planes[n].n.z = -a_planes[n*4+2];
    planes[n].d = a_planes[n*4+3];    
  }

  double start = omp_get_wtime();  
  r3d_poly poly;
  r3d_init_poly(&poly, verts, nvert, faces, verts_per_face, nface);
  double end = omp_get_wtime();
  a_times[0] = end - start;  

  start = omp_get_wtime();    
  r3d_clip(&poly, planes, a_number_of_planes);
  end = omp_get_wtime();
  a_times[1] = end-start;

  start = omp_get_wtime();
  r3d_reduce(&poly, a_volume, 0);  
  end = omp_get_wtime();
  a_times[2] = end-start;

  for (r3d_int f = 0; f < nface; ++f){
    free(faces[f]);
  }
  free(faces);  
}

void r3d_stelIcosahedronByPlanes(const double *a_stel_icosahedron_pts,
				 const int a_number_of_planes,
				 const double *a_planes, double *a_volume,
				 double *a_times) {
  // Full BREP for Stellated Icosahedron
  r3d_int nvert = 32;
  r3d_int nface = 60;
  
  r3d_int face_flat[180] = {
			    0, 8, 12,
			    8, 4, 12,
			    4, 0, 12,
			    0, 5, 13,
			    5, 10, 13,
			    10, 0, 13,
			    2, 4, 14,
			    4, 9, 14,
			    9, 2, 14,
			    2, 11, 15,
			    11, 5, 15,
			    5, 2, 15,
			    1, 6, 16,
			    6, 8, 16,
			    8, 1, 16,
			    1, 10, 17,
			    10, 7, 17,
			    7, 1, 17,
			    3, 9, 18,
			    9, 6, 18,
			    6, 3, 18,
			    3, 7, 19,
			    7, 11, 19,
			    11, 3, 19,
			    0, 10, 20,
			    10, 8, 20,
			    8, 0, 20,
			    1, 8, 21,
			    8, 10, 21,
			    10, 1, 21,
			    2, 9, 22,
			    9, 11, 22,
			    11, 2, 22,
			    3, 11, 23,
			    11, 9, 23,
			    9, 3, 23,
			    4, 2, 24,
			    2, 0, 24,
			    0, 4, 24,
			    5, 0, 25,
			    0, 2, 25,
			    2, 5, 25,
			    6, 1, 26,
			    1, 3, 26,
			    3, 6, 26,
			    7, 3, 27,
			    3, 1, 27,
			    1, 7, 27,
			    8, 6, 28,
			    6, 4, 28,
			    4, 8, 28,
			    9, 4, 29,
			    4, 6, 29,
			    6, 9, 29,
			    10, 5, 30,
			    5, 7, 30,
			    7, 10, 30,
			    11, 7, 31,
			    7, 5, 31,
			    5, 11, 31
  };
  r3d_int verts_per_face[60] = {
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3				
  };
  
  r3d_rvec3 verts[nvert];
  for(r3d_int p = 0; p < nvert; ++p){
    for(r3d_int d = 0; d < 3; ++d){
      verts[p].xyz[d] = a_stel_icosahedron_pts[p*3+d];
    }
  }
  
  r3d_int **faces = malloc(nface*sizeof(r3d_int*)); 
  for (r3d_int f = 0; f < nface; ++f){
    faces[f] = malloc(verts_per_face[f]*sizeof(r3d_int));
  }
  r3d_int ind = 0;
  for(r3d_int f = 0; f < nface; ++f){
    for(r3d_int v = 0; v < verts_per_face[f]; ++v){
      faces[f][v] = face_flat[ind];
      ++ind;
    }
  }
  
  
  // Multiply normal by -1.0 due to difference in convention.
  // Makes consistent with IRL convention     
  r3d_plane planes[a_number_of_planes];
  for(r3d_int n = 0; n < a_number_of_planes; ++n){
    planes[n].n.x = -a_planes[n*4+0];
    planes[n].n.y = -a_planes[n*4+1];
    planes[n].n.z = -a_planes[n*4+2];
    planes[n].d = a_planes[n*4+3];    
  }

  double start = omp_get_wtime();  
  r3d_poly poly;
  r3d_init_poly(&poly, verts, nvert, faces, verts_per_face, nface);
  double end = omp_get_wtime();
  a_times[0] = end - start;  

  start = omp_get_wtime();    
  r3d_clip(&poly, planes, a_number_of_planes);
  end = omp_get_wtime();
  a_times[1] = end-start;

  start = omp_get_wtime();
  r3d_reduce(&poly, a_volume, 0);  
  end = omp_get_wtime();
  a_times[2] = end-start;

  for (r3d_int f = 0; f < nface; ++f){
    free(faces[f]);
  }
  free(faces);  
}



///// Implementation of same functions from above but timing everything at once /////

void r3d_prismByPlanes_total(const double *a_prism_pts,
			     const int a_number_of_planes,
			     const double *a_planes, double *a_volume,
			     double *a_times) {

  // Full BREP for Triangular Prism
  r3d_int nvert = 6;
  r3d_int nface = 5;
  
  r3d_int face_flat[18] = {
			   0, 1, 2,
			   0, 3, 4, 1,
			   0, 2, 5, 3,
			   1, 4, 5, 2,
			   3, 5, 4
  };
  r3d_int verts_per_face[5] = {
			       3,4,4,4,3
  };
  
  r3d_rvec3 verts[nvert];
  for(r3d_int p = 0; p < nvert; ++p){
    for(r3d_int d = 0; d < 3; ++d){
      verts[p].xyz[d] = a_prism_pts[p*3+d];
    }
  }
  
  r3d_int **faces = malloc(nface*sizeof(r3d_int*)); 
  for (r3d_int f = 0; f < nface; ++f){
    faces[f] = malloc(verts_per_face[f]*sizeof(r3d_int));
  }
  r3d_int ind = 0;
  for(r3d_int f = 0; f < nface; ++f){
    for(r3d_int v = 0; v < verts_per_face[f]; ++v){
      faces[f][v] = face_flat[ind];
      ++ind;
    }
  }  

  // Multiply normal by -1.0 due to difference in convention.
  // Makes consistent with IRL convention     
  r3d_plane planes[a_number_of_planes];
  for(r3d_int n = 0; n < a_number_of_planes; ++n){
    planes[n].n.x = -a_planes[n*4+0];
    planes[n].n.y = -a_planes[n*4+1];
    planes[n].n.z = -a_planes[n*4+2];
    planes[n].d = a_planes[n*4+3];    
  }

  double start = omp_get_wtime();  
  r3d_poly poly;
  r3d_init_prism(&poly, verts);
  r3d_clip(&poly, planes, a_number_of_planes);
  r3d_reduce(&poly, a_volume, 0);  
  double end = omp_get_wtime();
  a_times[0] = end - start;  
}

void r3d_unitCubeByPlanes_total(const double *a_cube_pts,
				const int a_number_of_planes,
				const double *a_planes, double *a_volume,
				double *a_times) {
  r3d_rvec3 pts[2];
  pts[0].x = a_cube_pts[0];
  pts[0].y = a_cube_pts[1];
  pts[0].z = a_cube_pts[2];
  pts[1].x = a_cube_pts[3];
  pts[1].y = a_cube_pts[4];
  pts[1].z = a_cube_pts[5];  

  // Multiply normal by -1.0 due to difference in convention.
  // Makes consistent with IRL convention       
  r3d_plane planes[a_number_of_planes];
  for(r3d_int n = 0; n < a_number_of_planes; ++n){
    planes[n].n.x = -a_planes[n*4+0];
    planes[n].n.y = -a_planes[n*4+1];
    planes[n].n.z = -a_planes[n*4+2];
    planes[n].d = a_planes[n*4+3];    
  }

  double start = omp_get_wtime();  
  r3d_poly poly;  
  r3d_init_box(&poly, pts);   
  r3d_clip(&poly, planes, a_number_of_planes);
  r3d_reduce(&poly, a_volume, 0);
  double end = omp_get_wtime();
  a_times[0] = end - start;    
}

void r3d_triPrismByPlanes_total(const double *a_tri_prism_pts,
				const int a_number_of_planes,
				const double *a_planes, double *a_volume,
				double *a_times) {
  // Full BREP for Triangulated Triangular Prism
  r3d_int nvert = 6;
  r3d_int nface = 8;
  
  r3d_int face_flat[24] = {
			   0, 1, 2,
			   4, 3, 5,
			   4, 5, 2,
			   4, 2, 1,
			   4, 1, 0,
			   4, 0, 3,
			   0, 2, 5,
			   0, 5, 3			 
  };
  r3d_int verts_per_face[8] = {
				3,3,3,3,3,3,
				3,3
  };
  
  r3d_rvec3 verts[nvert];
  for(r3d_int p = 0; p < nvert; ++p){
    for(r3d_int d = 0; d < 3; ++d){
      verts[p].xyz[d] = a_tri_prism_pts[p*3+d];
    }
  }
  
  r3d_int **faces = malloc(nface*sizeof(r3d_int*)); 
  for (r3d_int f = 0; f < nface; ++f){
    faces[f] = malloc(verts_per_face[f]*sizeof(r3d_int));
  }
  r3d_int ind = 0;
  for(r3d_int f = 0; f < nface; ++f){
    for(r3d_int v = 0; v < verts_per_face[f]; ++v){
      faces[f][v] = face_flat[ind];
      ++ind;
    }
  }
  
  
  // Multiply normal by -1.0 due to difference in convention.
  // Makes consistent with IRL convention     
  r3d_plane planes[a_number_of_planes];
  for(r3d_int n = 0; n < a_number_of_planes; ++n){
    planes[n].n.x = -a_planes[n*4+0];
    planes[n].n.y = -a_planes[n*4+1];
    planes[n].n.z = -a_planes[n*4+2];
    planes[n].d = a_planes[n*4+3];    
  }

  double start = omp_get_wtime();  
  r3d_poly poly;
  r3d_init_poly(&poly, verts, nvert, faces, verts_per_face, nface);
  r3d_clip(&poly, planes, a_number_of_planes);
  r3d_reduce(&poly, a_volume, 0);  
  double end = omp_get_wtime();
  a_times[0] = end - start;  

  for (r3d_int f = 0; f < nface; ++f){
    free(faces[f]);
  }
  free(faces);
}


void r3d_triHexByPlanes_total(const double *a_tri_hex_pts,
			      const int a_number_of_planes,
			      const double *a_planes, double *a_volume,
			      double *a_times) {
  // Full BREP for Triangulated Hexahedron
  r3d_int nvert = 8;
  r3d_int nface = 12;

  
  r3d_int face_flat[36] = {
			   5, 7, 6,
			   5, 4, 7,
			   3, 0, 1,
			   3, 1, 2,
			   4, 3, 7,
			   4, 0, 3,
			   2, 5, 6,
			   2, 1, 5,
			   0, 5, 1,
			   0, 4, 5,
			   3, 6, 7,
			   3, 2, 6
  };
  r3d_int verts_per_face[12] = {
				3,3,3,3,3,3,
				3,3,3,3,3,3
  };
  
  r3d_rvec3 verts[nvert];
  for(r3d_int p = 0; p < nvert; ++p){
    for(r3d_int d = 0; d < 3; ++d){
      verts[p].xyz[d] = a_tri_hex_pts[p*3+d];
    }
  }
  
  r3d_int **faces = malloc(nface*sizeof(r3d_int*)); 
  for (r3d_int f = 0; f < nface; ++f){
    faces[f] = malloc(verts_per_face[f]*sizeof(r3d_int));
  }
  r3d_int ind = 0;
  for(r3d_int f = 0; f < nface; ++f){
    for(r3d_int v = 0; v < verts_per_face[f]; ++v){
      faces[f][v] = face_flat[ind];
      ++ind;
    }
  }
  
  
  // Multiply normal by -1.0 due to difference in convention.
  // Makes consistent with IRL convention     
  r3d_plane planes[a_number_of_planes];
  for(r3d_int n = 0; n < a_number_of_planes; ++n){
    planes[n].n.x = -a_planes[n*4+0];
    planes[n].n.y = -a_planes[n*4+1];
    planes[n].n.z = -a_planes[n*4+2];
    planes[n].d = a_planes[n*4+3];    
  }

  double start = omp_get_wtime();  
  r3d_poly poly;
  r3d_init_poly(&poly, verts, nvert, faces, verts_per_face, nface);
  r3d_clip(&poly, planes, a_number_of_planes);
  r3d_reduce(&poly, a_volume, 0);  
  double end = omp_get_wtime();
  a_times[0] = end - start;  

  for (r3d_int f = 0; f < nface; ++f){
    free(faces[f]);
  }
  free(faces);  
}

void r3d_symPrismByPlanes_total(const double *a_sym_prism_pts,
				const int a_number_of_planes,
				const double *a_planes, double *a_volume,
				double *a_times) {
  // Full BREP for Symmetric PRism
  r3d_int nvert = 11;
  r3d_int nface = 18;

  r3d_int face_flat[54] = {
			   6, 0, 1,
			   6, 1, 2,
			   6, 2, 0,
			   7, 1, 0,
			   7, 4, 1,
			   7, 3, 4,
			   7, 0, 3,
			   8, 2, 1,
			   8, 1, 4,
			   8, 4, 5,
			   8, 5, 2,
			   9, 0, 2,
			   9, 3, 0,
			   9, 5, 3,
			   9, 2, 5,
			   10, 4, 3,
			   10, 3, 5,
			   10, 5, 4
  };
  r3d_int verts_per_face[18] = {
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3							
  };
  
  r3d_rvec3 verts[nvert];
  for(r3d_int p = 0; p < nvert; ++p){
    for(r3d_int d = 0; d < 3; ++d){
      verts[p].xyz[d] = a_sym_prism_pts[p*3+d];
    }
  }
  
  r3d_int **faces = malloc(nface*sizeof(r3d_int*)); 
  for (r3d_int f = 0; f < nface; ++f){
    faces[f] = malloc(verts_per_face[f]*sizeof(r3d_int));
  }
  r3d_int ind = 0;
  for(r3d_int f = 0; f < nface; ++f){
    for(r3d_int v = 0; v < verts_per_face[f]; ++v){
      faces[f][v] = face_flat[ind];
      ++ind;
    }
  }
  
  
  // Multiply normal by -1.0 due to difference in convention.
  // Makes consistent with IRL convention     
  r3d_plane planes[a_number_of_planes];
  for(r3d_int n = 0; n < a_number_of_planes; ++n){
    planes[n].n.x = -a_planes[n*4+0];
    planes[n].n.y = -a_planes[n*4+1];
    planes[n].n.z = -a_planes[n*4+2];
    planes[n].d = a_planes[n*4+3];    
  }

  double start = omp_get_wtime();  
  r3d_poly poly;
  r3d_init_poly(&poly, verts, nvert, faces, verts_per_face, nface);
  r3d_clip(&poly, planes, a_number_of_planes);
  r3d_reduce(&poly, a_volume, 0);  
  double end = omp_get_wtime();
  a_times[0] = end - start;  

  for (r3d_int f = 0; f < nface; ++f){
    free(faces[f]);
  }
  free(faces);  
}


void r3d_symHexByPlanes_total(const double *a_sym_hex_pts,
			      const int a_number_of_planes,
			      const double *a_planes, double *a_volume,
			      double *a_times) {
  // Full BREP for Symmetric Hex
  r3d_int nvert = 14;
  r3d_int nface = 24;
  
  r3d_int face_flat[72] = {
			   0, 1, 8,
			   1, 2, 8,
			   2, 3, 8,
			   3, 0, 8,
			   5, 1, 9,
			   1, 0, 9,
			   0, 4, 9,
			   4, 5, 9,
			   1, 5, 10,
			   5, 6, 10,
			   6, 2, 10,
			   2, 1, 10,
			   2, 6, 11,
			   6, 7, 11,
			   7, 3, 11,
			   3, 2, 11,
			   0, 3, 12,
			   3, 7, 12,
			   7, 4, 12,
			   4, 0, 12,
			   5, 4, 13,
			   4, 7, 13,
			   7, 6, 13,
			   6, 5, 13
  };
  r3d_int verts_per_face[24] = {
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3   
  };
  
  r3d_rvec3 verts[nvert];
  for(r3d_int p = 0; p < nvert; ++p){
    for(r3d_int d = 0; d < 3; ++d){
      verts[p].xyz[d] = a_sym_hex_pts[p*3+d];
    }
  }
  
  r3d_int **faces = malloc(nface*sizeof(r3d_int*)); 
  for (r3d_int f = 0; f < nface; ++f){
    faces[f] = malloc(verts_per_face[f]*sizeof(r3d_int));
  }
  r3d_int ind = 0;
  for(r3d_int f = 0; f < nface; ++f){
    for(r3d_int v = 0; v < verts_per_face[f]; ++v){
      faces[f][v] = face_flat[ind];
      ++ind;
    }
  }
  
  
  // Multiply normal by -1.0 due to difference in convention.
  // Makes consistent with IRL convention     
  r3d_plane planes[a_number_of_planes];
  for(r3d_int n = 0; n < a_number_of_planes; ++n){
    planes[n].n.x = -a_planes[n*4+0];
    planes[n].n.y = -a_planes[n*4+1];
    planes[n].n.z = -a_planes[n*4+2];
    planes[n].d = a_planes[n*4+3];    
  }

  double start = omp_get_wtime();  
  r3d_poly poly;
  r3d_init_poly(&poly, verts, nvert, faces, verts_per_face, nface);
  r3d_clip(&poly, planes, a_number_of_planes);
  r3d_reduce(&poly, a_volume, 0);  
  double end = omp_get_wtime();
  a_times[0] = end - start;  

  for (r3d_int f = 0; f < nface; ++f){
    free(faces[f]);
  }
  free(faces);  
}

void r3d_stelDodecahedronByPlanes_total(const double *a_stel_dodecahedron_pts,
					const int a_number_of_planes,
					const double *a_planes, double *a_volume,
					double *a_times) {
  // Full BREP for Stellated Dodecahedron
  r3d_int nvert = 32;
  r3d_int nface = 60;
  
  r3d_int face_flat[180] = {
			    0, 8, 20,
			    8, 9, 20,
			    9, 4, 20,
			    4, 16, 20,
			    16, 0, 20,
			    0, 16, 21,
			    16, 17, 21,
			    17, 2, 21,
			    2, 12, 21,
			    12, 0, 21,
			    12, 2, 22,
			    2, 10, 22,
			    10, 3, 22,
			    3, 13, 22,
			    13, 12, 22,
			    9, 5, 23,
			    5, 15, 23,
			    15, 14, 23,
			    14, 4, 23,
			    4, 9, 23,
			    3, 19, 24,
			    19, 18, 24,
			    18, 1, 24,
			    1, 13, 24,
			    13, 3, 24,
			    7, 11, 25,
			    11, 6, 25,
			    6, 14, 25,
			    14, 15, 25,
			    15, 7, 25,
			    0, 12, 26,
			    12, 13, 26,
			    13, 1, 26,
			    1, 8, 26,
			    8, 0, 26,
			    8, 1, 27,
			    1, 18, 27,
			    18, 5, 27,
			    5, 9, 27,
			    9, 8, 27,
			    16, 4, 28,
			    4, 14, 28,
			    14, 6, 28,
			    6, 17, 28,
			    17, 16, 28,
			    6, 11, 29,
			    11, 10, 29,
			    10, 2, 29,
			    2, 17, 29,
			    17, 6, 29,
			    7, 15, 30,
			    15, 5, 30,
			    5, 18, 30,
			    18, 19, 30,
			    19, 7, 30,
			    7, 19, 31,
			    19, 3, 31,
			    3, 10, 31,
			    10, 11, 31,
			    11, 7, 31
  };
  r3d_int verts_per_face[60] = {
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3				
  };
  
  r3d_rvec3 verts[nvert];
  for(r3d_int p = 0; p < nvert; ++p){
    for(r3d_int d = 0; d < 3; ++d){
      verts[p].xyz[d] = a_stel_dodecahedron_pts[p*3+d];
    }
  }
  
  r3d_int **faces = malloc(nface*sizeof(r3d_int*)); 
  for (r3d_int f = 0; f < nface; ++f){
    faces[f] = malloc(verts_per_face[f]*sizeof(r3d_int));
  }
  r3d_int ind = 0;
  for(r3d_int f = 0; f < nface; ++f){
    for(r3d_int v = 0; v < verts_per_face[f]; ++v){
      faces[f][v] = face_flat[ind];
      ++ind;
    }
  }
  
  
  // Multiply normal by -1.0 due to difference in convention.
  // Makes consistent with IRL convention     
  r3d_plane planes[a_number_of_planes];
  for(r3d_int n = 0; n < a_number_of_planes; ++n){
    planes[n].n.x = -a_planes[n*4+0];
    planes[n].n.y = -a_planes[n*4+1];
    planes[n].n.z = -a_planes[n*4+2];
    planes[n].d = a_planes[n*4+3];    
  }

  double start = omp_get_wtime();  
  r3d_poly poly;
  r3d_init_poly(&poly, verts, nvert, faces, verts_per_face, nface);
  r3d_clip(&poly, planes, a_number_of_planes);
  r3d_reduce(&poly, a_volume, 0);  
  double end = omp_get_wtime();
  a_times[0] = end - start;  

  for (r3d_int f = 0; f < nface; ++f){
    free(faces[f]);
  }
  free(faces);  
}

void r3d_stelIcosahedronByPlanes_total(const double *a_stel_icosahedron_pts,
				       const int a_number_of_planes,
				       const double *a_planes, double *a_volume,
				       double *a_times) {
  // Full BREP for Stellated Icosahedron
  r3d_int nvert = 32;
  r3d_int nface = 60;
  
  r3d_int face_flat[180] = {
			    0, 8, 12,
			    8, 4, 12,
			    4, 0, 12,
			    0, 5, 13,
			    5, 10, 13,
			    10, 0, 13,
			    2, 4, 14,
			    4, 9, 14,
			    9, 2, 14,
			    2, 11, 15,
			    11, 5, 15,
			    5, 2, 15,
			    1, 6, 16,
			    6, 8, 16,
			    8, 1, 16,
			    1, 10, 17,
			    10, 7, 17,
			    7, 1, 17,
			    3, 9, 18,
			    9, 6, 18,
			    6, 3, 18,
			    3, 7, 19,
			    7, 11, 19,
			    11, 3, 19,
			    0, 10, 20,
			    10, 8, 20,
			    8, 0, 20,
			    1, 8, 21,
			    8, 10, 21,
			    10, 1, 21,
			    2, 9, 22,
			    9, 11, 22,
			    11, 2, 22,
			    3, 11, 23,
			    11, 9, 23,
			    9, 3, 23,
			    4, 2, 24,
			    2, 0, 24,
			    0, 4, 24,
			    5, 0, 25,
			    0, 2, 25,
			    2, 5, 25,
			    6, 1, 26,
			    1, 3, 26,
			    3, 6, 26,
			    7, 3, 27,
			    3, 1, 27,
			    1, 7, 27,
			    8, 6, 28,
			    6, 4, 28,
			    4, 8, 28,
			    9, 4, 29,
			    4, 6, 29,
			    6, 9, 29,
			    10, 5, 30,
			    5, 7, 30,
			    7, 10, 30,
			    11, 7, 31,
			    7, 5, 31,
			    5, 11, 31
  };
  r3d_int verts_per_face[60] = {
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3,
				3,3,3,3,3,3				
  };
  
  r3d_rvec3 verts[nvert];
  for(r3d_int p = 0; p < nvert; ++p){
    for(r3d_int d = 0; d < 3; ++d){
      verts[p].xyz[d] = a_stel_icosahedron_pts[p*3+d];
    }
  }
  
  r3d_int **faces = malloc(nface*sizeof(r3d_int*)); 
  for (r3d_int f = 0; f < nface; ++f){
    faces[f] = malloc(verts_per_face[f]*sizeof(r3d_int));
  }
  r3d_int ind = 0;
  for(r3d_int f = 0; f < nface; ++f){
    for(r3d_int v = 0; v < verts_per_face[f]; ++v){
      faces[f][v] = face_flat[ind];
      ++ind;
    }
  }
  
  
  // Multiply normal by -1.0 due to difference in convention.
  // Makes consistent with IRL convention     
  r3d_plane planes[a_number_of_planes];
  for(r3d_int n = 0; n < a_number_of_planes; ++n){
    planes[n].n.x = -a_planes[n*4+0];
    planes[n].n.y = -a_planes[n*4+1];
    planes[n].n.z = -a_planes[n*4+2];
    planes[n].d = a_planes[n*4+3];    
  }

  double start = omp_get_wtime();  
  r3d_poly poly;
  r3d_init_poly(&poly, verts, nvert, faces, verts_per_face, nface);
  r3d_clip(&poly, planes, a_number_of_planes);
  r3d_reduce(&poly, a_volume, 0);  
  double end = omp_get_wtime();
  a_times[0] = end - start;  

  for (r3d_int f = 0; f < nface; ++f){
    free(faces[f]);
  }
  free(faces);  
}



