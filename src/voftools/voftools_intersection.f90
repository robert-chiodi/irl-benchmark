module VOFtools_mod
  use, intrinsic :: iso_fortran_env, only: r8 => REAL64
  use VOFtools_wrapper
  use omp_lib
  use iso_c_binding  

contains

  subroutine voftools_prismByPlanes(a_prism_pts, &
       a_number_of_planes, a_planes, a_volume, a_times) &
       bind(C, name = "c_voftools_prismByPlanes")
    
    implicit none

    real(C_DOUBLE), intent(in) :: a_prism_pts(*)
    integer(C_INT), intent(in), value :: a_number_of_planes
    real(C_DOUBLE), intent(in) :: a_planes(*)
    real(C_DOUBLE), intent(out) :: a_volume
    real(C_DOUBLE), intent(out) :: a_times(*)

    real(r8) :: start, end
    type(polyhedron) :: poly
    integer(C_INT) :: p
    
    start = omp_get_wtime()
    call make_prism(poly,a_prism_pts)
    end = omp_get_wtime()    
    a_times(1) = end - start

    start = omp_get_wtime()
    do p = 0, a_number_of_planes-1
       ! Multiply normal by -1.0 due to difference in convention.
       ! Makes consistent with IRL convention       
       call VOFtools_INTE3D_wrapper(poly, [-a_planes(p*4+1:p*4+3),a_planes(p*4+4)])
       if(poly%NTS == 0) then
          exit
       end if        
    end do
    end = omp_get_wtime() 
    a_times(2) = end - start

    start = omp_get_wtime()
    a_volume = VOFtools_TOOLV3D_wrapper(poly)
    end = omp_get_wtime() 
    a_times(3) = end - start

  end subroutine voftools_prismByPlanes

  subroutine voftools_unitCubeByPlanes(a_cube_pts, &
       a_number_of_planes, a_planes, a_volume, a_times) &
       bind(C, name = "c_voftools_unitCubeByPlanes")
    
    implicit none

    real(C_DOUBLE), intent(in) :: a_cube_pts(*)
    integer(C_INT), intent(in), value :: a_number_of_planes
    real(C_DOUBLE), intent(in) :: a_planes(*)
    real(C_DOUBLE), intent(out) :: a_volume
    real(C_DOUBLE), intent(out) :: a_times(*)

    real(r8) :: start, end
    type(polyhedron) :: poly
    integer(C_INT) :: p
    
    start = omp_get_wtime()
    call make_cuboid(poly, a_cube_pts)
    end = omp_get_wtime()    
    a_times(1) = end - start

    start = omp_get_wtime()
    do p = 0, a_number_of_planes-1
       ! Multiply normal by -1.0 due to difference in convention.
       ! Makes consistent with IRL convention       
       call VOFtools_INTE3D_wrapper(poly, [-a_planes(p*4+1:p*4+3),a_planes(p*4+4)])
       if(poly%NTS == 0) then
          exit
       end if        
    end do
    end = omp_get_wtime() 
    a_times(2) = end - start

    start = omp_get_wtime()
    a_volume = VOFtools_TOOLV3D_wrapper(poly)
    end = omp_get_wtime() 
    a_times(3) = end - start

  end subroutine voftools_unitCubeByPlanes

  subroutine voftools_triPrismByPlanes(a_tri_prism_pts, &
       a_number_of_planes, a_planes, a_volume, a_times) &
       bind(C, name = "c_voftools_triPrismByPlanes")
    
    implicit none

    real(C_DOUBLE), intent(in) :: a_tri_prism_pts(*)
    integer(C_INT), intent(in), value :: a_number_of_planes
    real(C_DOUBLE), intent(in) :: a_planes(*)
    real(C_DOUBLE), intent(out) :: a_volume
    real(C_DOUBLE), intent(out) :: a_times(*)

    real(r8) :: start, end
    type(polyhedron) :: poly
    integer(C_INT) :: p
    
    start = omp_get_wtime()
    call make_tri_prism(poly, a_tri_prism_pts)
    end = omp_get_wtime()    
    a_times(1) = end - start

    start = omp_get_wtime()
    do p = 0, a_number_of_planes-1
       ! Multiply normal by -1.0 due to difference in convention.
       ! Makes consistent with IRL convention       
       call VOFtools_INTE3D_wrapper(poly, [-a_planes(p*4+1:p*4+3),a_planes(p*4+4)])
       if(poly%NTS == 0) then
          exit
       end if        
    end do
    end = omp_get_wtime() 
    a_times(2) = end - start

    start = omp_get_wtime()
    a_volume = VOFtools_TOOLV3D_wrapper(poly)
    end = omp_get_wtime() 
    a_times(3) = end - start

  end subroutine voftools_triPrismByPlanes

  subroutine voftools_triHexByPlanes(a_tri_hex_pts, &
       a_number_of_planes, a_planes, a_volume, a_times) &
       bind(C, name = "c_voftools_triHexByPlanes")
    
    implicit none

    real(C_DOUBLE), intent(in) :: a_tri_hex_pts(*)
    integer(C_INT), intent(in), value :: a_number_of_planes
    real(C_DOUBLE), intent(in) :: a_planes(*)
    real(C_DOUBLE), intent(out) :: a_volume
    real(C_DOUBLE), intent(out) :: a_times(*)

    real(r8) :: start, end
    type(polyhedron) :: poly
    integer(C_INT) :: p
    
    start = omp_get_wtime()
    call make_tri_hex(poly, a_tri_hex_pts)
    end = omp_get_wtime()    
    a_times(1) = end - start

    start = omp_get_wtime()
    do p = 0, a_number_of_planes-1
       ! Multiply normal by -1.0 due to difference in convention.
       ! Makes consistent with IRL convention       
       call VOFtools_INTE3D_wrapper(poly, [-a_planes(p*4+1:p*4+3),a_planes(p*4+4)])
       if(poly%NTS == 0) then
          exit
       end if        
    end do
    end = omp_get_wtime() 
    a_times(2) = end - start

    start = omp_get_wtime()
    a_volume = VOFtools_TOOLV3D_wrapper(poly)
    end = omp_get_wtime() 
    a_times(3) = end - start

  end subroutine voftools_triHexByPlanes

  subroutine voftools_symPrismByPlanes(a_sym_prism_pts, &
       a_number_of_planes, a_planes, a_volume, a_times) &
       bind(C, name = "c_voftools_symPrismByPlanes")
    
    implicit none

    real(C_DOUBLE), intent(in) :: a_sym_prism_pts(*)
    integer(C_INT), intent(in), value :: a_number_of_planes
    real(C_DOUBLE), intent(in) :: a_planes(*)
    real(C_DOUBLE), intent(out) :: a_volume
    real(C_DOUBLE), intent(out) :: a_times(*)

    real(r8) :: start, end
    type(polyhedron) :: poly
    integer(C_INT) :: p
    
    start = omp_get_wtime()
    call make_sym_prism(poly, a_sym_prism_pts)
    end = omp_get_wtime()    
    a_times(1) = end - start

    start = omp_get_wtime()
    do p = 0, a_number_of_planes-1
       ! Multiply normal by -1.0 due to difference in convention.
       ! Makes consistent with IRL convention       
       call VOFtools_INTE3D_wrapper(poly, [-a_planes(p*4+1:p*4+3),a_planes(p*4+4)])
       if(poly%NTS == 0) then
          exit
       end if        
    end do
    end = omp_get_wtime() 
    a_times(2) = end - start

    start = omp_get_wtime()
    a_volume = VOFtools_TOOLV3D_wrapper(poly)
    end = omp_get_wtime() 
    a_times(3) = end - start

  end subroutine voftools_symPrismByPlanes

  subroutine voftools_symHexByPlanes(a_sym_hex_pts, &
       a_number_of_planes, a_planes, a_volume, a_times) &
       bind(C, name = "c_voftools_symHexByPlanes")
    
    implicit none

    real(C_DOUBLE), intent(in) :: a_sym_hex_pts(*)
    integer(C_INT), intent(in), value :: a_number_of_planes
    real(C_DOUBLE), intent(in) :: a_planes(*)
    real(C_DOUBLE), intent(out) :: a_volume
    real(C_DOUBLE), intent(out) :: a_times(*)

    real(r8) :: start, end
    type(polyhedron) :: poly
    integer(C_INT) :: p
    
    start = omp_get_wtime()
    call make_sym_hex(poly, a_sym_hex_pts)
    end = omp_get_wtime()    
    a_times(1) = end - start

    start = omp_get_wtime()
    do p = 0, a_number_of_planes-1
       ! Multiply normal by -1.0 due to difference in convention.
       ! Makes consistent with IRL convention       
       call VOFtools_INTE3D_wrapper(poly, [-a_planes(p*4+1:p*4+3),a_planes(p*4+4)])
       if(poly%NTS == 0) then
          exit
       end if        
    end do
    end = omp_get_wtime() 
    a_times(2) = end - start

    start = omp_get_wtime()
    a_volume = VOFtools_TOOLV3D_wrapper(poly)
    end = omp_get_wtime() 
    a_times(3) = end - start

  end subroutine voftools_symHexByPlanes

  subroutine voftools_stelDodecahedronByPlanes(a_stel_dodecahedron_pts, &
       a_number_of_planes, a_planes, a_volume, a_times) &
       bind(C, name = "c_voftools_stelDodecahedronByPlanes")
    
    implicit none

    real(C_DOUBLE), intent(in) :: a_stel_dodecahedron_pts(*)
    integer(C_INT), intent(in), value :: a_number_of_planes
    real(C_DOUBLE), intent(in) :: a_planes(*)
    real(C_DOUBLE), intent(out) :: a_volume
    real(C_DOUBLE), intent(out) :: a_times(*)

    real(r8) :: start, end
    type(polyhedron) :: poly
    integer(C_INT) :: p
    
    start = omp_get_wtime()
    call make_stellated_dodecahedron(poly, a_stel_dodecahedron_pts)
    end = omp_get_wtime()    
    a_times(1) = end - start

    start = omp_get_wtime()
    do p = 0, a_number_of_planes-1
       ! Multiply normal by -1.0 due to difference in convention.
       ! Makes consistent with IRL convention       
       call VOFtools_INTE3D_wrapper(poly, [-a_planes(p*4+1:p*4+3),a_planes(p*4+4)])
       if(poly%NTS == 0) then
          exit
       end if        
    end do
    end = omp_get_wtime() 
    a_times(2) = end - start

    start = omp_get_wtime()
    a_volume = VOFtools_TOOLV3D_wrapper(poly)
    end = omp_get_wtime() 
    a_times(3) = end - start

  end subroutine voftools_stelDodecahedronByPlanes
  

  subroutine voftools_stelIcosahedronByPlanes(a_stel_icosahedron_pts, &
       a_number_of_planes, a_planes, a_volume, a_times) &
       bind(C, name = "c_voftools_stelIcosahedronByPlanes")
    
    implicit none

    real(C_DOUBLE), intent(in) :: a_stel_icosahedron_pts(*)
    integer(C_INT), intent(in), value :: a_number_of_planes
    real(C_DOUBLE), intent(in) :: a_planes(*)
    real(C_DOUBLE), intent(out) :: a_volume
    real(C_DOUBLE), intent(out) :: a_times(*)

    real(r8) :: start, end
    type(polyhedron) :: poly
    integer(C_INT) :: p
    
    start = omp_get_wtime()
    call make_stellated_icosahedron(poly, a_stel_icosahedron_pts)
    end = omp_get_wtime()    
    a_times(1) = end - start

    start = omp_get_wtime()
    do p = 0, a_number_of_planes-1
       ! Multiply normal by -1.0 due to difference in convention.
       ! Makes consistent with IRL convention       
       call VOFtools_INTE3D_wrapper(poly, [-a_planes(p*4+1:p*4+3),a_planes(p*4+4)])
       if(poly%NTS == 0) then
          exit
       end if        
    end do
    end = omp_get_wtime() 
    a_times(2) = end - start

    start = omp_get_wtime()
    a_volume = VOFtools_TOOLV3D_wrapper(poly)
    end = omp_get_wtime() 
    a_times(3) = end - start

  end subroutine voftools_stelIcosahedronByPlanes

  
!!!! Implementation of same functions from above but timing everything at once !!!!

  subroutine voftools_prismByPlanes_total(a_prism_pts, &
       a_number_of_planes, a_planes, a_volume, a_times) &
       bind(C, name = "c_voftools_prismByPlanes_total")
    
    implicit none

    real(C_DOUBLE), intent(in) :: a_prism_pts(*)
    integer(C_INT), intent(in), value :: a_number_of_planes
    real(C_DOUBLE), intent(in) :: a_planes(*)
    real(C_DOUBLE), intent(out) :: a_volume
    real(C_DOUBLE), intent(out) :: a_times(*)

    real(r8) :: start, end
    type(polyhedron) :: poly
    integer(C_INT) :: p
    
    start = omp_get_wtime()
    call make_prism(poly,a_prism_pts)

    do p = 0, a_number_of_planes-1
       ! Multiply normal by -1.0 due to difference in convention.
       ! Makes consistent with IRL convention       
       call VOFtools_INTE3D_wrapper(poly, [-a_planes(p*4+1:p*4+3),a_planes(p*4+4)])
       if(poly%NTS == 0) then
          exit
       end if       
    end do
    
    a_volume = VOFtools_TOOLV3D_wrapper(poly)
    end = omp_get_wtime()    
    a_times(1) = end - start
    
  end subroutine voftools_prismByPlanes_total

  subroutine voftools_unitCubeByPlanes_total(a_cube_pts, &
       a_number_of_planes, a_planes, a_volume, a_times) &
       bind(C, name = "c_voftools_unitCubeByPlanes_total")
    
    implicit none

    real(C_DOUBLE), intent(in) :: a_cube_pts(*)
    integer(C_INT), intent(in), value :: a_number_of_planes
    real(C_DOUBLE), intent(in) :: a_planes(*)
    real(C_DOUBLE), intent(out) :: a_volume
    real(C_DOUBLE), intent(out) :: a_times(*)

    real(r8) :: start, end
    type(polyhedron) :: poly
    integer(C_INT) :: p
    
    start = omp_get_wtime()
    call make_cuboid(poly, a_cube_pts)

    do p = 0, a_number_of_planes-1
       ! Multiply normal by -1.0 due to difference in convention.
       ! Makes consistent with IRL convention       
       call VOFtools_INTE3D_wrapper(poly, [-a_planes(p*4+1:p*4+3),a_planes(p*4+4)])
       if(poly%NTS == 0) then
          exit
       end if        
    end do

    a_volume = VOFtools_TOOLV3D_wrapper(poly)
    end = omp_get_wtime()    
    a_times(1) = end - start

  end subroutine voftools_unitCubeByPlanes_total

  subroutine voftools_triPrismByPlanes_total(a_tri_prism_pts, &
       a_number_of_planes, a_planes, a_volume, a_times) &
       bind(C, name = "c_voftools_triPrismByPlanes_total")
    
    implicit none

    real(C_DOUBLE), intent(in) :: a_tri_prism_pts(*)
    integer(C_INT), intent(in), value :: a_number_of_planes
    real(C_DOUBLE), intent(in) :: a_planes(*)
    real(C_DOUBLE), intent(out) :: a_volume
    real(C_DOUBLE), intent(out) :: a_times(*)

    real(r8) :: start, end
    type(polyhedron) :: poly
    integer(C_INT) :: p
    
    start = omp_get_wtime()
    call make_tri_prism(poly, a_tri_prism_pts)

    do p = 0, a_number_of_planes-1
       ! Multiply normal by -1.0 due to difference in convention.
       ! Makes consistent with IRL convention       
       call VOFtools_INTE3D_wrapper(poly, [-a_planes(p*4+1:p*4+3),a_planes(p*4+4)])
       if(poly%NTS == 0) then
          exit
       end if        
    end do

    a_volume = VOFtools_TOOLV3D_wrapper(poly)
    end = omp_get_wtime()    
    a_times(1) = end - start

  end subroutine voftools_triPrismByPlanes_total

  subroutine voftools_triHexByPlanes_total(a_tri_hex_pts, &
       a_number_of_planes, a_planes, a_volume, a_times) &
       bind(C, name = "c_voftools_triHexByPlanes_total")
    
    implicit none

    real(C_DOUBLE), intent(in) :: a_tri_hex_pts(*)
    integer(C_INT), intent(in), value :: a_number_of_planes
    real(C_DOUBLE), intent(in) :: a_planes(*)
    real(C_DOUBLE), intent(out) :: a_volume
    real(C_DOUBLE), intent(out) :: a_times(*)

    real(r8) :: start, end
    type(polyhedron) :: poly
    integer(C_INT) :: p
    
    start = omp_get_wtime()
    call make_tri_hex(poly, a_tri_hex_pts)

    do p = 0, a_number_of_planes-1
       ! Multiply normal by -1.0 due to difference in convention.
       ! Makes consistent with IRL convention       
       call VOFtools_INTE3D_wrapper(poly, [-a_planes(p*4+1:p*4+3),a_planes(p*4+4)])
       if(poly%NTS == 0) then
          exit
       end if        
    end do

    a_volume = VOFtools_TOOLV3D_wrapper(poly)
    end = omp_get_wtime()    
    a_times(1) = end - start
    
  end subroutine voftools_triHexByPlanes_total

  subroutine voftools_symPrismByPlanes_total(a_sym_prism_pts, &
       a_number_of_planes, a_planes, a_volume, a_times) &
       bind(C, name = "c_voftools_symPrismByPlanes_total")
    
    implicit none

    real(C_DOUBLE), intent(in) :: a_sym_prism_pts(*)
    integer(C_INT), intent(in), value :: a_number_of_planes
    real(C_DOUBLE), intent(in) :: a_planes(*)
    real(C_DOUBLE), intent(out) :: a_volume
    real(C_DOUBLE), intent(out) :: a_times(*)

    real(r8) :: start, end
    type(polyhedron) :: poly
    integer(C_INT) :: p
    
    start = omp_get_wtime()
    call make_sym_prism(poly, a_sym_prism_pts)

    do p = 0, a_number_of_planes-1
       ! Multiply normal by -1.0 due to difference in convention.
       ! Makes consistent with IRL convention       
       call VOFtools_INTE3D_wrapper(poly, [-a_planes(p*4+1:p*4+3),a_planes(p*4+4)])
       if(poly%NTS == 0) then
          exit
       end if        
    end do

    a_volume = VOFtools_TOOLV3D_wrapper(poly)
    end = omp_get_wtime()    
    a_times(1) = end - start

  end subroutine voftools_symPrismByPlanes_total

  subroutine voftools_symHexByPlanes_total(a_sym_hex_pts, &
       a_number_of_planes, a_planes, a_volume, a_times) &
       bind(C, name = "c_voftools_symHexByPlanes_total")
    
    implicit none

    real(C_DOUBLE), intent(in) :: a_sym_hex_pts(*)
    integer(C_INT), intent(in), value :: a_number_of_planes
    real(C_DOUBLE), intent(in) :: a_planes(*)
    real(C_DOUBLE), intent(out) :: a_volume
    real(C_DOUBLE), intent(out) :: a_times(*)

    real(r8) :: start, end
    type(polyhedron) :: poly
    integer(C_INT) :: p
    
    start = omp_get_wtime()
    call make_sym_hex(poly, a_sym_hex_pts)

    do p = 0, a_number_of_planes-1
       ! Multiply normal by -1.0 due to difference in convention.
       ! Makes consistent with IRL convention       
       call VOFtools_INTE3D_wrapper(poly, [-a_planes(p*4+1:p*4+3),a_planes(p*4+4)])
       if(poly%NTS == 0) then
          exit
       end if        
    end do

    a_volume = VOFtools_TOOLV3D_wrapper(poly)
    end = omp_get_wtime()    
    a_times(1) = end - start

  end subroutine voftools_symHexByPlanes_total

  subroutine voftools_stelDodecahedronByPlanes_total(a_stel_dodecahedron_pts, &
       a_number_of_planes, a_planes, a_volume, a_times) &
       bind(C, name = "c_voftools_stelDodecahedronByPlanes_total")
    
    implicit none

    real(C_DOUBLE), intent(in) :: a_stel_dodecahedron_pts(*)
    integer(C_INT), intent(in), value :: a_number_of_planes
    real(C_DOUBLE), intent(in) :: a_planes(*)
    real(C_DOUBLE), intent(out) :: a_volume
    real(C_DOUBLE), intent(out) :: a_times(*)

    real(r8) :: start, end
    type(polyhedron) :: poly
    integer(C_INT) :: p
    
    start = omp_get_wtime()
    call make_stellated_dodecahedron(poly, a_stel_dodecahedron_pts)

    do p = 0, a_number_of_planes-1
       ! Multiply normal by -1.0 due to difference in convention.
       ! Makes consistent with IRL convention       
       call VOFtools_INTE3D_wrapper(poly, [-a_planes(p*4+1:p*4+3),a_planes(p*4+4)])
       if(poly%NTS == 0) then
          exit
       end if        
    end do

    a_volume = VOFtools_TOOLV3D_wrapper(poly)
    end = omp_get_wtime()    
    a_times(1) = end - start

  end subroutine voftools_stelDodecahedronByPlanes_total
  

  subroutine voftools_stelIcosahedronByPlanes_total(a_stel_icosahedron_pts, &
       a_number_of_planes, a_planes, a_volume, a_times) &
       bind(C, name = "c_voftools_stelIcosahedronByPlanes_total")
    
    implicit none

    real(C_DOUBLE), intent(in) :: a_stel_icosahedron_pts(*)
    integer(C_INT), intent(in), value :: a_number_of_planes
    real(C_DOUBLE), intent(in) :: a_planes(*)
    real(C_DOUBLE), intent(out) :: a_volume
    real(C_DOUBLE), intent(out) :: a_times(*)

    real(r8) :: start, end
    type(polyhedron) :: poly
    integer(C_INT) :: p
    
    start = omp_get_wtime()
    call make_stellated_icosahedron(poly, a_stel_icosahedron_pts)
    
    do p = 0, a_number_of_planes-1
       ! Multiply normal by -1.0 due to difference in convention.
       ! Makes consistent with IRL convention       
       call VOFtools_INTE3D_wrapper(poly, [-a_planes(p*4+1:p*4+3),a_planes(p*4+4)])
       if(poly%NTS == 0) then
          exit
       end if        
    end do

    a_volume = VOFtools_TOOLV3D_wrapper(poly)
    end = omp_get_wtime()    
    a_times(1) = end - start

  end subroutine voftools_stelIcosahedronByPlanes_total
  


end module VOFtools_mod
