module VOFtools_wrapper
  use, intrinsic :: iso_fortran_env, only: r8 => REAL64
  use iso_c_binding

  integer, parameter :: NS = 100 ! Value should match that in VOFTools/dim.h
  integer, parameter :: NV = 140 ! Value should match that in VOFTools/dim.h
  
  type polyhedron
     integer, public :: NTS ! Number of faces
     integer, public :: NIPV(NS) ! Number of vertices on each face
     integer, public :: NTP ! Last global vertex  number
     integer, public :: NTV ! Total number of vertices     
     integer, public :: IPV(NS,NV) ! Definition of Faces by Vertex Index, (face, vertex)
     real(r8), public :: VERTP(NV,3) ! Vertex locations (vertex, dimension)
     real(r8), public :: NORMAL(NS,3) ! Normal for each face
  end type polyhedron

contains
  
  subroutine VOFtools_INTE3D_wrapper(poly, plane)
    type(polyhedron), intent(inout) :: poly
    real(r8), intent(in) :: plane(4)

    integer :: ICONTN, ICONTP

    call INTE3D(plane(4),ICONTN,ICONTP,poly%IPV,poly%NIPV,poly%NTP,poly%NTS,poly%NTV, &
         poly%VERTP,plane(1),poly%NORMAL(:,1),plane(2),poly%NORMAL(:,2),plane(3),poly%NORMAL(:,3))
    ! VOFTools does not set number of polyhedron faces/verts to 0
    ! if it is entirely above the clipping plane, so do that here.
    ! This will make sure subsequent volume calculations on this
    ! entirely polyhedron will be 0.
    if(ICONTP == 0) then
       poly%NTS = 0
       poly%NTV = 0
       poly%NTP = 0
    end if
       
    
  end subroutine VOFtools_INTE3D_wrapper

  function VOFtools_TOOLV3D_wrapper(poly) result(volume)
    type(polyhedron), intent(inout) :: poly
    real(r8) :: volume

    call TOOLV3D(poly%IPV,poly%NIPV,poly%NTS,poly%VERTP,volume,poly%NORMAL(:,1),poly%NORMAL(:,2),poly%NORMAL(:,3))
    
  end function VOFtools_TOOLV3D_wrapper

  pure function cross_product(vector1, vector2) result(cp_value)
    implicit none
    
    real(r8), intent(in) :: vector1(3), vector2(3)
    real(r8) :: cp_value(3)
   
    cp_value(1) = vector1(2) * vector2(3) - vector1(3) * vector2(2)
    cp_value(2) = vector1(3) * vector2(1) - vector1(1) * vector2(3)
    cp_value(3) = vector1(1) * vector2(2) - vector1(2) * vector2(1)
    
    return
    
  end function cross_product

  subroutine make_prism(poly, a_pts)
    implicit none

    type(polyhedron), intent(out) :: poly    
    real(r8), intent(in) :: a_pts(3,6)

    integer :: n, f
    integer :: tri(3)
    real(r8) :: face_normal(3)

    ! Definition of what a Prism is
    poly%NTS = 5
    poly%NIPV(1:5) = [3,4,4,4,3]
    poly%NTP = 6
    poly%NTV = 6

    poly%IPV(1,1:3) = [1, 2, 3] 
    poly%IPV(2,1:4) = [1, 4, 5, 2] 
    poly%IPV(3,1:4) = [1, 3, 6, 4] 
    poly%IPV(4,1:4) = [2, 5, 6, 3]
    poly%IPV(5,1:3) = [4, 6, 5]        

    ! Copy over vertex locations
    do n = 1,poly%NTV
       poly%VERTP(n,:) = a_pts(:,n)
    end do

    ! Calculate and set normal for each face
    do f = 1, poly%NTS
       face_normal = 0.0_r8
       do n = 1, poly%NIPV(f)-2 ! Triangulate and calculate normal
          tri(1) = poly%IPV(f,1)
          tri(2) = poly%IPV(f,n+1)
          tri(3) = poly%IPV(f,n+2)
          face_normal = face_normal + cross_product(poly%VERTP(tri(2),:)-poly%VERTP(tri(1),:), &
                                                    poly%VERTP(tri(3),:)-poly%VERTP(tri(1),:)  )  
       end do
       poly%NORMAL(f,1:3) = face_normal / norm2(face_normal)
    end do
    
  end subroutine make_prism

  subroutine make_cuboid(poly, a_cuboid_pts)
    implicit none

    type(polyhedron), intent(out) :: poly    
    real(C_DOUBLE), intent(in) :: a_cuboid_pts(6)

    integer :: n

    ! Definition of what a Cuboid is
    poly%NTS = 6
    poly%NIPV(1:6) = [4,4,4,4,4,4]
    poly%NTP = 8
    poly%NTV = 8

    poly%IPV(1,1:4) = [1,2,3,4]
    poly%IPV(2,1:4) = [2,1,5,6]
    poly%IPV(3,1:4) = [3,2,6,7]
    poly%IPV(4,1:4) = [4,3,7,8]
    poly%IPV(5,1:4) = [1,4,8,5]
    poly%IPV(6,1:4) = [8,7,6,5]        

    ! Copy over vertex locations
    poly%VERTP(1,:) = [a_cuboid_pts(4),a_cuboid_pts(2), a_cuboid_pts(3)]
    poly%VERTP(2,:) = [a_cuboid_pts(4),a_cuboid_pts(5), a_cuboid_pts(3)]
    poly%VERTP(3,:) = [a_cuboid_pts(4),a_cuboid_pts(5), a_cuboid_pts(6)]
    poly%VERTP(4,:) = [a_cuboid_pts(4),a_cuboid_pts(2), a_cuboid_pts(6)]
    poly%VERTP(5,:) = [a_cuboid_pts(1),a_cuboid_pts(2), a_cuboid_pts(3)]
    poly%VERTP(6,:) = [a_cuboid_pts(1),a_cuboid_pts(5), a_cuboid_pts(3)]
    poly%VERTP(7,:) = [a_cuboid_pts(1),a_cuboid_pts(5), a_cuboid_pts(6)]
    poly%VERTP(8,:) = [a_cuboid_pts(1),a_cuboid_pts(2), a_cuboid_pts(6)]        

    ! Calculate and set normal for each face
    poly%NORMAL(1,:) = [1.0, 0.0, 0.0]
    poly%NORMAL(2,:) = [0.0, 0.0, -1.0]
    poly%NORMAL(3,:) = [0.0, 1.0, 0.0]
    poly%NORMAL(4,:) = [0.0, 0.0, 1.0]
    poly%NORMAL(5,:) = [0.0, -1.0, 0.0]
    poly%NORMAL(6,:) = [-1.0, 0.0, 0.0]

  end subroutine make_cuboid

  subroutine make_tri_prism(poly, a_pts)

    type(polyhedron), intent(out) :: poly
    real(C_DOUBLE), intent(in) :: a_pts(3,6)
    
    integer :: n, f
    integer :: tri(3)
    real(r8) :: face_normal(3)
    
    ! Definition of what a Triangulated Prism is
    poly%NTS = 8
    poly%NIPV(1:8) = 3
    poly%NTP = 6
    poly%NTV = 6

    poly%IPV(1,1:3) = [1, 2, 3]
    poly%IPV(2,1:3) = [5, 4, 6]
    poly%IPV(3,1:3) = [5, 6, 3]
    poly%IPV(4,1:3) = [5, 3, 2]
    poly%IPV(5,1:3) = [5, 2, 1]
    poly%IPV(6,1:3) = [5, 1, 4]
    poly%IPV(7,1:3) = [1, 3, 6]
    poly%IPV(8,1:3) = [1, 6, 4]

    ! Copy over vertex locations
    do n = 1,poly%NTV
       poly%VERTP(n,:) = a_pts(:,n)
    end do

    ! Calculate and set normal for each face
    do f = 1, poly%NTS
       face_normal = 0.0_r8
       do n = 1, poly%NIPV(f)-2 ! Triangulate and calculate normal
          tri(1) = poly%IPV(f,1)
          tri(2) = poly%IPV(f,n+1)
          tri(3) = poly%IPV(f,n+2)
          face_normal = face_normal + cross_product(poly%VERTP(tri(2),:)-poly%VERTP(tri(1),:), &
               poly%VERTP(tri(3),:)-poly%VERTP(tri(1),:)  )  
       end do
       poly%NORMAL(f,1:3) = face_normal / norm2(face_normal)
    end do

  end subroutine make_tri_prism

  subroutine make_tri_hex(poly, a_pts)

    type(polyhedron), intent(out) :: poly
    real(C_DOUBLE), intent(in) :: a_pts(3,8)    

    integer :: n, f
    integer :: tri(3)
    real(r8) :: face_normal(3)

    
    ! Definition of what a Triangulated Hexahedron is
    poly%NTS = 12
    poly%NIPV(1:12) = 3
    poly%NTP = 8
    poly%NTV = 8
    
    poly%IPV(1,1:3) = [6, 8, 7]
    poly%IPV(2,1:3) = [6, 5, 8]
    poly%IPV(3,1:3) = [4, 1, 2]
    poly%IPV(4,1:3) = [4, 2, 3]
    poly%IPV(5,1:3) = [5, 4, 8]
    poly%IPV(6,1:3) = [5, 1, 4]
    poly%IPV(7,1:3) = [3, 6, 7]
    poly%IPV(8,1:3) = [3, 2, 6]
    poly%IPV(9,1:3) = [1, 6, 2]
    poly%IPV(10,1:3) = [1, 5, 6]
    poly%IPV(11,1:3) = [4, 7, 8]
    poly%IPV(12,1:3) = [4, 3, 7]

    ! Copy over vertex locations
    do n = 1,poly%NTV
       poly%VERTP(n,:) = a_pts(:,n)
    end do

    ! Calculate and set normal for each face
    do f = 1, poly%NTS
       face_normal = 0.0_r8
       do n = 1, poly%NIPV(f)-2 ! Triangulate and calculate normal
          tri(1) = poly%IPV(f,1)
          tri(2) = poly%IPV(f,n+1)
          tri(3) = poly%IPV(f,n+2)
          face_normal = face_normal + cross_product(poly%VERTP(tri(2),:)-poly%VERTP(tri(1),:), &
               poly%VERTP(tri(3),:)-poly%VERTP(tri(1),:)  )  
       end do
       poly%NORMAL(f,1:3) = face_normal / norm2(face_normal)
    end do

  end subroutine make_tri_hex

  subroutine make_sym_prism(poly, a_pts)

    type(polyhedron), intent(out) :: poly
    real(C_DOUBLE), intent(in) :: a_pts(3,11)
    
    integer :: n, f
    integer :: tri(3)
    real(r8) :: face_normal(3)
    
    ! Definition of what a Symmetric Prism is
    poly%NTS = 18
    poly%NIPV(1:18) = 3
    poly%NTP = 11
    poly%NTV = 11

    poly%IPV(1,1:3) = [7, 1, 2]
    poly%IPV(2,1:3) = [7, 2, 3]
    poly%IPV(3,1:3) = [7, 3, 1]
    poly%IPV(4,1:3) = [8, 2, 1]
    poly%IPV(5,1:3) = [8, 5, 2]
    poly%IPV(6,1:3) = [8, 4, 5]
    poly%IPV(7,1:3) = [8, 1, 4]
    poly%IPV(8,1:3) = [9, 3, 2]
    poly%IPV(9,1:3) = [9, 2, 5]
    poly%IPV(10,1:3) = [9, 5, 6]
    poly%IPV(11,1:3) = [9, 6, 3]
    poly%IPV(12,1:3) = [10, 1, 3]
    poly%IPV(13,1:3) = [10, 4, 1]
    poly%IPV(14,1:3) = [10, 6, 4]
    poly%IPV(15,1:3) = [10, 3, 6]
    poly%IPV(16,1:3) = [11, 5, 4]
    poly%IPV(17,1:3) = [11, 4, 6]
    poly%IPV(18,1:3) = [11, 6, 5]
    
    ! Copy over vertex locations
    do n = 1,poly%NTV
       poly%VERTP(n,:) = a_pts(:,n)
    end do

    ! Calculate and set normal for each face
    do f = 1, poly%NTS
       face_normal = 0.0_r8
       do n = 1, poly%NIPV(f)-2 ! Triangulate and calculate normal
          tri(1) = poly%IPV(f,1)
          tri(2) = poly%IPV(f,n+1)
          tri(3) = poly%IPV(f,n+2)
          face_normal = face_normal + cross_product(poly%VERTP(tri(2),:)-poly%VERTP(tri(1),:), &
               poly%VERTP(tri(3),:)-poly%VERTP(tri(1),:)  )  
       end do
       poly%NORMAL(f,1:3) = face_normal / norm2(face_normal)
    end do

  end subroutine make_sym_prism

  
  subroutine make_sym_hex(poly, a_pts)

    type(polyhedron), intent(out) :: poly
    real(C_DOUBLE), intent(in) :: a_pts(3,14)
    
    integer :: n, f
    integer :: tri(3)
    real(r8) :: face_normal(3)
    
    ! Definition of what a Symmetric Hexahedron is
    poly%NTS = 24
    poly%NIPV(1:24) = 3
    poly%NTP = 14
    poly%NTV = 14

    poly%IPV(1,1:3) = [1, 2, 9]
    poly%IPV(2,1:3) = [2, 3, 9]
    poly%IPV(3,1:3) = [3, 4, 9]
    poly%IPV(4,1:3) = [4, 1, 9]
    poly%IPV(5,1:3) = [6, 2, 10]
    poly%IPV(6,1:3) = [2, 1, 10]
    poly%IPV(7,1:3) = [1, 5, 10]
    poly%IPV(8,1:3) = [5, 6, 10]
    poly%IPV(9,1:3) = [2, 6, 11]
    poly%IPV(10,1:3) = [6, 7, 11]
    poly%IPV(11,1:3) = [7, 3, 11]
    poly%IPV(12,1:3) = [3, 2, 11]
    poly%IPV(13,1:3) = [3, 7, 12]
    poly%IPV(14,1:3) = [7, 8, 12]
    poly%IPV(15,1:3) = [8, 4, 12]
    poly%IPV(16,1:3) = [4, 3, 12]
    poly%IPV(17,1:3) = [1, 4, 13]
    poly%IPV(18,1:3) = [4, 8, 13]
    poly%IPV(19,1:3) = [8, 5, 13]
    poly%IPV(20,1:3) = [5, 1, 13]
    poly%IPV(21,1:3) = [6, 5, 14]
    poly%IPV(22,1:3) = [5, 8, 14]
    poly%IPV(23,1:3) = [8, 7, 14]
    poly%IPV(24,1:3) = [7, 6, 14]
    
    
    ! Copy over vertex locations
    do n = 1,poly%NTV
       poly%VERTP(n,:) = a_pts(:,n)
    end do

    ! Calculate and set normal for each face
    do f = 1, poly%NTS
       face_normal = 0.0_r8
       do n = 1, poly%NIPV(f)-2 ! Triangulate and calculate normal
          tri(1) = poly%IPV(f,1)
          tri(2) = poly%IPV(f,n+1)
          tri(3) = poly%IPV(f,n+2)
          face_normal = face_normal + cross_product(poly%VERTP(tri(2),:)-poly%VERTP(tri(1),:), &
               poly%VERTP(tri(3),:)-poly%VERTP(tri(1),:)  )  
       end do
       poly%NORMAL(f,1:3) = face_normal / norm2(face_normal)
    end do    
    
  end subroutine make_sym_hex

  subroutine make_stellated_dodecahedron(poly, a_pts)
    implicit none
    
    type(polyhedron), intent(out) :: poly
    real(C_DOUBLE), intent(in) :: a_pts(3,32)

    ! Call VOFtools function directly
    ! Vertex locations specifically chosen in timing driver to match
    ! the VOFTools description.
    call NCDODECAMESH(poly%IPV, poly%NIPV, poly%NTP, poly%NTS, poly%NTV, &
         poly%VERTP, poly%NORMAL(:,1), poly%NORMAL(:,2), poly%NORMAL(:,3))
    
  end subroutine make_stellated_dodecahedron
  

  subroutine make_stellated_icosahedron(poly, a_pts)
    implicit none
    
    type(polyhedron), intent(out) :: poly
    real(C_DOUBLE), intent(in) :: a_pts(3,32)

    ! Call VOFtools function directly
    ! Vertex locations specifically chosen in timing driver to match
    ! the VOFTools description.
    call NCICOSAMESH(poly%IPV, poly%NIPV, poly%NTP, poly%NTS, poly%NTV, &
         poly%VERTP, poly%NORMAL(:,1), poly%NORMAL(:,2), poly%NORMAL(:,3))
    
  end subroutine make_stellated_icosahedron
  
  
end module VOFtools_wrapper
