!-----------------------------------------------------------------------------
! This file is part of the CWIPI library.
!
! Copyright (C) 2011  ONERA
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 3 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library. If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------

program testf

  use pdm
  use mpi
  use pdm_mesh_location
  use iso_c_binding


  implicit none

  !
  ! About point cloud
  !

  integer, parameter :: n_point_cloud = 1
  integer, parameter :: n_part_cloud = 1
  integer, parameter :: i_point_cloud = 0
  integer, parameter :: n_points_into_cloud = 2

  double precision, pointer :: coords_cloud(:) ! pointer or allocatble, target
  integer (kind = pdm_g_num_s), pointer :: gnum_cloud(:) ! pointer or allocatble, target

  type(c_ptr)  :: cptr_coords_cloud
  type(c_ptr)  :: cptr_gnum_cloud

  !
  ! About mesh
  !

  integer, parameter :: n_part_mesh = 1

  integer, parameter :: n_cell = 2
  integer (c_int), pointer :: cell_face_idx
  integer (c_int), pointer :: cell_face
  integer (kind = pdm_g_num_s), pointer :: cell_ln_to_gn
  type(c_ptr)  :: cptr_cell_face_idx
  type(c_ptr)  :: cptr_cell_face
  type(c_ptr)  :: cptr_cell_ln_to_gn

  integer, parameter :: n_face = 11
  integer (c_int), pointer :: face_vtx_idx
  integer (c_int), pointer :: face_vtx
  integer (kind = pdm_g_num_s), pointer :: face_ln_to_gn
  type(c_ptr)  :: cptr_face_vtx_idx
  type(c_ptr)  :: cptr_face_vtx
  type(c_ptr)  :: cptr_face_ln_to_gn

  integer, parameter :: n_vtx = 12
  double precision, pointer :: coords_vtx(:) ! pointer or allocatble, target
  integer (kind = pdm_g_num_s), pointer :: gnum_vtx(:) ! pointer or allocatble, target
  type(c_ptr)  :: cptr_coords_vtx
  type(c_ptr)  :: cptr_gnum_vtx

  integer :: i

  integer :: id

  integer, parameter :: partial = 0 ! Put 1 to keep results when the subroutine closest_points_free is

  !
  ! Set point cloud
  !

  allocate(coords_cloud(3*n_points_into_cloud))
  allocate(gnum_cloud(n_points_into_cloud))

  coords_cloud(1) = 0.5
  coords_cloud(2) = 0.5
  coords_cloud(3) = 0.5

  coords_cloud(4) = 1.5
  coords_cloud(5) = 0.5
  coords_cloud(6) = 0.5

  do i = 1, n_points_into_cloud
    gnum_cloud(i) = i
  end do

  cptr_coords_cloud = c_loc (coords_cloud)
  cptr_gnum_cloud = c_loc (gnum_cloud)

  !
  ! Set mesh
  !

  allocate(coords_vtx(3*n_vtx))

  coords_vtx(1) = 0.d0
  coords_vtx(2) = 0.d0
  coords_vtx(3) = 0.d0

  coords_vtx(4) = 1.d0
  coords_vtx(5) = 0.d0
  coords_vtx(6) = 0.d0

  coords_vtx(7) = 1.d0
  coords_vtx(8) = 1.d0
  coords_vtx(9) = 0.d0

  coords_vtx(10) = 0.d0
  coords_vtx(11) = 1.d0
  coords_vtx(12) = 0.d0

  coords_vtx(13) = 0.d0
  coords_vtx(14) = 0.d0
  coords_vtx(15) = 1.d0

  coords_vtx(16) = 1.d0
  coords_vtx(17) = 0.d0
  coords_vtx(18) = 1.d0

  coords_vtx(19) = 1.d0
  coords_vtx(20) = 1.d0
  coords_vtx(21) = 1.d0

  coords_vtx(22) = 0.d0
  coords_vtx(23) = 1.d0
  coords_vtx(24) = 1.d0

  coords_vtx(25) = 2.d0
  coords_vtx(26) = 0.d0
  coords_vtx(27) = 0.d0

  coords_vtx(28) = 2.d0
  coords_vtx(29) = 1.d0
  coords_vtx(30) = 0.d0

  coords_vtx(31) = 2.d0
  coords_vtx(32) = 0.d0
  coords_vtx(33) = 1.d0

  coords_vtx(34) = 2.d0
  coords_vtx(35) = 1.d0
  coords_vtx(36) = 1.d0

  allocate(gnum_vtx(3*n_vtx))

  do i = 1, n_vtx
    gnum_vtx(i) = i
  end do

  cptr_coords_vtx = c_loc (coords_vtx)
  cptr_gnum_vtx = c_loc (gnum_vtx)


  allocate(face_vtx_idx(n_face+1))
  allocate(face_vtx(4*n_face))

  face_vtx_idx(0) = 0
  do i = 1, n_face
    face_vtx_idx(i+1) = face_vtx_idx(i) + 4
  end do

  face_vtx(1) = 1
  face_vtx(2) = 4
  face_vtx(3) = 3
  face_vtx(4) = 2

  face_vtx(5) = 5
  face_vtx(6) = 6
  face_vtx(7) = 7
  face_vtx(8) = 8

  face_vtx(9) = 1
  face_vtx(10) = 5
  face_vtx(11) = 8
  face_vtx(12) = 4

  face_vtx(13) = 6 ! face milieu
  face_vtx(14) = 2
  face_vtx(15) = 3
  face_vtx(16) = 7

  face_vtx(17) = 3
  face_vtx(18) = 4
  face_vtx(19) = 8
  face_vtx(20) = 7

  face_vtx(21) = 1
  face_vtx(22) = 2
  face_vtx(23) = 6
  face_vtx(24) = 5

  face_vtx(25) = 11
  face_vtx(26) = 12
  face_vtx(27) = 7
  face_vtx(28) = 6

  face_vtx(29) = 9
  face_vtx(30) = 2
  face_vtx(31) = 3
  face_vtx(32) = 10

  face_vtx(33) = 12
  face_vtx(34) = 10
  face_vtx(35) = 3
  face_vtx(36) = 7

  face_vtx(37) = 11
  face_vtx(38) = 6
  face_vtx(39) = 2
  face_vtx(40) = 9

  face_vtx(41) = 9
  face_vtx(42) = 10
  face_vtx(43) = 12
  face_vtx(44) = 11

  allocate(gnum_face(n_face))

  do i = 1, n_face
    gnum_face(i) = i
  end do

  cptr_face_vtx_idx = c_loc (face_vtx_idx)
  cptr_face_vtx = c_loc (face_vtx)
  cptr_gnum_face = c_loc (gnum_face)


  allocate(cell_face_idx(n_cell+1))
  allocate(cell_face(6*n_cell))

  cell_face_idx(0) = 0
  do i = 1, n_cell
    cell_face_idx(i+1) = cell_face_idx(i) + 6
  end do

  cell_face(1) = 1
  cell_face(2) = 2
  cell_face(3) = 3
  cell_face(4) = 4
  cell_face(5) = 5
  cell_face(6) = 6

  cell_face(7) = 7
  cell_face(8) = 8
  cell_face(9) = 9
  cell_face(10) = 10
  cell_face(11) = 11
  cell_face(12) = -4

  allocate(gnum_cell(n_cell))

  do i = 1, n_cell
    gnum_cell(i) = i
  end do

  cptr_cell_face_idx = c_loc (cell_face_idx)
  cptr_cell_face = c_loc (cell_face)
  cptr_gnum_cell = c_loc (gnum_cell)






  integer (kind = pdm_g_num_s), parameter :: n_g_points_src = 10
  integer (kind = pdm_g_num_s), parameter :: n_g_points_tgt = 10

  integer, parameter :: n_part_cloud_src = 1
  integer, parameter :: n_part_cloud_tgt = 1




  integer, parameter :: n_local_points_src = 5
  integer, parameter :: n_local_points_tgt = 5

  integer, parameter :: n_closest = 2

  double precision, pointer :: coords_src(:) ! pointer or allocatble, target
  double precision, pointer :: coords_tgt(:) ! pointer or allocatble, target

  type(c_ptr), pointer    :: cptr_coords_src(:)
  type(c_ptr), pointer    :: cptr_coords_tgt(:)

  integer :: code
  integer :: i_rank
  integer :: n_rank

  integer (kind = pdm_g_num_s), pointer :: gnum_src(:) ! pointer or allocatble, target
  integer (kind = pdm_g_num_s), pointer :: gnum_tgt(:) ! pointer or allocatble, target

  type(c_ptr), pointer    :: cptr_gnum_src(:)
  type(c_ptr), pointer    :: cptr_gnum_tgt(:)

  type(c_ptr)     :: cptr_closest_src_gnum
  type(c_ptr)     :: cptr_closest_src_distance

  integer (kind = pdm_g_num_s), pointer :: closest_src_gnum(:)
  double precision, pointer :: closest_src_distance(:)

  integer :: i

  integer :: id

  integer, parameter :: partial = 0 ! Put 1 to keep results when the subroutine closest_points_free is called


  call mpi_init(code)
  call mpi_comm_rank(mpi_comm_world, i_rank, code)
  call mpi_comm_size(mpi_comm_world, n_rank, code)

  if (n_rank .ne. 1) then
    print *,'Error : 1 MPI process is mandatory'
    call mpi_finalize(code)
    stop
  end if

  allocate (coords_src(3*n_local_points_src))
  allocate (coords_tgt(3*n_local_points_tgt))
  allocate (gnum_src(n_local_points_src))
  allocate (gnum_tgt(n_local_points_tgt))

  if (i_rank .eq. 0) then
    do i = 1, n_local_points_src
      coords_src(3*(i-1)+1) = 0. + i - 1
      coords_src(3*(i-1)+2) = 0. + i - 1
      coords_src(3*(i-1)+3) = 0. + i - 1
      gnum_src(i) = 2*(i-1) + 1
    end do
    do i = 1, n_local_points_tgt
      coords_tgt(3*(i-1)+1) = n_local_points_tgt + i - 1
      coords_tgt(3*(i-1)+2) = n_local_points_tgt + i - 1
      coords_tgt(3*(i-1)+3) = n_local_points_tgt + i - 1
      gnum_tgt(i) = 2*(i-1) + 2
    end do
  else
    do i = 1, n_local_points_src
      coords_src(3*(i-1)+1) = n_local_points_tgt + i - 1
      coords_src(3*(i-1)+2) = n_local_points_tgt + i - 1
      coords_src(3*(i-1)+3) = n_local_points_tgt + i - 1
      gnum_src(i) = 2*(i-1) + 2
    end do
    do i = 1, n_local_points_tgt
      coords_tgt(3*(i-1)+1) = 0. + i - 1
      coords_tgt(3*(i-1)+2) = 0. + i - 1
      coords_tgt(3*(i-1)+3) = 0. + i - 1
      gnum_tgt(i) = 2*(i-1) + 1
    end do
  endif


  do i = 1, n_local_points_tgt
    coords_tgt(3*(i-1)+1) = coords_tgt(3*(i-1)+1) / 10.
    coords_tgt(3*(i-1)+2) = coords_tgt(3*(i-1)+2) / 10.
    coords_tgt(3*(i-1)+3) = coords_tgt(3*(i-1)+3) / 10.
  enddo

  do i = 1, n_local_points_src
    coords_src(3*(i-1)+1) = coords_src(3*(i-1)+1) / 10.
    coords_src(3*(i-1)+2) = coords_src(3*(i-1)+2) / 10.
    coords_src(3*(i-1)+3) = coords_src(3*(i-1)+3) / 10.
  enddo

  !
  ! Create a new PDM_mesh_location structure
  !   The MPI communicator and the number of point cloud are setted
  !


  call PDM_mesh_location_create (PDM_MESH_NATURE_MESH_SETTED, &
                                 n_point_cloud, &
                                 MPI_COMM_WORLD, &
                                 id)

  !
  ! Set the local number partition for any point cloud
  !

  call PDM_mesh_location_n_part_cloud_set (id, &
                                           i_point_cloud, &
                                           n_part_cloud) &


  !
  ! Set point properties for any partition of point cloud
  !

  call PDM_mesh_location_cloud_set (id, &
                                    i_point_cloud, &
                                    i_part, &
                                    n_points, &
                                    coords, &
                                    gnum) &


  !
  ! Set the point coordinates and the global numbering for any partition of the source
  !

  allocate (cptr_gnum_src(n_part_cloud_src))
  allocate (cptr_coords_src(n_part_cloud_src))

  if (n_part_cloud_src .ne. 1) then
    print *, "For this test, n_part_cloud_src must be equal to 1"
    stop
  end if

  do i = 1, n_part_cloud_src
    cptr_gnum_src = c_loc(gnum_src)
    cptr_coords_src = c_loc(coords_src)
    call PDM_closest_points_src_cloud_set (id, &
                                           i-1, & !!! ipart : 0 -> n_part-1 !!!
                                           n_local_points_src, &
                                           cptr_coords_src(i), &
                                           cptr_gnum_src(i))
  end do

  !
  ! Set the point coordinates and the global numbering for any partition of the target
  !

  allocate (cptr_gnum_tgt(n_part_cloud_tgt))
  allocate (cptr_coords_tgt(n_part_cloud_tgt))

  if (n_part_cloud_tgt .ne. 1) then
    print *, "For this test, n_part_cloud_tgt must be equal to 1"
    stop
  end if

  do i = 1, n_part_cloud_tgt
    cptr_gnum_tgt(i) = c_loc(gnum_tgt)
    cptr_coords_tgt(i) = c_loc(coords_tgt)
    call PDM_closest_points_tgt_cloud_set (id, &
                                           i-1, &  !!! ipart : 0 -> n_part-1 !!!
                                           n_local_points_tgt, &
                                           cptr_coords_tgt(i), &
                                           cptr_gnum_tgt(i))
  end do


  !
  ! Compute the 'n' closest neighbors into the source point cloud for any taget point
  !

  call PDM_closest_points_compute (id)

  !
  ! Dump the time used to compute
  !

  call PDM_closest_points_dump_times (id)

  !
  ! Get the 'n' closest neighbors into the source point cloud for any taget point
  !


  do i = 1, n_part_cloud_tgt
    call PDM_closest_points_get (id, &
                                 i-1, & !!! ipart : 0 -> n_part-1 !!!
                                 cptr_closest_src_gnum, &
                                 cptr_closest_src_distance)
    call c_f_pointer(cptr_closest_src_gnum, closest_src_gnum, [n_closest * n_local_points_tgt])
    call c_f_pointer(cptr_closest_src_distance, closest_src_distance, [n_closest * n_local_points_tgt])
  end do


  do i = 1, n_local_points_tgt
    print *, "Closest points for ", gnum_tgt(i), ":", closest_src_gnum(2*(i-1)+1),"/", closest_src_distance(2*(i-1)+1), " and ", closest_src_gnum(2*(i-1)+2),"/", closest_src_distance(2*(i-1)+2)

  end do

  !
  ! Free the current cloest_point structure
  !

  call PDM_closest_points_free (id, partial)

  deallocate (coords_src)
  deallocate (coords_tgt)
  deallocate (gnum_src)
  deallocate (gnum_tgt)
  deallocate (cptr_gnum_src)
  deallocate (cptr_coords_src)
  deallocate (cptr_gnum_tgt)
  deallocate (cptr_coords_tgt)


  call mpi_finalize(code)


end program testf
