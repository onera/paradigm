!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2020  ONERA
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

#include "pdm_configf.h"


module pdm_generate_mesh

  use pdm

  implicit none

  interface PDM_generate_mesh_rectangle_simplified ; module procedure &
  PDM_generate_mesh_rectangle_simplified_
  end interface

  ! interface PDM_generate_mesh_sphere_simplified ; module procedure &
  ! PDM_generate_mesh_sphere_simplified_
  ! end interface

  ! interface PDM_generate_mesh_ball_simplified ; module procedure &
  ! PDM_generate_mesh_ball_simplified_
  ! end interface

  ! interface PDM_generate_mesh_parallelepiped_simplified ; module procedure &
  ! PDM_generate_mesh_parallelepiped_simplified_
  ! end interface

  ! interface PDM_generate_mesh_sphere ; module procedure &
  ! PDM_generate_mesh_sphere_
  ! end interface

  ! interface PDM_generate_mesh_rectangle ; module procedure &
  ! PDM_generate_mesh_rectangle_
  ! end interface

  ! interface PDM_generate_mesh_ball ; module procedure &
  ! PDM_generate_mesh_ball_
  ! end interface

  ! interface PDM_generate_mesh_parallelepiped ; module procedure &
  ! PDM_generate_mesh_parallelepiped_
  ! end interface

  ! interface PDM_generate_mesh_rectangle_ngon ; module procedure &
  ! PDM_generate_mesh_rectangle_ngon_
  ! end interface

  ! interface PDM_generate_mesh_sphere_ngon ; module procedure &
  ! PDM_generate_mesh_sphere_ngon_
  ! end interface

  ! interface PDM_generate_mesh_ball_ngon ; module procedure &
  ! PDM_generate_mesh_ball_ngon_
  ! end interface

  ! interface PDM_generate_mesh_parallelepiped_ngon ; module procedure &
  ! PDM_generate_mesh_parallelepiped_ngon_
  ! end interface


  private :: PDM_generate_mesh_rectangle_simplified_
  ! private :: PDM_generate_mesh_sphere_simplified
  ! private :: PDM_generate_mesh_ball_simplified
  ! private :: PDM_generate_mesh_parallelepiped_simplified
  ! private :: PDM_generate_mesh_sphere
  ! private :: PDM_generate_mesh_rectangle
  ! private :: PDM_generate_mesh_ball
  ! private :: PDM_generate_mesh_parallelepiped
  ! private :: PDM_generate_mesh_rectangle_ngon
  ! private :: PDM_generate_mesh_sphere_ngon
  ! private :: PDM_generate_mesh_ball_ngon
  ! private :: PDM_generate_mesh_parallelepiped_ngon

  interface

  subroutine PDM_generate_mesh_rectangle_simplified_cf(comm,        &
                                                       n_vtx_seg,   &
                                                       n_vtx,       &
                                                       n_elt,       &
                                                       coords,      &
                                                       elt_vtx_idx, &
                                                       elt_vtx)     &

      bind (c, name = 'PDM_generate_mesh_rectangle_simplified')

      use iso_c_binding
      implicit none

      integer(c_int),  value :: comm
#ifdef PDM_LONG_G_NUM
      integer(c_long), value :: n_vtx_seg
#else
      integer(c_int),  value :: n_vtx_seg
#endif
      integer(c_int)         :: n_vtx
      integer(c_int)         :: n_elt
      type(c_ptr)            :: coords
      type(c_ptr)            :: elt_vtx_idx
      type(c_ptr)            :: elt_vtx

  end subroutine PDM_generate_mesh_rectangle_simplified_cf

  subroutine PDM_generate_mesh_sphere_simplified_cf(comm,        &
                                                    n_vtx,       &
                                                    n_elt,       &
                                                    coords,      &
                                                    elt_vtx_idx, &
                                                    elt_vtx)     &

    bind (c, name = 'PDM_generate_mesh_sphere_simplified')

    use iso_c_binding
    implicit none

    integer(c_int), value  :: comm 
    integer(c_int)         :: n_vtx
    integer(c_int)         :: n_elt
    type(c_ptr)            :: coords
    type(c_ptr)            :: elt_vtx_idx
    type(c_ptr)            :: elt_vtx
  end subroutine PDM_generate_mesh_sphere_simplified_cf

  subroutine PDM_generate_mesh_ball_simplified_cf(comm,        &
                                                  n_vtx,       &
                                                  n_elt,       &
                                                  coords,      &
                                                  elt_vtx_idx, &
                                                  elt_vtx)     &

    bind (c, name = 'PDM_generate_mesh_ball_simplified')

    use iso_c_binding
    implicit none

    integer(c_int), value :: comm 
    integer(c_int)        :: n_vtx
    integer(c_int)        :: n_elt
    type(c_ptr)           :: coords
    type(c_ptr)           :: elt_vtx_idx
    type(c_ptr)           :: elt_vtx
  end subroutine PDM_generate_mesh_ball_simplified_cf

  subroutine PDM_generate_mesh_parallelepiped_simplified_cf(comm,        &
                                                            n_vtx_seg,   &
                                                            n_vtx,       &
                                                            n_elt,       &
                                                            coords,      &
                                                            elt_vtx_idx, &
                                                            elt_vtx)     &

    bind (c, name = 'PDM_generate_mesh_parallelepiped_simplified')

    use iso_c_binding
    implicit none

    integer(c_int), value :: comm 
#ifdef PDM_LONG_G_NUM
    integer(c_long), value :: n_vtx_seg
#else
    integer(c_int),  value :: n_vtx_seg
#endif
    integer(c_int)         :: n_vtx
    integer(c_int)         :: n_elt
    type(c_ptr)            :: coords
    type(c_ptr)            :: elt_vtx_idx
    type(c_ptr)            :: elt_vtx
  end subroutine PDM_generate_mesh_parallelepiped_simplified_cf

  function PDM_generate_mesh_sphere_cf(comm,        &
                                       elt_type,    &
                                       order,       &
                                       ho_ordering, &
                                       radius,      &
                                       center_x,    &
                                       center_y,    &
                                       center_z,    &
                                       n_u,         &
                                       n_v,         &
                                       n_part,      &
                                       part_method) &
    result(mesh_nodal) &
    bind (c, name = 'PDM_generate_mesh_sphere')

    use iso_c_binding
    implicit none

    integer(c_int), value :: comm 
    integer(c_int), value :: elt_type
    integer(c_int), value :: order
    type(c_ptr), value    :: ho_ordering
    real(c_double), value :: radius
    real(c_double), value :: center_x
    real(c_double), value :: center_y
    real(c_double), value :: center_z
#ifdef PDM_LONG_G_NUM
    integer(c_long), value :: n_u
    integer(c_long), value :: n_v
#else
    integer(c_int), value :: n_u
    integer(c_int), value :: n_v
#endif
    integer(c_int), value :: n_part
    integer(c_int), value :: part_method

    type(c_ptr)            :: mesh_nodal
  end function PDM_generate_mesh_sphere_cf

  function PDM_generate_mesh_rectangle_cf (comm,        &
                                           elt_type,    &
                                           order,       &
                                           ho_ordering, &
                                           xmin,        &
                                           ymin,        &
                                           zmin,        &
                                           lengthx,     &
                                           lengthy,     &
                                           n_x,         &
                                           n_y,         &
                                           n_part,      &
                                           part_method) &

    result(mesh_nodal) &
    bind (c, name = 'PDM_generate_mesh_rectangle')

    use iso_c_binding
    implicit none

    integer(c_int), value :: comm
    integer(c_int), value :: elt_type
    integer(c_int), value :: order
    type(c_ptr), value    :: ho_ordering
    real(c_double), value :: xmin
    real(c_double), value :: ymin
    real(c_double), value :: zmin
    real(c_double), value :: lengthx
    real(c_double), value :: lengthy
#ifdef PDM_LONG_G_NUM
    integer(c_long), value :: n_x
    integer(c_long), value :: n_y
#else
    integer(c_int), value :: n_x
    integer(c_int), value :: n_y
#endif
    integer(c_int), value :: n_part
    integer(c_int), value :: part_method

    type(c_ptr)            :: mesh_nodal

  end function PDM_generate_mesh_rectangle_cf

  function PDM_generate_mesh_ball_cf (comm,            &
                                        elt_type,        &
                                        order,           &
                                        ho_ordering,     &
                                        radius,          &
                                        hole_radius,     &
                                        center_x,        &
                                        center_y,        &
                                        center_z,        &
                                        n_x,             &
                                        n_y,             &
                                        n_z,             &
                                        n_layer,         &
                                        geometric_ratio, &
                                        n_part,          &
                                        part_method)     &

    result(mesh_nodal) &
    bind (c, name = 'PDM_generate_mesh_ball')

    use iso_c_binding
    implicit none

    integer(c_int), value :: comm
    integer(c_int), value :: elt_type
    integer(c_int), value :: order
    type(c_ptr), value    :: ho_ordering
    real(c_double), value :: radius
    real(c_double), value :: hole_radius
    real(c_double), value :: center_x
    real(c_double), value :: center_y
    real(c_double), value :: center_z
#ifdef PDM_LONG_G_NUM
    integer(c_long), value :: n_x
    integer(c_long), value :: n_y
    integer(c_long), value :: n_z
    integer(c_long), value :: n_layer
#else
    integer(c_int), value :: n_x
    integer(c_int), value :: n_y
    integer(c_int), value :: n_z
    integer(c_int), value :: n_layer
#endif
    real(c_double), value :: geometric_ratio
    integer(c_int), value :: n_part
    integer(c_int), value :: part_method

    type(c_ptr)            :: mesh_nodal

  end function PDM_generate_mesh_ball_cf

  ! subroutine PDM_generate_mesh_parallelepiped_cf
  !   bind (c, name = 'PDM_generate_mesh_parallelepiped')

  !   use iso_c_binding
  !   implicit none

  ! end subroutine PDM_generate_mesh_parallelepiped_cf

  ! subroutine PDM_generate_mesh_rectangle_ngon_cf
  !   bind (c, name = 'PDM_generate_mesh_rectangle_ngon')

  !   use iso_c_binding
  !   implicit none

  !   integer(c_int), value :: comm 
  !   integer(c_int)         :: n_vtx
  !   integer(c_int)         :: n_elt
  !   type(c_ptr)            :: coords
  !   type(c_ptr)            :: elt_vtx_idx
  !   type(c_ptr)            :: elt_vtx
  ! end subroutine PDM_generate_mesh_rectangle_ngon_cf

  ! subroutine PDM_generate_mesh_sphere_ngon_cf
  !   bind (c, name = 'PDM_generate_mesh_sphere_ngon')

  !   use iso_c_binding
  !   implicit none

  !   integer(c_int), value :: comm 
  ! end subroutine PDM_generate_mesh_sphere_ngon_cf

  ! subroutine PDM_generate_mesh_ball_ngon_cf
  !   bind (c, name = 'PDM_generate_mesh_ball_ngon')

  !   use iso_c_binding
  !   implicit none

  !   integer(c_int), value :: comm 
  !   integer(c_int)         :: n_vtx
  !   integer(c_int)         :: n_elt
  !   type(c_ptr)            :: coords
  !   type(c_ptr)            :: elt_vtx_idx
  !   type(c_ptr)            :: elt_vtx
  ! end subroutine PDM_generate_mesh_ball_ngon_cf

  ! subroutine PDM_generate_mesh_parallelepiped_ngon_cf
  !   bind (c, name = 'PDM_generate_mesh_parallelepiped_ngon')

  !   use iso_c_binding
  !   implicit none

  !   integer(c_int), value :: comm 
  !   integer(c_int)         :: n_vtx
  !   integer(c_int)         :: n_elt
  !   type(c_ptr)            :: coords
  !   type(c_ptr)            :: elt_vtx_idx
  !   type(c_ptr)            :: elt_vtx
  ! end subroutine PDM_generate_mesh_parallelepiped_ngon_cf


  end interface

  contains

  !>
  !!
  !! \brief Create a simple partitionned rectangle mesh (2D).
  !!
  !! \param [in]   comm        MPI communicator
  !! \param [in]   n_vtx_seg   Number of vertices along each side of the rectangle
  !! \param [out]  n_vtx       Number of vertices
  !! \param [out]  n_elt       Number of elements
  !! \param [out]  coords      Array of vertex coordinates
  !! \param [out]  elt_vtx_idx Index array of the element vertex connectivity
  !! \param [out]  elt_vtx     Array of the element vertex connectivity
  !!
  !!

  subroutine PDM_generate_mesh_rectangle_simplified_(comm,        &
                                                     n_vtx_seg,   &
                                                     n_vtx,       &
                                                     n_elt,       &
                                                     coords,      &
                                                     elt_vtx_idx, &
                                                     elt_vtx)

      use iso_c_binding
      implicit none

      integer,                     intent(in) :: comm
      integer(kind=pdm_g_num_s),   intent(in) :: n_vtx_seg
      integer,                    intent(out) :: n_vtx
      integer,                    intent(out) :: n_elt
      double precision,               pointer :: coords(:,:)
      integer(kind=pdm_l_num_s),      pointer :: elt_vtx_idx(:)
      integer(kind=pdm_l_num_s),      pointer :: elt_vtx(:)

      integer(c_int)                          :: c_comm
      integer(c_int)                          :: c_n_vtx
      integer(c_int)                          :: c_n_elt
      type(c_ptr)                             :: c_coords
      type(c_ptr)                             :: c_elt_vtx_idx
      type(c_ptr)                             :: c_elt_vtx

      c_comm = PDM_MPI_Comm_f2c(comm)

      call PDM_generate_mesh_rectangle_simplified_cf(c_comm,        &
                                                     n_vtx_seg,     &
                                                     c_n_vtx,       &
                                                     c_n_elt,       &
                                                     c_coords,      &
                                                     c_elt_vtx_idx, &
                                                     c_elt_vtx)

      n_vtx = c_n_vtx
      n_elt = c_n_elt

      call c_f_pointer(c_coords, &
                       coords,   &
                       [3, n_vtx])

      call c_f_pointer(c_elt_vtx_idx, &
                       elt_vtx_idx,   &
                       [n_elt + 1])

      call c_f_pointer(c_elt_vtx, &
                       elt_vtx,   &
                       [elt_vtx_idx(n_elt + 1)])

  end subroutine PDM_generate_mesh_rectangle_simplified_

end module pdm_generate_mesh
