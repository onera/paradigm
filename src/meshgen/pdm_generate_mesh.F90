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
  use pdm_pointer_array

  implicit none

  interface PDM_generate_mesh_rectangle_simplified ; module procedure &
  PDM_generate_mesh_rectangle_simplified_
  end interface

  interface PDM_generate_mesh_sphere_simplified ; module procedure &
  PDM_generate_mesh_sphere_simplified_
  end interface

  interface PDM_generate_mesh_ball_simplified ; module procedure &
  PDM_generate_mesh_ball_simplified_
  end interface

  interface PDM_generate_mesh_parallelepiped_simplified ; module procedure &
  PDM_generate_mesh_parallelepiped_simplified_
  end interface

  interface PDM_generate_mesh_sphere ; module procedure &
  PDM_generate_mesh_sphere_
  end interface

  interface PDM_generate_mesh_rectangle ; module procedure &
  PDM_generate_mesh_rectangle_
  end interface

  interface PDM_generate_mesh_ball ; module procedure &
  PDM_generate_mesh_ball_
  end interface

  interface PDM_generate_mesh_parallelepiped ; module procedure &
  PDM_generate_mesh_parallelepiped_
  end interface

  interface PDM_generate_mesh_rectangle_ngon ; module procedure &
  PDM_generate_mesh_rectangle_ngon_
  end interface

  interface PDM_generate_mesh_sphere_ngon ; module procedure &
  PDM_generate_mesh_sphere_ngon_
  end interface

  interface PDM_generate_mesh_ball_ngon ; module procedure &
  PDM_generate_mesh_ball_ngon_
  end interface

  interface PDM_generate_mesh_parallelepiped_ngon ; module procedure &
  PDM_generate_mesh_parallelepiped_ngon_
  end interface


  private :: PDM_generate_mesh_rectangle_simplified_
  private :: PDM_generate_mesh_sphere_simplified_
  private :: PDM_generate_mesh_ball_simplified_
  private :: PDM_generate_mesh_parallelepiped_simplified_
  private :: PDM_generate_mesh_sphere_
  private :: PDM_generate_mesh_rectangle_
  private :: PDM_generate_mesh_ball_
  private :: PDM_generate_mesh_parallelepiped_
  private :: PDM_generate_mesh_rectangle_ngon_
  private :: PDM_generate_mesh_sphere_ngon_
  private :: PDM_generate_mesh_ball_ngon_
  private :: PDM_generate_mesh_parallelepiped_ngon_

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

  function PDM_generate_mesh_parallelepiped_cf (comm,         &
                                                elt_type,     & 
                                                order,        &
                                                ho_ordering,  &
                                                xmin,         &
                                                ymin,         &
                                                zmin,         &
                                                lengthx,      &
                                                lengthy,      &
                                                lengthz,      &
                                                n_x,          &
                                                n_y,          &
                                                n_z,          &
                                                n_part,       &
                                                part_method)  &

    result(mesh_nodal) &
    bind (c, name = 'PDM_generate_mesh_parallelepiped')
 
    use iso_c_binding
    implicit none
! 
    integer(c_int), value :: comm
    integer(c_int), value :: elt_type
    integer(c_int), value :: order
    type(c_ptr), value    :: ho_ordering
    real(c_double), value :: xmin
    real(c_double), value :: ymin
    real(c_double), value :: zmin
    real(c_double), value :: lengthx
    real(c_double), value :: lengthy
    real(c_double), value :: lengthz
#ifdef PDM_LONG_G_NUM
    integer(c_long), value :: n_x
    integer(c_long), value :: n_y
    integer(c_long), value :: n_z
#else
    integer(c_int), value :: n_x
    integer(c_int), value :: n_y
    integer(c_int), value :: n_z
#endif
    integer(c_int), value :: n_part
    integer(c_int), value :: part_method

    type(c_ptr)           :: mesh_nodal

  end function PDM_generate_mesh_parallelepiped_cf
! 
  subroutine PDM_generate_mesh_rectangle_ngon_cf (comm, &
                                                  elt_type, &
                                                  xmin, &
                                                  ymin, &
                                                  zmin, &
                                                  lengthx, &
                                                  lengthy, &
                                                  n_x, &
                                                  n_y, &
                                                  n_part, &
                                                  part_method, &
                                                  pn_vtx, &
                                                  pn_edge, &
                                                  pn_face, &
                                                  pvtx_coord, &
                                                  pedge_vtx, &
                                                  pface_edge_idx, &
                                                  pface_edge, &
                                                  pface_vtx, &
                                                  pvtx_ln_to_gn, &
                                                  pedge_ln_to_gn, &
                                                  pface_ln_to_gn) &

    bind (c, name = 'PDM_generate_mesh_rectangle_ngon')
! 
    use iso_c_binding
    implicit none
! 
    integer(c_int), value :: comm
    integer(c_int), value :: elt_type
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

    type(c_ptr) :: pn_vtx
    type(c_ptr) :: pn_edge
    type(c_ptr) :: pn_face
    type(c_ptr) :: pvtx_coord
    type(c_ptr) :: pedge_vtx
    type(c_ptr) :: pface_edge_idx
    type(c_ptr) :: pface_edge
    type(c_ptr) :: pface_vtx
    type(c_ptr) :: pvtx_ln_to_gn
    type(c_ptr) :: pedge_ln_to_gn
    type(c_ptr) :: pface_ln_to_gn


  end subroutine PDM_generate_mesh_rectangle_ngon_cf

  subroutine PDM_generate_mesh_sphere_ngon_cf (comm, &
                                               elt_type, &
                                               order, &
                                               ho_ordering, &
                                               radius, &
                                               center_x, &
                                               center_y, &
                                               center_z, &
                                               n_u, &
                                               n_v, &
                                               n_part, &
                                               part_method, &
                                               pn_vtx, &
                                               pn_edge, &
                                               pn_face, &
                                               pvtx_coord, &
                                               pedge_vtx, &
                                               pface_edge_idx, &
                                               pface_edge, &
                                               pface_vtx, &
                                               pvtx_ln_to_gn, &
                                               pedge_ln_to_gn, &
                                               pface_ln_to_gn) &

    bind (c, name = 'PDM_generate_mesh_sphere_ngon')
! 
    use iso_c_binding
    implicit none
! 
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

    type(c_ptr) :: pn_vtx
    type(c_ptr) :: pn_edge
    type(c_ptr) :: pn_face
    type(c_ptr) :: pvtx_coord
    type(c_ptr) :: pedge_vtx
    type(c_ptr) :: pface_vtx
    type(c_ptr) :: pface_edge_idx
    type(c_ptr) :: pface_edge
    type(c_ptr) :: pvtx_ln_to_gn
    type(c_ptr) :: pedge_ln_to_gn
    type(c_ptr) :: pface_ln_to_gn

  end subroutine PDM_generate_mesh_sphere_ngon_cf

  subroutine PDM_generate_mesh_ball_ngon_cf (comm, &
                                             elt_type, &
                                             order, &
                                             ho_ordering, &
                                             radius, &
                                             hole_radius, &
                                             center_x, &
                                             center_y, &
                                             center_z, &
                                             n_x, &
                                             n_y, &
                                             n_z, &
                                             n_layer, &
                                             geometric_ratio, &
                                             n_part, &
                                             part_method, &
                                             pn_vtx, &
                                             pn_edge, &
                                             pn_face, &
                                             pn_cell, &
                                             pvtx_coord, &
                                             pedge_vtx, &
                                             pface_edge_idx, &
                                             pface_edge, &
                                             pface_vtx, &
                                             pcell_face_idx, &
                                             pcell_face, &
                                             pvtx_ln_to_gn, &
                                             pedge_ln_to_gn, &
                                             pface_ln_to_gn, &
                                             pcell_ln_to_gn, &
                                             pn_surface, &
                                             psurface_face_idx, &
                                             psurface_face, &
                                             psurface_face_ln_to_gn) &

    bind (c, name = 'PDM_generate_mesh_ball_ngon')

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

    type(c_ptr) :: pn_vtx
    type(c_ptr) :: pn_edge
    type(c_ptr) :: pn_face
    type(c_ptr) :: pn_cell
    type(c_ptr) :: pvtx_coord
    type(c_ptr) :: pedge_vtx
    type(c_ptr) :: pface_edge_idx
    type(c_ptr) :: pface_edge
    type(c_ptr) :: pface_vtx
    type(c_ptr) :: pcell_face_idx
    type(c_ptr) :: pcell_face
    type(c_ptr) :: pvtx_ln_to_gn
    type(c_ptr) :: pedge_ln_to_gn
    type(c_ptr) :: pface_ln_to_gn
    type(c_ptr) :: pcell_ln_to_gn
    type(c_ptr) :: pn_surface
    type(c_ptr) :: psurface_face_idx
    type(c_ptr) :: psurface_face
    type(c_ptr) :: psurface_face_ln_to_gn
  end subroutine PDM_generate_mesh_ball_ngon_cf

  subroutine PDM_generate_mesh_parallelepiped_ngon_cf (comm, &
                                                       elt_type, &
                                                       order, &
                                                       ho_ordering, &
                                                       xmin, &
                                                       ymin, &
                                                       zmin, &
                                                       lengthx, &
                                                       lengthy, &
                                                       lengthz, &
                                                       n_x, &
                                                       n_y, &
                                                       n_z, &
                                                       n_part, &
                                                       part_method, &
                                                       pn_vtx, &
                                                       pn_edge, &
                                                       pn_face, &
                                                       pn_cell, &
                                                       pvtx_coord, &
                                                       pedge_vtx, &
                                                       pface_edge_idx, &
                                                       pface_edge, &
                                                       pface_vtx, &
                                                       pcell_face_idx, &
                                                       pcell_face, &
                                                       pvtx_ln_to_gn, &
                                                       pedge_ln_to_gn, &
                                                       pface_ln_to_gn, &
                                                       pcell_ln_to_gn, &
                                                       pn_surface, &
                                                       psurface_face_idx, &
                                                       psurface_face, &
                                                       psurface_face_ln_to_gn, &
                                                       pn_ridge, &
                                                       pridge_edge_idx, &
                                                       pridge_edge, &
                                                       pridge_edge_ln_to_gn) &

    bind (c, name = 'PDM_generate_mesh_parallelepiped_ngon')

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
    real(c_double), value :: lengthz
#ifdef PDM_LONG_G_NUM
    integer(c_long), value :: n_x
    integer(c_long), value :: n_y
    integer(c_long), value :: n_z
#else
    integer(c_int), value :: n_x
    integer(c_int), value :: n_y
    integer(c_int), value :: n_z
#endif
    integer(c_int), value :: n_part
    integer(c_int), value :: part_method

    type(c_ptr) :: pn_vtx
    type(c_ptr) :: pn_edge
    type(c_ptr) :: pn_face
    type(c_ptr) :: pn_cell
    type(c_ptr) :: pvtx_coord
    type(c_ptr) :: pedge_vtx
    type(c_ptr) :: pface_edge_idx
    type(c_ptr) :: pface_edge
    type(c_ptr) :: pface_vtx
    type(c_ptr) :: pcell_face_idx
    type(c_ptr) :: pcell_face
    type(c_ptr) :: pvtx_ln_to_gn
    type(c_ptr) :: pedge_ln_to_gn
    type(c_ptr) :: pface_ln_to_gn
    type(c_ptr) :: pcell_ln_to_gn
    type(c_ptr) :: pn_surface
    type(c_ptr) :: psurface_face_idx
    type(c_ptr) :: psurface_face
    type(c_ptr) :: psurface_face_ln_to_gn
    type(c_ptr) :: pn_ridge
    type(c_ptr) :: pridge_edge_idx
    type(c_ptr) :: pridge_edge
    type(c_ptr) :: pridge_edge_ln_to_gn

  end subroutine PDM_generate_mesh_parallelepiped_ngon_cf


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
    integer(c_int)                          :: c_n_vtx       = 0  
    integer(c_int)                          :: c_n_elt       = 0
    type(c_ptr)                             :: c_coords      = C_NULL_PTR
    type(c_ptr)                             :: c_elt_vtx_idx = C_NULL_PTR
    type(c_ptr)                             :: c_elt_vtx     = C_NULL_PTR

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

!>
!!
!! \brief Create a simple partitionned sphere mesh (2D).
!!
!! \param [in]   comm        MPI communicator
!! \param [out]  n_vtx       Number of vertices
!! \param [out]  n_elt       Number of elements
!! \param [out]  coords      Array of vertex coordinates
!! \param [out]  elt_vtx_idx Index array of the element vertex connectivity
!! \param [out]  elt_vtx     Array of the element vertex connectivity
!!
!!

  subroutine PDM_generate_mesh_sphere_simplified_(comm,        &
                                                  n_vtx,       &
                                                  n_elt,       &
                                                  coords,      &
                                                  elt_vtx_idx, &
                                                  elt_vtx)     
    use iso_c_binding
    implicit none

    integer, intent(in)                     :: comm 
    integer, intent(out)                    :: n_vtx
    integer, intent(out)                    :: n_elt
    double precision,               pointer :: coords(:,:)
    integer(kind=pdm_l_num_s),      pointer :: elt_vtx_idx(:)
    integer(kind=pdm_l_num_s),      pointer :: elt_vtx(:)

    integer(c_int)         :: c_comm 
    integer(c_int)         :: c_n_vtx       = 0
    integer(c_int)         :: c_n_elt       = 0
    type(c_ptr)            :: c_coords      = C_NULL_PTR
    type(c_ptr)            :: c_elt_vtx_idx = C_NULL_PTR
    type(c_ptr)            :: c_elt_vtx     = C_NULL_PTR

    c_comm = PDM_MPI_Comm_f2c(comm)

    call PDM_generate_mesh_sphere_simplified_cf(c_comm,        &
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


  end subroutine PDM_generate_mesh_sphere_simplified_

 !>
 !!
 !! \brief Create a simple partitionned ball mesh (3D).
 !!
 !! \param [in]   comm        MPI communicator
 !! \param [out]  n_vtx       Number of vertices
 !! \param [out]  n_elt       Number of elements
 !! \param [out]  coords      Array of vertex coordinates
 !! \param [out]  elt_vtx_idx Index array of the element vertex connectivity
 !! \param [out]  elt_vtx     Array of the element vertex connectivity
 !!
 !!

  subroutine PDM_generate_mesh_ball_simplified_(comm,        &
                                                n_vtx,       &
                                                n_elt,       &
                                                coords,      &
                                                elt_vtx_idx, &
                                                elt_vtx)     
    use iso_c_binding
    implicit none

    integer, intent(in)                     :: comm 
    integer, intent(out)                    :: n_vtx
    integer, intent(out)                    :: n_elt
    double precision,               pointer :: coords(:,:)
    integer(kind=pdm_l_num_s),      pointer :: elt_vtx_idx(:)
    integer(kind=pdm_l_num_s),      pointer :: elt_vtx(:)

    integer(c_int)         :: c_comm 
    integer(c_int)         :: c_n_vtx       = 0 
    integer(c_int)         :: c_n_elt       = 0
    type(c_ptr)            :: c_coords      = C_NULL_PTR
    type(c_ptr)            :: c_elt_vtx_idx = C_NULL_PTR
    type(c_ptr)            :: c_elt_vtx     = C_NULL_PTR

    c_comm = PDM_MPI_Comm_f2c(comm)

    call PDM_generate_mesh_ball_simplified_cf(c_comm,        &
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


  end subroutine PDM_generate_mesh_ball_simplified_

!>
!!
!! \brief Create a simple partitionned parallelepiped mesh (3D).
!!
!! \param [in]   comm        MPI communicator
!! \param [in]   n_vtx_seg   Number of vertices along each side of the parallelepiped
!! \param [out]  n_vtx       Number of vertices
!! \param [out]  n_elt       Number of elements
!! \param [out]  coords      Array of vertex coordinates
!! \param [out]  elt_vtx_idx Index array of the element vertex connectivity
!! \param [out]  elt_vtx     Array of the element vertex connectivity
!!
!!

  subroutine PDM_generate_mesh_parallelepiped_simplified_(comm,        &
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

    call PDM_generate_mesh_parallelepiped_simplified_cf(c_comm,        &
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

  end subroutine PDM_generate_mesh_parallelepiped_simplified_

!>
!!
!! \brief Create a partitionned sphere mesh (2D).
!!
!! \param [in]  comm        MPI communicator
!! \param [in]  elt_type    Mesh element type
!! \param [in]  order       Mesh element order
!! \param [in]  ho_ordering High order nodes ordering type
!! \param [in]  radius      Radius of the sphere
!! \param [in]  center_x    x-coordinate of the sphere center
!! \param [in]  center_y    y-coordinate of the sphere center
!! \param [in]  center_z    z-coordinate of the sphere center
!! \param [in]  n_u         Number of points in longitude
!! \param [in]  n_v         Number of points in latitude
!! \param [in]  n_part      Number of mesh partitions
!! \param [in]  part_method Mesh partitionning method
!!
!! \return PDM_part_mesh_nodal_t
!!
!!

  function PDM_generate_mesh_sphere_(comm,        &
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
    result(mesh_nodal) 

    use iso_c_binding
    implicit none

    integer, intent(in)                   :: comm 
    integer, intent(in)                   :: elt_type
    integer, intent(in)                   :: order
    type(c_ptr), intent(in)               :: ho_ordering
    double precision, intent(in)          :: radius
    double precision, intent(in)          :: center_x
    double precision, intent(in)          :: center_y
    double precision, intent(in)          :: center_z
    integer(kind=pdm_g_num_s), intent(in) :: n_u
    integer(kind=pdm_g_num_s), intent(in) :: n_v
    integer, intent(in)                   :: n_part
    integer, intent(in)                   :: part_method

    integer(c_int)        :: c_comm 
    integer(c_int)        :: c_elt_type
    integer(c_int)        :: c_order
    type(c_ptr)           :: c_ho_ordering
    real(c_double)        :: c_radius
    real(c_double)        :: c_center_x
    real(c_double)        :: c_center_y
    real(c_double)        :: c_center_z
#ifdef PDM_LONG_G_NUM
    integer(c_long)       :: c_n_u
    integer(c_long)       :: c_n_v
#else
    integer(c_int)        :: c_n_u
    integer(c_int)        :: c_n_v
#endif
    integer(c_int)        :: c_n_part
    integer(c_int)        :: c_part_method

    type(c_ptr)            :: mesh_nodal

    c_comm = PDM_MPI_Comm_f2c(comm)

    c_elt_type    = elt_type
    c_order       = order
    c_ho_ordering = ho_ordering
    c_radius      = radius
    c_center_x    = center_x
    c_center_y    = center_y
    c_center_z    = center_z
    c_n_u         = n_u
    c_n_v         = n_v
    c_n_part      = n_part
    c_part_method = part_method

    mesh_nodal = PDM_generate_mesh_sphere_cf(c_comm,        &
                                             c_elt_type,    &
                                             c_order,       &
                                             c_ho_ordering, &
                                             c_radius,      &
                                             c_center_x,    &
                                             c_center_y,    &
                                             c_center_z,    &
                                             c_n_u,         &
                                             c_n_v,         &
                                             c_n_part,      &
                                             c_part_method) 
 
  end function PDM_generate_mesh_sphere_

!>
!!
!! \brief Create a partitionned rectangle mesh (2D).
!!
!! \param [in]  comm        MPI communicator
!! \param [in]  elt_type    Mesh element type
!! \param [in]  order       Mesh element order
!! \param [in]  ho_ordering High order nodes ordering type
!! \param [in]  xmin        x-coordinate of the rctangle minimum corner
!! \param [in]  ymin        y-coordinate of the rctangle minimum corner
!! \param [in]  zmin        z-coordinate of the rctangle minimum corner
!! \param [in]  lengthx     Length of the rectangle in the x-direction
!! \param [in]  lengthy     Length of the rectangle in the y-direction
!! \param [in]  n_x         Number of points in the x-direction
!! \param [in]  n_y         Number of points in the y-direction
!! \param [in]  n_part      Number of mesh partitions
!! \param [in]  part_method Mesh partitionning method
!!
!! \return PDM_part_mesh_t or PDM_part_mesh_nodal_t
!!
!!

  function PDM_generate_mesh_rectangle_ (comm,        &
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

    result(mesh_nodal)

    use iso_c_binding
    implicit none

    integer, intent(in)                   :: comm
    integer, intent(in)                   :: elt_type
    integer, intent(in)                   :: order
    type(c_ptr), intent(in)               :: ho_ordering
    double precision, intent(in)          :: xmin
    double precision, intent(in)          :: ymin
    double precision, intent(in)          :: zmin
    double precision, intent(in)          :: lengthx
    double precision, intent(in)          :: lengthy
    integer(kind=pdm_g_num_s), intent(in) :: n_x
    integer(kind=pdm_g_num_s), intent(in) :: n_y
    integer, intent(in)                   :: n_part
    integer, intent(in)                   :: part_method

    type(c_ptr)            :: mesh_nodal

    integer(c_int)        :: c_comm
    integer(c_int)        :: c_elt_type
    integer(c_int)        :: c_order
    type(c_ptr)           :: c_ho_ordering
    real(c_double)        :: c_xmin
    real(c_double)        :: c_ymin
    real(c_double)        :: c_zmin
    real(c_double)        :: c_lengthx
    real(c_double)        :: c_lengthy
#ifdef PDM_LONG_G_NUM
    integer(c_long)       :: c_n_x
    integer(c_long)       :: c_n_y
#else
    integer(c_int)        :: c_n_x
    integer(c_int)        :: c_n_y
#endif
    integer(c_int)        :: c_n_part
    integer(c_int)        :: c_part_method

    c_comm = PDM_MPI_Comm_f2c(comm)

    c_elt_type = elt_type
    c_order = order
    c_ho_ordering = ho_ordering
    c_xmin = xmin
    c_ymin = ymin
    c_zmin = zmin
    c_lengthx = lengthx
    c_lengthy = lengthy
    c_n_x = n_x
    c_n_y = n_y
    c_n_part = n_part
    c_part_method = part_method

    mesh_nodal = PDM_generate_mesh_rectangle_cf (c_comm,        &
                                                 c_elt_type,    &
                                                 c_order,       &
                                                 c_ho_ordering, &
                                                 c_xmin,        &
                                                 c_ymin,        &
                                                 c_zmin,        &
                                                 c_lengthx,     &
                                                 c_lengthy,     &
                                                 c_n_x,         &
                                                 c_n_y,         &
                                                 c_n_part,      &
                                                 c_part_method) 

  end function PDM_generate_mesh_rectangle_


!>
!!
!! \brief Create a partitionned ball mesh (3D).
!!
!! \param [in]  comm            MPI communicator
!! \param [in]  elt_type        Mesh element type
!! \param [in]  order           Mesh element order
!! \param [in]  ho_ordering     High order nodes ordering type
!! \param [in]  radius          Radius of the ball
!! \param [in]  hole_radius     Radius of the hole of the ball
!! \param [in]  center_x        x-coordinate of the ball center
!! \param [in]  center_y        y-coordinate of the ball center
!! \param [in]  center_z        z-coordinate of the ball center
!! \param [in]  n_x             Number of vertices on segments in x-direction
!! \param [in]  n_y             Number of vertices on segments in y-direction
!! \param [in]  n_z             Number of vertices on segments in z-direction
!! \param [in]  n_layer         Number of extrusion layers
!! \param [in]  geometric_ratio Geometric ratio for layer thickness
!! \param [in]  n_part          Number of mesh partitions
!! \param [in]  part_method     Mesh partitionning method
!!
!! \return PDM_part_mesh_t or PDM_part_mesh_nodal_t
!!
!!

  function PDM_generate_mesh_ball_ (comm,            &
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

    result(mesh_nodal)

    use iso_c_binding
    implicit none

    integer, intent(in)                   :: comm
    integer, intent(in)                   :: elt_type
    integer, intent(in)                   :: order
    type(c_ptr), intent(in)               :: ho_ordering
    double precision, intent(in)          :: radius
    double precision, intent(in)          :: hole_radius
    double precision, intent(in)          :: center_x
    double precision, intent(in)          :: center_y
    double precision, intent(in)          :: center_z
    integer(kind=pdm_g_num_s), intent(in) :: n_x
    integer(kind=pdm_g_num_s), intent(in) :: n_y
    integer(kind=pdm_g_num_s), intent(in) :: n_z
    integer(kind=pdm_g_num_s), intent(in) :: n_layer
    double precision, intent(in)          :: geometric_ratio
    integer, intent(in)                   :: n_part
    integer, intent(in)                   :: part_method

    type(c_ptr)            :: mesh_nodal

    integer(c_int) :: c_comm
    integer(c_int) :: c_elt_type
    integer(c_int) :: c_order
    type(c_ptr)    :: c_ho_ordering
    real(c_double) :: c_radius
    real(c_double) :: c_hole_radius
    real(c_double) :: c_center_x
    real(c_double) :: c_center_y
    real(c_double) :: c_center_z
#ifdef PDM_LONG_G_NUM
    integer(c_long) :: c_n_x
    integer(c_long) :: c_n_y
    integer(c_long) :: c_n_z
    integer(c_long) :: c_n_layer
#else
    integer(c_int) :: c_n_x
    integer(c_int) :: c_n_y
    integer(c_int) :: c_n_z
    integer(c_int) :: c_n_layer
#endif
    real(c_double) :: c_geometric_ratio
    integer(c_int) :: c_n_part
    integer(c_int) :: c_part_method

    c_comm = PDM_MPI_Comm_f2c(comm)

    c_elt_type = elt_type
    c_order = order
    c_ho_ordering = ho_ordering
    c_radius = radius
    c_hole_radius = hole_radius
    c_center_x = center_x
    c_center_y = center_y
    c_center_z = center_z
    c_n_x = n_x
    c_n_y = n_y
    c_n_z = n_z
    c_n_layer = n_layer
    c_geometric_ratio = geometric_ratio
    c_n_part = n_part
    c_part_method = part_method

    mesh_nodal = PDM_generate_mesh_ball_cf(c_comm,            &
                                            c_elt_type,        &
                                            c_order,           &
                                            c_ho_ordering,     &
                                            c_radius,          &
                                            c_hole_radius,     &
                                            c_center_x,        &
                                            c_center_y,        &
                                            c_center_z,        &
                                            c_n_x,             &
                                            c_n_y,             &
                                            c_n_z,             &
                                            c_n_layer,         &
                                            c_geometric_ratio, &
                                            c_n_part,          &
                                            c_part_method)     

  end function PDM_generate_mesh_ball_

!>
!!
!! \brief Create a partitionned parallelepiped mesh (3D).
!!
!! \param [in]  comm        MPI communicator
!! \param [in]  elt_type    Mesh element type
!! \param [in]  order       Mesh element order
!! \param [in]  ho_ordering High order nodes ordering type
!! \param [in]  xmin        x-coordinate of the rctangle minimum corner
!! \param [in]  ymin        y-coordinate of the rctangle minimum corner
!! \param [in]  zmin        z-coordinate of the rctangle minimum corner
!! \param [in]  lengthx     Length of the rectangle in the x-direction
!! \param [in]  lengthy     Length of the rectangle in the y-direction
!! \param [in]  lengthz     Length of the rectangle in the z-direction
!! \param [in]  n_x         Number of points in the x-direction
!! \param [in]  n_y         Number of points in the y-direction
!! \param [in]  n_z         Number of points in the z-direction
!! \param [in]  n_part      Number of mesh partitions
!! \param [in]  part_method Mesh partitionning method
!!
!! \return PDM_part_mesh_t or PDM_part_mesh_nodal_t
!!
!!

  function PDM_generate_mesh_parallelepiped_ (comm,         &
                                              elt_type,     & 
                                              order,        &
                                              ho_ordering,  &
                                              xmin,         &
                                              ymin,         &
                                              zmin,         &
                                              lengthx,      &
                                              lengthy,      &
                                              lengthz,      &
                                              n_x,          &
                                              n_y,          &
                                              n_z,          &
                                              n_part,       &
                                              part_method)  &

    result(mesh_nodal) 

    use iso_c_binding
    implicit none
! 
    integer, intent(in)                   :: comm
    integer, intent(in)                   :: elt_type
    integer, intent(in)                   :: order
    type(c_ptr), intent(in)               :: ho_ordering
    double precision, intent(in)          :: xmin
    double precision, intent(in)          :: ymin
    double precision, intent(in)          :: zmin
    double precision, intent(in)          :: lengthx
    double precision, intent(in)          :: lengthy
    double precision, intent(in)          :: lengthz
    integer(kind=pdm_g_num_s), intent(in) :: n_x
    integer(kind=pdm_g_num_s), intent(in) :: n_y
    integer(kind=pdm_g_num_s), intent(in) :: n_z
    integer, intent(in)                   :: n_part
    integer, intent(in)                   :: part_method

    type(c_ptr)           :: mesh_nodal

    integer(c_int)   :: c_comm
    integer(c_int)   :: c_elt_type
    integer(c_int)   :: c_order
    type(c_ptr)      :: c_ho_ordering
    real(c_double)   :: c_xmin
    real(c_double)   :: c_ymin
    real(c_double)   :: c_zmin
    real(c_double)   :: c_lengthx
    real(c_double)   :: c_lengthy
    real(c_double)   :: c_lengthz
#ifdef PDM_LONG_G_NUM
    integer(c_long)  :: c_n_x
    integer(c_long)  :: c_n_y
    integer(c_long)  :: c_n_z
#else
    integer(c_int)   :: c_n_x
    integer(c_int)   :: c_n_y
    integer(c_int)   :: c_n_z
#endif
    integer(c_int)   :: c_n_part
    integer(c_int)   :: c_part_method

    c_comm = PDM_MPI_Comm_f2c(comm)

    c_elt_type = elt_type
    c_order = order
    c_ho_ordering = ho_ordering
    c_xmin = xmin
    c_ymin = ymin
    c_zmin = zmin
    c_lengthx = lengthx
    c_lengthy = lengthy
    c_lengthz = lengthz
    c_n_x = n_x
    c_n_y = n_y
    c_n_z = n_z
    c_n_part = n_part
    c_part_method = part_method

    mesh_nodal =  PDM_generate_mesh_parallelepiped_cf (c_comm,         &
                                                       c_elt_type,     & 
                                                       c_order,        &
                                                       c_ho_ordering,  &
                                                       c_xmin,         &
                                                       c_ymin,         &
                                                       c_zmin,         &
                                                       c_lengthx,      &
                                                       c_lengthy,      &
                                                       c_lengthz,      &
                                                       c_n_x,          &
                                                       c_n_y,          &
                                                       c_n_z,          &
                                                       c_n_part,       &
                                                       c_part_method)  

  end function PDM_generate_mesh_parallelepiped_

!>
!!
!! \brief Create a partitionned rectangle mesh (2D) with descending connectivities.
!!
!! \param [in]   comm           MPI communicator
!! \param [in]   elt_type       Element type
!! \param [in]   xmin           Minimal x-coordinate
!! \param [in]   ymin           Minimal y-coordinate
!! \param [in]   zmin           Minimal z-coordinate
!! \param [in]   lengthx        Length of the rectangle in the x-direction
!! \param [in]   lengthy        Length of the rectangle in the y-direction
!! \param [in]   n_x            Number of points in the x-direction
!! \param [in]   n_y            Number of points in the y-direction
!! \param [in]   n_part         Number of partitions
!! \param [in]   part_method    Paritioning method
!! \param [in]   pn_vtx         Number of vertices
!! \param [in]   pn_edge        Number of edges
!! \param [in]   pn_face        Number of faces
!! \param [in]   pvtx_coord     Vertex coordinates
!! \param [in]   pedge_vtx      edge->vertex connectivity
!! \param [in]   pface_edge_idx Index of face->edge connectivity
!! \param [in]   pface_edge     face->edge connectivity
!! \param [in]   pvtx_ln_to_gn  Vertex global number
!! \param [in]   pedge_ln_to_gn Edge global number
!! \param [in]   pface_ln_to_gn Face global number
!!
!!

  subroutine PDM_generate_mesh_rectangle_ngon_ (comm, &
                                                elt_type, &
                                                xmin, &
                                                ymin, &
                                                zmin, &
                                                lengthx, &
                                                lengthy, &
                                                n_x, &
                                                n_y, &
                                                n_part, &
                                                part_method, &
                                                pn_vtx, &
                                                pn_edge, &
                                                pn_face, &
                                                pvtx_coord, &
                                                pedge_vtx, &
                                                pface_edge_idx, &
                                                pface_edge, &
                                                pface_vtx, &
                                                pvtx_ln_to_gn, &
                                                pedge_ln_to_gn, &
                                                pface_ln_to_gn)

    use iso_c_binding
    implicit none
 
    integer, intent(in)                   :: comm
    integer, intent(in)                   :: elt_type
    double precision, intent(in)          :: xmin
    double precision, intent(in)          :: ymin
    double precision, intent(in)          :: zmin
    double precision, intent(in)          :: lengthx
    double precision, intent(in)          :: lengthy
    integer(kind=pdm_g_num_s), intent(in) :: n_x
    integer(kind=pdm_g_num_s), intent(in) :: n_y
    integer, intent(in)                   :: n_part
    integer, intent(in)                   :: part_method

    integer(kind=pdm_l_num_s), pointer    :: pn_vtx(:)
    integer(kind=pdm_l_num_s), pointer    :: pn_edge(:)
    integer(kind=pdm_l_num_s), pointer    :: pn_face(:)
    type(PDM_pointer_array_t), pointer    :: pvtx_coord
    type(PDM_pointer_array_t), pointer    :: pedge_vtx
    type(PDM_pointer_array_t), pointer    :: pface_edge_idx
    type(PDM_pointer_array_t), pointer    :: pface_edge
    type(PDM_pointer_array_t), pointer    :: pface_vtx
    type(PDM_pointer_array_t), pointer    :: pvtx_ln_to_gn
    type(PDM_pointer_array_t), pointer    :: pedge_ln_to_gn
    type(PDM_pointer_array_t), pointer    :: pface_ln_to_gn

    integer(c_int)        :: c_comm
    integer(c_int)        :: c_elt_type
    real(c_double)        :: c_xmin
    real(c_double)        :: c_ymin
    real(c_double)        :: c_zmin
    real(c_double)        :: c_lengthx
    real(c_double)        :: c_lengthy
#ifdef PDM_LONG_G_NUM
    integer(c_long)       :: c_n_x
    integer(c_long)       :: c_n_y
#else
    integer(c_int)        :: c_n_x
    integer(c_int)        :: c_n_y
#endif
    integer(c_int)        :: c_n_part
    integer(c_int)        :: c_part_method

    type(c_ptr)           :: c_pn_vtx
    type(c_ptr)           :: c_pn_edge
    type(c_ptr)           :: c_pn_face
    type(c_ptr)           :: c_pvtx_coord
    type(c_ptr)           :: c_pedge_vtx
    type(c_ptr)           :: c_pface_edge_idx
    type(c_ptr)           :: c_pface_edge
    type(c_ptr)           :: c_pface_vtx
    type(c_ptr)           :: c_pvtx_ln_to_gn
    type(c_ptr)           :: c_pedge_ln_to_gn
    type(c_ptr)           :: c_pface_ln_to_gn

    integer, allocatable  :: s_array(:)

    integer :: i  
    integer, pointer      :: ipart_pface_edge_idx(:)

    c_comm = PDM_MPI_Comm_f2c(comm)

    c_elt_type = elt_type
    c_xmin = xmin
    c_ymin = ymin
    c_zmin = zmin
    c_lengthx = lengthx
    c_lengthy = lengthy
    c_n_x = n_x
    c_n_y = n_y
    c_n_part = n_part
    c_part_method = part_method
   
    call PDM_generate_mesh_rectangle_ngon_cf (c_comm, &
                                              c_elt_type, &
                                              c_xmin, &
                                              c_ymin, &
                                              c_zmin, &
                                              c_lengthx, &
                                              c_lengthy, &
                                              c_n_x, &
                                              c_n_y, &
                                              c_n_part, &
                                              c_part_method, &
                                              c_pn_vtx, &
                                              c_pn_edge, &
                                              c_pn_face, &
                                              c_pvtx_coord, &
                                              c_pedge_vtx, &
                                              c_pface_edge_idx, &
                                              c_pface_edge, &
                                              c_pface_vtx, &
                                              c_pvtx_ln_to_gn, &
                                              c_pedge_ln_to_gn, &
                                              c_pface_ln_to_gn)
 
    call c_f_pointer(c_pn_vtx, &
                     pn_vtx,   &
                     [n_part])

    call c_f_pointer(c_pn_edge, &
                     pn_edge,   &
                     [n_part])

    call c_f_pointer(c_pn_face, &
                     pn_face,   &
                     [n_part])

    allocate(s_array(n_part)) 

    do i = 1, n_part
      s_array(i) = 3 * pn_vtx(i)
    enddo

    call PDM_pointer_array_create (pvtx_coord,         &
                                   n_part,             &
                                   PDM_TYPE_DOUBLE,    &
                                   c_pvtx_coord,       &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = 2 * pn_edge(i)
    enddo

    call PDM_pointer_array_create (pedge_vtx,          &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_pedge_vtx,        &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = pn_face(i) + 1
    enddo

    call PDM_pointer_array_create (pface_edge_idx,     &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_pface_edge_idx,   &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      call PDM_pointer_array_part_get (pface_edge_idx, &
                                       i-1,            &
                                       ipart_pface_edge_idx)      
      s_array(i) = ipart_pface_edge_idx(pn_face(i) + 1)
    enddo

    call PDM_pointer_array_create (pface_edge,         &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_pface_edge,       &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    call PDM_pointer_array_create (pface_vtx,          &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_pface_vtx,        &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = pn_vtx(i)
    enddo

    call PDM_pointer_array_create (pvtx_ln_to_gn,      &
                                   n_part,             &
                                   PDM_TYPE_G_NUM,     &
                                   c_pvtx_ln_to_gn,    &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = pn_edge(i)
    enddo

    call PDM_pointer_array_create (pedge_ln_to_gn,     &
                                   n_part,             &
                                   PDM_TYPE_G_NUM,     &
                                   c_pedge_ln_to_gn,   &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = pn_face(i)
    enddo

    call PDM_pointer_array_create (pface_ln_to_gn,     &
                                   n_part,             &
                                   PDM_TYPE_G_NUM,     &
                                   c_pface_ln_to_gn,   &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    deallocate(s_array)

  end subroutine PDM_generate_mesh_rectangle_ngon_


!>
!!
!! \brief Create a partitionned sphere mesh (2D) with descending connectivities.
!!
!! \param [in]   comm           MPI communicator
!! \param [in]   elt_type       Element type
!! \param [in]   order          Element order
!! \param [in]   ho_ordering    Ordering of nodes of the HO element
!! \param [in]   radius         Radius of the sphere
!! \param [in]   center_x       x-coordinate of the sphere center
!! \param [in]   center_y       y-coordinate of the sphere center
!! \param [in]   center_z       z-coordinate of the sphere center
!! \param [in]   n_u            Number of vertices in the u-direction
!! \param [in]   n_v            Number of vertices in the v-direction
!! \param [in]   n_part         Number of partitions
!! \param [in]   part_method    Paritioning method
!! \param [in]   pn_vtx         Number of vertices
!! \param [in]   pn_edge        Number of edges
!! \param [in]   pn_face        Number of faces
!! \param [in]   pvtx_coord     Vertex coordinates
!! \param [in]   pedge_vtx      edge->vertex connectivity
!! \param [in]   pface_edge_idx Index of face->edge connectivity
!! \param [in]   pface_edge     face->edge connectivity
!! \param [in]   pface_vtx      face->vtx connectivity
!! \param [in]   pvtx_ln_to_gn  Vertex global number
!! \param [in]   pedge_ln_to_gn Edge global number
!! \param [in]   pface_ln_to_gn Face global number
!!
!!

  subroutine PDM_generate_mesh_sphere_ngon_ (comm,           &
                                             elt_type,       &
                                             order,          &
                                             ho_ordering,    &
                                             radius,         &
                                             center_x,       &
                                             center_y,       &
                                             center_z,       &
                                             n_u,            &
                                             n_v,            &
                                             n_part,         &
                                             part_method,    &
                                             pn_vtx,         &
                                             pn_edge,        &
                                             pn_face,        &
                                             pvtx_coord,     & 
                                             pedge_vtx,      &
                                             pface_edge_idx, &
                                             pface_edge,     &
                                             pface_vtx,      &
                                             pvtx_ln_to_gn,  &
                                             pedge_ln_to_gn, &
                                             pface_ln_to_gn)


    use iso_c_binding
    implicit none

    integer, intent(in)                   :: comm 
    integer, intent(in)                   :: elt_type
    integer, intent(in)                   :: order
    type(c_ptr), intent(in)               :: ho_ordering
    double precision, intent(in)          :: radius
    double precision, intent(in)          :: center_x
    double precision, intent(in)          :: center_y
    double precision, intent(in)          :: center_z
    integer(kind=pdm_g_num_s), intent(in) :: n_u
    integer(kind=pdm_g_num_s), intent(in) :: n_v
    integer, intent(in)                   :: n_part
    integer, intent(in)                   :: part_method

    integer(c_int)        :: c_comm 
    integer(c_int)        :: c_elt_type
    integer(c_int)        :: c_order
    type(c_ptr)           :: c_ho_ordering
    real(c_double)        :: c_radius
    real(c_double)        :: c_center_x
    real(c_double)        :: c_center_y
    real(c_double)        :: c_center_z
#ifdef PDM_LONG_G_NUM
    integer(c_long)       :: c_n_u
    integer(c_long)       :: c_n_v
#else
    integer(c_int)        :: c_n_u
    integer(c_int)        :: c_n_v
#endif
    integer(c_int)        :: c_n_part
    integer(c_int)        :: c_part_method

    integer(kind=pdm_l_num_s), pointer    :: pn_vtx(:)
    integer(kind=pdm_l_num_s), pointer    :: pn_edge(:)
    integer(kind=pdm_l_num_s), pointer    :: pn_face(:)

    type(PDM_pointer_array_t), pointer    :: pvtx_coord
    type(PDM_pointer_array_t), pointer    :: pedge_vtx
    type(PDM_pointer_array_t), pointer    :: pface_edge_idx
    type(PDM_pointer_array_t), pointer    :: pface_edge
    type(PDM_pointer_array_t), pointer    :: pface_vtx
    type(PDM_pointer_array_t), pointer    :: pvtx_ln_to_gn
    type(PDM_pointer_array_t), pointer    :: pedge_ln_to_gn
    type(PDM_pointer_array_t), pointer    :: pface_ln_to_gn

    type(c_ptr)           :: c_pn_vtx
    type(c_ptr)           :: c_pn_edge
    type(c_ptr)           :: c_pn_face
    type(c_ptr)           :: c_pvtx_coord
    type(c_ptr)           :: c_pedge_vtx
    type(c_ptr)           :: c_pface_edge_idx
    type(c_ptr)           :: c_pface_edge
    type(c_ptr)           :: c_pface_vtx
    type(c_ptr)           :: c_pvtx_ln_to_gn
    type(c_ptr)           :: c_pedge_ln_to_gn
    type(c_ptr)           :: c_pface_ln_to_gn

    integer, allocatable  :: s_array(:)

    integer :: i  
    integer, pointer      :: ipart_pface_edge_idx(:)

    c_comm = PDM_MPI_Comm_f2c(comm)

    c_elt_type    = elt_type
    c_order       = order
    c_ho_ordering = ho_ordering
    c_radius      = radius
    c_center_x    = center_x
    c_center_y    = center_y
    c_center_z    = center_z
    c_n_u         = n_u
    c_n_v         = n_v
    c_n_part      = n_part
    c_part_method = part_method

    call PDM_generate_mesh_sphere_ngon_cf (c_comm, &
                                           c_elt_type, &
                                           c_order, &
                                           c_ho_ordering, &
                                           c_radius, &
                                           c_center_x, &
                                           c_center_y, &
                                           c_center_z, &
                                           c_n_u, &
                                           c_n_v, &
                                           c_n_part, &
                                           c_part_method, &
                                           c_pn_vtx, &
                                           c_pn_edge, &
                                           c_pn_face, &
                                           c_pvtx_coord, &
                                           c_pedge_vtx, &
                                           c_pface_edge_idx, &
                                           c_pface_edge, &
                                           c_pface_vtx, &
                                           c_pvtx_ln_to_gn, &
                                           c_pedge_ln_to_gn, &
                                           c_pface_ln_to_gn)
 
    call c_f_pointer(c_pn_vtx, &
                     pn_vtx,   &
                     [n_part])

    call c_f_pointer(c_pn_edge, &
                     pn_edge,   &
                     [n_part])

    call c_f_pointer(c_pn_face, &
                     pn_face,   &
                     [n_part])

    allocate(s_array(n_part)) 

    do i = 1, n_part
      s_array(i) = 3 * pn_vtx(i)
    enddo

    call PDM_pointer_array_create (pvtx_coord,         &
                                   n_part,             &
                                   PDM_TYPE_DOUBLE,    &
                                   c_pvtx_coord,       &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = 2 * pn_edge(i)
    enddo

    call PDM_pointer_array_create (pedge_vtx,          &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_pedge_vtx,        &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = pn_face(i) + 1
    enddo

    call PDM_pointer_array_create (pface_edge_idx,     &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_pface_edge_idx,   &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      call PDM_pointer_array_part_get (pface_edge_idx, &
                                       i-1,            &
                                       ipart_pface_edge_idx)      
      s_array(i) = ipart_pface_edge_idx(pn_face(i) + 1)
    enddo

    call PDM_pointer_array_create (pface_edge,         &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_pface_edge,       &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    call PDM_pointer_array_create (pface_vtx,          &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_pface_vtx,        &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = pn_vtx(i)
    enddo

    call PDM_pointer_array_create (pvtx_ln_to_gn,      &
                                   n_part,             &
                                   PDM_TYPE_G_NUM,     &
                                   c_pvtx_ln_to_gn,    &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = pn_edge(i)
    enddo

    call PDM_pointer_array_create (pedge_ln_to_gn,     &
                                   n_part,             &
                                   PDM_TYPE_G_NUM,     &
                                   c_pedge_ln_to_gn,   &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = pn_face(i)
    enddo

    call PDM_pointer_array_create (pface_ln_to_gn,     &
                                   n_part,             &
                                   PDM_TYPE_G_NUM,     &
                                   c_pface_ln_to_gn,   &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    deallocate(s_array)

  end subroutine PDM_generate_mesh_sphere_ngon_ 


!>
!!
!! \brief Create a partitionned ball mesh (3D) with descending connectivities.
!!
!! \param [in]  comm                      MPI communicator
!! \param [in]  elt_type                  Mesh element type
!! \param [in]  order                     Mesh element order
!! \param [in]  ho_ordering               High order nodes ordering type
!! \param [in]  radius                    Radius of the ball
!! \param [in]  hole_radius               Radius of the hole of the ball
!! \param [in]  center_x                  x-coordinate of the ball center
!! \param [in]  center_y                  y-coordinate of the ball center
!! \param [in]  center_z                  z-coordinate of the ball center
!! \param [in]  n_x                       Number of vertices on segments in x-direction
!! \param [in]  n_y                       Number of vertices on segments in y-direction
!! \param [in]  n_z                       Number of vertices on segments in z-direction
!! \param [in]  n_layer                   Number of extrusion layers
!! \param [in]  geometric_ratio           Geometric ratio for layer thickness
!! \param [in]  n_part                    Number of mesh partitions
!! \param [in]  part_method               Mesh partitionning method
!! \param [out] pn_vtx                    Number of vertices
!! \param [out] pn_edge                   Number of edges
!! \param [out] pn_face                   Number of faces
!! \param [out] pvtx_coord                Vertex coordinates
!! \param [out] pedge_vtx                 edge->vertex connectivity
!! \param [out] pface_edge_idx            Index of face->edge connectivity
!! \param [out] pface_edge                face->edge connectivity
!! \param [out] pface_vtx                face->vtx connectivity
!! \param [out] pvtx_ln_to_gn             Vertex global number
!! \param [out] pedge_ln_to_gn            Edge global number
!! \param [out] pface_ln_to_gn            Face global number
!! \param [out] pn_surface                Number of surfaces
!! \param [out] psurface_face_idx         surface->face connectivity index
!! \param [out] psurface_face             surface->face connectivity
!! \param [out] psurface_face_ln_to_gn    surface->face connectivity with global numbers
!!
!!

  subroutine PDM_generate_mesh_ball_ngon_ (comm,            &
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
                                           part_method,     &
                                           pn_vtx,          &
                                           pn_edge,         &
                                           pn_face,         &
                                           pn_cell,         &
                                           pvtx_coord,      &
                                           pedge_vtx,       & 
                                           pface_edge_idx,  &
                                           pface_edge,      &
                                           pface_vtx,       &
                                           pcell_face_idx,  & 
                                           pcell_face,      & 
                                           pvtx_ln_to_gn,   &
                                           pedge_ln_to_gn,  &
                                           pface_ln_to_gn,  &
                                           pcell_ln_to_gn,  &
                                           pn_surface,      & 
                                           psurface_face_idx,&
                                           psurface_face,    &
                                           psurface_face_ln_to_gn)


    use iso_c_binding
    implicit none

    integer, intent(in)                   :: comm
    integer, intent(in)                   :: elt_type
    integer, intent(in)                   :: order
    type(c_ptr), intent(in)               :: ho_ordering
    double precision, intent(in)          :: radius
    double precision, intent(in)          :: hole_radius
    double precision, intent(in)          :: center_x
    double precision, intent(in)          :: center_y
    double precision, intent(in)          :: center_z
    integer(kind=pdm_g_num_s), intent(in) :: n_x
    integer(kind=pdm_g_num_s), intent(in) :: n_y
    integer(kind=pdm_g_num_s), intent(in) :: n_z
    integer(kind=pdm_g_num_s), intent(in) :: n_layer
    double precision, intent(in)          :: geometric_ratio
    integer, intent(in)                   :: n_part
    integer, intent(in)                   :: part_method

    integer(kind=pdm_l_num_s), pointer    :: pn_vtx(:)
    integer(kind=pdm_l_num_s), pointer    :: pn_edge(:)
    integer(kind=pdm_l_num_s), pointer    :: pn_face(:)
    integer(kind=pdm_l_num_s), pointer    :: pn_cell(:)
    integer(kind=pdm_l_num_s), pointer    :: pn_surface(:)
    type(PDM_pointer_array_t), pointer    :: pvtx_coord
    type(PDM_pointer_array_t), pointer    :: pedge_vtx
    type(PDM_pointer_array_t), pointer    :: pface_edge_idx
    type(PDM_pointer_array_t), pointer    :: pface_edge
    type(PDM_pointer_array_t), pointer    :: pface_vtx
    type(PDM_pointer_array_t), pointer    :: pcell_face_idx
    type(PDM_pointer_array_t), pointer    :: pcell_face
    type(PDM_pointer_array_t), pointer    :: pvtx_ln_to_gn
    type(PDM_pointer_array_t), pointer    :: pedge_ln_to_gn
    type(PDM_pointer_array_t), pointer    :: pface_ln_to_gn
    type(PDM_pointer_array_t), pointer    :: pcell_ln_to_gn
    type(PDM_pointer_array_t), pointer    :: psurface_face_idx
    type(PDM_pointer_array_t), pointer    :: psurface_face
    type(PDM_pointer_array_t), pointer    :: psurface_face_ln_to_gn

    type(c_ptr) :: c_pn_vtx
    type(c_ptr) :: c_pn_edge
    type(c_ptr) :: c_pn_face
    type(c_ptr) :: c_pn_cell
    type(c_ptr) :: c_pn_surface
    type(c_ptr) :: c_pvtx_coord
    type(c_ptr) :: c_pedge_vtx
    type(c_ptr) :: c_pface_edge_idx
    type(c_ptr) :: c_pface_edge
    type(c_ptr) :: c_pface_vtx
    type(c_ptr) :: c_pcell_face_idx
    type(c_ptr) :: c_pcell_face
    type(c_ptr) :: c_pvtx_ln_to_gn
    type(c_ptr) :: c_pedge_ln_to_gn
    type(c_ptr) :: c_pface_ln_to_gn
    type(c_ptr) :: c_pcell_ln_to_gn
    type(c_ptr) :: c_psurface_face_idx
    type(c_ptr) :: c_psurface_face
    type(c_ptr) :: c_psurface_face_ln_to_gn

    integer(c_int) :: c_comm
    integer(c_int) :: c_elt_type
    integer(c_int) :: c_order
    type(c_ptr)    :: c_ho_ordering
    real(c_double) :: c_radius
    real(c_double) :: c_hole_radius
    real(c_double) :: c_center_x
    real(c_double) :: c_center_y
    real(c_double) :: c_center_z
#ifdef PDM_LONG_G_NUM
    integer(c_long) :: c_n_x
    integer(c_long) :: c_n_y
    integer(c_long) :: c_n_z
    integer(c_long) :: c_n_layer
#else
    integer(c_int) :: c_n_x
    integer(c_int) :: c_n_y
    integer(c_int) :: c_n_z
    integer(c_int) :: c_n_layer
#endif
    real(c_double) :: c_geometric_ratio
    integer(c_int) :: c_n_part
    integer(c_int) :: c_part_method

    integer, allocatable  :: s_array(:)

    integer          :: i  
    integer, pointer :: ipart_pface_edge_idx(:)
    integer, pointer :: ipart_pcell_face_idx(:)
    integer, pointer :: ipsurface_face_idx(:)

    c_comm = PDM_MPI_Comm_f2c(comm)

    c_elt_type        = elt_type
    c_order           = order
    c_ho_ordering     = ho_ordering
    c_radius          = radius
    c_hole_radius     = hole_radius
    c_center_x        = center_x
    c_center_y        = center_y
    c_center_z        = center_z
    c_n_layer         = n_layer
    c_geometric_ratio = geometric_ratio
    c_n_x             = n_x
    c_n_y             = n_y
    c_n_z             = n_z
    c_n_part          = n_part
    c_part_method     = part_method

    call PDM_generate_mesh_ball_ngon_cf (c_comm, &
                                         c_elt_type, &
                                         c_order, &
                                         c_ho_ordering, &
                                         c_radius, &
                                         c_hole_radius, &
                                         c_center_x, &
                                         c_center_y, &
                                         c_center_z, &
                                         c_n_x, &
                                         c_n_y, &
                                         c_n_z, &
                                         c_n_layer, &
                                         c_geometric_ratio, &
                                         c_n_part, &
                                         c_part_method, &
                                         c_pn_vtx, &
                                         c_pn_edge, &
                                         c_pn_face, &
                                         c_pn_cell, &
                                         c_pvtx_coord, &
                                         c_pedge_vtx, &
                                         c_pface_edge_idx, &
                                         c_pface_edge, &
                                         c_pface_vtx, &
                                         c_pcell_face_idx, &
                                         c_pcell_face, &
                                         c_pvtx_ln_to_gn, &
                                         c_pedge_ln_to_gn, &
                                         c_pface_ln_to_gn, &
                                         c_pcell_ln_to_gn, &
                                         c_pn_surface, &
                                         c_psurface_face_idx, &
                                         c_psurface_face, &
                                         c_psurface_face_ln_to_gn)

    call c_f_pointer(c_pn_vtx, &
                     pn_vtx,   &
                     [n_part])

    call c_f_pointer(c_pn_edge, &
                     pn_edge,   &
                     [n_part])

    call c_f_pointer(c_pn_face, &
                     pn_face,   &
                     [n_part])

    call c_f_pointer(c_pn_cell, &
                     pn_cell,   &
                     [n_part])

    call c_f_pointer(c_pn_surface, &
                     pn_surface,   &
                     [n_part])

    allocate(s_array(n_part)) 

    do i = 1, n_part
      s_array(i) = 3 * pn_vtx(i)
    enddo

    call PDM_pointer_array_create (pvtx_coord,         &
                                   n_part,             &
                                   PDM_TYPE_DOUBLE,    &
                                   c_pvtx_coord,       &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = 2 * pn_edge(i)
    enddo

    call PDM_pointer_array_create (pedge_vtx,          &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_pedge_vtx,        &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = pn_face(i) + 1
    enddo

    call PDM_pointer_array_create (pface_edge_idx,     &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_pface_edge_idx,   &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      call PDM_pointer_array_part_get (pface_edge_idx, &
                                       i-1,            &
                                       ipart_pface_edge_idx)      
      s_array(i) = ipart_pface_edge_idx(pn_face(i) + 1)
    enddo

    call PDM_pointer_array_create (pface_edge,         &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_pface_edge,       &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    call PDM_pointer_array_create (pface_vtx,          &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_pface_vtx,        &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = pn_cell(i) + 1
    enddo

    call PDM_pointer_array_create (pcell_face_idx,     &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_pcell_face_idx,   &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      call PDM_pointer_array_part_get (pcell_face_idx, &
                                       i-1,            &
                                       ipart_pcell_face_idx)      
      s_array(i) = ipart_pcell_face_idx(pn_cell(i) + 1)
    enddo

    call PDM_pointer_array_create (pcell_face,         &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_pcell_face,       &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)
    do i = 1, n_part
      s_array(i) = pn_vtx(i)
    enddo

    call PDM_pointer_array_create (pvtx_ln_to_gn,      &
                                   n_part,             &
                                   PDM_TYPE_G_NUM,     &
                                   c_pvtx_ln_to_gn,    &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = pn_edge(i)
    enddo

    call PDM_pointer_array_create (pedge_ln_to_gn,     &
                                   n_part,             &
                                   PDM_TYPE_G_NUM,     &
                                   c_pedge_ln_to_gn,   &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = pn_face(i)
    enddo

    call PDM_pointer_array_create (pface_ln_to_gn,     &
                                   n_part,             &
                                   PDM_TYPE_G_NUM,     &
                                   c_pface_ln_to_gn,   &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = pn_cell(i)
    enddo

    call PDM_pointer_array_create (pcell_ln_to_gn,     &
                                   n_part,             &
                                   PDM_TYPE_G_NUM,     &
                                   c_pcell_ln_to_gn,   &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = pn_surface(i) + 1
    enddo

    call PDM_pointer_array_create (psurface_face_idx,  &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_psurface_face_idx,&
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      call PDM_pointer_array_part_get (psurface_face_idx, &
                                       i-1,               &
                                       ipsurface_face_idx)      
      s_array(i) = ipsurface_face_idx(pn_surface(i) + 1)
    enddo

    call PDM_pointer_array_create (psurface_face,      &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_psurface_face,    &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    call PDM_pointer_array_create (psurface_face_ln_to_gn,   &
                                   n_part,                   &
                                   PDM_TYPE_G_NUM,           &
                                   c_psurface_face_ln_to_gn, &
                                   s_array,                  &
                                   PDM_OWNERSHIP_KEEP)

    deallocate(s_array)

  end subroutine PDM_generate_mesh_ball_ngon_

!>
!!
!! \brief Create a partitionned parallelepiped mesh (3D) with descending connectivities.
!!
!! \param [in]  comm                      MPI communicator
!! \param [in]  elt_type                  Mesh element type
!! \param [in]  order                     Mesh element order
!! \param [in]  ho_ordering               High order nodes ordering type
!! \param [in]  radius                    Radius of the ball
!! \param [in]  hole_radius               Radius of the hole of the ball
!! \param [in]  center_x                  x-coordinate of the ball center
!! \param [in]  center_y                  y-coordinate of the ball center
!! \param [in]  center_z                  z-coordinate of the ball center
!! \param [in]  n_x                       Number of vertices on segments in x-direction
!! \param [in]  n_y                       Number of vertices on segments in y-direction
!! \param [in]  n_z                       Number of vertices on segments in z-direction
!! \param [in]  n_layer                   Number of extrusion layers
!! \param [in]  geometric_ratio           Geometric ratio for layer thickness
!! \param [in]  n_part                    Number of mesh partitions
!! \param [in]  part_method               Mesh partitionning method
!! \param [out] pn_vtx                    Number of vertices
!! \param [out] pn_edge                   Number of edges
!! \param [out] pn_face                   Number of faces
!! \param [out] pvtx_coord                Vertex coordinates
!! \param [out] pedge_vtx                 edge->vertex connectivity
!! \param [out] pface_edge_idx            Index of face->edge connectivity
!! \param [out] pface_edge                face->edge connectivity
!! \param [out] pface_vtx                 face->vtx connectivity
!! \param [out] pvtx_ln_to_gn             Vertex global number
!! \param [out] pedge_ln_to_gn            Edge global number
!! \param [out] pface_ln_to_gn            Face global number
!! \param [out] pn_surface                Number of surfaces
!! \param [out] psurface_face_idx         surface->face connectivity index
!! \param [out] psurface_face             surface->face connectivity
!! \param [out] psurface_face_ln_to_gn    surface->face connectivity with global numbers
!! \param [out] pn_ridge                  Number of ridges
!! \param [out] pridge_edge_idx           ridge->edge connectivity index
!! \param [out] pridge_edge               ridge->edge connectivity
!! \param [out] pridge_edge_ln_to_gn      ridge->edge connectivity with global numbers
!!
!!

  subroutine PDM_generate_mesh_parallelepiped_ngon_ (comm,                   &
                                                     elt_type,               & 
                                                     order,                  &
                                                     ho_ordering,            &
                                                     xmin,                   &
                                                     ymin,                   &
                                                     zmin,                   &
                                                     lengthx,                &
                                                     lengthy,                &
                                                     lengthz,                &
                                                     n_x,                    &
                                                     n_y,                    &
                                                     n_z,                    &
                                                     n_part,                 &
                                                     part_method,            &
                                                     pn_vtx,                 &
                                                     pn_edge,                &
                                                     pn_face,                &
                                                     pn_cell,                &
                                                     pvtx_coord,             &
                                                     pedge_vtx,              &
                                                     pface_edge_idx,         &
                                                     pface_edge,             &
                                                     pface_vtx,              &
                                                     pcell_face_idx,         &
                                                     pcell_face,             &
                                                     pvtx_ln_to_gn,          &
                                                     pedge_ln_to_gn,         &
                                                     pface_ln_to_gn,         &
                                                     pcell_ln_to_gn,         &
                                                     pn_surface,             &
                                                     psurface_face_idx,      &
                                                     psurface_face,          &
                                                     psurface_face_ln_to_gn, &
                                                     pn_ridge,               &
                                                     pridge_edge_idx,        &
                                                     pridge_edge,            &
                                                     pridge_edge_ln_to_gn)
      
    use iso_c_binding
    implicit none
! 
    integer, intent(in)                   :: comm
    integer, intent(in)                   :: elt_type
    integer, intent(in)                   :: order
    type(c_ptr), intent(in)               :: ho_ordering
    double precision, intent(in)          :: xmin
    double precision, intent(in)          :: ymin
    double precision, intent(in)          :: zmin
    double precision, intent(in)          :: lengthx
    double precision, intent(in)          :: lengthy
    double precision, intent(in)          :: lengthz
    integer(kind=pdm_g_num_s), intent(in) :: n_x
    integer(kind=pdm_g_num_s), intent(in) :: n_y
    integer(kind=pdm_g_num_s), intent(in) :: n_z
    integer, intent(in)                   :: n_part
    integer, intent(in)                   :: part_method

    integer(kind=pdm_l_num_s), pointer    :: pn_vtx(:)
    integer(kind=pdm_l_num_s), pointer    :: pn_edge(:)
    integer(kind=pdm_l_num_s), pointer    :: pn_face(:)
    integer(kind=pdm_l_num_s), pointer    :: pn_cell(:)
    integer(kind=pdm_l_num_s), pointer    :: pn_surface(:)
    integer(kind=pdm_l_num_s), pointer    :: pn_ridge(:)
    type(PDM_pointer_array_t), pointer    :: pvtx_coord
    type(PDM_pointer_array_t), pointer    :: pedge_vtx
    type(PDM_pointer_array_t), pointer    :: pface_edge_idx
    type(PDM_pointer_array_t), pointer    :: pface_edge
    type(PDM_pointer_array_t), pointer    :: pface_vtx
    type(PDM_pointer_array_t), pointer    :: pcell_face_idx
    type(PDM_pointer_array_t), pointer    :: pcell_face
    type(PDM_pointer_array_t), pointer    :: pvtx_ln_to_gn
    type(PDM_pointer_array_t), pointer    :: pedge_ln_to_gn
    type(PDM_pointer_array_t), pointer    :: pface_ln_to_gn
    type(PDM_pointer_array_t), pointer    :: pcell_ln_to_gn
    type(PDM_pointer_array_t), pointer    :: psurface_face_idx
    type(PDM_pointer_array_t), pointer    :: psurface_face
    type(PDM_pointer_array_t), pointer    :: psurface_face_ln_to_gn
    type(PDM_pointer_array_t), pointer    :: pridge_edge_idx
    type(PDM_pointer_array_t), pointer    :: pridge_edge
    type(PDM_pointer_array_t), pointer    :: pridge_edge_ln_to_gn

    integer(c_int)   :: c_comm
    integer(c_int)   :: c_elt_type
    integer(c_int)   :: c_order
    type(c_ptr)      :: c_ho_ordering
    real(c_double)   :: c_xmin
    real(c_double)   :: c_ymin
    real(c_double)   :: c_zmin
    real(c_double)   :: c_lengthx
    real(c_double)   :: c_lengthy
    real(c_double)   :: c_lengthz
#ifdef PDM_LONG_G_NUM
    integer(c_long)  :: c_n_x
    integer(c_long)  :: c_n_y
    integer(c_long)  :: c_n_z
#else
    integer(c_int)   :: c_n_x
    integer(c_int)   :: c_n_y
    integer(c_int)   :: c_n_z
#endif
    integer(c_int)   :: c_n_part
    integer(c_int)   :: c_part_method

    type(c_ptr) :: c_pn_vtx
    type(c_ptr) :: c_pn_edge
    type(c_ptr) :: c_pn_face
    type(c_ptr) :: c_pn_cell
    type(c_ptr) :: c_pvtx_coord
    type(c_ptr) :: c_pedge_vtx
    type(c_ptr) :: c_pface_edge_idx
    type(c_ptr) :: c_pface_edge
    type(c_ptr) :: c_pface_vtx
    type(c_ptr) :: c_pcell_face_idx
    type(c_ptr) :: c_pcell_face
    type(c_ptr) :: c_pvtx_ln_to_gn
    type(c_ptr) :: c_pedge_ln_to_gn
    type(c_ptr) :: c_pface_ln_to_gn
    type(c_ptr) :: c_pcell_ln_to_gn
    type(c_ptr) :: c_pn_surface
    type(c_ptr) :: c_psurface_face_idx
    type(c_ptr) :: c_psurface_face
    type(c_ptr) :: c_psurface_face_ln_to_gn
    type(c_ptr) :: c_pn_ridge
    type(c_ptr) :: c_pridge_edge_idx
    type(c_ptr) :: c_pridge_edge
    type(c_ptr) :: c_pridge_edge_ln_to_gn

    integer, allocatable  :: s_array(:)

    integer          :: i  
    integer, pointer :: ipart_pface_edge_idx(:)
    integer, pointer :: ipart_pcell_face_idx(:)
    integer, pointer :: ipsurface_face_idx(:)
    integer, pointer :: ipridge_edge_idx(:)

    c_comm = PDM_MPI_Comm_f2c(comm)

    c_elt_type = elt_type
    c_order = order
    c_ho_ordering = ho_ordering
    c_xmin = xmin
    c_ymin = ymin
    c_zmin = zmin
    c_lengthx = lengthx
    c_lengthy = lengthy
    c_lengthz = lengthz
    c_n_x = n_x
    c_n_y = n_y
    c_n_z = n_z
    c_n_part = n_part
    c_part_method = part_method
 
    call PDM_generate_mesh_parallelepiped_ngon_cf (c_comm,                   &
                                                   c_elt_type,               & 
                                                   c_order,                  &
                                                   c_ho_ordering,            &
                                                   c_xmin,                   &
                                                   c_ymin,                   &
                                                   c_zmin,                   &
                                                   c_lengthx,                &
                                                   c_lengthy,                &
                                                   c_lengthz,                &
                                                   c_n_x,                    &
                                                   c_n_y,                    &
                                                   c_n_z,                    &
                                                   c_n_part,                 &
                                                   c_part_method,            &
                                                   c_pn_vtx,                 &
                                                   c_pn_edge,                &
                                                   c_pn_face,                &
                                                   c_pn_cell,                &
                                                   c_pvtx_coord,             &
                                                   c_pedge_vtx,              &
                                                   c_pface_edge_idx,         &
                                                   c_pface_edge,             &
                                                   c_pface_vtx,              &
                                                   c_pcell_face_idx,         &
                                                   c_pcell_face,             &
                                                   c_pvtx_ln_to_gn,          &
                                                   c_pedge_ln_to_gn,         &
                                                   c_pface_ln_to_gn,         &
                                                   c_pcell_ln_to_gn,         &
                                                   c_pn_surface,             &
                                                   c_psurface_face_idx,      &
                                                   c_psurface_face,          &
                                                   c_psurface_face_ln_to_gn, &
                                                   c_pn_ridge,               &
                                                   c_pridge_edge_idx,        &
                                                   c_pridge_edge,            &
                                                   c_pridge_edge_ln_to_gn)

    call c_f_pointer(c_pn_vtx, &
                     pn_vtx,   &
                     [n_part])

    call c_f_pointer(c_pn_edge, &
                     pn_edge,   &
                     [n_part])

    call c_f_pointer(c_pn_face, &
                     pn_face,   &
                     [n_part])

    call c_f_pointer(c_pn_cell, &
                     pn_cell,   &
                     [n_part])

    call c_f_pointer(c_pn_surface, &
                     pn_surface,   &
                     [n_part])

    call c_f_pointer(c_pn_ridge, &
                     pn_ridge,   &
                     [n_part])

    allocate(s_array(n_part)) 

    do i = 1, n_part
      s_array(i) = 3 * pn_vtx(i)
    enddo

    call PDM_pointer_array_create (pvtx_coord,         &
                                   n_part,             &
                                   PDM_TYPE_DOUBLE,    &
                                   c_pvtx_coord,       &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = 2 * pn_edge(i)
    enddo

    call PDM_pointer_array_create (pedge_vtx,          &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_pedge_vtx,        &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = pn_face(i) + 1
    enddo

    call PDM_pointer_array_create (pface_edge_idx,     &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_pface_edge_idx,   &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      call PDM_pointer_array_part_get (pface_edge_idx, &
                                       i-1,            &
                                       ipart_pface_edge_idx)      
      s_array(i) = ipart_pface_edge_idx(pn_face(i) + 1)
    enddo

    call PDM_pointer_array_create (pface_edge,         &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_pface_edge,       &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    call PDM_pointer_array_create (pface_vtx,          &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_pface_vtx,        &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = pn_cell(i) + 1
    enddo

    call PDM_pointer_array_create (pcell_face_idx,     &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_pcell_face_idx,   &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      call PDM_pointer_array_part_get (pcell_face_idx, &
                                       i-1,            &
                                       ipart_pcell_face_idx)      
      s_array(i) = ipart_pcell_face_idx(pn_cell(i) + 1)
    enddo

    call PDM_pointer_array_create (pcell_face,         &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_pcell_face,       &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)
    do i = 1, n_part
      s_array(i) = pn_vtx(i)
    enddo

    call PDM_pointer_array_create (pvtx_ln_to_gn,      &
                                   n_part,             &
                                   PDM_TYPE_G_NUM,     &
                                   c_pvtx_ln_to_gn,    &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = pn_edge(i)
    enddo

    call PDM_pointer_array_create (pedge_ln_to_gn,     &
                                   n_part,             &
                                   PDM_TYPE_G_NUM,     &
                                   c_pedge_ln_to_gn,   &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = pn_face(i)
    enddo

    call PDM_pointer_array_create (pface_ln_to_gn,     &
                                   n_part,             &
                                   PDM_TYPE_G_NUM,     &
                                   c_pface_ln_to_gn,   &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = pn_cell(i)
    enddo

    call PDM_pointer_array_create (pcell_ln_to_gn,     &
                                   n_part,             &
                                   PDM_TYPE_G_NUM,     &
                                   c_pcell_ln_to_gn,   &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = pn_surface(i) + 1
    enddo

    call PDM_pointer_array_create (psurface_face_idx,  &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_psurface_face_idx,&
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      call PDM_pointer_array_part_get (psurface_face_idx, &
                                       i-1,                 &
                                       ipsurface_face_idx)      
      s_array(i) = ipsurface_face_idx(pn_surface(i) + 1)
    enddo

    call PDM_pointer_array_create (psurface_face,      &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_psurface_face,    &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    call PDM_pointer_array_create (psurface_face_ln_to_gn,   &
                                   n_part,                   &
                                   PDM_TYPE_G_NUM,           &
                                   c_psurface_face_ln_to_gn, &
                                   s_array,                  &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      s_array(i) = pn_ridge(i) + 1
    enddo

    call PDM_pointer_array_create (pridge_edge_idx,  &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_pridge_edge_idx,&
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    do i = 1, n_part
      call PDM_pointer_array_part_get (pridge_edge_idx, &
                                       i-1,         &
                                       ipridge_edge_idx)      
      s_array(i) = ipsurface_face_idx(pn_surface(i) + 1)
    enddo

    call PDM_pointer_array_create (pridge_edge,      &
                                   n_part,             &
                                   PDM_TYPE_INT,       &
                                   c_pridge_edge,    &
                                   s_array,            &
                                   PDM_OWNERSHIP_KEEP)

    call PDM_pointer_array_create (pridge_edge_ln_to_gn,   &
                                   n_part,                   &
                                   PDM_TYPE_G_NUM,           &
                                   c_pridge_edge_ln_to_gn, &
                                   s_array,                  &
                                   PDM_OWNERSHIP_KEEP)

    deallocate(s_array)

  end subroutine PDM_generate_mesh_parallelepiped_ngon_

end module pdm_generate_mesh
