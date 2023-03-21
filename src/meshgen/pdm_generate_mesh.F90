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

  private :: PDM_generate_mesh_rectangle_simplified_

  interface

  subroutine PDM_generate_mesh_rectangle_simplified_cf(comm,        &
                                                       n_vtx,       &
                                                       n_elt,       &
                                                       coords,      &
                                                       elt_vtx_idx, &
                                                       elt_vtx)     &

      bind (c, name = 'PDM_generate_mesh_rectangle_simplified')

      use iso_c_binding
      implicit none

      integer(c_int),  value :: comm
      integer(c_int)         :: n_vtx
      integer(c_int)         :: n_elt
      type(c_ptr)            :: coords
      type(c_ptr)            :: elt_vtx_idx
      type(c_ptr)            :: elt_vtx

  end subroutine PDM_generate_mesh_rectangle_simplified_cf

  end interface

  contains

  !>
  !!
  !! \brief Create a simple partitionned rectangle mesh (2D).
  !!
  !! \param [in]   comm        MPI communicator
  !! \param [out]  n_vtx       Number of vertices
  !! \param [out]  n_elt       Number of elements
  !! \param [out]  coords      Array of vertex coordinates
  !! \param [out]  elt_vtx_idx Index array of the element vertex connectivity
  !! \param [out]  elt_vtx     Array of the element vertex connectivity
  !!
  !!

  subroutine PDM_generate_mesh_rectangle_simplified_(comm,        &
                                                     n_vtx,       &
                                                     n_elt,       &
                                                     coords,      &
                                                     elt_vtx_idx, &
                                                     elt_vtx)

      use iso_c_binding
      implicit none

      integer,                     intent(in) :: comm
      integer,                    intent(out) :: n_vtx
      integer,                    intent(out) :: n_elt
      double precision, dimension(:), pointer :: coords
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
                                                     c_n_vtx,       &
                                                     c_n_elt,       &
                                                     c_coords,      &
                                                     c_elt_vtx_idx, &
                                                     c_elt_vtx)

      n_vtx = c_n_vtx
      n_elt = c_n_elt

      call c_f_pointer(c_coords, &
                       coords,   &
                       [3 * n_vtx])

      call c_f_pointer(c_elt_vtx_idx, &
                       elt_vtx_idx,   &
                       [n_elt + 1])

      call c_f_pointer(c_elt_vtx, &
                       elt_vtx,   &
                       [elt_vtx_idx(n_elt + 1)])

  end subroutine PDM_generate_mesh_rectangle_simplified_

end module pdm_generate_mesh
