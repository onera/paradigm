!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2023  ONERA
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

module pdm_mesh_nodal

  use pdm
  use iso_c_binding

  implicit none

  interface PDM_Mesh_nodal_n_vtx_elt_get ; module procedure &
  PDM_Mesh_nodal_n_vtx_elt_get_
  end interface

  private :: PDM_Mesh_nodal_n_vtx_elt_get_

  interface

    !>
    !!
    !! \brief Get the number of vertices of an element type
    !!
    !! \param [in]   type     Element type
    !! \param [in]   comm     Element order
    !!
    !! \return       Number of vertices
    !!

    function PDM_Mesh_nodal_n_vtx_elt_get_cf(elt_t, &
                                             order) &

      result(n_vtx_per_elt) &
      bind (c, name = 'PDM_Mesh_nodal_n_vtx_elt_get')

      use iso_c_binding
      implicit none

      integer(c_int), value :: elt_t, order
      integer(c_int)        :: n_vtx_per_elt

    end function PDM_Mesh_nodal_n_vtx_elt_get_cf

  end interface

  contains

    !>
    !!
    !! \brief Get the number of vertices of an element type
    !!
    !! \param [in]   type           Element type
    !! \param [in]   comm           Element order
    !!
    !! \param [out]  n_vtx_per_elt  Number of vertices per element
    !!

    subroutine PDM_Mesh_nodal_n_vtx_elt_get_(elt_t, &
                                             order, &
                                             n_vtx_per_elt)

      use iso_c_binding
      implicit none

      integer, intent(in)  :: elt_t, order
      integer, intent(out) :: n_vtx_per_elt

      integer :: c_elt_t, c_order, c_n_vtx_per_elt

      c_elt_t = elt_t
      c_order = order

      c_n_vtx_per_elt = PDM_Mesh_nodal_n_vtx_elt_get_cf(c_elt_t, &
                                                        c_order)

      n_vtx_per_elt = c_n_vtx_per_elt

    end subroutine PDM_Mesh_nodal_n_vtx_elt_get_


end module PDM_mesh_nodal
