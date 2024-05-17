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

#include "pdm_configf.h"

module pdm_dmesh_nodal

  use pdm
  use iso_c_binding

  implicit none

  integer(c_int), parameter :: PDM_GEOMETRY_KIND_VOLUMIC  = 0
  integer(c_int), parameter :: PDM_GEOMETRY_KIND_SURFACIC = 1
  integer(c_int), parameter :: PDM_GEOMETRY_KIND_RIDGE    = 2
  integer(c_int), parameter :: PDM_GEOMETRY_KIND_CORNER   = 3

  interface PDM_DMesh_nodal_create ; module procedure &
  PDM_DMesh_nodal_create_
  end interface

  interface PDM_DMesh_nodal_section_add ; module procedure &
  PDM_DMesh_nodal_section_add_
  end interface

  interface PDM_DMesh_nodal_section_poly2d_set ; module procedure &
  PDM_DMesh_nodal_section_poly2d_set_
  end interface

  interface PDM_DMesh_nodal_section_group_elmt_set ; module procedure &
  PDM_DMesh_nodal_section_group_elmt_set_
  end interface

  interface PDM_DMesh_nodal_coord_set ; module procedure &
  PDM_DMesh_nodal_coord_set_
  end interface

  private :: PDM_DMesh_nodal_create_
  private :: PDM_DMesh_nodal_section_add_
  private :: PDM_DMesh_nodal_section_poly2d_set_
  private :: PDM_DMesh_nodal_section_group_elmt_set_
  private :: PDM_DMesh_nodal_coord_set_

  interface

    !>
    !!
    !! \brief Setup global distribution of all elements register in current structure
    !!
    !! \param [inout]   dmn   Pointer to \ref PDM_dmesh_nodal_t object
    !!

    subroutine PDM_dmesh_nodal_generate_distribution (dmn) &
    bind(c, name='PDM_dmesh_nodal_generate_distribution')

      use iso_c_binding
      implicit none

      type (c_ptr), value :: dmn
    end subroutine PDM_dmesh_nodal_generate_distribution

    !>
    !!
    !! \brief Create a Mesh nodal structure
    !!
    !! \param [in]   n_part   Number of partition on the current process
    !! \param [in]   comm     MPI communicator
    !!
    !! \return       New mesh nodal handle
    !!

    function PDM_DMesh_nodal_create_cf (comm, mesh_dimension, n_vtx, n_cell, n_face, n_edge) result(dmn) &
    bind(c, name='PDM_DMesh_nodal_create')

      use iso_c_binding
      implicit none

      integer(c_int),  value :: comm
      integer(c_int),  value :: mesh_dimension
#ifdef PDM_LONG_G_NUM
      integer(c_long), value :: n_vtx
      integer(c_long), value :: n_cell
      integer(c_long), value :: n_face
      integer(c_long), value :: n_edge
#else
      integer(c_int),  value :: n_vtx
      integer(c_int),  value :: n_cell
      integer(c_int),  value :: n_face
      integer(c_int),  value :: n_edge
#endif
      type (c_ptr)           :: dmn
    end function PDM_DMesh_nodal_create_cf

    function PDM_DMesh_nodal_section_add_cf (dmn, geom_kind, t_elt) result(id_section) &
    bind(c, name='PDM_DMesh_nodal_section_add')

      use iso_c_binding
      implicit none

      type (c_ptr),   value :: dmn
      integer(c_int), value :: geom_kind
      integer(c_int), value :: t_elt
      integer(c_int)        :: id_section
    end function PDM_DMesh_nodal_section_add_cf

    !>
    !!
    !! \brief Define a polygon section
    !!
    !! \param [in]  hdl            Distributed nodal mesh handle
    !! \param [in]  id_section     Block identifier
    !! \param [in]  n_elt          Number of elements
    !! \param [in]  connect_idx    Connectivity index (size = \ref n_elt + 1)
    !! \param [in]  connect        Connectivity (size = \ref connect_idx[\ref n_elt])
    !!

    subroutine PDM_DMesh_nodal_section_poly2d_set_cf (dmn,         &
                                                      geom_kind,   &
                                                      id_section,  &
                                                      n_elt,       &
                                                      connec_idx,  &
                                                      connec,      &
                                                      owner)       &
    bind(c, name = 'PDM_DMesh_nodal_section_poly2d_set')

      use iso_c_binding
      implicit none

      type (c_ptr),   value :: dmn
      integer(c_int), value :: geom_kind
      integer(c_int), value :: id_section
      integer(c_int), value :: n_elt
      type(c_ptr),    value :: connec_idx
      type(c_ptr),    value :: connec
      integer(c_int), value :: owner

    end subroutine PDM_DMesh_nodal_section_poly2d_set_cf

    subroutine PDM_DMesh_nodal_section_group_elmt_set_cf (dmn,             &
                                                          geom_kind,       &
                                                          n_group_elmt,    &
                                                          dgroup_elmt_idx, &
                                                          dgroup_elmt,     &
                                                          owner)           &
    bind(c, name = 'PDM_DMesh_nodal_section_group_elmt_set')

      use iso_c_binding
      implicit none

      type (c_ptr),   value :: dmn
      integer(c_int), value :: geom_kind
      integer(c_int), value :: n_group_elmt
      type(c_ptr),    value :: dgroup_elmt_idx
      type(c_ptr),    value :: dgroup_elmt
      integer(c_int), value :: owner

    end subroutine PDM_DMesh_nodal_section_group_elmt_set_cf

    !>
    !!
    !! \brief Define partition vertices
    !!
    !! \param [in]  hdl       Distributed nodal mesh handle
    !! \param [in]  n_vtx     Number of vertices
    !! \param [in]  coords    Interlaced coordinates (size = 3 * \ref n_vtx)
    !!

    subroutine PDM_DMesh_nodal_coord_set_cf (dmn,    &
                                             n_vtx,  &
                                             coords, &
                                             owner)  &
    bind(c, name = 'PDM_DMesh_nodal_coord_set')

      use iso_c_binding
      implicit none

      type (c_ptr),   value :: dmn
      integer(c_int), value :: n_vtx
      type(c_ptr),    value :: coords
      integer(c_int), value :: owner

    end subroutine PDM_DMesh_nodal_coord_set_cf

    !>
    !!
    !! \brief Free \ref PDM_dmesh_nodal_t object
    !!
    !! \param [inout]   dmn   Pointer to \ref PDM_dmesh_nodal_t object
    !!

    subroutine PDM_DMesh_nodal_free (dmn) &
    bind(c, name='PDM_DMesh_nodal_free')

      use iso_c_binding
      implicit none

      type (c_ptr), value :: dmn
    end subroutine PDM_DMesh_nodal_free

  end interface

  contains

  !>
  !!
  !! \brief Create a Mesh nodal structure
  !!
  !! \param [in]   n_part   Number of partition on the current process
  !! \param [in]   comm     MPI communicator
  !!
  !! \return       New mesh nodal handle
  !!

  subroutine PDM_DMesh_nodal_create_ (dmn,            &
                                      f_comm,         &
                                      mesh_dimension, &
                                      n_vtx,          &
                                      n_cell,         &
                                      n_face,         &
                                      n_edge)

    use iso_c_binding
    implicit none

    type (c_ptr) :: dmn
    integer, intent(in) :: f_comm
    integer, intent(in) :: mesh_dimension
#ifdef PDM_LONG_G_NUM
    integer(PDM_g_num_s), intent(in) :: n_vtx
    integer(PDM_g_num_s), intent(in) :: n_cell
    integer(PDM_g_num_s), intent(in) :: n_face
    integer(PDM_g_num_s), intent(in) :: n_edge
#else
    integer(PDM_l_num_s), intent(in) :: n_vtx
    integer(PDM_l_num_s), intent(in) :: n_cell
    integer(PDM_l_num_s), intent(in) :: n_face
    integer(PDM_l_num_s), intent(in) :: n_edge
#endif

  integer(c_int) :: c_comm

  c_comm = PDM_MPI_Comm_f2c(f_comm)

  dmn = PDM_DMesh_nodal_create_cf(c_comm,         &
                                  mesh_dimension, &
                                  n_vtx,          &
                                  n_cell,         &
                                  n_face,         &
                                  n_edge)

  end subroutine PDM_DMesh_nodal_create_

  subroutine PDM_DMesh_nodal_section_add_ (dmn,       &
                                           geom_kind, &
                                           t_elt,     &
                                           id_section)

    use iso_c_binding
    implicit none

    type (c_ptr), value :: dmn
    integer, intent(in) :: geom_kind
    integer, intent(in) :: t_elt
    integer             :: id_section

    id_section = PDM_DMesh_nodal_section_add_cf (dmn,       &
                                                 geom_kind, &
                                                 t_elt)

  end subroutine PDM_DMesh_nodal_section_add_

  !>
  !!
  !! \brief Define a polygon section
  !!
  !! \param [in]  hdl            Distributed nodal mesh handle
  !! \param [in]  id_section     Block identifier
  !! \param [in]  n_elt          Number of elements
  !! \param [in]  connect_idx    Connectivity index (size = \ref n_elt + 1)
  !! \param [in]  connect        Connectivity (size = \ref connect_idx[\ref n_elt])
  !!

  subroutine PDM_DMesh_nodal_section_poly2d_set_ (dmn,        &
                                                  geom_kind,  &
                                                  id_section, &
                                                  n_elt,      &
                                                  connec_idx, &
                                                  connec,     &
                                                  owner)

    use iso_c_binding
    implicit none

    type (c_ptr), value           :: dmn
    integer, intent(in)           :: geom_kind
    integer, intent(in)           :: id_section
    integer, intent(in)           :: n_elt
    integer, intent(in)           :: owner
    integer(PDM_l_num_s), pointer :: connec_idx(:)
    integer(PDM_g_num_s), pointer :: connec(:)

    integer(c_int) :: c_n_elt
    type(c_ptr)    :: c_connec_idx
    type(c_ptr)    :: c_connec

    c_n_elt = n_elt

    c_connec_idx = C_NULL_PTR
    if (associated(connec_idx)) then
      c_connec_idx = c_loc(connec_idx)
    end if

    c_connec = C_NULL_PTR
    if (associated(connec)) then
      c_connec = c_loc(connec)
    end if

    call PDM_DMesh_nodal_section_poly2d_set_cf(dmn,          &
                                               geom_kind,    &
                                               id_section,   &
                                               c_n_elt,      &
                                               c_connec_idx, &
                                               c_connec,     &
                                               owner)

  end subroutine PDM_DMesh_nodal_section_poly2d_set_

  subroutine PDM_DMesh_nodal_section_group_elmt_set_ (dmn,             &
                                                      geom_kind,       &
                                                      n_group_elmt,    &
                                                      dgroup_elmt_idx, &
                                                      dgroup_elmt,     &
                                                      owner)

    use iso_c_binding
    implicit none

    type (c_ptr), value           :: dmn
    integer, intent(in)           :: geom_kind
    integer, intent(in)           :: n_group_elmt
    integer, intent(in)           :: owner
    integer(PDM_l_num_s), pointer :: dgroup_elmt_idx(:)
    integer(PDM_g_num_s), pointer :: dgroup_elmt(:)

    integer(c_int) :: c_n_group_elmt
    type(c_ptr)    :: c_dgroup_elmt_idx
    type(c_ptr)    :: c_dgroup_elmt

    c_n_group_elmt = n_group_elmt

    c_dgroup_elmt_idx = C_NULL_PTR
    if (associated(dgroup_elmt_idx)) then
      c_dgroup_elmt_idx = c_loc(dgroup_elmt_idx)
    end if

    c_dgroup_elmt = C_NULL_PTR
    if (associated(dgroup_elmt)) then
      c_dgroup_elmt = c_loc(dgroup_elmt)
    end if

    call PDM_DMesh_nodal_section_group_elmt_set_cf(dmn,               &
                                                   geom_kind,         &
                                                   c_n_group_elmt,    &
                                                   c_dgroup_elmt_idx, &
                                                   c_dgroup_elmt,     &
                                                   owner)

  end subroutine PDM_DMesh_nodal_section_group_elmt_set_

  !>
  !!
  !! \brief Define partition vertices
  !!
  !! \param [in]  hdl       Distributed nodal mesh handle
  !! \param [in]  n_vtx     Number of vertices
  !! \param [in]  coords    Interlaced coordinates (size = 3 * \ref n_vtx)
  !!

  subroutine PDM_DMesh_nodal_coord_set_ (dmn,    &
                                         n_vtx,  &
                                         coords, &
                                         owner)

    use iso_c_binding
    implicit none

    type (c_ptr), value       :: dmn
    integer, intent(in)       :: n_vtx
    integer, intent(in)       :: owner
    double precision, pointer :: coords(:,:)

    integer(c_int) :: c_n_vtx
    type(c_ptr)    :: c_coords

    c_n_vtx = n_vtx

    c_coords = C_NULL_PTR
    if (associated(coords)) then
      c_coords = c_loc(coords)
    end if

    call PDM_DMesh_nodal_coord_set_cf(dmn,      &
                                      c_n_vtx,  &
                                      c_coords, &
                                      owner)

  end subroutine PDM_DMesh_nodal_coord_set_

end module pdm_dmesh_nodal
