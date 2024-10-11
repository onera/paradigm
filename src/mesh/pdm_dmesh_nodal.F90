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

  interface PDM_dmesh_nodal_sections_id_get ; module procedure &
  PDM_dmesh_nodal_sections_id_get_
  end interface

  interface PDM_DMesh_nodal_section_g_dims_get ; module procedure &
  PDM_DMesh_nodal_section_g_dims_get_
  end interface

  interface PDM_DMesh_nodal_section_add ; module procedure &
  PDM_DMesh_nodal_section_add_
  end interface

  interface PDM_DMesh_nodal_section_std_set ; module procedure &
  PDM_DMesh_nodal_section_std_set_
  end interface

  interface PDM_DMesh_nodal_section_std_get ; module procedure &
  PDM_DMesh_nodal_section_std_get_
  end interface

  interface PDM_DMesh_nodal_section_poly2d_set ; module procedure &
  PDM_DMesh_nodal_section_poly2d_set_
  end interface

  interface PDM_DMesh_nodal_section_poly2d_get ; module procedure &
  PDM_DMesh_nodal_section_poly2d_get_
  end interface

  interface PDM_DMesh_nodal_section_poly3d_set ; module procedure &
  PDM_DMesh_nodal_section_poly3d_set_
  end interface

  interface PDM_DMesh_nodal_section_poly3d_get ; module procedure &
  PDM_DMesh_nodal_section_poly3d_get_
  end interface

  interface PDM_DMesh_nodal_section_group_elmt_set ; module procedure &
  PDM_DMesh_nodal_section_group_elmt_set_
  end interface

  interface PDM_DMesh_nodal_section_group_elmt_get ; module procedure &
  PDM_DMesh_nodal_section_group_elmt_get_
  end interface

  interface PDM_DMesh_nodal_coord_set ; module procedure &
  PDM_DMesh_nodal_coord_set_
  end interface

  interface PDM_DMesh_nodal_vtx_get ; module procedure &
  PDM_DMesh_nodal_vtx_get_
  end interface

  private :: PDM_DMesh_nodal_create_
  private :: PDM_DMesh_nodal_section_g_dims_get_
  private :: PDM_DMesh_nodal_section_add_
  private :: PDM_DMesh_nodal_section_std_set_
  private :: PDM_DMesh_nodal_section_std_get_
  private :: PDM_DMesh_nodal_section_poly2d_set_
  private :: PDM_DMesh_nodal_section_poly2d_get_
  private :: PDM_DMesh_nodal_section_poly3d_set_
  private :: PDM_DMesh_nodal_section_poly3d_get_
  private :: PDM_DMesh_nodal_section_group_elmt_set_
  private :: PDM_DMesh_nodal_section_group_elmt_get_
  private :: PDM_DMesh_nodal_coord_set_
  private :: PDM_DMesh_nodal_vtx_get_

  interface

    function PDM_dmesh_nodal_n_vtx_get (dmn)          &
                                        result(n_vtx) &
    bind(c, name='PDM_DMesh_nodal_n_vtx_get')

      use iso_c_binding
      implicit none

      type (c_ptr), value :: dmn
      integer(c_int)      :: n_vtx
    end function PDM_dmesh_nodal_n_vtx_get

    function PDM_dmesh_nodal_n_section_get (dmn,              &
                                            geom_kind)        &
                                            result(n_section) &
    bind(c, name='PDM_DMesh_nodal_n_section_get')

      use iso_c_binding
      implicit none

      type (c_ptr),   value :: dmn
      integer(c_int), value :: geom_kind
      integer(c_int)        :: n_section
    end function PDM_dmesh_nodal_n_section_get

    function PDM_dmesh_nodal_section_n_elt_get (dmn,          &
                                                geom_kind,    &
                                                id_section)   &
                                                result(n_elt) &
    bind(c, name='PDM_DMesh_nodal_section_n_elt_get')

      use iso_c_binding
      implicit none

      type (c_ptr),   value :: dmn
      integer(c_int), value :: geom_kind
      integer(c_int), value :: id_section
      integer(c_int)        :: n_elt
    end function PDM_dmesh_nodal_section_n_elt_get

    function PDM_DMesh_nodal_section_elt_type_get (dmn,          &
                                                   geom_kind,    &
                                                   id_section)   &
                                                   result(elt_type) &
    bind(c, name='PDM_DMesh_nodal_section_elt_type_get')

      use iso_c_binding
      implicit none

      type (c_ptr),   value :: dmn
      integer(c_int), value :: geom_kind
      integer(c_int), value :: id_section
      integer(c_int)        :: elt_type
    end function PDM_DMesh_nodal_section_elt_type_get

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

    function PDM_dmesh_nodal_sections_id_get_cf (dmn,               &
                                                 geom_kind)         &
                                                 result(id_section) &
    bind(c, name='PDM_DMesh_nodal_sections_id_get')

      use iso_c_binding
      implicit none

      type (c_ptr),   value :: dmn
      integer(c_int), value :: geom_kind
      type (c_ptr)          :: id_section
    end function PDM_dmesh_nodal_sections_id_get_cf

    subroutine PDM_DMesh_nodal_section_g_dims_get_cf (dmn,    &
                                                      n_cell, &
                                                      n_face, &
                                                      n_edge, &
                                                      n_vtx)  &
    bind(c, name='PDM_DMesh_nodal_section_g_dims_get')

      use iso_c_binding
      implicit none

      type (c_ptr), value :: dmn
#ifdef PDM_LONG_G_NUM
      integer(c_long)     :: n_cell
      integer(c_long)     :: n_face
      integer(c_long)     :: n_edge
      integer(c_long)     :: n_vtx
#else
      integer(c_int)      :: n_cell
      integer(c_int)      :: n_face
      integer(c_int)      :: n_edge
      integer(c_int)      :: n_vtx
#endif

    end subroutine PDM_DMesh_nodal_section_g_dims_get_cf

    function PDM_DMesh_nodal_section_add_cf (dmn, geom_kind, t_elt) result(id_section) &
    bind(c, name='PDM_DMesh_nodal_section_add')

      use iso_c_binding
      implicit none

      type (c_ptr),   value :: dmn
      integer(c_int), value :: geom_kind
      integer(c_int), value :: t_elt
      integer(c_int)        :: id_section
    end function PDM_DMesh_nodal_section_add_cf

    subroutine PDM_DMesh_nodal_section_std_set_cf (dmn,         &
                                                   geom_kind,   &
                                                   id_section,  &
                                                   n_elt,       &
                                                   connec,      &
                                                   owner)       &
    bind(c, name = 'PDM_DMesh_nodal_section_std_set')

      use iso_c_binding
      implicit none

      type (c_ptr),   value :: dmn
      integer(c_int), value :: geom_kind
      integer(c_int), value :: id_section
      integer(c_int), value :: n_elt
      type(c_ptr),    value :: connec
      integer(c_int), value :: owner

    end subroutine PDM_DMesh_nodal_section_std_set_cf

    function PDM_DMesh_nodal_section_std_get_cf (dmn,           &
                                                 geom_kind,     &
                                                 id_section,    &
                                                 owner)         &
                                                 result(connec) &
    bind(c, name = 'PDM_DMesh_nodal_section_std_get')

      use iso_c_binding
      implicit none

      type (c_ptr),   value :: dmn
      integer(c_int), value :: geom_kind
      integer(c_int), value :: id_section
      type(c_ptr)           :: connec
      integer(c_int), value :: owner

    end function PDM_DMesh_nodal_section_std_get_cf

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

    subroutine PDM_DMesh_nodal_section_poly2d_get_cf (dmn,         &
                                                      geom_kind,   &
                                                      id_section,  &
                                                      connec_idx,  &
                                                      connec,      &
                                                      owner)       &
    bind(c, name = 'PDM_DMesh_nodal_section_poly2d_get')

      use iso_c_binding
      implicit none

      type (c_ptr),   value :: dmn
      integer(c_int), value :: geom_kind
      integer(c_int), value :: id_section
      type(c_ptr)           :: connec_idx
      type(c_ptr)           :: connec
      integer(c_int), value :: owner

    end subroutine PDM_DMesh_nodal_section_poly2d_get_cf

    !>
    !!
    !! \brief Define a polyhedra section
    !!
    !! \param [in]  hdl            Distributed nodal mesh handle
    !! \param [in]  id_section       Block identifier
    !! \param [in]  n_elt          Number of polyhedra
    !! \param [in]  n_face         Number of faces used to describe polyhedra
    !! \param [in]  facvtx_idx     Index of face vertex connectivity
    !! \param [in]  facvtx         Face vertex connectivity
    !! \param [in]  cellfac_idx    Index of cell face connectivity
    !! \param [in]  cellfac        Cell face connectivity
    !!

    subroutine PDM_DMesh_nodal_section_poly3d_set_cf (dmn,           &
                                                      geom_kind,     &
                                                      id_section,    &
                                                      n_elt,         &
                                                      n_face,        &
                                                      face_vtx_idx,  &
                                                      face_vtx,      &
                                                      cell_face_idx, &
                                                      cell_face,     &
                                                      owner)         &
    bind(c, name = 'PDM_DMesh_nodal_section_poly3d_set')

      use iso_c_binding
      implicit none

      type (c_ptr),   value :: dmn
      integer(c_int), value :: geom_kind
      integer(c_int), value :: id_section
      integer(c_int), value :: n_elt
      integer(c_int), value :: n_face
      type(c_ptr),    value :: face_vtx_idx
      type(c_ptr),    value :: face_vtx
      type(c_ptr),    value :: cell_face_idx
      type(c_ptr),    value :: cell_face
      integer(c_int), value :: owner

    end subroutine PDM_DMesh_nodal_section_poly3d_set_cf

    subroutine PDM_DMesh_nodal_section_poly3d_get_cf (dmn,           &
                                                      geom_kind,     &
                                                      id_section,    &
                                                      n_face,        &
                                                      face_vtx_idx,  &
                                                      face_vtx,      &
                                                      cell_face_idx, &
                                                      cell_face,     &
                                                      owner)         &
    bind(c, name = 'PDM_DMesh_nodal_section_poly3d_get')

      use iso_c_binding
      implicit none

      type (c_ptr),   value :: dmn
      integer(c_int), value :: geom_kind
      integer(c_int), value :: id_section
      integer(c_int)        :: n_face
      type(c_ptr)           :: face_vtx_idx
      type(c_ptr)           :: face_vtx
      type(c_ptr)           :: cell_face_idx
      type(c_ptr)           :: cell_face
      integer(c_int), value :: owner

    end subroutine PDM_DMesh_nodal_section_poly3d_get_cf

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

    subroutine PDM_DMesh_nodal_section_group_elmt_get_cf (dmn,             &
                                                          geom_kind,       &
                                                          n_group_elmt,    &
                                                          dgroup_elmt_idx, &
                                                          dgroup_elmt,     &
                                                          owner)           &
    bind(c, name = 'PDM_DMesh_nodal_section_group_elmt_get')

      use iso_c_binding
      implicit none

      type (c_ptr),   value :: dmn
      integer(c_int), value :: geom_kind
      integer(c_int)        :: n_group_elmt
      type(c_ptr)           :: dgroup_elmt_idx
      type(c_ptr)           :: dgroup_elmt
      integer(c_int), value :: owner

    end subroutine PDM_DMesh_nodal_section_group_elmt_get_cf

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

    function PDM_DMesh_nodal_vtx_get_cf (dmn,           &
                                         owner)         &
                                         result(coords) &
    bind(c, name = 'PDM_DMesh_nodal_vtx_get')

      use iso_c_binding
      implicit none

      type (c_ptr),   value :: dmn
      integer(c_int), value :: owner
      type(c_ptr)           :: coords

    end function PDM_DMesh_nodal_vtx_get_cf

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

  subroutine PDM_dmesh_nodal_sections_id_get_ (dmn,       &
                                               geom_kind, &
                                               id_section)

    use iso_c_binding
    implicit none

    type (c_ptr), value           :: dmn
    integer, intent(in)           :: geom_kind
    integer(PDM_l_num_s), pointer :: id_section(:)

    integer(c_int) :: c_n_section
    type(c_ptr)    :: c_id_section

    c_n_section = PDM_DMesh_nodal_n_section_get(dmn,       &
                                                geom_kind)

    c_id_section = PDM_dmesh_nodal_sections_id_get_cf(dmn,       &
                                                      geom_kind)

    call c_f_pointer(c_id_section, &
                     id_section,   &
                     [c_n_section])

  end subroutine PDM_dmesh_nodal_sections_id_get_

  subroutine PDM_DMesh_nodal_section_g_dims_get_ (dmn,    &
                                                  n_cell, &
                                                  n_face, &
                                                  n_edge, &
                                                  n_vtx)

    use iso_c_binding
    implicit none

    type (c_ptr), value  :: dmn
#ifdef PDM_LONG_G_NUM
    integer(PDM_g_num_s) :: n_cell
    integer(PDM_g_num_s) :: n_face
    integer(PDM_g_num_s) :: n_edge
    integer(PDM_g_num_s) :: n_vtx
#else
    integer(PDM_l_num_s) :: n_cell
    integer(PDM_l_num_s) :: n_face
    integer(PDM_l_num_s) :: n_edge
    integer(PDM_l_num_s) :: n_vtx
#endif

    call PDM_DMesh_nodal_section_g_dims_get_cf(dmn,    &
                                               n_cell, &
                                               n_face, &
                                               n_edge, &
                                               n_vtx)

  end subroutine PDM_DMesh_nodal_section_g_dims_get_

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

  subroutine PDM_DMesh_nodal_section_std_set_ (dmn,        &
                                               geom_kind,  &
                                               id_section, &
                                               n_elt,      &
                                               connec,     &
                                               owner)

    use iso_c_binding
    implicit none

    type (c_ptr), value           :: dmn
    integer, intent(in)           :: geom_kind
    integer, intent(in)           :: id_section
    integer, intent(in)           :: n_elt
    integer, intent(in)           :: owner
    integer(PDM_g_num_s), pointer :: connec(:)

    integer(c_int) :: c_n_elt
    type(c_ptr)    :: c_connec

    c_n_elt = n_elt

    c_connec = C_NULL_PTR
    if (associated(connec)) then
      c_connec = c_loc(connec)
    end if

    call PDM_DMesh_nodal_section_std_set_cf(dmn,        &
                                            geom_kind,  &
                                            id_section, &
                                            c_n_elt,    &
                                            c_connec,   &
                                            owner)

  end subroutine PDM_DMesh_nodal_section_std_set_

  subroutine PDM_DMesh_nodal_section_std_get_ (dmn,        &
                                               geom_kind,  &
                                               id_section, &
                                               connec,     &
                                               owner)

    use iso_c_binding
    implicit none

    type (c_ptr), value           :: dmn
    integer, intent(in)           :: geom_kind
    integer, intent(in)           :: id_section
    integer, intent(in)           :: owner
    integer(PDM_g_num_s), pointer :: connec(:)

    integer(c_int) :: c_n_elt
    integer(c_int) :: c_n_elt_size
    type(c_ptr)    :: c_connec

    c_n_elt = PDM_DMesh_nodal_section_n_elt_get(dmn,       &
                                                geom_kind, &
                                                id_section)

    c_n_elt_size = PDM_DMesh_nodal_section_elt_type_get(dmn,       &
                                                        geom_kind, &
                                                        id_section)

    c_connec = PDM_DMesh_nodal_section_std_get_cf(dmn,        &
                                                  geom_kind,  &
                                                  id_section, &
                                                  owner)

    call c_f_pointer(c_connec, &
                     connec,   &
                     [c_n_elt_size*c_n_elt])

  end subroutine PDM_DMesh_nodal_section_std_get_

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

  subroutine PDM_DMesh_nodal_section_poly2d_get_ (dmn,        &
                                                  geom_kind,  &
                                                  id_section, &
                                                  connec_idx, &
                                                  connec,     &
                                                  owner)

    use iso_c_binding
    implicit none

    type (c_ptr), value           :: dmn
    integer, intent(in)           :: geom_kind
    integer, intent(in)           :: id_section
    integer, intent(in)           :: owner
    integer(PDM_l_num_s), pointer :: connec_idx(:)
    integer(PDM_g_num_s), pointer :: connec(:)

    integer(c_int) :: c_n_elt
    type(c_ptr)    :: c_connec_idx
    type(c_ptr)    :: c_connec

    c_n_elt = PDM_DMesh_nodal_section_n_elt_get(dmn,       &
                                                geom_kind, &
                                                id_section)

    call PDM_DMesh_nodal_section_poly2d_get_cf(dmn,          &
                                               geom_kind,    &
                                               id_section,   &
                                               c_connec_idx, &
                                               c_connec,     &
                                               owner)

    call c_f_pointer(c_connec_idx, &
                     connec_idx,   &
                     [c_n_elt+1])

    call c_f_pointer(c_connec, &
                     connec,   &
                     [connec_idx(c_n_elt+1)])

  end subroutine PDM_DMesh_nodal_section_poly2d_get_

  !>
  !!
  !! \brief Define a polyhedra section
  !!
  !! \param [in]  hdl            Distributed nodal mesh handle
  !! \param [in]  id_section       Block identifier
  !! \param [in]  n_elt          Number of polyhedra
  !! \param [in]  n_face         Number of faces used to describe polyhedra
  !! \param [in]  facvtx_idx     Index of face vertex connectivity
  !! \param [in]  facvtx         Face vertex connectivity
  !! \param [in]  cellfac_idx    Index of cell face connectivity
  !! \param [in]  cellfac        Cell face connectivity
  !!

  subroutine PDM_DMesh_nodal_section_poly3d_set_ (dmn,           &
                                                  geom_kind,     &
                                                  id_section,    &
                                                  n_elt,         &
                                                  n_face,        &
                                                  face_vtx_idx,  &
                                                  face_vtx,      &
                                                  cell_face_idx, &
                                                  cell_face,     &
                                                  owner)
    use iso_c_binding
    implicit none

    type (c_ptr), value           :: dmn
    integer, intent(in)           :: geom_kind
    integer, intent(in)           :: id_section
    integer, intent(in)           :: n_elt
    integer, intent(in)           :: n_face
    integer, intent(in)           :: owner
    integer(PDM_l_num_s), pointer :: face_vtx_idx(:)
    integer(PDM_g_num_s), pointer :: face_vtx(:)
    integer(PDM_l_num_s), pointer :: cell_face_idx(:)
    integer(PDM_g_num_s), pointer :: cell_face(:)

    integer(c_int) :: c_n_elt
    integer(c_int) :: c_n_face
    type(c_ptr)    :: c_face_vtx_idx
    type(c_ptr)    :: c_face_vtx
    type(c_ptr)    :: c_cell_face_idx
    type(c_ptr)    :: c_cell_face

    c_n_elt  = n_elt
    c_n_face = n_face

    c_face_vtx_idx = C_NULL_PTR
    if (associated(face_vtx_idx)) then
      c_face_vtx_idx = c_loc(face_vtx_idx)
    end if

    c_face_vtx = C_NULL_PTR
    if (associated(face_vtx)) then
      c_face_vtx = c_loc(face_vtx)
    end if

    c_cell_face_idx = C_NULL_PTR
    if (associated(cell_face_idx)) then
      c_cell_face_idx = c_loc(cell_face_idx)
    end if

    c_cell_face = C_NULL_PTR
    if (associated(cell_face)) then
      c_cell_face = c_loc(cell_face)
    end if

    call PDM_DMesh_nodal_section_poly3d_set_cf(dmn,             &
                                               geom_kind,       &
                                               id_section,      &
                                               c_n_elt,         &
                                               c_n_face,        &
                                               c_face_vtx_idx,  &
                                               c_face_vtx,      &
                                               c_cell_face_idx, &
                                               c_cell_face,     &
                                               owner)

  end subroutine PDM_DMesh_nodal_section_poly3d_set_

  subroutine PDM_DMesh_nodal_section_poly3d_get_ (dmn,           &
                                                  geom_kind,     &
                                                  id_section,    &
                                                  n_face,        &
                                                  face_vtx_idx,  &
                                                  face_vtx,      &
                                                  cell_face_idx, &
                                                  cell_face,     &
                                                  owner)
    use iso_c_binding
    implicit none

    type (c_ptr), value           :: dmn
    integer, intent(in)           :: geom_kind
    integer, intent(in)           :: id_section
    integer                       :: n_face
    integer, intent(in)           :: owner
    integer(PDM_l_num_s), pointer :: face_vtx_idx(:)
    integer(PDM_g_num_s), pointer :: face_vtx(:)
    integer(PDM_l_num_s), pointer :: cell_face_idx(:)
    integer(PDM_g_num_s), pointer :: cell_face(:)

    integer(c_int) :: c_n_elt
    integer(c_int) :: c_n_face
    type(c_ptr)    :: c_face_vtx_idx
    type(c_ptr)    :: c_face_vtx
    type(c_ptr)    :: c_cell_face_idx
    type(c_ptr)    :: c_cell_face

    c_n_elt = PDM_DMesh_nodal_section_n_elt_get(dmn,       &
                                                geom_kind, &
                                                id_section)

    call PDM_DMesh_nodal_section_poly3d_get_cf(dmn,             &
                                               geom_kind,       &
                                               id_section,      &
                                               c_n_face,        &
                                               c_face_vtx_idx,  &
                                               c_face_vtx,      &
                                               c_cell_face_idx, &
                                               c_cell_face,     &
                                               owner)

    n_face = c_n_face

    call c_f_pointer(c_face_vtx_idx, &
                     face_vtx_idx,   &
                     [c_n_face+1])

    call c_f_pointer(c_face_vtx, &
                     face_vtx,   &
                     [face_vtx_idx(c_n_face+1)])

    call c_f_pointer(c_cell_face_idx, &
                     cell_face_idx,   &
                     [c_n_elt+1])

    call c_f_pointer(c_cell_face, &
                     cell_face,   &
                     [cell_face_idx(c_n_elt+1)])

  end subroutine PDM_DMesh_nodal_section_poly3d_get_

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

  subroutine PDM_DMesh_nodal_section_group_elmt_get_ (dmn,             &
                                                      geom_kind,       &
                                                      n_group_elmt,    &
                                                      dgroup_elmt_idx, &
                                                      dgroup_elmt,     &
                                                      owner)

    use iso_c_binding
    implicit none

    type (c_ptr), value           :: dmn
    integer, intent(in)           :: geom_kind
    integer                       :: n_group_elmt
    integer, intent(in)           :: owner
    integer(PDM_l_num_s), pointer :: dgroup_elmt_idx(:)
    integer(PDM_g_num_s), pointer :: dgroup_elmt(:)

    integer(c_int) :: c_n_group_elmt
    type(c_ptr)    :: c_dgroup_elmt_idx
    type(c_ptr)    :: c_dgroup_elmt

    call PDM_DMesh_nodal_section_group_elmt_get_cf(dmn,               &
                                                   geom_kind,         &
                                                   c_n_group_elmt,    &
                                                   c_dgroup_elmt_idx, &
                                                   c_dgroup_elmt,     &
                                                   owner)

    n_group_elmt = c_n_group_elmt

    call c_f_pointer(c_dgroup_elmt_idx, &
                     dgroup_elmt_idx,   &
                     [n_group_elmt+1])

    call c_f_pointer(c_dgroup_elmt, &
                     dgroup_elmt,   &
                     [dgroup_elmt_idx(n_group_elmt+1)])

  end subroutine PDM_DMesh_nodal_section_group_elmt_get_

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

  subroutine PDM_DMesh_nodal_vtx_get_ (dmn,    &
                                       coords, &
                                       owner)

    use iso_c_binding
    implicit none

    type (c_ptr), value       :: dmn
    integer(c_int)            :: c_n_vtx
    integer, intent(in)       :: owner
    double precision, pointer :: coords(:,:)

    type(c_ptr) :: c_coords

    c_n_vtx = PDM_DMesh_nodal_n_vtx_get(dmn)

    c_coords = PDM_DMesh_nodal_vtx_get_cf(dmn, &
                                          owner)

    call c_f_pointer(c_coords, &
                     coords,   &
                     [3,c_n_vtx])

  end subroutine PDM_DMesh_nodal_vtx_get_

end module pdm_dmesh_nodal
