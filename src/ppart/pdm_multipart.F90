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

! WARNING: not tested yet / need to make fortran interface of dmesh for that

#include "pdm_configf.h"

module pdm_multipart

  use pdm

  implicit none

  integer, parameter :: PDM_PART_SIZE_HOMOGENEOUS   = 1
  integer, parameter :: PDM_PART_SIZE_HETEROGENEOUS = 2

  interface PDM_multipart_create ; module procedure  &
    PDM_multipart_create_
  end interface

  interface PDM_multipart_register_joins ; module procedure  &
    PDM_multipart_register_joins_
  end interface

  interface PDM_multipart_set_reordering_options ; module procedure  &
    PDM_multipart_set_reordering_options_
  end interface

  interface PDM_multipart_set_reordering_options_vtx ; module procedure  &
    PDM_multipart_set_reordering_options_vtx_
  end interface

  interface PDM_multipart_get_part_mesh_nodal ; module procedure  &
    PDM_multipart_get_part_mesh_nodal_
  end interface

  interface PDM_multipart_dconnectivity_set ; module procedure  &
    PDM_multipart_dconnectivity_set_
  end interface

  interface PDM_multipart_dvtx_coord_set ; module procedure  &
    PDM_multipart_dvtx_coord_set_
  end interface

  interface PDM_multipart_part_dim_get ; module procedure  &
    PDM_multipart_part_dim_get_
  end interface

  interface PDM_multipart_part_graph_comm_vtx_dim_get ; module procedure  &
    PDM_multipart_part_graph_comm_vtx_dim_get_
  end interface

  interface PDM_multipart_part_val_get ; module procedure  &
    PDM_multipart_part_val_get_
  end interface

  interface PDM_multipart_part_connectivity_get ; module procedure  &
    PDM_multipart_part_connectivity_get_
  end interface

  interface PDM_multipart_part_ln_to_gn_get ; module procedure  &
    PDM_multipart_part_ln_to_gn_get_
  end interface

  interface PDM_multipart_partition_color_get ; module procedure  &
    PDM_multipart_partition_color_get_
  end interface

  interface PDM_multipart_part_color_get ; module procedure  &
    PDM_multipart_part_color_get_
  end interface

  interface PDM_multipart_part_ghost_infomation_get ; module procedure  &
    PDM_multipart_part_ghost_infomation_get_
  end interface

  interface PDM_multipart_part_vtx_coord_get ; module procedure  &
    PDM_multipart_part_vtx_coord_get_
  end interface

  interface PDM_multipart_bound_get ; module procedure  &
    PDM_multipart_bound_get_
  end interface

interface

  !>
  !!
  !! \brief Build a multipart structure
  !!
  !! \param [in]   n_zone           Number of zones in the original mesh
  !! \param [in]   n_part           Number of partition per proc in each zone
  !! \param [in]   merge_blocks     Merge or not the zones before splitting
  !! \param [in]   split_method     Choice of library used to split the mesh
  !! \param [in]   part_size_method Choice of homogeneous or heterogeneous partitions
  !! \param [in]   part_fraction    Weight (in %) of each partition in heterogeneous case
  !! \param [in]   comm             PDM_MPI communicator
  !! \param [in]   owner            Data ownership
  !!
  !! \return   multipart   Pointer to a new \ref PDM_multipart_t object
  !!

  function PDM_multipart_create_c (n_zone, &
                                   n_part, &
                                   merge_blocks, &
                                   split_method, &
                                   part_size_method, &
                                   part_fraction, &
                                   comm, &
                                   owner) &

  result(multipart) &
  bind (c, name='PDM_multipart_create')

    use iso_c_binding
    implicit none

    type(c_ptr)            :: multipart
    integer(c_int), value  :: n_zone
    type(c_ptr),    value  :: n_part
    integer(c_int), value  :: merge_blocks
    integer(c_int), value  :: split_method
    integer(c_int), value  :: part_size_method
    type(c_ptr),    value  :: part_fraction
    integer(c_int), value  :: comm
    integer(c_int), value  :: owner

  end function PDM_multipart_create_c

  !>
  !!
  !! \brief Set connecting data between all the zones
  !!
  !! \param [in]   multipart        Pointer to \ref PDM_multipart_t object
  !! \param [in]   n_total_joins    Total number of interfaces
  !! \param [in]   join_to_opposite For each global join id, give the global id
  !!                                of the opposite join (size = n_total_joins)
  !!

  subroutine PDM_multipart_register_joins_c (multipart, &
                                             n_total_joins, &
                                             join_to_opposite) &
  bind (c, name='PDM_multipart_register_joins')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: n_total_joins
    type(c_ptr),    value  :: join_to_opposite

  end subroutine PDM_multipart_register_joins_c

  !>
  !!
  !! \brief Set the reordering methods to be used after partitioning
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   renum_cell_method     Choice of renumbering method for cells
  !! \param [in]   renum_cell_properties Parameters used by cacheblocking method :
  !!                                     [n_cell_per_cache_wanted, is_asynchrone, is_vectorisation,
  !!                                     n_vect_face, split_method]
  !! \param [in]   renum_face_method     Choice of renumbering method for faces
  !!

  subroutine PDM_multipart_set_reordering_options_c (multipart, &
                                                     i_zone, &
                                                     renum_cell_method, &
                                                     renum_cell_properties, &
                                                     renum_face_method) &
  bind (c, name='PDM_multipart_set_reordering_options')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_zone
    character(c_char)      :: renum_cell_method(*)
    type(c_ptr),    value  :: renum_cell_properties
    character(c_char)      :: renum_face_method(*)

  end subroutine PDM_multipart_set_reordering_options_c

  !>
  !!
  !! \brief Set the reordering methods to be used after partitioning
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   renum_vtx_method      Choice of renumbering method for vertices

  subroutine PDM_multipart_set_reordering_options_vtx_c (multipart, &
                                                         i_zone, &
                                                         renum_vtx_method) &
  bind (c, name='PDM_multipart_set_reordering_options_vtx')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_zone
    character(c_char)      :: renum_vtx_method(*)

  end subroutine PDM_multipart_set_reordering_options_vtx_c

  !>
  !!
  !! \brief ???
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   pmesh_nodal           Partitionned nodal mesh
  !! \param [in]   ownership             Data ownership
  !!

  subroutine PDM_multipart_get_part_mesh_nodal_c (multipart, &
                                                  i_zone, &
                                                  pmesh_nodal, &
                                                  ownership) &
  bind (c, name='PDM_multipart_get_part_mesh_nodal')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_zone
    type(c_ptr),    value  :: pmesh_nodal
    integer(c_int), value  :: ownership

  end subroutine PDM_multipart_get_part_mesh_nodal_c

  !>
  !!
  !! \brief Set number connectivity for current block
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   connectivity_type     Type of connectivity
  !! \param [in]   connect               Connectivity (size = connect_idx[dn_entity])
  !! \param [in]   connect_idx           Index of connectivity or NULL if face_cell for example  (size = dn_entity)
  !!

  subroutine PDM_multipart_dconnectivity_set_c (multipart, &
                                                i_zone, &
                                                connectivity_type, &
                                                dconnect, &
                                                dconnect_idx) &
  bind (c, name='PDM_multipart_dconnectivity_set')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: connectivity_type
    type(c_ptr),    value  :: dconnect
    type(c_ptr),    value  :: dconnect_idx

  end subroutine PDM_multipart_dconnectivity_set_c

  !>
  !!
  !! \brief Set group connectivity by kind
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   bound_type            Type of boundary
  !! \param [in]   connect               Connectivity (size = connect_idx[dn_entity])
  !! \param [in]   connect_idx           Index of connectivity or NULL if face_cell for example  (size = dn_entity)
  !!

  subroutine PDM_multipart_dgroup_set_c (multipart, &
                                         i_zone, &
                                         bound_type, &
                                         dconnect, &
                                         dconnect_idx) &
  bind (c, name='PDM_multipart_dgroup_set')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: bound_type
    type(c_ptr),    value  :: dconnect
    type(c_ptr),    value  :: dconnect_idx

  end subroutine PDM_multipart_dgroup_set_c

  !>
  !!
  !! \brief Set coordinates
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   dvtx_coord            Mesh coordinates (size = 3 * dn_vtx)
  !!

  subroutine PDM_multipart_dvtx_coord_set_c (multipart, &
                                             i_zone, &
                                             dvtx_coord) &
  bind (c, name='PDM_multipart_dvtx_coord_set')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_zone
    type(c_ptr),    value  :: dvtx_coord

  end subroutine PDM_multipart_dvtx_coord_set_c

  !>
  !!
  !! \brief Returns the dimensions of a given partition
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                Partition index
  !! \param [out]  n_elt                 Number of elements
  !! \param [out]  n_cell                Number of cells
  !! \param [out]  n_face                Number of faces
  !! \param [out]  n_face_part_bound     Number of boundary faces in partition
  !! \param [out]  n_vtx                 Number of vertices
  !! \param [out]  n_proc                Number of processes
  !! \param [out]  n_total_part          Total number of partitions
  !! \param [out]  s_cell_face           Size of cell->face connectivity
  !! \param [out]  s_face_vtx            Size of face->vtx connectivity
  !! \param [out]  s_face_bound          Size of face->boundary connectivity
  !! \param [out]  n_bound_groups        Number of boundary groups
  !! \param [out]  s_face_join           ??
  !! \param [out]  n_join_groups         ??
  !!
  !!

  subroutine PDM_multipart_part_dim_get_c (multipart, &
                                           i_zone, &
                                           i_part, &
                                           n_section, &
                                           n_elt, &
                                           n_cell, &
                                           n_face, &
                                           n_face_part_bound, &
                                           n_vtx, &
                                           n_proc, &
                                           n_total_part, &
                                           s_cell_face, &
                                           s_face_vtx, &
                                           s_face_bound, &
                                           n_bound_groups, &
                                           s_face_join, &
                                           n_join_groups) &
  bind (c, name='PDM_multipart_part_dim_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: i_part
    integer(c_int)         :: n_section
    type(c_ptr)            :: n_elt
    integer(c_int)         :: n_cell
    integer(c_int)         :: n_face
    integer(c_int)         :: n_face_part_bound
    integer(c_int)         :: n_vtx
    integer(c_int)         :: n_proc
    integer(c_int)         :: n_total_part
    integer(c_int)         :: s_cell_face
    integer(c_int)         :: s_face_vtx
    integer(c_int)         :: s_face_bound
    integer(c_int)         :: n_bound_groups
    integer(c_int)         :: s_face_join
    integer(c_int)         :: n_join_groups

  end subroutine PDM_multipart_part_dim_get_c

  !>
  !!
  !! \brief Returns the data arrays of a given partition
  !!
  !! \param [in]   multipart                Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                   Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                   Partition index
  !! \param [out]  elt_vtx_idx              Index of element->vertex connecvitiy
  !! \param [out]  elt_vtx                  Element->vertex connecvitiy
  !! \param [out]  elt_section_ln_to_gn     Element local number to global number
  !! \param [out]  cell_tag                 Cell->tag connecvitiy
  !! \param [out]  cell_face_idx            Index of cell->face connectivity
  !! \param [out]  cell_face                Cell->face connectivity
  !! \param [out]  cell_ln_to_gn            Cell local number to global number
  !! \param [out]  face_tag                 Face->tag connecvitiy
  !! \param [out]  face_cell                Face->cell connectivity
  !! \param [out]  face_vtx_idx             Index of face->vertex connectivity
  !! \param [out]  face_vtx                 Face->vertex connectivity
  !! \param [out]  face_ln_to_gn            Face local number to global number
  !! \param [out]  face_part_bound_proc_idx ??
  !! \param [out]  face_part_bound_part_idx ??
  !! \param [out]  face_part_bound          ??
  !! \param [out]  vtx_tag                  Vertex->tag connectivity
  !! \param [out]  vtx                      Vertices coordinates
  !! \param [out]  vtx_ln_to_gn             Vertex local number to global number
  !! \param [out]  face_bound_idx           Index of face->boundary connectivity
  !! \param [out]  face_bound               Boundary faces
  !! \param [out]  face_bound_ln_to_gn      Boundary face local number to global number
  !! \param [out]  face_join_idx            Index of face->join connectivity
  !! \param [out]  face_join                Face->join connectivity
  !! \param [out]  face_join_ln_to_gn       Join face local number to global number
  !!

  subroutine PDM_multipart_part_val_get_c (multipart, &
                                           i_zone, &
                                           i_part, &
                                           elt_vtx_idx, &
                                           elt_vtx, &
                                           elt_section_ln_to_gn, &
                                           cell_tag, &
                                           cell_face_idx, &
                                           cell_face, &
                                           cell_ln_to_gn, &
                                           face_tag, &
                                           face_cell, &
                                           face_vtx_idx, &
                                           face_vtx, &
                                           face_ln_to_gn, &
                                           face_part_bound_proc_idx, &
                                           face_part_bound_part_idx, &
                                           face_part_bound, &
                                           vtx_tag, &
                                           vtx, &
                                           vtx_ln_to_gn, &
                                           face_bound_idx, &
                                           face_bound, &
                                           face_bound_ln_to_gn, &
                                           face_join_idx, &
                                           face_join, &
                                           face_join_ln_to_gn) &
  bind (c, name='PDM_multipart_part_val_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: i_part
    type(c_ptr)            :: elt_vtx_idx
    type(c_ptr)            :: elt_vtx
    type(c_ptr)            :: elt_section_ln_to_gn
    type(c_ptr)            :: cell_tag
    type(c_ptr)            :: cell_face_idx
    type(c_ptr)            :: cell_face
    type(c_ptr)            :: cell_ln_to_gn
    type(c_ptr)            :: face_tag
    type(c_ptr)            :: face_cell
    type(c_ptr)            :: face_vtx_idx
    type(c_ptr)            :: face_vtx
    type(c_ptr)            :: face_ln_to_gn
    type(c_ptr)            :: face_part_bound_proc_idx
    type(c_ptr)            :: face_part_bound_part_idx
    type(c_ptr)            :: face_part_bound
    type(c_ptr)            :: vtx_tag
    type(c_ptr)            :: vtx
    type(c_ptr)            :: vtx_ln_to_gn
    type(c_ptr)            :: face_bound_idx
    type(c_ptr)            :: face_bound
    type(c_ptr)            :: face_bound_ln_to_gn
    type(c_ptr)            :: face_join_idx
    type(c_ptr)            :: face_join
    type(c_ptr)            :: face_join_ln_to_gn

  end subroutine PDM_multipart_part_val_get_c

  !>
  !!
  !! \brief Returns the dimensions of a given partition
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                Partition index
  !! \param [out]  connectivity_type     Type of connectivity
  !! \param [out]  connect               Connectivity
  !! \param [out]  connect_idx           Connectivity index
  !! \param [out]  ownership             Data ownership
  !!

  function PDM_multipart_part_connectivity_get_c (multipart, &
                                                  i_zone, &
                                                  i_part, &
                                                  connectivity_type, &
                                                  connect, &
                                                  connect_idx, &
                                                  ownership) &
  result (pn_entity) &
  bind (c, name='PDM_multipart_part_connectivity_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value :: multipart
    integer(c_int), value :: i_zone
    integer(c_int), value :: i_part
    integer(c_int), value :: connectivity_type
    type(c_ptr)           :: connect
    type(c_ptr)           :: connect_idx
    integer(c_int), value :: ownership

    integer(c_int)        :: pn_entity

  end function PDM_multipart_part_connectivity_get_c

  !>
  !!
  !! \brief Returns the dimensions of a given partition
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                Partition index
  !! \param [out]  entity_type           Type of entity
  !! \param [out]  entity_ln_to_gn       Entity local number to global number
  !! \param [out]  ownership             Data ownership
  !!

  function PDM_multipart_part_ln_to_gn_get_c (multipart, &
                                              i_zone, &
                                              i_part, &
                                              entity_type, &
                                              entity_ln_to_gn, &
                                              ownership) &
  result (pn_entity) &
  bind (c, name='PDM_multipart_part_ln_to_gn_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: i_part
    integer(c_int), value  :: entity_type
    type(c_ptr)            :: entity_ln_to_gn
    integer(c_int), value  :: ownership

    integer(c_int)         :: pn_entity

  end function PDM_multipart_part_ln_to_gn_get_c

  !>
  !!
  !! \brief Returns the partitions color
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                Partition index
  !! \param [out]  entity_type           Type of entity
  !! \param [out]  entity_color          Entity color
  !! \param [out]  ownership             Data ownership
  !!

  function PDM_multipart_partition_color_get_c (multipart, &
                                                i_zone, &
                                                i_part, &
                                                entity_type, &
                                                entity_color, &
                                                ownership) &
  result (pn_entity) &
  bind (c, name='PDM_multipart_partition_color_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: i_part
    integer(c_int), value  :: entity_type
    type(c_ptr)            :: entity_color
    integer(c_int), value  :: ownership

    integer(c_int)         :: pn_entity

  end function PDM_multipart_partition_color_get_c

  !>
  !!
  !! \brief Get cache blocking data ??
  !!
  !! \param [in]   multipart               Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                  Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                  Partition index
  !! \param [out]  cell_color              Cell->color connectivity
  !! \param [out]  face_color              Face->color connectivity
  !! \param [out]  face_hp_color           ??
  !! \param [out]  thread_color            ??
  !! \param [out]  hyperplane_color        ??
  !!

  subroutine PDM_multipart_part_color_get_c (multipart, &
                                             i_zone, &
                                             i_part, &
                                             cell_color, &
                                             face_color, &
                                             face_hp_color, &
                                             thread_color, &
                                             hyperplane_color) &
  bind (c, name='PDM_multipart_part_color_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: i_part
    type(c_ptr)            :: cell_color
    type(c_ptr)            :: face_color
    type(c_ptr)            :: face_hp_color
    type(c_ptr)            :: thread_color
    type(c_ptr)            :: hyperplane_color

  end subroutine PDM_multipart_part_color_get_c

  !>
  !!
  !! \brief Get information of ghost vertices
  !!
  !! \param [in]   multipart               Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                  Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                  Partition index
  !! \param [out]  vtx_ghost_information   Information of ghost vertices
  !!

  subroutine PDM_multipart_part_ghost_infomation_get_c (multipart, &
                                                        i_zone, &
                                                        i_part, &
                                                        vtx_ghost_information) &
  bind (c, name='PDM_multipart_part_ghost_infomation_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: i_part
    type(c_ptr)            :: vtx_ghost_information

  end subroutine PDM_multipart_part_ghost_infomation_get_c

  !>
  !!
  !! \brief Get coordinates of vertices
  !!
  !! \param [in]   multipart               Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                  Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                  Partition index
  !! \param [out]  vtx_coord               Coordinates of vertices
  !! \param [out]  ownership               Data ownership
  !!
  !! \return Number of vertices
  !!

  function PDM_multipart_part_vtx_coord_get_c (multipart, &
                                               i_zone, &
                                               i_part, &
                                               vtx_coord, &
                                               ownership) &
  result (n_vtx) &
  bind (c, name='PDM_multipart_part_vtx_coord_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: i_part
    type(c_ptr)            :: vtx_coord
    integer(c_int)         :: ownership
    integer(c_int)         :: n_vtx

  end function PDM_multipart_part_vtx_coord_get_c

  !>
  !!
  !! \brief Get boundary data
  !!
  !! \param [in]   multipart               Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                  Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                  Partition index
  !! \param [out]  bound_type              Boundary type
  !! \param [out]  n_bound                 Number of boundaries
  !! \param [out]  bound_idx               Boundary index
  !! \param [out]  bound                   Boundaries
  !! \param [out]  bound_ln_to_gn          Boundary local number to global number
  !!

  subroutine PDM_multipart_bound_get_c (multipart, &
                                        i_zone, &
                                        i_part, &
                                        bound_type, &
                                        n_bound, &
                                        bound_idx, &
                                        bound, &
                                        bound_ln_to_gn) &
  bind (c, name='PDM_multipart_bound_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: i_part
    integer(c_int), value  :: bound_type
    integer(c_int)         :: n_bound
    type(c_ptr)            :: bound_idx
    type(c_ptr)            :: bound
    type(c_ptr)            :: bound_ln_to_gn

  end subroutine PDM_multipart_bound_get_c

  !>
  !!
  !! \brief Set distributed mesh data for the input zone
  !!
  !! \param [in]   multipart      Pointer to \ref PDM_multipart_t object
  !! \param [in]   zone_id        Global zone id
  !! \param [in]   dmesh          Distributed mesh structure
  !!

  subroutine PDM_multipart_register_block (multipart, &
                                           zone_id, &
                                           dmesh) &
  bind (c, name='PDM_multipart_register_block')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: zone_id
    type(c_ptr),    value  :: dmesh

  end subroutine PDM_multipart_register_block

  !>
  !!
  !! \brief Set distributed mesh data for the input zone
  !!
  !! \param [in]   multipart      Pointer to \ref PDM_multipart_t object
  !! \param [in]   zone_id        Global zone id
  !! \param [in]   dmesh_nodal    Distributed nodal mesh structure
  !!

  subroutine PDM_multipart_register_dmesh_nodal (multipart, &
                                                 zone_id, &
                                                 dmesh_nodal) &
  bind (c, name='PDM_multipart_register_dmesh_nodal')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: zone_id
    type(c_ptr),    value  :: dmesh_nodal

  end subroutine PDM_multipart_register_dmesh_nodal

  !>
  !!
  !! \brief Construct the partitioned meshes on every zones
  !!
  !! \param [in]   multipart   Pointer to \ref PDM_multipart_t object
  !!

  subroutine PDM_multipart_run_ppart (multipart) &
  bind (c, name='PDM_multipart_run_ppart')

    use iso_c_binding
    implicit none

    type(c_ptr), value  :: multipart

  end subroutine PDM_multipart_run_ppart

  !>
  !!
  !! \brief Set number of element in the block entity
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   entity_type           Type of entity (can be cell/face/edge/vtx)
  !! \param [in]   dn_entity             Distributed number of entity in current process
  !!

  subroutine PDM_multipart_dn_entity_set (multipart, &
                                          i_zone, &
                                          entity_type, &
                                          dn_entity) &
  bind (c, name='PDM_multipart_dn_entity_set')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: entity_type
    integer(c_int), value  :: dn_entity

  end subroutine PDM_multipart_dn_entity_set

  !>
  !!
  !! \brief Set the domain interface
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   ditrf                 Domain interface
  !!

  subroutine PDM_multipart_domain_interface_shared_set (multipart, &
                                                        ditrf) &
  bind (c, name='PDM_multipart_domain_interface_shared_set')

    use iso_c_binding
    implicit none

    type(c_ptr), value :: multipart
    type(c_ptr), value :: ditrf

  end subroutine PDM_multipart_domain_interface_shared_set

  !>
  !!
  !! \brief Free the structure
  !!
  !! \param [in]   multipart   Pointer to \ref PDM_multipart_t object
  !!

  subroutine PDM_multipart_free (multipart) &
  bind (c, name='PDM_multipart_free')

    use iso_c_binding
    implicit none

    type(c_ptr), value     :: multipart

  end subroutine PDM_multipart_free

end interface

private :: PDM_multipart_create_,&
           PDM_multipart_register_joins_,&
           PDM_multipart_set_reordering_options_,&
           PDM_multipart_set_reordering_options_vtx_,&
           PDM_multipart_get_part_mesh_nodal_,&
           PDM_multipart_dconnectivity_set_,&
           PDM_multipart_dvtx_coord_set_,&
           PDM_multipart_part_dim_get_,&
           PDM_multipart_part_graph_comm_vtx_dim_get_,&
           PDM_multipart_part_val_get_,&
           PDM_multipart_part_connectivity_get_,&
           PDM_multipart_part_ln_to_gn_get_,&
           PDM_multipart_partition_color_get_,&
           PDM_multipart_part_color_get_,&
           PDM_multipart_part_ghost_infomation_get_,&
           PDM_multipart_part_vtx_coord_get_,&
           PDM_multipart_bound_get_

contains

  !>
  !!
  !! \brief Build a multipart structure
  !!
  !! \param [out]  multipart        Pointer to a new \ref PDM_multipart_t object
  !! \param [in]   n_zone           Number of zones in the original mesh
  !! \param [in]   n_part           Number of partition per proc in each zone
  !! \param [in]   merge_blocks     Merge or not the zones before splitting
  !! \param [in]   split_method     Choice of library used to split the mesh
  !! \param [in]   part_size_method Choice of homogeneous or heterogeneous partitions
  !! \param [in]   part_fraction    Weight (in %) of each partition in heterogeneous case
  !! \param [in]   comm             PDM_MPI communicator
  !! \param [in]   owner            Data ownership
  !!

  subroutine PDM_multipart_create_ (multipart, &
                                    n_zone, &
                                    n_part, &
                                    merge_blocks, &
                                    split_method, &
                                    part_size_method, &
                                    part_fraction, &
                                    comm, &
                                    owner)

    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr)                        :: multipart
    integer(c_int),            value   :: n_zone
    integer(kind=PDM_l_num_s), pointer :: n_part(:)
    type(c_ptr)                        :: c_n_part         = C_NULL_PTR
    integer(c_int),            value   :: merge_blocks
    integer(c_int),            value   :: split_method
    integer(c_int),            value   :: part_size_method
    double precision,          pointer :: part_fraction(:)
    type(c_ptr)                        :: c_part_fraction  = C_NULL_PTR
    integer(c_int),            value   :: comm
    integer(c_int),            value   :: owner

    c_n_part        = c_loc(n_part)
    c_part_fraction = c_loc(part_fraction)

    multipart = PDM_multipart_create_c(n_zone, &
                                       c_n_part, &
                                       merge_blocks, &
                                       split_method, &
                                       part_size_method, &
                                       c_part_fraction, &
                                       comm, &
                                       owner)

  end subroutine PDM_multipart_create_

  !>
  !!
  !! \brief Set connecting data between all the zones
  !!
  !! \param [in]   multipart        Pointer to \ref PDM_multipart_t object
  !! \param [in]   n_total_joins    Total number of interfaces
  !! \param [in]   join_to_opposite For each global join id, give the global id
  !!                                of the opposite join (size = n_total_joins)
  !!

  subroutine PDM_multipart_register_joins_ (multipart, &
                                            n_total_joins, &
                                            join_to_opposite)

    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr),               value   :: multipart
    integer(c_int),            value   :: n_total_joins
    integer(kind=PDM_l_num_s), pointer :: join_to_opposite(:)
    type(c_ptr)                        :: c_join_to_opposite = C_NULL_PTR

    c_join_to_opposite = c_loc(join_to_opposite)

    call PDM_multipart_register_joins_c(multipart, &
                                        n_total_joins, &
                                        c_join_to_opposite)

  end subroutine PDM_multipart_register_joins_

  !>
  !!
  !! \brief Set the reordering methods to be used after partitioning
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   renum_cell_method     Choice of renumbering method for cells
  !! \param [in]   renum_cell_properties Parameters used by cacheblocking method :
  !!                                     [n_cell_per_cache_wanted, is_asynchrone, is_vectorisation,
  !!                                     n_vect_face, split_method]
  !! \param [in]   renum_face_method     Choice of renumbering method for faces
  !!

  subroutine PDM_multipart_set_reordering_options_ (multipart, &
                                                    i_zone, &
                                                    renum_cell_method, &
                                                    renum_cell_properties, &
                                                    renum_face_method)

    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr),               value   :: multipart
    integer(c_int),            value   :: i_zone
    character (len=*)                  :: renum_cell_method
    integer(kind=PDM_l_num_s), pointer :: renum_cell_properties(:)
    type(c_ptr)                        :: c_renum_cell_properties = C_NULL_PTR
    character (len=*)                  :: renum_face_method

    c_renum_cell_properties = c_loc(renum_cell_properties)

    call PDM_multipart_set_reordering_options_c(multipart, &
                                                i_zone, &
                                                trim(renum_cell_method)//C_NULL_CHAR, &
                                                c_renum_cell_properties, &
                                                trim(renum_face_method)//C_NULL_CHAR)

  end subroutine PDM_multipart_set_reordering_options_

  !>
  !!
  !! \brief Set the reordering methods to be used after partitioning
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   renum_vtx_method      Choice of renumbering method for vertices

  subroutine PDM_multipart_set_reordering_options_vtx_ (multipart, &
                                                        i_zone, &
                                                        renum_vtx_method)

    use iso_c_binding
    implicit none

    type(c_ptr),      value  :: multipart
    integer(c_int),   value  :: i_zone
    character (len=*)        :: renum_vtx_method

    call PDM_multipart_set_reordering_options_vtx_c(multipart, &
                                                    i_zone, &
                                                    trim(renum_vtx_method)//C_NULL_CHAR)

  end subroutine PDM_multipart_set_reordering_options_vtx_

  !>
  !!
  !! \brief ???
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   pmesh_nodal           Partitionned nodal mesh
  !! \param [in]   ownership             Data ownership
  !!

  subroutine PDM_multipart_get_part_mesh_nodal_ (multipart, &
                                                 i_zone, &
                                                 pmesh_nodal, &
                                                 ownership)

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_zone
    type(c_ptr),    value  :: pmesh_nodal
    integer(c_int), value  :: ownership

    call PDM_multipart_get_part_mesh_nodal_c(multipart, &
                                             i_zone, &
                                             pmesh_nodal, &
                                             ownership)

  end subroutine PDM_multipart_get_part_mesh_nodal_

  !>
  !!
  !! \brief Set number connectivity for current block
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   connectivity_type     Type of connectivity
  !! \param [in]   connect               Connectivity (size = connect_idx[dn_entity])
  !! \param [in]   connect_idx           Index of connectivity or NULL if face_cell for example  (size = dn_entity)
  !!

  subroutine PDM_multipart_dconnectivity_set_ (multipart, &
                                               i_zone, &
                                               connectivity_type, &
                                               dconnect, &
                                               dconnect_idx)

    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr),               value   :: multipart
    integer(c_int),            value   :: i_zone
    integer(c_int),            value   :: connectivity_type
    integer(kind=PDM_l_num_s), pointer :: dconnect(:)
    type(c_ptr)                        :: c_dconnect      = C_NULL_PTR
    integer(kind=PDM_g_num_s), pointer :: dconnect_idx(:)
    type(c_ptr)                        :: c_dconnect_idx  = C_NULL_PTR

    c_dconnect_idx = c_loc(dconnect_idx)

    if (associated(dconnect)) then
      c_dconnect = c_loc(dconnect)
    endif

    call PDM_multipart_dconnectivity_set_c(multipart, &
                                           i_zone, &
                                           connectivity_type, &
                                           c_dconnect, &
                                           c_dconnect_idx)

  end subroutine PDM_multipart_dconnectivity_set_

  !>
  !!
  !! \brief Set group connectivity by kind
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   bound_type            Type of bound
  !! \param [in]   connect               Connectivity (size = connect_idx[dn_entity])
  !! \param [in]   connect_idx           Index of connectivity or NULL if face_cell for example  (size = dn_entity)
  !!

  subroutine PDM_multipart_dgroup_set_ (multipart, &
                                        i_zone, &
                                        bound_type, &
                                        dconnect, &
                                        dconnect_idx)

    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr),               value   :: multipart
    integer(c_int),            value   :: i_zone
    integer(c_int),            value   :: bound_type
    integer(kind=PDM_l_num_s), pointer :: dconnect(:)
    type(c_ptr)                        :: c_dconnect     = C_NULL_PTR
    integer(kind=PDM_g_num_s), pointer :: dconnect_idx(:)
    type(c_ptr)                        :: c_dconnect_idx = C_NULL_PTR

    c_dconnect_idx = c_loc(dconnect_idx)

    if (associated(dconnect)) then
      c_dconnect = c_loc(dconnect)
    endif

    call PDM_multipart_dgroup_set_c(multipart, &
                                    i_zone, &
                                    bound_type, &
                                    c_dconnect, &
                                    c_dconnect_idx)

  end subroutine PDM_multipart_dgroup_set_

  !>
  !!
  !! \brief Set coordinates
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   dvtx_coord            Mesh coordinates (size = 3 * dn_vtx)
  !!

  subroutine PDM_multipart_dvtx_coord_set_ (multipart, &
                                            i_zone, &
                                            dvtx_coord)

    use iso_c_binding
    implicit none

    type(c_ptr),      value   :: multipart
    integer(c_int),   value   :: i_zone
    double precision, pointer :: dvtx_coord(:)
    type(c_ptr)               :: c_dvtx_coord = C_NULL_PTR

    if (associated(dvtx_coord)) then
      c_dvtx_coord = c_loc(dvtx_coord)
    endif

    call PDM_multipart_dvtx_coord_set_c(multipart, &
                                        i_zone, &
                                        c_dvtx_coord)

  end subroutine PDM_multipart_dvtx_coord_set_

  !>
  !!
  !! \brief Returns the dimensions of a given partition
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                Partition index
  !! \param [out]  n_elt                 Number of elements
  !! \param [out]  n_cell                Number of cells
  !! \param [out]  n_face                Number of faces
  !! \param [out]  n_face_part_bound     Number of boundary faces in partition
  !! \param [out]  n_vtx                 Number of vertices
  !! \param [out]  n_proc                Number of processes
  !! \param [out]  n_total_part          Total number of partitions
  !! \param [out]  s_cell_face           Size of cell->face connectivity
  !! \param [out]  s_face_vtx            Size of face->vtx connectivity
  !! \param [out]  s_face_bound          Size of boundary faces
  !! \param [out]  n_bound_groups        Number of boundary groups
  !! \param [out]  s_face_join           ??
  !! \param [out]  n_join_groups         ??
  !!
  !!

  subroutine PDM_multipart_part_dim_get_ (multipart, &
                                          i_zone, &
                                          i_part, &
                                          n_section, &
                                          n_elt, &
                                          n_cell, &
                                          n_face, &
                                          n_face_part_bound, &
                                          n_vtx, &
                                          n_proc, &
                                          n_total_part, &
                                          s_cell_face, &
                                          s_face_vtx, &
                                          s_face_bound, &
                                          n_bound_groups, &
                                          s_face_join, &
                                          n_join_groups)

    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr), value       :: multipart
    integer,     intent(in)  :: i_zone
    integer,     intent(in)  :: i_part
    integer,     intent(out) :: n_section
    integer,     intent(out) :: n_cell
    integer,     intent(out) :: n_face
    integer,     intent(out) :: n_face_part_bound
    integer,     intent(out) :: n_vtx
    integer,     intent(out) :: n_proc
    integer,     intent(out) :: n_total_part
    integer,     intent(out) :: s_cell_face
    integer,     intent(out) :: s_face_vtx
    integer,     intent(out) :: s_face_bound
    integer,     intent(out) :: n_bound_groups
    integer,     intent(out) :: s_face_join
    integer,     intent(out) :: n_join_groups

    integer(c_int) :: c_n_section
    integer(c_int) :: c_n_cell
    integer(c_int) :: c_n_face
    integer(c_int) :: c_n_face_part_bound
    integer(c_int) :: c_n_vtx
    integer(c_int) :: c_n_proc
    integer(c_int) :: c_n_total_part
    integer(c_int) :: c_s_cell_face
    integer(c_int) :: c_s_face_vtx
    integer(c_int) :: c_s_face_bound
    integer(c_int) :: c_n_bound_groups
    integer(c_int) :: c_s_face_join
    integer(c_int) :: c_n_join_groups

    integer (kind = PDM_l_num_s), pointer :: n_elt(:)
    type(c_ptr)                           :: c_n_elt = C_NULL_PTR

    c_n_elt = c_loc(n_elt)

    call PDM_multipart_part_dim_get_c(multipart, &
                                      i_zone, &
                                      i_part, &
                                      c_n_section, &
                                      c_n_elt, &
                                      c_n_cell, &
                                      c_n_face, &
                                      c_n_face_part_bound, &
                                      c_n_vtx, &
                                      c_n_proc, &
                                      c_n_total_part, &
                                      c_s_cell_face, &
                                      c_s_face_vtx, &
                                      c_s_face_bound, &
                                      c_n_bound_groups, &
                                      c_s_face_join, &
                                      c_n_join_groups)

    n_section         = c_n_section
    n_cell            = c_n_cell
    n_face            = c_n_face
    n_face_part_bound = c_n_face_part_bound
    n_vtx             = c_n_vtx
    n_proc            = c_n_proc
    n_total_part      = c_n_total_part
    s_cell_face       = c_s_cell_face
    s_face_vtx        = c_s_face_vtx
    s_face_bound      = c_s_face_bound
    n_bound_groups    = c_n_bound_groups
    s_face_join       = c_s_face_join
    n_join_groups     = c_n_join_groups

    call c_f_pointer(c_n_elt, &
                     n_elt,   &
                     [n_section])

  end subroutine PDM_multipart_part_dim_get_

  !>
  !!
  !! \brief Returns the data arrays of a given partition
  !!
  !! \param [in]   multipart                Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                   Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                   Partition index
  !! \param [out]  elt_vtx_idx              ??
  !! \param [out]  elt_vtx                  ??
  !! \param [out]  elt_section_ln_to_gn     ??
  !! \param [out]  cell_tag                 ??
  !! \param [out]  cell_face_idx            ??
  !! \param [out]  cell_face                ??
  !! \param [out]  cell_ln_to_gn            ??
  !! \param [out]  face_tag                 ??
  !! \param [out]  face_cell                ??
  !! \param [out]  face_vtx_idx             ??
  !! \param [out]  face_vtx                 ??
  !! \param [out]  face_ln_to_gn            ??
  !! \param [out]  face_part_bound_proc_idx ??
  !! \param [out]  face_part_bound_part_idx ??
  !! \param [out]  face_part_bound          ??
  !! \param [out]  vtx_tag                  ??
  !! \param [out]  vtx                      ??
  !! \param [out]  vtx_ln_to_gn             ??
  !! \param [out]  face_bound_idx           ??
  !! \param [out]  face_bound               ??
  !! \param [out]  face_bound_ln_to_gn      ??
  !! \param [out]  face_join_idx            ??
  !! \param [out]  face_join                ??
  !! \param [out]  face_join_ln_to_gn       ??
  !!

  subroutine PDM_multipart_part_val_get_(multipart, &
                                         i_zone, &
                                         i_part, &
                                         elt_vtx_idx, &
                                         elt_vtx, &
                                         elt_section_ln_to_gn, &
                                         cell_tag, &
                                         cell_face_idx, &
                                         cell_face, &
                                         cell_ln_to_gn, &
                                         face_tag, &
                                         face_cell, &
                                         face_vtx_idx, &
                                         face_vtx, &
                                         face_ln_to_gn, &
                                         face_part_bound_proc_idx, &
                                         face_part_bound_part_idx, &
                                         face_part_bound, &
                                         vtx_tag, &
                                         vtx, &
                                         vtx_ln_to_gn, &
                                         face_bound_idx, &
                                         face_bound, &
                                         face_bound_ln_to_gn, &
                                         face_join_idx, &
                                         face_join, &
                                         face_join_ln_to_gn)

    use pdm
    use iso_c_binding
    use pdm_pointer_array
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: i_part

    type(PDM_pointer_array_t), target :: elt_vtx_idx
    type(PDM_pointer_array_t), target :: elt_vtx
    type(PDM_pointer_array_t), target :: elt_section_ln_to_gn

    type(c_ptr) :: c_elt_vtx_idx          = C_NULL_PTR
    type(c_ptr) :: c_elt_vtx              = C_NULL_PTR
    type(c_ptr) :: c_elt_section_ln_to_gn = C_NULL_PTR

    integer (kind = PDM_l_num_s), pointer :: cell_tag(:)
    integer (kind = PDM_l_num_s), pointer :: cell_face_idx(:)
    integer (kind = PDM_l_num_s), pointer :: cell_face(:)
    integer (kind = PDM_g_num_s), pointer :: cell_ln_to_gn(:)
    integer (kind = PDM_l_num_s), pointer :: face_tag(:)
    integer (kind = PDM_l_num_s), pointer :: face_cell(:)
    integer (kind = PDM_l_num_s), pointer :: face_vtx_idx(:)
    integer (kind = PDM_l_num_s), pointer :: face_vtx(:)
    integer (kind = PDM_g_num_s), pointer :: face_ln_to_gn(:)
    integer (kind = PDM_l_num_s), pointer :: face_part_bound_proc_idx(:)
    integer (kind = PDM_l_num_s), pointer :: face_part_bound_part_idx(:)
    integer (kind = PDM_l_num_s), pointer :: face_part_bound(:)
    integer (kind = PDM_l_num_s), pointer :: vtx_tag(:)
    integer (kind = PDM_l_num_s), pointer :: vtx(:)
    integer (kind = PDM_g_num_s), pointer :: vtx_ln_to_gn(:)
    integer (kind = PDM_l_num_s), pointer :: face_bound_idx(:)
    integer (kind = PDM_l_num_s), pointer :: face_bound(:)
    integer (kind = PDM_g_num_s), pointer :: face_bound_ln_to_gn(:)
    integer (kind = PDM_l_num_s), pointer :: face_join_idx(:)
    integer (kind = PDM_l_num_s), pointer :: face_join(:)
    integer (kind = PDM_g_num_s), pointer :: face_join_ln_to_gn(:)

    type(c_ptr)  :: c_cell_tag
    type(c_ptr)  :: c_cell_face_idx
    type(c_ptr)  :: c_cell_face
    type(c_ptr)  :: c_cell_ln_to_gn
    type(c_ptr)  :: c_face_tag
    type(c_ptr)  :: c_face_cell
    type(c_ptr)  :: c_face_vtx_idx
    type(c_ptr)  :: c_face_vtx
    type(c_ptr)  :: c_face_ln_to_gn
    type(c_ptr)  :: c_face_part_bound_proc_idx
    type(c_ptr)  :: c_face_part_bound_part_idx
    type(c_ptr)  :: c_face_part_bound
    type(c_ptr)  :: c_vtx_tag
    type(c_ptr)  :: c_vtx
    type(c_ptr)  :: c_vtx_ln_to_gn
    type(c_ptr)  :: c_face_bound_idx
    type(c_ptr)  :: c_face_bound
    type(c_ptr)  :: c_face_bound_ln_to_gn
    type(c_ptr)  :: c_face_join_idx
    type(c_ptr)  :: c_face_join
    type(c_ptr)  :: c_face_join_ln_to_gn

    integer(c_int) :: c_n_section
    integer(c_int) :: c_n_cell
    integer(c_int) :: c_n_face
    integer(c_int) :: c_n_face_part_bound
    integer(c_int) :: c_n_vtx
    integer(c_int) :: c_n_proc
    integer(c_int) :: c_n_total_part
    integer(c_int) :: c_s_cell_face
    integer(c_int) :: c_s_face_vtx
    integer(c_int) :: c_s_face_bound
    integer(c_int) :: c_n_bound_groups
    integer(c_int) :: c_s_face_join
    integer(c_int) :: c_n_join_groups

    integer :: s_face_vtx
    integer :: n_section
    integer :: n_cell
    integer :: n_face
    integer :: n_face_part_bound
    integer :: n_vtx
    integer :: n_bound_groups
    integer :: n_join_groups

    integer :: taille, i

    integer (kind = PDM_l_num_s), pointer :: n_elt(:) => null()
    type(c_ptr)                           :: c_n_elt

    c_n_elt = c_loc(n_elt)

    call PDM_multipart_part_dim_get_c(multipart, &
                                      i_zone, &
                                      i_part, &
                                      c_n_section, &
                                      c_n_elt, &
                                      c_n_cell, &
                                      c_n_face, &
                                      c_n_face_part_bound, &
                                      c_n_vtx, &
                                      c_n_proc, &
                                      c_n_total_part, &
                                      c_s_cell_face, &
                                      c_s_face_vtx, &
                                      c_s_face_bound, &
                                      c_n_bound_groups, &
                                      c_s_face_join, &
                                      c_n_join_groups)

    call c_f_pointer(c_n_elt, &
                     n_elt,   &
                     [c_n_section])

    s_face_vtx         = c_s_face_vtx
    n_section          = c_n_section
    n_cell             = c_n_cell
    n_face             = c_n_face
    n_face_part_bound  = c_n_face_part_bound
    n_vtx              = c_n_vtx
    n_bound_groups     = c_n_bound_groups
    n_join_groups      = c_n_join_groups

    c_cell_tag                 = c_loc(cell_tag)
    c_cell_face_idx            = c_loc(cell_face_idx)
    c_cell_face                = c_loc(cell_face)
    c_cell_ln_to_gn            = c_loc(cell_ln_to_gn)
    c_face_tag                 = c_loc(face_tag)
    c_face_cell                = c_loc(face_cell)
    c_face_vtx_idx             = c_loc(face_vtx_idx)
    c_face_vtx                 = c_loc(face_vtx)
    c_face_ln_to_gn            = c_loc(face_ln_to_gn)
    c_face_part_bound_proc_idx = c_loc(face_part_bound_proc_idx)
    c_face_part_bound_part_idx = c_loc(face_part_bound_part_idx)
    c_face_part_bound          = c_loc(face_part_bound)
    c_vtx_tag                  = c_loc(vtx_tag)
    c_vtx                      = c_loc(vtx)
    c_vtx_ln_to_gn             = c_loc(vtx_ln_to_gn)
    c_face_bound_idx           = c_loc(face_bound_idx)
    c_face_bound               = c_loc(face_bound)
    c_face_bound_ln_to_gn      = c_loc(face_bound_ln_to_gn)
    c_face_join_idx            = c_loc(face_join_idx)
    c_face_join                = c_loc(face_join)
    c_face_join_ln_to_gn       = c_loc(face_join_ln_to_gn)

    call PDM_multipart_part_val_get_c(multipart, &
                                      i_zone, &
                                      i_part, &
                                      c_elt_vtx_idx, &
                                      c_elt_vtx, &
                                      c_elt_section_ln_to_gn, &
                                      c_cell_tag, &
                                      c_cell_face_idx, &
                                      c_cell_face, &
                                      c_cell_ln_to_gn, &
                                      c_face_tag, &
                                      c_face_cell, &
                                      c_face_vtx_idx, &
                                      c_face_vtx, &
                                      c_face_ln_to_gn, &
                                      c_face_part_bound_proc_idx, &
                                      c_face_part_bound_part_idx, &
                                      c_face_part_bound, &
                                      c_vtx_tag, &
                                      c_vtx, &
                                      c_vtx_ln_to_gn, &
                                      c_face_bound_idx, &
                                      c_face_bound, &
                                      c_face_bound_ln_to_gn, &
                                      c_face_join_idx, &
                                      c_face_join, &
                                      c_face_join_ln_to_gn)

    call c_f_pointer(c_elt_vtx_idx, &
                     elt_vtx_idx%cptr, &
                     [n_section])

    call  PDM_pointer_array_create_type (elt_vtx_idx, &
                                         n_section, &
                                         PDM_TYPE_INT)

    do i = 1, n_section
      elt_vtx_idx%length(i) = n_elt(i) + 1
    end do

    call c_f_pointer(c_elt_vtx,    &
                     elt_vtx%cptr, &
                     [n_section])

    call  PDM_pointer_array_create_type (elt_vtx, &
                                         n_section, &
                                         PDM_TYPE_INT)

    ! do i = 1, n_section
    !   elt_vtx%length(i) = elt_vtx_idx(i)(n_elt(i))
    ! end do

    call c_f_pointer(c_elt_section_ln_to_gn,    &
                     elt_section_ln_to_gn%cptr, &
                     [n_section])

    call  PDM_pointer_array_create_type (elt_section_ln_to_gn, &
                                         n_section, &
                                         PDM_TYPE_G_NUM)

    do i = 1, n_section
      elt_section_ln_to_gn%length(i) = n_elt(i)
    end do

    call c_f_pointer(c_cell_tag, &
                     cell_tag,   &
                     [n_cell])

    call c_f_pointer(c_cell_face_idx, &
                     cell_face_idx,   &
                     [n_cell + 1])

    call c_f_pointer(c_cell_face, &
                     cell_face,   &
                     [cell_face_idx(n_cell)])

    call c_f_pointer(c_cell_ln_to_gn, &
                     cell_ln_to_gn,   &
                     [n_cell])

    call c_f_pointer(c_face_tag, &
                     face_tag,   &
                     [n_face])

    call c_f_pointer(c_face_cell, &
                     face_cell,   &
                     [2 * n_face])

    call c_f_pointer(c_face_vtx_idx, &
                     face_vtx_idx,   &
                     [n_face + 1])

    taille = face_vtx_idx(n_face)

    call c_f_pointer(c_face_vtx, &
                     face_vtx,   &
                     [taille])

    call c_f_pointer(c_face_ln_to_gn, &
                     face_ln_to_gn,   &
                     [n_face])

    call c_f_pointer(c_face_part_bound_proc_idx, &
                     face_part_bound_proc_idx,   &
                     [n_face_part_bound + 1])

    call c_f_pointer(c_face_part_bound_part_idx, &
                     face_part_bound_part_idx,   &
                     [n_face_part_bound + 1])

    taille = face_part_bound_part_idx(i_part + 1) - face_part_bound_part_idx(i_part)

    call c_f_pointer(c_face_part_bound, &
                     face_part_bound,   &
                     [taille])

    call c_f_pointer(c_vtx_tag, &
                     vtx_tag,   &
                     [n_vtx])

    taille = s_face_vtx * n_vtx

    call c_f_pointer(c_vtx, &
                     vtx,   &
                     [s_face_vtx * n_vtx])

    call c_f_pointer(c_vtx_ln_to_gn, &
                     vtx_ln_to_gn,   &
                     [n_vtx])

    call c_f_pointer(c_face_bound_idx, &
                     face_bound_idx,   &
                     [n_bound_groups])

    taille = face_bound_idx(n_bound_groups)

    call c_f_pointer(c_face_bound, &
                     face_bound,   &
                     [taille])

    call c_f_pointer(c_face_bound_ln_to_gn, &
                     face_bound_ln_to_gn,   &
                     [n_bound_groups])

    call c_f_pointer(c_face_join_idx, &
                     face_join_idx,   &
                     [n_join_groups + 1])

    taille = face_join(n_join_groups)

    call c_f_pointer(c_face_join, &
                     face_join,   &
                     [taille])

    call c_f_pointer(c_face_join_ln_to_gn, &
                     face_join_ln_to_gn,   &
                     [n_join_groups])

  end subroutine PDM_multipart_part_val_get_

  !>
  !!
  !! \brief Returns the dimensions of a given partition
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                Partition index
  !! \param [out]  connectivity_type     Connectivity type
  !! \param [out]  connect               Connectivity
  !! \param [out]  connect_idx           Connectivity index
  !! \param [out]  ownership             Data ownership
  !! \param [out]  pn_entity             Number of entities
  !!

  subroutine PDM_multipart_part_connectivity_get_ (multipart, &
                                                   i_zone, &
                                                   i_part, &
                                                   connectivity_type, &
                                                   connect, &
                                                   connect_idx, &
                                                   ownership, &
                                                   pn_entity)

    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr),               value   :: multipart
    integer(c_int),            value   :: i_zone
    integer(c_int),            value   :: i_part
    integer(c_int),            value   :: connectivity_type
    integer(kind=PDM_l_num_s), pointer :: connect(:)
    type(c_ptr)                        :: c_connect = C_NULL_PTR
    integer(kind=PDM_l_num_s), pointer :: connect_idx(:)
    type(c_ptr)                        :: c_connect_idx = C_NULL_PTR
    integer(c_int),            value   :: ownership
    integer(c_int),            value   :: pn_entity

    integer :: taille

    c_connect     = c_loc(connect)
    c_connect_idx = c_loc(connect_idx)

    pn_entity = PDM_multipart_part_connectivity_get_c(multipart, &
                                                      i_zone, &
                                                      i_part, &
                                                      connectivity_type, &
                                                      c_connect, &
                                                      c_connect_idx, &
                                                      ownership)

    call c_f_pointer(c_connect_idx, &
                     connect_idx,   &
                     [pn_entity])

    taille = connect_idx(pn_entity)

    call c_f_pointer(c_connect, &
                     connect,   &
                     [taille])

  end subroutine PDM_multipart_part_connectivity_get_

  !>
  !!
  !! \brief Returns the dimensions of a given partition
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                Partition index
  !! \param [out]  entity_type           Entity type
  !! \param [out]  entity_ln_to_gn       Entity local number to global number
  !! \param [out]  ownership             Data ownership
  !! \param [out]  pn_entity             Number of entities
  !!

  subroutine PDM_multipart_part_ln_to_gn_get_ (multipart, &
                                               i_zone, &
                                               i_part, &
                                               entity_type, &
                                               entity_ln_to_gn, &
                                               ownership, &
                                               pn_entity)

    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr),               value   :: multipart
    integer(c_int),            value   :: i_zone
    integer(c_int),            value   :: i_part
    integer(c_int),            value   :: entity_type
    integer(kind=PDM_g_num_s), pointer :: entity_ln_to_gn(:)
    type(c_ptr)                        :: c_entity_ln_to_gn = C_NULL_PTR
    integer(c_int),            value   :: ownership
    integer(c_int),            value   :: pn_entity

    c_entity_ln_to_gn = c_loc(entity_ln_to_gn)

    pn_entity = PDM_multipart_part_ln_to_gn_get_c(multipart, &
                                                  i_zone, &
                                                  i_part, &
                                                  entity_type, &
                                                  c_entity_ln_to_gn, &
                                                  ownership)

    call c_f_pointer(c_entity_ln_to_gn, &
                     entity_ln_to_gn,   &
                     [pn_entity])

  end subroutine PDM_multipart_part_ln_to_gn_get_

  !>
  !!
  !! \brief Returns the partitions color
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                Partition index
  !! \param [out]  entity_type           Entity type
  !! \param [out]  entity_color          Entity->color connectivity
  !! \param [out]  ownership             Data ownership
  !! \param [out]  pn_entity             Number of entities
  !!

  subroutine PDM_multipart_partition_color_get_(multipart, &
                                                i_zone, &
                                                i_part, &
                                                entity_type, &
                                                entity_color, &
                                                ownership, &
                                                pn_entity)

    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr),               value   :: multipart
    integer(c_int),            value   :: i_zone
    integer(c_int),            value   :: i_part
    integer(c_int),            value   :: entity_type
    integer(kind=PDM_l_num_s), pointer :: entity_color(:)
    type(c_ptr)                        :: c_entity_color = C_NULL_PTR
    integer(c_int),            value   :: ownership
    integer(c_int),            value   :: pn_entity

    c_entity_color = c_loc(entity_color)

    pn_entity = PDM_multipart_partition_color_get_c(multipart, &
                                                    i_zone, &
                                                    i_part, &
                                                    entity_type, &
                                                    c_entity_color, &
                                                    ownership)

    call c_f_pointer(c_entity_color, &
                     entity_color,   &
                     [pn_entity])

  end subroutine PDM_multipart_partition_color_get_

  !>
  !!
  !! \brief Get cache blocking data ??
  !!
  !! \param [in]   multipart               Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                  Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                  Partition index
  !! \param [out]  cell_color              Cell->color connectivity
  !! \param [out]  face_color              Face->color connectivity
  !! \param [out]  face_hp_color           ??
  !! \param [out]  thread_color            ??
  !! \param [out]  hyperplane_color        ??
  !!

  subroutine PDM_multipart_part_color_get_(multipart, &
                                           i_zone, &
                                           i_part, &
                                           cell_color, &
                                           face_color, &
                                           face_hp_color, &
                                           thread_color, &
                                           hyperplane_color)

    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr),               value   :: multipart
    integer(c_int),            value   :: i_zone
    integer(c_int),            value   :: i_part
    integer(kind=PDM_l_num_s), pointer :: cell_color(:)
    type(c_ptr)                        :: c_cell_color = C_NULL_PTR
    integer(kind=PDM_l_num_s), pointer :: face_color(:)
    type(c_ptr)                        :: c_face_color = C_NULL_PTR
    integer(kind=PDM_l_num_s), pointer :: face_hp_color(:)
    type(c_ptr)                        :: c_face_hp_color = C_NULL_PTR
    integer(kind=PDM_l_num_s), pointer :: thread_color(:)
    type(c_ptr)                        :: c_thread_color = C_NULL_PTR
    integer(kind=PDM_l_num_s), pointer :: hyperplane_color(:)
    type(c_ptr)                        :: c_hyperplane_color = C_NULL_PTR

    integer(c_int) :: c_n_section
    integer(c_int) :: c_n_cell
    integer(c_int) :: c_n_face
    integer(c_int) :: c_n_face_part_bound
    integer(c_int) :: c_n_vtx
    integer(c_int) :: c_n_proc
    integer(c_int) :: c_n_total_part
    integer(c_int) :: c_s_cell_face
    integer(c_int) :: c_s_face_vtx
    integer(c_int) :: c_s_face_bound
    integer(c_int) :: c_n_bound_groups
    integer(c_int) :: c_s_face_join
    integer(c_int) :: c_n_join_groups

    integer (kind = PDM_l_num_s), pointer :: n_elt(:) => null()
    type(c_ptr)                           :: c_n_elt

    c_n_elt = c_loc(n_elt)

    call PDM_multipart_part_dim_get_c(multipart, &
                                      i_zone, &
                                      i_part, &
                                      c_n_section, &
                                      c_n_elt, &
                                      c_n_cell, &
                                      c_n_face, &
                                      c_n_face_part_bound, &
                                      c_n_vtx, &
                                      c_n_proc, &
                                      c_n_total_part, &
                                      c_s_cell_face, &
                                      c_s_face_vtx, &
                                      c_s_face_bound, &
                                      c_n_bound_groups, &
                                      c_s_face_join, &
                                      c_n_join_groups)

    c_cell_color       = c_loc(cell_color)
    c_face_color       = c_loc(face_color)
    c_face_hp_color    = c_loc(face_hp_color)
    c_thread_color     = c_loc(thread_color)
    c_hyperplane_color = c_loc(hyperplane_color)

    call PDM_multipart_part_color_get_c(multipart, &
                                        i_zone, &
                                        i_part, &
                                        c_cell_color, &
                                        c_face_color, &
                                        c_face_hp_color, &
                                        c_thread_color, &
                                        c_hyperplane_color)

    call c_f_pointer(c_cell_color, &
                     cell_color,   &
                     [c_n_cell])

    call c_f_pointer(c_face_color, &
                     face_color,   &
                     [c_n_face])

    call c_f_pointer(c_face_hp_color, &
                     face_hp_color,   &
                     [c_n_face])

    call c_f_pointer(c_thread_color, &
                     thread_color,   &
                     [c_n_cell])

    call c_f_pointer(c_hyperplane_color, &
                     hyperplane_color,   &
                     [c_n_cell])

  end subroutine PDM_multipart_part_color_get_

  !>
  !!
  !! \brief Get ghost vertex information
  !!
  !! \param [in]   multipart               Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                  Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                  Partition index
  !! \param [out]  vtx_ghost_information   An integer for each vertex to describe it's kind
  !!

  subroutine PDM_multipart_part_ghost_infomation_get_(multipart, &
                                                      i_zone, &
                                                      i_part, &
                                                      vtx_ghost_information)

    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr),    value              :: multipart
    integer(c_int), value              :: i_zone
    integer(c_int), value              :: i_part
    integer(kind=PDM_l_num_s), pointer :: vtx_ghost_information(:)
    type(c_ptr)                        :: c_vtx_ghost_information = C_NULL_PTR

    integer(c_int) :: c_n_section
    integer(c_int) :: c_n_cell
    integer(c_int) :: c_n_face
    integer(c_int) :: c_n_face_part_bound
    integer(c_int) :: c_n_vtx
    integer(c_int) :: c_n_proc
    integer(c_int) :: c_n_total_part
    integer(c_int) :: c_s_cell_face
    integer(c_int) :: c_s_face_vtx
    integer(c_int) :: c_s_face_bound
    integer(c_int) :: c_n_bound_groups
    integer(c_int) :: c_s_face_join
    integer(c_int) :: c_n_join_groups

    integer (kind = PDM_l_num_s), pointer :: n_elt(:) => null()
    type(c_ptr)                           :: c_n_elt

    c_n_elt = c_loc(n_elt)

    call PDM_multipart_part_dim_get_c(multipart, &
                                      i_zone, &
                                      i_part, &
                                      c_n_section, &
                                      c_n_elt, &
                                      c_n_cell, &
                                      c_n_face, &
                                      c_n_face_part_bound, &
                                      c_n_vtx, &
                                      c_n_proc, &
                                      c_n_total_part, &
                                      c_s_cell_face, &
                                      c_s_face_vtx, &
                                      c_s_face_bound, &
                                      c_n_bound_groups, &
                                      c_s_face_join, &
                                      c_n_join_groups)

    c_vtx_ghost_information = c_loc(vtx_ghost_information)

    call PDM_multipart_part_ghost_infomation_get_c(multipart, &
                                                   i_zone, &
                                                   i_part, &
                                                   c_vtx_ghost_information)

    call c_f_pointer(c_vtx_ghost_information, &
                     vtx_ghost_information,   &
                     [c_n_vtx])

  end subroutine PDM_multipart_part_ghost_infomation_get_

  !>
  !!
  !! \brief ??
  !!
  !! \param [in]   multipart               Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                  Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                  Partition index
  !! \param [out]  vtx_coord               Coordinates of vertices
  !! \param [out]  ownership               Data ownership
  !! \param [out]  n_vtx                   Number of vertices
  !!

  subroutine PDM_multipart_part_vtx_coord_get_(multipart, &
                                               i_zone, &
                                               i_part, &
                                               vtx_coord, &
                                               ownership, &
                                               n_vtx)

    use iso_c_binding
    implicit none

    type(c_ptr),      value     :: multipart
    integer(c_int),   value     :: i_zone
    integer(c_int),   value     :: i_part
    double precision, pointer   :: vtx_coord(:)
    type(c_ptr)                 :: c_vtx_coord = C_NULL_PTR
    integer(c_int),   value     :: ownership
    integer(c_int),   value     :: n_vtx

    c_vtx_coord = c_loc(vtx_coord)

    n_vtx = PDM_multipart_part_vtx_coord_get_c(multipart, &
                                               i_zone, &
                                               i_part, &
                                               c_vtx_coord, &
                                               ownership)

    call c_f_pointer(c_vtx_coord, &
                     vtx_coord,   &
                     [n_vtx])

  end subroutine PDM_multipart_part_vtx_coord_get_

  !>
  !!
  !! \brief Get boundary information
  !!
  !! \param [in]   multipart               Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                  Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                  Partition index
  !! \param [out]  bound_type              Boundary type
  !! \param [out]  n_bound                 Number of boundaries
  !! \param [out]  bound_idx               Boundary index
  !! \param [out]  bound                   Boundaries
  !! \param [out]  bound_ln_to_gn          Boundary local number to global number
  !!

  subroutine PDM_multipart_bound_get_(multipart, &
                                      i_zone, &
                                      i_part, &
                                      bound_type, &
                                      n_bound, &
                                      bound_idx, &
                                      bound, &
                                      bound_ln_to_gn)

    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr),               value   :: multipart
    integer(c_int),            value   :: i_zone
    integer(c_int),            value   :: i_part
    integer(c_int)                     :: bound_type
    integer(c_int)                     :: n_bound
    integer(kind=PDM_l_num_s), pointer :: bound_idx(:)
    type(c_ptr)                        :: c_bound_idx = C_NULL_PTR
    integer(kind=PDM_l_num_s), pointer :: bound(:)
    type(c_ptr)                        :: c_bound = C_NULL_PTR
    integer(kind=PDM_g_num_s), pointer :: bound_ln_to_gn(:)
    type(c_ptr)                        :: c_bound_ln_to_gn = C_NULL_PTR

    integer :: taille

    c_bound_idx      = c_loc(bound_idx)
    c_bound          = c_loc(bound)
    c_bound_ln_to_gn = c_loc(bound_ln_to_gn)

    call PDM_multipart_bound_get_c(multipart, &
                                   i_zone, &
                                   i_part, &
                                   bound_type, &
                                   n_bound, &
                                   c_bound_idx, &
                                   c_bound, &
                                   c_bound_ln_to_gn)

    call c_f_pointer(c_bound_idx, &
                     bound_idx,   &
                     [n_bound + 1])

    taille = bound_idx(n_bound)

    call c_f_pointer(c_bound, &
                     bound,   &
                     [taille])

    call c_f_pointer(c_bound_ln_to_gn, &
                     bound_ln_to_gn,   &
                     [n_bound])

  end subroutine PDM_multipart_bound_get_

end module pdm_multipart
