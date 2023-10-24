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

  interface PDM_multipart_block_set ; module procedure  &
    PDM_multipart_block_set_
  end interface

  interface PDM_multipart_part_dim_get ; module procedure  &
    PDM_multipart_part_dim_get_
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

  interface PDM_multipart_part_ghost_infomation_get ; module procedure  &
    PDM_multipart_part_ghost_infomation_get_
  end interface

  interface PDM_multipart_part_vtx_coord_get ; module procedure  &
    PDM_multipart_part_vtx_coord_get_
  end interface

  interface PDM_multipart_bound_get ; module procedure  &
    PDM_multipart_bound_get_
  end interface

  interface PDM_multipart_part_graph_comm_get ; module procedure  &
    PDM_multipart_part_graph_comm_get_
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
    type(c_ptr)            :: pmesh_nodal
    integer(c_int), value  :: ownership

  end subroutine PDM_multipart_get_part_mesh_nodal_c

  !>
  !! \brief Set block
  !!
  !! \param [in]   multipart              Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                 Id of zone
  !! \param [in]   dn_cell                Number of distributed cells
  !! \param [in]   dn_face                Number of distributed faces
  !! \param [in]   dn_vtx                 Number of distributed vertices
  !! \param [in]   n_face_group           Number of face groups
  !! \param [in]   dcell_face_idx         Distributed cell face connectivity index or NULL
  !!                                      (size : dn_cell + 1, numbering : 0 to n-1)
  !! \param [in]   dcell_face             Distributed cell face connectivity or NULL
  !!                                      (size : dface_vtx_idx[dn_cell], numbering : 1 to n)
  !! \param [in]   dface_cell             Distributed face cell connectivity or NULL
  !!                                      (size : 2 * dn_face, numbering : 1 to n)
  !! \param [in]   dface_vtx_idx          Distributed face to vertex connectivity index
  !!                                      (size : dn_face + 1, numbering : 0 to n-1)
  !! \param [in]   dface_vtx              Distributed face to vertex connectivity
  !!                                      (size : dface_vtx_idx[dn_face], numbering : 1 to n)
  !! \param [in]   dvtx_coord             Distributed vertex coordinates
  !!                                      (size : 3*dn_vtx)
  !! \param [in]   dface_group_idx        Index of distributed faces list of each group
  !!                                      (size = n_face_group + 1) or NULL
  !! \param [in]   dface_group            Distributed faces list of each group
  !!                                      (size = dface_group[dface_group_idx[n_face_group]], numbering : 1 to n)
  !!                                      or NULL
  !!

  subroutine PDM_multipart_block_set_c (multipart, &
                                        i_zone, &
                                        dn_cell, &
                                        dn_face, &
                                        dn_vtx, &
                                        n_face_group, &
                                        dcell_face_idx, &
                                        dcell_face, &
                                        dface_cell, &
                                        dface_vtx_idx, &
                                        dface_vtx, &
                                        dvtx_coord, &
                                        dface_group_idx, &
                                        dface_group) &
  bind (c, name='PDM_multipart_block_set')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: dn_cell
    integer(c_int), value  :: dn_face
    integer(c_int), value  :: dn_vtx
    integer(c_int), value  :: n_face_group
    type(c_ptr),    value  :: dcell_face_idx
    type(c_ptr),    value  :: dcell_face
    type(c_ptr),    value  :: dface_cell
    type(c_ptr),    value  :: dface_vtx_idx
    type(c_ptr),    value  :: dface_vtx
    type(c_ptr),    value  :: dvtx_coord
    type(c_ptr),    value  :: dface_group_idx
    type(c_ptr),    value  :: dface_group

  end subroutine PDM_multipart_block_set_c

  !>
  !!
  !! \brief Returns the dimensions of a given partition
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                Partition index
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
  !!
  !!

  subroutine PDM_multipart_part_dim_get_c (multipart, &
                                           i_zone, &
                                           i_part, &
                                           n_cell, &
                                           n_face, &
                                           n_face_part_bound, &
                                           n_vtx, &
                                           n_proc, &
                                           n_total_part, &
                                           s_cell_face, &
                                           s_face_vtx, &
                                           s_face_bound, &
                                           n_bound_groups) &
  bind (c, name='PDM_multipart_part_dim_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: i_part
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

  end subroutine PDM_multipart_part_dim_get_c

  !>
  !!
  !! \brief Returns the data arrays of a given partition
  !!
  !! \param [in]   multipart                Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                   Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                   Partition index
  !! \param [out]  cell_face_idx            Index of cell->face connectivity
  !! \param [out]  cell_face                Cell->face connectivity
  !! \param [out]  cell_ln_to_gn            Cell local number to global number
  !! \param [out]  face_cell                Face->cell connectivity
  !! \param [out]  face_vtx_idx             Index of face->vertex connectivity
  !! \param [out]  face_vtx                 Face->vertex connectivity
  !! \param [out]  face_ln_to_gn            Face local number to global number
  !! \param [out]  face_part_bound_proc_idx ??
  !! \param [out]  face_part_bound_part_idx ??
  !! \param [out]  face_part_bound          ??
  !! \param [out]  vtx                      Vertices coordinates
  !! \param [out]  vtx_ln_to_gn             Vertex local number to global number
  !! \param [out]  face_bound_idx           Index of face->boundary connectivity
  !! \param [out]  face_bound               Boundary faces
  !! \param [out]  face_bound_ln_to_gn      Boundary face local number to global number
  !!

  subroutine PDM_multipart_part_val_get_c (multipart, &
                                           i_zone, &
                                           i_part, &
                                           cell_face_idx, &
                                           cell_face, &
                                           cell_ln_to_gn, &
                                           face_cell, &
                                           face_vtx_idx, &
                                           face_vtx, &
                                           face_ln_to_gn, &
                                           face_part_bound_proc_idx, &
                                           face_part_bound_part_idx, &
                                           face_part_bound, &
                                           vtx, &
                                           vtx_ln_to_gn, &
                                           face_bound_idx, &
                                           face_bound, &
                                           face_bound_ln_to_gn) &
  bind (c, name='PDM_multipart_part_val_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: i_part
    type(c_ptr)            :: cell_face_idx
    type(c_ptr)            :: cell_face
    type(c_ptr)            :: cell_ln_to_gn
    type(c_ptr)            :: face_cell
    type(c_ptr)            :: face_vtx_idx
    type(c_ptr)            :: face_vtx
    type(c_ptr)            :: face_ln_to_gn
    type(c_ptr)            :: face_part_bound_proc_idx
    type(c_ptr)            :: face_part_bound_part_idx
    type(c_ptr)            :: face_part_bound
    type(c_ptr)            :: vtx
    type(c_ptr)            :: vtx_ln_to_gn
    type(c_ptr)            :: face_bound_idx
    type(c_ptr)            :: face_bound
    type(c_ptr)            :: face_bound_ln_to_gn

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
  !! \brief Get the connection graph between partition for the requested bound type
  !!

  subroutine PDM_multipart_part_graph_comm_get_c (multipart,            &
                                                  i_zone,               &
                                                  i_part,               &
                                                  bound_type,           &
                                                  ppart_bound_proc_idx, &
                                                  ppart_bound_part_idx, &
                                                  ppart_bound,          &
                                                  ownership)            &
  bind (c, name='PDM_multipart_part_graph_comm_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: i_part
    integer(c_int), value  :: bound_type
    type(c_ptr)            :: ppart_bound_proc_idx
    type(c_ptr)            :: ppart_bound_part_idx
    type(c_ptr)            :: ppart_bound
    integer(c_int), value  :: ownership

  end subroutine PDM_multipart_part_graph_comm_get_c

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

  subroutine PDM_multipart_compute (multipart) &
  bind (c, name='PDM_multipart_compute')

    use iso_c_binding
    implicit none

    type(c_ptr), value  :: multipart

  end subroutine PDM_multipart_compute

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
           PDM_multipart_block_set_,&
           PDM_multipart_part_dim_get_,&
           PDM_multipart_part_val_get_,&
           PDM_multipart_part_connectivity_get_,&
           PDM_multipart_part_ln_to_gn_get_,&
           PDM_multipart_partition_color_get_,&
           PDM_multipart_part_ghost_infomation_get_,&
           PDM_multipart_part_vtx_coord_get_,&
           PDM_multipart_bound_get_,&
           PDM_multipart_part_graph_comm_get_

contains

  ! Build a multipart structure instance

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

    type(c_ptr)                        :: multipart                     ! Pointer to a new \ref PDM_multipart_t object
    integer(c_int),            value   :: n_zone                        ! Number of zones in the original mesh
    integer(kind=PDM_l_num_s), pointer :: n_part(:)                     ! Number of partition per proc in each zone
    type(c_ptr)                        :: c_n_part         = C_NULL_PTR
    integer(c_int),            value   :: merge_blocks                  ! Merge or not the zones before splitting
    integer(c_int),            value   :: split_method                  ! Choice of method used to split the mesh
    integer(c_int),            value   :: part_size_method              ! Choice of homogeneous or heterogeneous partitions
    double precision,          pointer :: part_fraction(:)              ! Weight (in %) of each partition in heterogeneous case
    type(c_ptr)                        :: c_part_fraction  = C_NULL_PTR
    integer(c_int),            value   :: comm                          ! PDM_MPI communicator
    integer(c_int),            value   :: owner                         ! Data ownership
    integer(c_int)                     :: c_comm

    c_n_part        = c_loc(n_part)
    c_part_fraction = c_loc(part_fraction)
    c_comm = PDM_MPI_Comm_f2c(comm)

    multipart = PDM_multipart_create_c(n_zone, &
                                       c_n_part, &
                                       merge_blocks, &
                                       split_method, &
                                       part_size_method, &
                                       c_part_fraction, &
                                       c_comm, &
                                       owner)

  end subroutine PDM_multipart_create_

  ! Set connecting data between all the zones

  subroutine PDM_multipart_register_joins_ (multipart, &
                                            n_total_joins, &
                                            join_to_opposite)

    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr),               value   :: multipart                       ! Pointer to \ref PDM_multipart_t object
    integer(c_int),            value   :: n_total_joins                   ! Total number of interfaces
    integer(kind=PDM_l_num_s), pointer :: join_to_opposite(:)             ! For each global join id, give the global id of the opposite join (size = n_total_joins)
    type(c_ptr)                        :: c_join_to_opposite = C_NULL_PTR

    c_join_to_opposite = c_loc(join_to_opposite)

    call PDM_multipart_register_joins_c(multipart, &
                                        n_total_joins, &
                                        c_join_to_opposite)

  end subroutine PDM_multipart_register_joins_

  ! Set the reordering methods to be used after partitioning

  subroutine PDM_multipart_set_reordering_options_ (multipart, &
                                                    i_zone, &
                                                    renum_cell_method, &
                                                    renum_cell_properties, &
                                                    renum_face_method)

    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr),               value   :: multipart                            ! Pointer to \ref PDM_multipart_t object
    integer(c_int),            value   :: i_zone                               ! Id of zone which parameters apply (or -1 for all zones)
    character (len=*)                  :: renum_cell_method                    ! Choice of renumbering method for cells
    integer(kind=PDM_l_num_s), pointer :: renum_cell_properties(:)             ! Parameters used by cacheblocking method : [n_cell_per_cache_wanted, is_asynchrone, is_vectorisation, n_vect_face, split_method]
    type(c_ptr)                        :: c_renum_cell_properties = C_NULL_PTR
    character (len=*)                  :: renum_face_method                    ! Choice of renumbering method for faces

    c_renum_cell_properties = c_loc(renum_cell_properties)

    call PDM_multipart_set_reordering_options_c(multipart, &
                                                i_zone, &
                                                trim(renum_cell_method)//C_NULL_CHAR, &
                                                c_renum_cell_properties, &
                                                trim(renum_face_method)//C_NULL_CHAR)

  end subroutine PDM_multipart_set_reordering_options_

  ! Set the reordering methods to be used after partitioning

  subroutine PDM_multipart_set_reordering_options_vtx_ (multipart, &
                                                        i_zone, &
                                                        renum_vtx_method)

    use iso_c_binding
    implicit none

    type(c_ptr),      value  :: multipart        ! Pointer to \ref PDM_multipart_t object
    integer(c_int),   value  :: i_zone           ! Id of zone which parameters apply (or -1 for all zones)
    character (len=*)        :: renum_vtx_method ! Choice of renumbering method for vertices

    call PDM_multipart_set_reordering_options_vtx_c(multipart, &
                                                    i_zone, &
                                                    trim(renum_vtx_method)//C_NULL_CHAR)

  end subroutine PDM_multipart_set_reordering_options_vtx_

  ! Get the \ref PDM_part_mesh_nodal_t object

  subroutine PDM_multipart_get_part_mesh_nodal_ (multipart, &
                                                 i_zone, &
                                                 pmesh_nodal, &
                                                 ownership)

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart   ! Pointer to \ref PDM_multipart_t object
    integer(c_int), value  :: i_zone      ! Id of zone which parameters apply (or -1 for all zones)
    type(c_ptr)            :: pmesh_nodal ! Partitionned nodal mesh
    integer(c_int), value  :: ownership   ! Data ownership

    call PDM_multipart_get_part_mesh_nodal_c(multipart, &
                                             i_zone, &
                                             pmesh_nodal, &
                                             ownership)

  end subroutine PDM_multipart_get_part_mesh_nodal_

  ! Set block data

  subroutine PDM_multipart_block_set_ (multipart, &
                                       i_zone, &
                                       dn_cell, &
                                       dn_face, &
                                       dn_vtx, &
                                       n_face_group, &
                                       dcell_face_idx, &
                                       dcell_face, &
                                       dface_cell, &
                                       dface_vtx_idx, &
                                       dface_vtx, &
                                       dvtx_coord, &
                                       dface_group_idx, &
                                       dface_group)

    use iso_c_binding
    implicit none

    type(c_ptr)                        :: multipart          ! Pointer to \ref PDM_multipart_t object
    integer, intent(in)                :: i_zone             ! Id of zone
    integer, intent(in)                :: dn_cell            ! Number of distributed cells
    integer, intent(in)                :: dn_face            ! Number of distributed faces
    integer, intent(in)                :: dn_vtx             ! Number of distributed vertices
    integer, intent(in)                :: n_face_group       ! Number of face groups
    integer(kind=PDM_l_num_s), pointer :: dcell_face_idx(:)  ! Distributed cell face connectivity index or NULL (size : dn_cell + 1, numbering : 0 to n-1)
    integer(kind=PDM_g_num_s), pointer :: dcell_face(:)      ! Distributed cell face connectivity or NULL (size : dface_vtx_idx(dn_cell+1), numbering : 1 to n)
    integer(kind=PDM_g_num_s), pointer :: dface_cell(:)      ! Distributed face cell connectivity
    integer(kind=PDM_l_num_s), pointer :: dface_vtx_idx(:)   ! Distributed face to vertex connectivity index (size : dn_face + 1, numbering : 0 to n-1)
    integer(kind=PDM_g_num_s), pointer :: dface_vtx(:)       ! Distributed face to vertex connectivity (size : dface_vtx_idx(dn_face+1), numbering : 1 to n)
    double precision,          pointer :: dvtx_coord(:,:)    ! Distributed vertex coordinates (shape = [3, n_vtx])
    integer(kind=PDM_l_num_s), pointer :: dface_group_idx(:) ! Index of distributed faces list of each group (size = n_face_group + 1) or NULL
    integer(kind=PDM_g_num_s), pointer :: dface_group(:)     ! Distributed faces list of each group or NULL (size = dface_group(dface_group_idx(n_face_group+1)+1), numbering : 1 to n)

    integer(kind=c_int)                :: c_i_zone
    integer(kind=c_int)                :: c_dn_cell
    integer(kind=c_int)                :: c_dn_face
    integer(kind=c_int)                :: c_dn_vtx
    integer(kind=c_int)                :: c_n_face_group
    type(c_ptr)                        :: c_dcell_face_idx  = C_NULL_PTR
    type(c_ptr)                        :: c_dcell_face      = C_NULL_PTR
    type(c_ptr)                        :: c_dface_cell      = C_NULL_PTR
    type(c_ptr)                        :: c_dface_vtx_idx   = C_NULL_PTR
    type(c_ptr)                        :: c_dface_vtx       = C_NULL_PTR
    type(c_ptr)                        :: c_dvtx_coord      = C_NULL_PTR
    type(c_ptr)                        :: c_dface_group_idx = C_NULL_PTR
    type(c_ptr)                        :: c_dface_group     = C_NULL_PTR

    c_i_zone          = i_zone
    c_dn_cell         = dn_cell
    c_dn_face         = dn_face
    c_dn_vtx          = dn_vtx
    c_n_face_group    = n_face_group

    c_dcell_face_idx  = c_loc(dcell_face_idx )
    c_dcell_face      = c_loc(dcell_face     )
    c_dface_vtx_idx   = c_loc(dface_vtx_idx  )
    c_dface_vtx       = c_loc(dface_vtx      )
    c_dvtx_coord      = c_loc(dvtx_coord     )
    c_dface_group_idx = c_loc(dface_group_idx)
    c_dface_group     = c_loc(dface_group    )

    if (associated(dface_cell)) then
      c_dface_cell = c_loc(dface_cell)
    endif

    call PDM_multipart_block_set_c (multipart, &
                                    c_i_zone, &
                                    c_dn_cell, &
                                    c_dn_face, &
                                    c_dn_vtx, &
                                    c_n_face_group, &
                                    c_dcell_face_idx, &
                                    c_dcell_face, &
                                    c_dface_cell, &
                                    c_dface_vtx_idx, &
                                    c_dface_vtx, &
                                    c_dvtx_coord, &
                                    c_dface_group_idx, &
                                    c_dface_group)

  end subroutine PDM_multipart_block_set_

  ! Returns the dimensions of a given partition

  subroutine PDM_multipart_part_dim_get_ (multipart, &
                                          i_zone, &
                                          i_part, &
                                          n_cell, &
                                          n_face, &
                                          n_face_part_bound, &
                                          n_vtx, &
                                          n_proc, &
                                          n_total_part, &
                                          s_cell_face, &
                                          s_face_vtx, &
                                          s_face_bound, &
                                          n_bound_groups)

    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr), value       :: multipart         ! Pointer to \ref PDM_multipart_t object
    integer,     intent(in)  :: i_zone            ! Id of zone which parameters apply (or -1 for all zones)
    integer,     intent(in)  :: i_part            ! Partition index
    integer,     intent(out) :: n_cell            ! Number of cells
    integer,     intent(out) :: n_face            ! Number of faces
    integer,     intent(out) :: n_face_part_bound ! Number of boundary faces in partition
    integer,     intent(out) :: n_vtx             ! Number of vertices
    integer,     intent(out) :: n_proc            ! Number of processes
    integer,     intent(out) :: n_total_part      ! Total number of partitions
    integer,     intent(out) :: s_cell_face       ! Size of cell->face connectivity
    integer,     intent(out) :: s_face_vtx        ! Size of face->vtx connectivity
    integer,     intent(out) :: s_face_bound      ! Size of boundary faces
    integer,     intent(out) :: n_bound_groups    ! Number of boundary groups

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

    call PDM_multipart_part_dim_get_c(multipart, &
                                      i_zone, &
                                      i_part, &
                                      c_n_cell, &
                                      c_n_face, &
                                      c_n_face_part_bound, &
                                      c_n_vtx, &
                                      c_n_proc, &
                                      c_n_total_part, &
                                      c_s_cell_face, &
                                      c_s_face_vtx, &
                                      c_s_face_bound, &
                                      c_n_bound_groups)

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

  end subroutine PDM_multipart_part_dim_get_

  ! Returns the data arrays of a given partition (Deprecated)

  subroutine PDM_multipart_part_val_get_(multipart, &
                                         i_zone, &
                                         i_part, &
                                         cell_face_idx, &
                                         cell_face, &
                                         cell_ln_to_gn, &
                                         face_cell, &
                                         face_vtx_idx, &
                                         face_vtx, &
                                         face_ln_to_gn, &
                                         face_part_bound_proc_idx, &
                                         face_part_bound_part_idx, &
                                         face_part_bound, &
                                         vtx, &
                                         vtx_ln_to_gn, &
                                         face_bound_idx, &
                                         face_bound, &
                                         face_bound_ln_to_gn)

    use pdm
    use iso_c_binding
    use pdm_pointer_array
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: i_part

    integer (kind = PDM_l_num_s), pointer :: cell_face_idx(:)
    integer (kind = PDM_l_num_s), pointer :: cell_face(:)
    integer (kind = PDM_g_num_s), pointer :: cell_ln_to_gn(:)
    integer (kind = PDM_l_num_s), pointer :: face_cell(:)
    integer (kind = PDM_l_num_s), pointer :: face_vtx_idx(:)
    integer (kind = PDM_l_num_s), pointer :: face_vtx(:)
    integer (kind = PDM_g_num_s), pointer :: face_ln_to_gn(:)
    integer (kind = PDM_l_num_s), pointer :: face_part_bound_proc_idx(:)
    integer (kind = PDM_l_num_s), pointer :: face_part_bound_part_idx(:)
    integer (kind = PDM_l_num_s), pointer :: face_part_bound(:)
    double precision,             pointer :: vtx(:,:)
    integer (kind = PDM_g_num_s), pointer :: vtx_ln_to_gn(:)
    integer (kind = PDM_l_num_s), pointer :: face_bound_idx(:)
    integer (kind = PDM_l_num_s), pointer :: face_bound(:)
    integer (kind = PDM_g_num_s), pointer :: face_bound_ln_to_gn(:)

    type(c_ptr)  :: c_cell_face_idx            = C_NULL_PTR
    type(c_ptr)  :: c_cell_face                = C_NULL_PTR
    type(c_ptr)  :: c_cell_ln_to_gn            = C_NULL_PTR
    type(c_ptr)  :: c_face_cell                = C_NULL_PTR
    type(c_ptr)  :: c_face_vtx_idx             = C_NULL_PTR
    type(c_ptr)  :: c_face_vtx                 = C_NULL_PTR
    type(c_ptr)  :: c_face_ln_to_gn            = C_NULL_PTR
    type(c_ptr)  :: c_face_part_bound_proc_idx = C_NULL_PTR
    type(c_ptr)  :: c_face_part_bound_part_idx = C_NULL_PTR
    type(c_ptr)  :: c_face_part_bound          = C_NULL_PTR
    type(c_ptr)  :: c_vtx                      = C_NULL_PTR
    type(c_ptr)  :: c_vtx_ln_to_gn             = C_NULL_PTR
    type(c_ptr)  :: c_face_bound_idx           = C_NULL_PTR
    type(c_ptr)  :: c_face_bound               = C_NULL_PTR
    type(c_ptr)  :: c_face_bound_ln_to_gn      = C_NULL_PTR

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

    integer :: n_cell
    integer :: n_face
    integer :: n_face_part_bound
    integer :: n_vtx
    integer :: n_bound_groups

    call PDM_multipart_part_dim_get_c(multipart, &
                                      i_zone, &
                                      i_part, &
                                      c_n_cell, &
                                      c_n_face, &
                                      c_n_face_part_bound, &
                                      c_n_vtx, &
                                      c_n_proc, &
                                      c_n_total_part, &
                                      c_s_cell_face, &
                                      c_s_face_vtx, &
                                      c_s_face_bound, &
                                      c_n_bound_groups)

    n_cell             = c_n_cell
    n_face             = c_n_face
    n_face_part_bound  = c_n_face_part_bound
    n_vtx              = c_n_vtx
    n_bound_groups     = c_n_bound_groups

    call PDM_multipart_part_val_get_c(multipart, &
                                      i_zone, &
                                      i_part, &
                                      c_cell_face_idx, &
                                      c_cell_face, &
                                      c_cell_ln_to_gn, &
                                      c_face_cell, &
                                      c_face_vtx_idx, &
                                      c_face_vtx, &
                                      c_face_ln_to_gn, &
                                      c_face_part_bound_proc_idx, &
                                      c_face_part_bound_part_idx, &
                                      c_face_part_bound, &
                                      c_vtx, &
                                      c_vtx_ln_to_gn, &
                                      c_face_bound_idx, &
                                      c_face_bound, &
                                      c_face_bound_ln_to_gn)

    call c_f_pointer(c_cell_face_idx, &
                     cell_face_idx,   &
                     [n_cell + 1])

    call c_f_pointer(c_cell_face, &
                     cell_face,   &
                     [cell_face_idx(n_cell+1)])

    call c_f_pointer(c_cell_ln_to_gn, &
                     cell_ln_to_gn,   &
                     [n_cell])

    call c_f_pointer(c_face_cell, &
                     face_cell,   &
                     [2 * n_face])

    call c_f_pointer(c_face_vtx_idx, &
                     face_vtx_idx,   &
                     [n_face + 1])

    call c_f_pointer(c_face_vtx, &
                     face_vtx,   &
                     [face_vtx_idx(n_face+1)])

    call c_f_pointer(c_face_ln_to_gn, &
                     face_ln_to_gn,   &
                     [n_face])

    call c_f_pointer(c_face_part_bound_proc_idx, &
                     face_part_bound_proc_idx,   &
                     [n_face_part_bound + 1])

    call c_f_pointer(c_face_part_bound_part_idx, &
                     face_part_bound_part_idx,   &
                     [n_face_part_bound + 1])

    call c_f_pointer(c_face_part_bound, &
                     face_part_bound,   &
                     [face_part_bound_part_idx(n_face_part_bound+1)])

    call c_f_pointer(c_vtx, &
                     vtx,   &
                     [3, n_vtx])

    call c_f_pointer(c_vtx_ln_to_gn, &
                     vtx_ln_to_gn,   &
                     [n_vtx])

    call c_f_pointer(c_face_bound_idx, &
                     face_bound_idx,   &
                     [n_bound_groups + 1])

    call c_f_pointer(c_face_bound, &
                     face_bound,   &
                     [face_bound_idx(n_bound_groups+1)])

    call c_f_pointer(c_face_bound_ln_to_gn, &
                     face_bound_ln_to_gn,   &
                     [face_bound_idx(n_bound_groups+1)])

  end subroutine PDM_multipart_part_val_get_

  ! Returns the dimensions of a given partition

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

    type(c_ptr),               value   :: multipart                  ! Pointer to \ref PDM_multipart_t object
    integer(c_int),            value   :: i_zone                     ! Id of zone which parameters apply (or -1 for all zones)
    integer(c_int),            value   :: i_part                     ! Partition index
    integer(c_int),            value   :: connectivity_type          ! Type of connectivity to be getted (enumerated type)
    integer(kind=PDM_l_num_s), pointer :: connect(:)                 ! Connectivity
    type(c_ptr)                        :: c_connect = C_NULL_PTR
    integer(kind=PDM_l_num_s), pointer :: connect_idx(:)             ! Connectivity index
    type(c_ptr)                        :: c_connect_idx = C_NULL_PTR
    integer(c_int),            value   :: ownership                  ! Data ownership
    integer(c_int)                     :: pn_entity                  ! Number of entities
    integer(c_int)                     :: connec_size

    pn_entity = PDM_multipart_part_connectivity_get_c(multipart,         &
                                                      i_zone,            &
                                                      i_part,            &
                                                      connectivity_type, &
                                                      c_connect,         &
                                                      c_connect_idx,     &
                                                      ownership)

    connect_idx => null()
    if ( c_associated(c_connect_idx) ) then
      call c_f_pointer(c_connect_idx, &
                       connect_idx,   &
                       [pn_entity+1])
    end if

    connec_size = 2 * pn_entity ! TO DO is c_connect_idx a C_NULL_PTR in other cases than edge_vtx ?
    if ( c_associated(c_connect_idx) ) then
      connec_size = connect_idx(pn_entity+1)
    end if

    call c_f_pointer(c_connect, &
                     connect,   &
                     [connec_size])

  end subroutine PDM_multipart_part_connectivity_get_

  ! Returns the dimensions of a given partition

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

    type(c_ptr),               value   :: multipart                      ! Pointer to \ref PDM_multipart_t object
    integer(c_int),            value   :: i_zone                         ! Id of zone which parameters apply (or -1 for all zones)
    integer(c_int),            value   :: i_part                         ! Partition index
    integer(c_int),            value   :: entity_type                    ! Type of entity to be getted
    integer(kind=PDM_g_num_s), pointer :: entity_ln_to_gn(:)             ! Entity local number to global number
    type(c_ptr)                        :: c_entity_ln_to_gn = C_NULL_PTR
    integer(c_int),            value   :: ownership                      ! Data ownership
    integer(c_int)                     :: pn_entity                      ! Number of entities

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

  ! Returns the partitions color

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

    type(c_ptr),               value   :: multipart                   ! Pointer to \ref PDM_multipart_t object
    integer(c_int),            value   :: i_zone                      ! Id of zone which parameters apply (or -1 for all zones)
    integer(c_int),            value   :: i_part                      ! Partition index
    integer(c_int),            value   :: entity_type                 ! Type of entity to be getted
    integer(kind=PDM_l_num_s), pointer :: entity_color(:)             ! Entity->color connectivity
    type(c_ptr)                        :: c_entity_color = C_NULL_PTR
    integer(c_int),            value   :: ownership                   ! Data ownership
    integer(c_int)                     :: pn_entity                   ! Number of entities

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

  ! Get ghost vertex information

  subroutine PDM_multipart_part_ghost_infomation_get_(multipart, &
                                                      i_zone, &
                                                      i_part, &
                                                      vtx_ghost_information)

    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr),    value              :: multipart                            ! Pointer to \ref PDM_multipart_t object
    integer(c_int), value              :: i_zone                               ! Id of zone which parameters apply (or -1 for all zones)
    integer(c_int), value              :: i_part                               ! Partition index
    integer(kind=PDM_l_num_s), pointer :: vtx_ghost_information(:)             ! An integer for each vertex to describe it's kind
    type(c_ptr)                        :: c_vtx_ghost_information = C_NULL_PTR

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

    call PDM_multipart_part_dim_get_c(multipart, &
                                      i_zone, &
                                      i_part, &
                                      c_n_cell, &
                                      c_n_face, &
                                      c_n_face_part_bound, &
                                      c_n_vtx, &
                                      c_n_proc, &
                                      c_n_total_part, &
                                      c_s_cell_face, &
                                      c_s_face_vtx, &
                                      c_s_face_bound, &
                                      c_n_bound_groups)

    call PDM_multipart_part_ghost_infomation_get_c(multipart, &
                                                   i_zone, &
                                                   i_part, &
                                                   c_vtx_ghost_information)

    call c_f_pointer(c_vtx_ghost_information, &
                     vtx_ghost_information,   &
                     [c_n_vtx])

  end subroutine PDM_multipart_part_ghost_infomation_get_

  ! Get partitionned mesh vertices coordiantes

  subroutine PDM_multipart_part_vtx_coord_get_(multipart, &
                                               i_zone, &
                                               i_part, &
                                               vtx_coord, &
                                               ownership, &
                                               n_vtx)

    use iso_c_binding
    implicit none

    type(c_ptr),      value     :: multipart                ! Pointer to \ref PDM_multipart_t object
    integer(c_int),   value     :: i_zone                   ! Id of zone which parameters apply (or -1 for all zones)
    integer(c_int),   value     :: i_part                   ! Partition index
    double precision, pointer   :: vtx_coord(:,:)           ! Coordinates of vertices
    type(c_ptr)                 :: c_vtx_coord = C_NULL_PTR
    integer(c_int),   value     :: ownership                ! Data ownership
    integer(c_int)              :: n_vtx                    ! Number of vertices

    n_vtx = PDM_multipart_part_vtx_coord_get_c(multipart, &
                                               i_zone, &
                                               i_part, &
                                               c_vtx_coord, &
                                               ownership)

    call c_f_pointer(c_vtx_coord, &
                     vtx_coord,   &
                     [3, n_vtx])

  end subroutine PDM_multipart_part_vtx_coord_get_

  ! Get boundary information

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

    type(c_ptr),               value   :: multipart                     ! Pointer to \ref PDM_multipart_t object
    integer(c_int),            value   :: i_zone                        ! Id of zone which parameters apply (or -1 for all zones)
    integer(c_int),            value   :: i_part                        ! Partition index
    integer(c_int)                     :: bound_type                    ! Boundary type
    integer(c_int)                     :: n_bound                       ! Number of boundaries
    integer(kind=PDM_l_num_s), pointer :: bound_idx(:)                  ! Boundary index
    type(c_ptr)                        :: c_bound_idx = C_NULL_PTR
    integer(kind=PDM_l_num_s), pointer :: bound(:)                      ! Boundaries
    type(c_ptr)                        :: c_bound = C_NULL_PTR
    integer(kind=PDM_g_num_s), pointer :: bound_ln_to_gn(:)             ! Boundary local number to global number
    type(c_ptr)                        :: c_bound_ln_to_gn = C_NULL_PTR

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

    call c_f_pointer(c_bound, &
                     bound,   &
                     [bound_idx(n_bound+1)])

    call c_f_pointer(c_bound_ln_to_gn, &
                     bound_ln_to_gn,   &
                     [bound_idx(n_bound+1)])

  end subroutine PDM_multipart_bound_get_

  ! Get the connection graph between partition for the requested bound type

  subroutine PDM_multipart_part_graph_comm_get_(multipart,            &
                                                i_zone,               &
                                                i_part,               &
                                                bound_type,           &
                                                ppart_bound_proc_idx, &
                                                ppart_bound_part_idx, &
                                                ppart_bound,          &
                                                ownership)

    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr),               value   :: multipart                     ! Pointer to \ref PDM_multipart_t object
    integer(c_int),            value   :: i_zone                        ! Id of zone which parameters apply (or -1 for all zones)
    integer(c_int),            value   :: i_part                        ! Partition index
    integer(c_int)                     :: bound_type                    ! Boundary type
    integer(kind=PDM_l_num_s), pointer :: ppart_bound_proc_idx(:)       ! Partitioning boundary entities index from process (size = n_proc + 1)
    integer(kind=PDM_l_num_s), pointer :: ppart_bound_part_idx(:)       ! Partitioning boundary entities index from partition (size = n_total_part + 1)
    integer(kind=PDM_l_num_s), pointer :: ppart_bound(:)                ! Partitioning boundary entities (size = 4 * n_entity_part_bound)
    integer(c_int)                     :: ownership                     ! Data ownership

    type(c_ptr)                        :: c_ppart_bound_proc_idx = C_NULL_PTR
    type(c_ptr)                        :: c_ppart_bound_part_idx = C_NULL_PTR
    type(c_ptr)                        :: c_ppart_bound          = C_NULL_PTR

    integer(c_int)                     :: n_cell
    integer(c_int)                     :: n_face
    integer(c_int)                     :: n_face_part_bound
    integer(c_int)                     :: n_vtx
    integer(c_int)                     :: n_proc
    integer(c_int)                     :: n_total_part
    integer(c_int)                     :: s_cell_face
    integer(c_int)                     :: s_face_vtx
    integer(c_int)                     :: s_face_bound
    integer(c_int)                     :: n_bound_groups

    call PDM_multipart_part_graph_comm_get_c(multipart,              &
                                             i_zone,                 &
                                             i_part,                 &
                                             bound_type,             &
                                             c_ppart_bound_proc_idx, &
                                             c_ppart_bound_part_idx, &
                                             c_ppart_bound,          &
                                             ownership)

    call PDM_multipart_part_dim_get_ (multipart, &
                                      i_zone, &
                                      i_part, &
                                      n_cell, &
                                      n_face, &
                                      n_face_part_bound, &
                                      n_vtx, &
                                      n_proc, &
                                      n_total_part, &
                                      s_cell_face, &
                                      s_face_vtx, &
                                      s_face_bound, &
                                      n_bound_groups)

    call c_f_pointer(c_ppart_bound_proc_idx, &
                     ppart_bound_proc_idx,   &
                     [n_proc+1])

    call c_f_pointer(c_ppart_bound_part_idx  , &
                     ppart_bound_part_idx  ,   &
                     [n_total_part+1])

    call c_f_pointer(c_ppart_bound, &
                     ppart_bound,   &
                     [4 * ppart_bound_part_idx(n_total_part+1)])

  end subroutine PDM_multipart_part_graph_comm_get_

end module pdm_multipart
