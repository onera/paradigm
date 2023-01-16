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

module pdm_part

  use pdm

  implicit none

  integer, parameter :: PDM_PART_SIZE_HOMOGENEOUS   = 1
  integer, parameter :: PDM_PART_SIZE_HETEROGENEOUS = 2

  interface PDM_multipart_create ; module procedure  &
    PDM_multipart_create_
  end interface

  interface PDM_multipart_register_block ; module procedure  &
    PDM_multipart_register_block_
  end interface

  interface PDM_multipart_register_dmesh_nodal ; module procedure  &
    PDM_multipart_register_dmesh_nodal_
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

  interface PDM_multipart_run_ppart ; module procedure  &
    PDM_multipart_run_ppart_
  end interface

  interface PDM_multipart_get_part_mesh_nodal ; module procedure  &
    PDM_multipart_get_part_mesh_nodal_
  end interface

  interface PDM_multipart_dn_entity_set ; module procedure  &
    PDM_multipart_dn_entity_set_
  end interface

  interface PDM_multipart_dconnectivity_set ; module procedure  &
    PDM_multipart_dconnectivity_set_
  end interface

  interface PDM_multipart_dvtx_coord_set ; module procedure  &
    PDM_multipart_dvtx_coord_set_
  end interface

  interface PDM_multipart_domain_interface_shared_set ; module procedure  &
    PDM_multipart_domain_interface_shared_set_
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

  interface PDM_multipart_part_graph_comm_vtx_data_get ; module procedure  &
    PDM_multipart_part_graph_comm_vtx_data_get_
  end interface

  interface PDM_multipart_part_color_get ; module procedure  &
    PDM_multipart_part_color_get_
  end interface

  interface PDM_multipart_part_ghost_infomation_get ; module procedure  &
    PDM_multipart_part_ghost_infomation_get_
  end interface

  interface PDM_multipart_time_get ; module procedure  &
    PDM_multipart_time_get_
  end interface

  interface PDM_multipart_free ; module procedure  &
    PDM_multipart_free_
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
  !! \brief Set distributed mesh data for the input zone
  !!
  !! \param [in]   multipart      Pointer to \ref PDM_multipart_t object
  !! \param [in]   zone_id        Global zone id
  !! \param [in]   dmesh          Distributed mesh structure
  !!

  subroutine PDM_multipart_register_block_c (multipart, &
                                             zone_id, &
                                             dmesh) &
  bind (c, name='PDM_multipart_register_block')

    use iso_c_binding
    implicit none

    type(c_ptr), value     :: multipart
    integer(c_int), value  :: zone_id
    type(c_ptr)            :: dmesh

  end subroutine PDM_multipart_register_block_c

  !>
  !!
  !! \brief Set distributed mesh data for the input zone
  !!
  !! \param [in]   multipart      Pointer to \ref PDM_multipart_t object
  !! \param [in]   zone_id        Global zone id
  !! \param [in]   dmesh_nodal    Distributed nodal mesh structure
  !!

  subroutine PDM_multipart_register_dmesh_nodal_c (multipart, &
                                                   zone_id, &
                                                   dmesh_nodal) &
  bind (c, name='PDM_multipart_register_dmesh_nodal')

    use iso_c_binding
    implicit none

    type(c_ptr), value     :: multipart
    integer(c_int), value  :: zone_id
    type(c_ptr)            :: dmesh_nodal

  end subroutine PDM_multipart_register_dmesh_nodal_c

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

    type(c_ptr), value     :: multipart
    integer(c_int), value  :: n_total_joins
    type(c_ptr)            :: join_to_opposite

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

    type(c_ptr), value     :: multipart
    integer(c_int), value  :: i_zone
    character(c_char)      :: renum_cell_method(*)
    type(c_ptr)            :: renum_cell_properties(*)
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

    type(c_ptr), value     :: multipart
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

    type(c_ptr), value     :: multipart
    integer(c_int), value  :: i_zone
    type(c_ptr)            :: pmesh_nodal
    type(c_ptr)            :: c_pmesh_nodal       = C_NULL_PTR
    integer(c_int), value  :: ownership

    ! TO DO: c_pmesh_nodal handle right
    if (associated(pmesh_nodal)) then
      c_pmesh_nodal = c_loc(pmesh_nodal)
    endif

  end subroutine PDM_multipart_get_part_mesh_nodal_c

  !>
  !!
  !! \brief Set number of element in the block entity
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   entity_type           Type of entity (can be cell/face/edge/vtx)
  !! \param [in]   dn_entity             Distributed number of entity in current process
  !!

  subroutine PDM_multipart_dn_entity_set_c (multipart, &
                                            i_zone, &
                                            entity_type, &
                                            dn_entity) &
  bind (c, name='PDM_multipart_dn_entity_set')

    use iso_c_binding
    implicit none

    type(c_ptr), value     :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: entity_type
    integer(c_int), value  :: dn_entity

  end subroutine PDM_multipart_dn_entity_set_c

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

    type(c_ptr), value     :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: connectivity_type
    type(c_ptr), value     :: dconnect
    type(c_ptr), value     :: dconnect_idx

  end subroutine PDM_multipart_dconnectivity_set_c

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

  subroutine PDM_multipart_dgroup_set_c (multipart, &
                                         i_zone, &
                                         bound_type, &
                                         dconnect, &
                                         dconnect_idx) &
  bind (c, name='PDM_multipart_dgroup_set')

    use iso_c_binding
    implicit none

    type(c_ptr), value     :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: bound_type
    type(c_ptr), value     :: dconnect
    type(c_ptr), value     :: dconnect_idx

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

    type(c_ptr), value     :: multipart
    integer(c_int), value  :: i_zone
    type(c_ptr), value     :: dvtx_coord

  end subroutine PDM_multipart_dvtx_coord_set_c

  !>
  !!
  !! \brief ??
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   ditrf                 ???
  !!

  subroutine PDM_multipart_domain_interface_shared_set_c (multipart, &
                                                          ditrf) &
  bind (c, name='PDM_multipart_domain_interface_shared_set')

    use iso_c_binding
    implicit none

    type(c_ptr), value     :: multipart
    type(c_ptr), value     :: ditrf

  end subroutine PDM_multipart_domain_interface_shared_set_c

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
  !! \param [out]  n_face_part_bound     ??
  !! \param [out]  n_vtx                 Number of vertices
  !! \param [out]  n_proc                Number of processes
  !! \param [out]  n_total_part          ??
  !! \param [out]  s_cell_face           ??
  !! \param [out]  s_face_vtx            ??
  !! \param [out]  s_face_bound          ??
  !! \param [out]  n_bound_groups        ??
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

    type(c_ptr), value     :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: i_part
    type(c_ptr)            :: n_section
    type(c_ptr)            :: n_elt
    type(c_ptr)            :: c_n_elt       = C_NULL_PTR
    type(c_ptr)            :: n_cell
    type(c_ptr)            :: n_face
    type(c_ptr)            :: n_face_part_bound
    type(c_ptr)            :: n_vtx
    type(c_ptr)            :: n_proc
    type(c_ptr)            :: n_total_part
    type(c_ptr)            :: s_cell_face
    type(c_ptr)            :: s_face_vtx
    type(c_ptr)            :: s_face_bound
    type(c_ptr)            :: n_bound_groups
    type(c_ptr)            :: s_face_join
    type(c_ptr)            :: n_join_groups

    ! TO DO: c_n_elt handle right
    if (associated(n_elt)) then
      c_n_elt = c_loc(n_elt)
    endif

  end subroutine PDM_multipart_part_dim_get_c

  !>
  !!
  !! \brief Returns the dimensions of a given partition
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                Partition index
  !! \param [out]  n_vtx_part_bound      ??
  !!

  subroutine PDM_multipart_part_graph_comm_vtx_dim_get_c (multipart, &
                                                          i_zone, &
                                                          i_part, &
                                                          n_vtx_part_bound) &
  bind (c, name='PDM_multipart_part_graph_comm_vtx_dim_get')

    use iso_c_binding
    implicit none

    type(c_ptr), value     :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: i_part
    type(c_ptr)            :: n_vtx_part_bound

  end subroutine PDM_multipart_part_graph_comm_vtx_dim_get_c

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

    type(c_ptr), value     :: multipart
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
  !! \param [out]  connectivity_type     ??
  !! \param [out]  connect               ??
  !! \param [out]  connect_idx           ??
  !! \param [out]  ownership             ??
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

    type(c_ptr), value     :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: i_part
    integer(c_int), value  :: connectivity_type
    type(c_ptr)            :: connect
    type(c_ptr)            :: connect_idx
    integer(c_int), value  :: ownership
    integer(c_int)         :: pn_entity

  end function PDM_multipart_part_connectivity_get_c

  !>
  !!
  !! \brief Returns the dimensions of a given partition
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                Partition index
  !! \param [out]  entity_type           ??
  !! \param [out]  entity_ln_to_gn       ??
  !! \param [out]  ownership             ??
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

    type(c_ptr), value     :: multipart
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
  !! \param [out]  entity_type           ??
  !! \param [out]  entity_color          ??
  !! \param [out]  ownership             ??
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

    type(c_ptr), value     :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: i_part
    integer(c_int), value  :: entity_type
    type(c_ptr)            :: entity_color
    integer(c_int), value  :: ownership
    integer(c_int)         :: pn_entity

  end function PDM_multipart_partition_color_get_c

  !>
  !!
  !! \brief ??
  !!
  !! \param [in]   multipart               Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                  Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                  Partition index
  !! \param [out]  vtx_part_bound_proc_idx ??
  !! \param [out]  vtx_part_bound_part_idx ??
  !! \param [out]  vtx_part_bound          ??
  !!

  subroutine PDM_multipart_part_graph_comm_vtx_data_get_c (multipart, &
                                                         i_zone, &
                                                         i_part, &
                                                         vtx_part_bound_proc_idx, &
                                                         vtx_part_bound_part_idx, &
                                                         vtx_part_bound) &
  bind (c, name='PDM_multipart_part_graph_comm_vtx_data_get')

    use iso_c_binding
    implicit none

    type(c_ptr), value     :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: i_part
    type(c_ptr)            :: vtx_part_bound_proc_idx
    type(c_ptr)            :: vtx_part_bound_part_idx
    type(c_ptr)            :: vtx_part_bound

  end subroutine PDM_multipart_part_graph_comm_vtx_data_get_c

  !>
  !!
  !! \brief ??
  !!
  !! \param [in]   multipart               Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                  Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                  Partition index
  !! \param [out]  cell_color              ??
  !! \param [out]  face_color              ??
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

    type(c_ptr), value     :: multipart
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
  !! \brief ??
  !!
  !! \param [in]   multipart               Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                  Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                  Partition index
  !! \param [out]  vtx_ghost_information   ??
  !!

  subroutine PDM_multipart_part_ghost_infomation_get_c (multipart, &
                                                        i_zone, &
                                                        i_part, &
                                                        vtx_ghost_information) &
  bind (c, name='PDM_multipart_part_ghost_infomation_get')

    use iso_c_binding
    implicit none

    type(c_ptr), value     :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: i_part
    type(c_ptr)            :: vtx_ghost_information

  end subroutine PDM_multipart_part_ghost_infomation_get_c

  !>
  !!
  !! \brief ??
  !!
  !! \param [in]   multipart  Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone     Id of zone which parameters apply (or -1 for all zones)
  !! \param [out]  elapsed    ??
  !! \param [out]  cpu        ??
  !! \param [out]  cpu_user   ??
  !! \param [out]  cpu_sys    ??
  !!

  subroutine PDM_multipart_time_get_c (multipart, &
                                       i_zone, &
                                       elapsed, &
                                       cpu, &
                                       cpu_user, &
                                       cpu_sys) &
  bind (c, name='PDM_multipart_time_get')

    use iso_c_binding
    implicit none

    type(c_ptr), value     :: multipart
    integer(c_int), value  :: i_zone
    type(c_ptr)            :: elapsed
    type(c_ptr)            :: cpu
    type(c_ptr)            :: cpu_user
    type(c_ptr)            :: cpu_sys

  end subroutine PDM_multipart_time_get_c

  !>
  !!
  !! \brief ??
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

    type(c_ptr), value     :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: i_part
    type(c_ptr)            :: vtx_coord
    integer(c_int)         :: ownership
    integer(c_int)         :: n_vtx

  end function PDM_multipart_part_vtx_coord_get_c

  !>
  !!
  !! \brief ??
  !!
  !! \param [in]   multipart               Pointer to \ref PDM_multipart_t object
  !! \param [in]   i_zone                  Id of zone which parameters apply (or -1 for all zones)
  !! \param [in]   i_part                  Partition index
  !! \param [out]  bound_type              ??
  !! \param [out]  n_bound                 ??
  !! \param [out]  bound_idx               ??
  !! \param [out]  bound                   ??
  !! \param [out]  bound_ln_to_gn          ??
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

    type(c_ptr), value     :: multipart
    integer(c_int), value  :: i_zone
    integer(c_int), value  :: i_part
    integer(c_int)         :: bound_type
    type(c_ptr)            :: n_bound
    type(c_ptr)            :: bound_idx
    type(c_ptr)            :: bound
    type(c_ptr)            :: bound_ln_to_gn

  end subroutine PDM_multipart_bound_get_c

end interface

private :: PDM_multipart_create_,&
           PDM_multipart_register_block_,&
           PDM_multipart_register_dmesh_nodal_,&
           PDM_multipart_register_joins_,&
           PDM_multipart_set_reordering_options_,&
           PDM_multipart_set_reordering_options_vtx_,&
           PDM_multipart_run_ppart_,&
           PDM_multipart_get_part_mesh_nodal_,&
           PDM_multipart_dn_entity_set_,&
           PDM_multipart_dconnectivity_set_,&
           PDM_multipart_dvtx_coord_set_,&
           PDM_multipart_domain_interface_shared_set_,&
           PDM_multipart_part_dim_get_,&
           PDM_multipart_part_graph_comm_vtx_dim_get_,&
           PDM_multipart_part_val_get_,&
           PDM_multipart_part_connectivity_get_,&
           PDM_multipart_part_ln_to_gn_get_,&
           PDM_multipart_partition_color_get_,&
           PDM_multipart_part_graph_comm_vtx_data_get_,&
           PDM_multipart_part_color_get_,&
           PDM_multipart_part_ghost_infomation_get_,&
           PDM_multipart_time_get_,&
           PDM_multipart_free_,&
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
    integer(c_int), value              :: n_zone
    integer(kind=PDM_l_num_s), pointer :: n_part(:)
    type(c_ptr)                        :: c_n_part
    integer(c_int), value              :: merge_blocks
    integer(c_int), value              :: split_method
    integer(c_int), value              :: part_size_method
    double precision, pointer          :: part_fraction(:)
    type(c_ptr)                        :: c_part_fraction
    integer(c_int), value              :: comm
    integer(c_int), value              :: owner

    c_n_part        = c_loc(n_part)
    c_part_fraction = c_loc(part_fraction)

    multipart = PDM_multipart_create_c(n_zone,
                                       c_n_part,
                                       merge_blocks,
                                       split_method,
                                       part_size_method,
                                       c_part_fraction,
                                       comm,
                                       owner)

  end subroutine PDM_multipart_create_

  !>
  !!
  !! \brief Construct the partitioned meshes on every zones
  !!
  !! \param [in]   multipart   Pointer to \ref PDM_multipart_t object
  !!

  subroutine PDM_multipart_run_ppart_ (multipart) &
  bind (c, name='PDM_multipart_run_ppart')

    use iso_c_binding
    implicit none

    type(c_ptr), value     :: multipart

  end subroutine PDM_multipart_run_ppart_

  !>
  !!
  !! \brief Free the structure
  !!
  !! \param [in]   multipart   Pointer to \ref PDM_multipart_t object
  !!

  subroutine PDM_multipart_free_ (multipart) &
  bind (c, name='PDM_multipart_free')

    use iso_c_binding
    implicit none

    type(c_ptr), value     :: multipart

  end subroutine PDM_multipart_free_

end module pdm_part
