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
  !! \brief Construct the partitioned meshes on every zones
  !!
  !! \param [in]   multipart             Pointer to \ref PDM_multipart_t object
  !!

  subroutine PDM_multipart_run_ppart_c (multipart) &
  bind (c, name='PDM_multipart_run_ppart')

    use iso_c_binding
    implicit none

    type(c_ptr), value     :: multipart

  end subroutine PDM_multipart_run_ppart_c

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

 ! TO DO

end module pdm_part
