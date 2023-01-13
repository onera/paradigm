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
