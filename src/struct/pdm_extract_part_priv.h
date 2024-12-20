#ifndef __PDM_EXTRACT_PART_PRIV_H__
#define __PDM_EXTRACT_PART_PRIV_H__

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2017       ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_part_to_part.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_part_mesh_nodal.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _pdm_extract_part_t
 * \brief  Define a partition mesh. Arrays are shared
 *
 */

struct _pdm_extract_part_t
{

  int                     dim;
  int                     n_part_in;
  int                     n_part_out;
  PDM_extract_part_kind_t extract_kind;
  PDM_bool_t              compute_child_gnum;

  PDM_split_dual_t        split_dual_method;
  PDM_ownership_t         ownership;
  PDM_MPI_Comm            comm;

  int                     renum_method           [PDM_MESH_ENTITY_MAX];
  const int              *renum_method_properties[PDM_MESH_ENTITY_MAX];

  /* Partitioned view - To do with extract for selected gnum + part_to_part */
  int                 *n_cell;
  int                 *n_face;
  int                 *n_edge;
  int                 *n_vtx;
  int                **pcell_face;
  int                **pcell_face_idx;
  int                **pface_edge;
  int                **pface_edge_idx;
  int                **pedge_vtx;
  int                **pface_vtx;
  int                **pface_vtx_idx;

  PDM_g_num_t        **cell_ln_to_gn;
  PDM_g_num_t        **face_ln_to_gn;
  PDM_g_num_t        **edge_ln_to_gn;
  PDM_g_num_t        **vtx_ln_to_gn;

  double             **pvtx_coord;

  int                  have_user_entity_center;
  double             **entity_center;

  int                  n_group              [PDM_BOUND_TYPE_MAX];
  int                **n_group_entity       [PDM_BOUND_TYPE_MAX];
  int               ***group_entity         [PDM_BOUND_TYPE_MAX];
  PDM_g_num_t       ***group_entity_ln_to_gn[PDM_BOUND_TYPE_MAX];

  int                  have_connectivity[PDM_CONNECTIVITY_TYPE_MAX];


  /* If partition is described by elements */
  int is_nodal;
  // PDM_part_mesh_nodal_elmts_t *pmne;
  PDM_part_mesh_nodal_t       *pmn;


  /* Which cell or face is selected */
  int                 *n_extract;
  int                **extract_lnum;
  PDM_ownership_t      owner_extract_lnum;

  /* Each rank specify the target AND order they want */
  int                  from_target;
  int                 *n_target;
  PDM_g_num_t        **target_gnum;
  PDM_ownership_t      target_ownership;
  PDM_ownership_t      owner_target_gnum;
  int                **target_location;
  PDM_mesh_entities_t  master_entity;

  /* Extracted part */
  double             **pextract_vtx_coord;

  PDM_ownership_t     owner_connectivity   [PDM_CONNECTIVITY_TYPE_MAX];
  PDM_ownership_t     owner_ln_to_gn       [PDM_MESH_ENTITY_MAX];
  PDM_ownership_t     owner_parent_ln_to_gn[PDM_MESH_ENTITY_MAX];
  PDM_ownership_t     owner_parent_lnum    [PDM_MESH_ENTITY_MAX];
  PDM_ownership_t     owner_vtx_coord;
  PDM_ownership_t     owner_init_location  [PDM_MESH_ENTITY_MAX];

  /* Only for mapping and clear API */
  int                *pextract_n_entity              [PDM_MESH_ENTITY_MAX];
  int               **pextract_connectivity          [PDM_CONNECTIVITY_TYPE_MAX];
  int               **pextract_connectivity_idx      [PDM_CONNECTIVITY_TYPE_MAX];
  PDM_g_num_t       **pextract_entity_ln_to_gn       [PDM_MESH_ENTITY_MAX];
  PDM_g_num_t       **pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_MAX];
  int               **pextract_entity_parent_lnum    [PDM_MESH_ENTITY_MAX];
  int               **pextract_entity_init_location  [PDM_MESH_ENTITY_MAX];

  // For renumbering
  PDM_ownership_t     owner_color                    [PDM_MESH_ENTITY_MAX];
  int               **pextract_entity_color          [PDM_MESH_ENTITY_MAX];
  int               **pextract_entity_order          [PDM_MESH_ENTITY_MAX];

  /* If partition is described by elements */
  PDM_ownership_t        owner_extract_pmn;
  PDM_part_mesh_nodal_t *extract_pmn;

  /* Part-to-part objects */
  PDM_part_to_part_t *ptp_entity   [PDM_MESH_ENTITY_MAX];
  PDM_ownership_t     ptp_ownership[PDM_MESH_ENTITY_MAX];

  /* Part-to-part objects */
  PDM_part_to_part_t  **ptp_group_entity                     [PDM_BOUND_TYPE_MAX];
  PDM_ownership_t      *ptp_group_ownership                  [PDM_BOUND_TYPE_MAX];

  int                 **pn_extract_group_entity              [PDM_BOUND_TYPE_MAX];
  int                ***pextract_group_entity                [PDM_BOUND_TYPE_MAX];
  PDM_g_num_t        ***pextract_group_entity_ln_to_gn       [PDM_BOUND_TYPE_MAX];
  PDM_g_num_t        ***pextract_group_entity_parent_ln_to_gn[PDM_BOUND_TYPE_MAX];
  PDM_ownership_t      *group_array_ownership                [PDM_BOUND_TYPE_MAX];
  PDM_ownership_t      *group_parent_ownership               [PDM_BOUND_TYPE_MAX];

};


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_EXTRACT_PART_PRIV_H__ */
