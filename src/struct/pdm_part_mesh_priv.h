#ifndef __PDM_PART_MESH_PRIV_H__
#define __PDM_PART_MESH_PRIV_H__

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
 * \struct _pdm_part_mesh_t
 * \brief  Define a partition mesh. Arrays are shared
 *
 */

struct _pdm_part_mesh_t
{
  int                  n_part;
  int                  tn_part;
  PDM_MPI_Comm         comm;

  int                 *pn_entity[PDM_MESH_ENTITY_MAX];                       /* Size for each entity (size = PDM_MESH_ENTITY_MAX)            */

  int                **pconnectivity    [PDM_CONNECTIVITY_TYPE_MAX];                   /* Array of connectivty (size = PDM_CONNECTIVITY_TYPE_MAX)            */
  int                **pconnectivity_idx[PDM_CONNECTIVITY_TYPE_MAX];               /* Array of connectivty_idx if any (size = PDM_CONNECTIVITY_TYPE_MAX) */

  PDM_g_num_t        **pentity_ln_to_gn[PDM_MESH_ENTITY_MAX];                /* Array of connectivty (size = PDM_MESH_ENTITY_MAX)            */
  int                **pentity_color   [PDM_MESH_ENTITY_MAX];

  double             **vtx_coords;

  PDM_bool_t          is_owner_connectivity[PDM_CONNECTIVITY_TYPE_MAX];
  PDM_bool_t          is_owner_ln_to_gn    [PDM_MESH_ENTITY_MAX];
  PDM_bool_t          is_owner_color       [PDM_MESH_ENTITY_MAX];
  PDM_bool_t          is_owner_vtx_coord;

  int                  n_group_bnd    [PDM_BOUND_TYPE_MAX];     /*!< Number of group by elememnt type */
  int                **pn_bound       [PDM_BOUND_TYPE_MAX];
  PDM_g_num_t       ***pbound_ln_to_gn[PDM_BOUND_TYPE_MAX];
  int               ***pbound         [PDM_BOUND_TYPE_MAX];
  PDM_bool_t           is_owner_bound [PDM_BOUND_TYPE_MAX];

  int                **pconcat_bound_idx      [PDM_BOUND_TYPE_MAX];
  PDM_g_num_t        **pconcat_bound_ln_to_gn [PDM_BOUND_TYPE_MAX];
  int                **pconcat_bound          [PDM_BOUND_TYPE_MAX];
  PDM_bool_t           is_owner_concat_bound  [PDM_BOUND_TYPE_MAX];
  PDM_bool_t           is_compute_concat_bound[PDM_BOUND_TYPE_MAX];

  /* Comm graph */
  int                **ppart_bound_proc_idx[PDM_BOUND_TYPE_MAX];
  int                **ppart_bound_part_idx[PDM_BOUND_TYPE_MAX];
  int                **ppart_bound         [PDM_BOUND_TYPE_MAX];
  PDM_bool_t           is_owner_part_bound [PDM_BOUND_TYPE_MAX];

};


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DMESH_H__ */
