#ifndef __PDM_PART_MESH_NODAL_TO_PART_MESH_PRIV_H__
#define __PDM_PART_MESH_NODAL_TO_PART_MESH_PRIV_H__

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2024       ONERA

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
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh.h"

#include "pdm_part_mesh_nodal_to_part_mesh.h"

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

typedef enum {

  PMN_TO_PM_STATE_NOT_REQ      = 0, /*!< Not requested and not done        */
  PMN_TO_PM_STATE_REQ          = 1, /*!< Requested but not yet done        */
  PMN_TO_PM_STATE_NOT_REQ_DONE = 2, /*!< Not requested but done (to clean) */
  PMN_TO_PM_STATE_REQ_DONE     = 3  /*!< Requested and done                */

} _pmn_to_m_state_t;


struct _pdm_part_mesh_nodal_to_part_mesh_t
{
  /* Options */
  _pmn_to_m_state_t      connectivity_state   [PDM_CONNECTIVITY_TYPE_MAX]; /*!< Part Mesh connectivities to build */
  _pmn_to_m_state_t      g_num_state          [PDM_MESH_ENTITY_MAX];       /*!< Global IDs to generate */
  _pmn_to_m_state_t      group_state          [PDM_BOUND_TYPE_MAX];        /*!< Groups to transfer */
  _pmn_to_m_state_t      part_comm_graph_state[PDM_BOUND_TYPE_MAX];        /*!< Inter-partition communication graphs to generate */
  PDM_bool_t             keep_link_elmt_to_entity;                         /*!< Keep link between nodal elements and their corresponding entities in the Part Mesh */
  PDM_ownership_t        vtx_ownership_pmesh;                              /*!< Part Mesh's ownership for vertices */

  /* Input */
  PDM_part_mesh_nodal_t *pmesh_nodal; /*!< Input Part Mesh Nodal */

  /* Output */
  PDM_part_mesh_t       *pmesh;       /*!< Output Part Mesh */

  PDM_ownership_t        owner_pmesh; /*!< Part Mesh ownership */

  /* Internal */
  PDM_MPI_Comm comm;
  int          n_part;
};

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_MESH_NODAL_TO_PART_MESH_PRIV_H__ */
