/*
 * \file
 */

#ifndef __PDM_PART_MESH_NODAL_TO_PART_MESH_H__
#define __PDM_PART_MESH_NODAL_TO_PART_MESH_H__

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

typedef struct _pdm_part_mesh_nodal_to_part_mesh_t PDM_part_mesh_nodal_to_part_mesh_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


PDM_part_mesh_nodal_to_part_mesh_t *
PDM_part_mesh_nodal_to_part_mesh_create
(
  PDM_part_mesh_nodal_t *pmesh_nodal,
  PDM_bool_t             keep_link_elmt_to_entity
);


void
PDM_part_mesh_nodal_to_part_mesh_enable_connectivity
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm,
  PDM_connectivity_type_t             connectivity_type
);


void
PDM_part_mesh_nodal_to_part_mesh_enable_g_nums
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm,
  PDM_mesh_entities_t                 entity_type
);


void
PDM_part_mesh_nodal_to_part_mesh_enable_groups
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm,
  PDM_bound_type_t                    bound_type
);


void
PDM_part_mesh_nodal_to_part_mesh_enable_part_comm_graph
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm,
  PDM_bound_type_t                    bound_type
);


void
PDM_part_mesh_nodal_to_part_mesh_compute
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm
);


void
PDM_part_mesh_nodal_to_part_mesh_part_mesh_get
(
  PDM_part_mesh_nodal_to_part_mesh_t  *pmn_to_pm,
  PDM_part_mesh_t                    **pmesh,
  PDM_ownership_t                      ownership
);


void
PDM_part_mesh_nodal_to_part_mesh_free
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm
);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_MESH_NODAL_TO_PART_MESH_H__ */
