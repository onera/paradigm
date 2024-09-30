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
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_priv.h"

#include "pdm_error.h"
#include "pdm_logging.h"

#include "pdm_array.h"
#include "pdm_mesh_nodal.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_vtk.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_part_mesh_nodal_elmts_utils.h"

#include "pdm_part_to_block.h"
#include "pdm_multipart.h"

#include "pdm_part_mesh_nodal_to_part_mesh.h"
#include "pdm_part_mesh_nodal_to_part_mesh_priv.h"

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

#define CHECK_PMN_TO_PM(pmn_to_pm)                                       \
  if ((pmn_to_pm) == NULL) {                                             \
    PDM_error(__FILE__, __LINE__, 0,                                     \
              "Invalid PDM_part_mesh_nodal_to_part_mesh_t instance.\n"); \
  }

/*=============================================================================
 * Local structure definitions
 *============================================================================*/


/*============================================================================
 * Global variable
 *============================================================================*/


/*=============================================================================
 * Private function definitions
 *============================================================================*/


/*=============================================================================
 * Public function prototypes
 *============================================================================*/

PDM_part_mesh_nodal_to_part_mesh_t *
PDM_part_mesh_nodal_to_part_mesh_create
(
  PDM_part_mesh_nodal_t *pmesh_nodal,
  PDM_bool_t             keep_link_elmt_to_entity
)
{
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm = NULL;

  PDM_malloc(pmn_to_pm, 1, PDM_part_mesh_nodal_to_part_mesh_t);

  memset(pmn_to_pm, 0, sizeof(PDM_part_mesh_nodal_to_part_mesh_t));

  pmn_to_pm->keep_link_elmt_to_entity = keep_link_elmt_to_entity;
  pmn_to_pm->owner_pmesh              = PDM_OWNERSHIP_BAD_VALUE;

  return pmn_to_pm;
}


void
PDM_part_mesh_nodal_to_part_mesh_enable_connectivity
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm,
  PDM_connectivity_type_t             connectivity_type
)
{
  CHECK_PMN_TO_PM(pmn_to_pm);
  // TODO: error/warning if requested connectivity is not handled

  pmn_to_pm->build_connectivity[connectivity_type] = PDM_TRUE;
}


void
PDM_part_mesh_nodal_to_part_mesh_enable_g_nums
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm,
  PDM_mesh_entities_t                 entity_type
)
{
  CHECK_PMN_TO_PM(pmn_to_pm);

  pmn_to_pm->compute_g_nums[entity_type] = PDM_TRUE;
}


void
PDM_part_mesh_nodal_to_part_mesh_enable_groups
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm,
  PDM_bound_type_t                    bound_type
)
{
  CHECK_PMN_TO_PM(pmn_to_pm);

  pmn_to_pm->transfer_groups[bound_type] = PDM_TRUE;
}


void
PDM_part_mesh_nodal_to_part_mesh_enable_part_comm_graph
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm,
  PDM_bound_type_t                    bound_type
)
{
  CHECK_PMN_TO_PM(pmn_to_pm);

  pmn_to_pm->build_part_comm_graph[bound_type] = PDM_TRUE;
}


void
PDM_part_mesh_nodal_to_part_mesh_compute
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm
)
{
  CHECK_PMN_TO_PM(pmn_to_pm);

  PDM_error(__FILE__, __LINE__, 0, "Not implemented yet ¯\\_(ツ)_/¯\n");
}


void
PDM_part_mesh_nodal_to_part_mesh_part_mesh_get
(
  PDM_part_mesh_nodal_to_part_mesh_t  *pmn_to_pm,
  PDM_part_mesh_t                    **pmesh,
  PDM_ownership_t                      ownership
)
{
  CHECK_PMN_TO_PM(pmn_to_pm);

  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    pmn_to_pm->owner_pmesh = ownership;
  }

  *pmesh = pmn_to_pm->pmesh;
}


void
PDM_part_mesh_nodal_to_part_mesh_free
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm
)
{
  CHECK_PMN_TO_PM(pmn_to_pm);

  if (pmn_to_pm->owner_pmesh == PDM_OWNERSHIP_KEEP) {
    PDM_part_mesh_free(pmn_to_pm->pmesh);
  }

  PDM_free(pmn_to_pm);
}


#undef CHECK_PMN_TO_PM

#ifdef  __cplusplus
}
#endif
