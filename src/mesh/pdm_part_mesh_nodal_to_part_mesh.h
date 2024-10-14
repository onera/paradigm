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

/**
 * \brief Create a structure for generating a \ref PDM_part_mesh_t instance from a \ref PDM_part_mesh_nodal_t instance
 *
 * \param [in] pmesh_nodal               Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in] keep_link_elmt_to_entity  Preserve link element (Part Mesh Nodal) -> entity (Part Mesh)
 * \param [in] vtx_ownership_pmesh       Part Mesh's ownership for vertices
 *
 * \return Pointer to a new \ref PDM_part_mesh_nodal_to_part_mesh_t instance.
 *
 * \note The vertices are identical in both data structures. In order to avoid ownership conflicts, the vertices are treated as follows :
 *   - If \p pmn is PDM_OWNERSHIP_KEEP :
 *     - if \p vtx_ownership_pmesh is PDM_OWNERSHIP_KEEP : the vertices are deep-copied (each struct holds its own memory)
 *     - if \p vtx_ownership_pmesh is PDM_OWNERSHIP_USER : only the Part Mesh Nodal owns the vertices
 *   - If \p pmn is PDM_OWNERSHIP_USER :
 *     - if \p vtx_ownership_pmesh is PDM_OWNERSHIP_KEEP : only the Part Mesh owns the vertices
 *     - if \p vtx_ownership_pmesh is PDM_OWNERSHIP_USER : neither the Part Mesh Nodal nor the Part Mesh own the vertices
 *
 */
PDM_part_mesh_nodal_to_part_mesh_t *
PDM_part_mesh_nodal_to_part_mesh_create
(
  PDM_part_mesh_nodal_t *pmesh_nodal,
  PDM_bool_t             keep_link_elmt_to_entity,
  PDM_ownership_t        vtx_ownership_pmesh
);


/**
 * \brief Enable the construction of a specific connectivity
 *
 * \param [in] pmn_to_pm          Pointer to \ref PDM_part_mesh_nodal_to_part_mesh_t instance
 * \param [in] connectivity_type  Connectivity to enable
 *
 */
void
PDM_part_mesh_nodal_to_part_mesh_connectivity_enable
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm,
  PDM_connectivity_type_t             connectivity_type
);


/**
 * \brief Enable the generation of global IDs for a specific entity type
 *
 * \param [in] pmn_to_pm    Pointer to \ref PDM_part_mesh_nodal_to_part_mesh_t instance
 * \param [in] entity_type  Entity type
 *
 */
void
PDM_part_mesh_nodal_to_part_mesh_g_nums_enable
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm,
  PDM_mesh_entities_t                 entity_type
);


/**
 * \brief Enable the transfer of groups a specific bound type
 *
 * \param [in] pmn_to_pm   Pointer to \ref PDM_part_mesh_nodal_to_part_mesh_t instance
 * \param [in] bound_type  Bound type
 *
 */
void
PDM_part_mesh_nodal_to_part_mesh_groups_enable
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm,
  PDM_bound_type_t                    bound_type
);


/**
 * \brief Enable the construction of inter-partition communication graph for a specific bound type
 *
 * \param [in] pmn_to_pm   Pointer to \ref PDM_part_mesh_nodal_to_part_mesh_t instance
 * \param [in] bound_type  Bound type
 *
 */
void
PDM_part_mesh_nodal_to_part_mesh_part_comm_graph_enable
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm,
  PDM_bound_type_t                    bound_type
);


/**
 * \brief Generate the Part Mesh
 *
 * \param [in] pmn_to_pm  Pointer to \ref PDM_part_mesh_nodal_to_part_mesh_t instance
 *
 */
void
PDM_part_mesh_nodal_to_part_mesh_compute
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm
);


/**
 * \brief Get the constructed \ref PDM_part_mesh_t instance
 *
 * \param [in]  pmn_to_pm  Pointer to \ref PDM_part_mesh_nodal_to_part_mesh_t instance
 * \param [out] pmesh      Pointer to \ref PDM_part_mesh_t instance
 * \param [in]  ownership  Ownership for \p pmesh
 *
 */
void
PDM_part_mesh_nodal_to_part_mesh_part_mesh_get
(
  PDM_part_mesh_nodal_to_part_mesh_t  *pmn_to_pm,
  PDM_part_mesh_t                    **pmesh,
  PDM_ownership_t                      ownership
);


/**
 * \brief Free a \ref PDM_part_mesh_nodal_to_part_mesh_t instance
 *
 * \param [inout] pmn_to_pm  Pointer to \ref PDM_part_mesh_nodal_to_part_mesh_t instance
 *
 */
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
