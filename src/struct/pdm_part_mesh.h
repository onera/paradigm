/*
 * \file
 */

#ifndef __PDM_PART_MESH_H__
#define __PDM_PART_MESH_H__

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
#include "pdm_part_comm_graph.h"

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

typedef struct _pdm_part_mesh_t PDM_part_mesh_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief Create a \ref PDM_part_mesh_t structure
 *
 * \param [in]   n_part           Number of partition on the current process
 * \param [in]   comm             MPI communicator
 *
 * \return       Pointer to \ref PDM_part_mesh_t object
 *
 */
PDM_part_mesh_t*
PDM_part_mesh_create
(
 const int             n_part,
       PDM_MPI_Comm    comm
);

/**
 * \brief Free a part mesh structure
 *
 * \param [in]  pmesh          Pointer to \ref PDM_part_mesh_t object
 *
 */
void
PDM_part_mesh_free
(
 PDM_part_mesh_t        *pmesh
);

/**
 * \brief Define partition vertices
 *
 * \param [in]   pmesh               PDM_part_mesh_t structure
 * \param [in]   i_part              Part identifier
 * \param [in]   entity_type         Kind of entity required \ref PDM_mesh_entities_t
 * \param [in]   pn_entity           Number of entity of entity_type
 *
 */
void
PDM_part_mesh_n_entity_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 int                       pn_entity
);

/**
 *
 * \brief Get number of entity on current \ref PDM_part_mesh_t
 *
 * \param [in]   pmesh               PDM_part_mesh_t structure
 * \param [in]   i_part              part identifier
 * \param [in]   PDM_mesh_entities_t Kind of entity required \ref PDM_mesh_entities_t
 *
 * \return Number of entity
 */
int
PDM_part_mesh_n_entity_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type
);

/**
 *
 * \brief Set connectivity
 *
 * \param [in]  pmesh              \ref PDM_part_mesh_t instance
 * \param [in]  i_part             Partition identifier
 * \param [in]  connectivity_type  Type of connectivity
 * \param [in]  connect            Connectivity (1-based)
 * \param [in]  connect_idx        Index for connectivity (can be \p NULL for \p PDM_CONNECTIVITY_TYPE_EDGE_VTX)
 * \param [in]  ownership          Ownership
 *
 */
void
PDM_part_mesh_connectivity_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_connectivity_type_t   connectivity_type,
 int                      *connect,
 int                      *connect_idx,
 PDM_ownership_t           ownership
);

/**
 * \brief Define partition vertices
 *
 * \param [in]  pmesh     Pointer to \ref PDM_part_mesh_t object
 * \param [in]  id_part   Partition identifier
 * \param [in]  coords    Interlaced coordinates (size = 3 * \ref n_vtx)
 * \param [in]  owner     Vertices ownship
 *
 */
void
PDM_part_mesh_vtx_coord_set
(
 PDM_part_mesh_t   *pmesh,
 int                i_part,
 double            *vtx_coord,
 PDM_ownership_t    ownership
);

/**
 *
 * \brief Get coordinates of part mesh vertices
 *
 * \param [in]  pmesh          \ref PDM_part_mesh_t instance
 * \param [in]  i_part         Partition identifier
 * \param [out] vtx_coord      Vertex coordinates (size = 3 * *n_vtx*)
 * \param [in]  ownership      Ownership
 *
 */
void
PDM_part_mesh_vtx_coord_get
(
 PDM_part_mesh_t   *pmesh,
 int                i_part,
 double           **vtx_coord,
 PDM_ownership_t    ownership
);

/**
 *
 * \brief Get part mesh mesh connectivity
 *
 * \param [in]  pmesh              \ref PDM_part_mesh_t instance
 * \param [in]  i_part             Partition identifier
 * \param [in]  connectivity_type  Connectivity type
 * \param [out] connect            Connectivity
 * \param [out] connect_idx        Connectivity index
 * \param [in]  ownership          Ownership
 *
 */
void
PDM_part_mesh_connectivity_get
(
 PDM_part_mesh_t           *pmesh,
 int                        i_part,
 PDM_connectivity_type_t    connectivity_type,
 int                      **connect,
 int                      **connect_idx,
 PDM_ownership_t           ownership
);


/**
 *
 * \brief Set global ids
 *
 * \param [in]  pmesh            \ref PDM_part_mesh_t instance
 * \param [in]  i_part           Partition identifier
 * \param [in]  entity_type      Type of mesh entity
 * \param [in]  pentity_ln_to_gn Global ids
 * \param [in]  ownership        Ownership
 *
 */
void
PDM_part_mesh_entity_ln_to_gn_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 PDM_g_num_t              *pentity_ln_to_gn,
 PDM_ownership_t           ownership
);

/**
 *
 * \brief Get global ids of part mesh entities
 *
 * \param [in]  pmesh            \ref PDM_part_mesh_t instance
 * \param [in]  i_part           Partition identifier
 * \param [in]  entity_type      Entity type
 * \param [out] pentity_ln_to_gn Global ids
 * \param [in]  ownership        Ownership
 *
 */
void
PDM_part_mesh_entity_ln_to_gn_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 PDM_g_num_t             **pentity_ln_to_gn,
 PDM_ownership_t           ownership
);

/**
 *
 * \brief Set color for current entity
 *
 * \param [in]  pmesh            \ref PDM_part_mesh_t instance
 * \param [in]  i_part           Partition identifier
 * \param [in]  entity_type      Type of mesh entity
 * \param [in]  pentity_color    Global ids
 * \param [in]  ownership        Ownership
 *
 */
void
PDM_part_mesh_entity_color_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 int                      *pentity_color,
 PDM_ownership_t           ownership
);

/**
 *
 * \brief Get color for current entity
 *
 * \param [in]  pmesh         \ref PDM_part_mesh_t instance
 * \param [in]  i_part        Partition identifier
 * \param [in]  entity_type   Entity type
 * \param [out] pentity_color Color for current entity
 * \param [in]  ownership     Ownership
 *
 */
void
PDM_part_mesh_entity_color_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 int                     **pentity_color,
 PDM_ownership_t           ownership
);

/**
 * \brief  Set number of group for a current geometry kind
 *
 * \param [in]  pmesh        Pointer to \ref PDM_part_mesh_t object
 * \param [in]  bound_type   Bound type \ref PDM_bound_type_t
 * \param [in]  n_bound      Number of group
 */
void
PDM_part_mesh_n_bound_set
(
 PDM_part_mesh_t          *pmesh,
 PDM_bound_type_t          bound_type,
 int                       n_bound
);

/**
 * \brief  Get number of group for a current geometry kind
 *
 * \param [in]  pmesh        Pointer to \ref PDM_part_mesh_t object
 * \param [in]  bound_type   Bound type \ref PDM_bound_type_t
 *
 * \return Number of group
 */
int
PDM_part_mesh_n_bound_get
(
 PDM_part_mesh_t          *pmesh,
 PDM_bound_type_t          bound_type
);

/**
 * \brief  Get the total number of partition among all process
 *
 * \param [in]  pmesh        Pointer to \ref PDM_part_mesh_t object
 *
 * \return Number of total partition number
 */
int
PDM_part_mesh_tn_part_get
(
 PDM_part_mesh_t          *pmesh
);

/**
 * \brief  Get the number of partition on current process
 *
 * \param [in]  pmesh        Pointer to \ref PDM_part_mesh_t object
 *
 * \return Number of partition on current process
 */
int
PDM_part_mesh_n_part_get
(
 PDM_part_mesh_t          *pmesh
);

/**
 *
 * \brief Set partition group
 *
 * \param [in]   pmesh                 PDM_part_mesh_t
 * \param [in]   i_part                part identifier
 * \param [in]   i_group               group identifier
 * \param [in]   bound_type            Kind of group
 * \param [in]   pn_bound              Number of entity in current group
 * \param [in]   pbound                List of entity in group (size = n_group_entity)
 * \param [in]   pbound_ln_to_gn       List of global identifier in group (size = n_group_entity)
 * \param [in]   ownership             Ownership
 *
 */
void
PDM_part_mesh_bound_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 int                       i_group,
 PDM_bound_type_t          bound_type,
 int                       pn_bound,
 int                      *pbound,
 PDM_g_num_t              *pbound_ln_to_gn,
 PDM_ownership_t           ownership
);

/**
 *
 * \brief Get partition group
 *
 * \param [in]   pmesh                 PDM_part_mesh_t
 * \param [in]   i_part                part identifier
 * \param [in]   i_group               group identifier
 * \param [in]   bound_type            Kind of group
 * \param [out]  pn_bound              Number of entity in current group
 * \param [out]  pbound                List of entity in group (size = n_group_entity)
 * \param [out]  pbound_ln_to_gn       List of global identifier in group (size = n_group_entity)
 * \param [in]   ownership             Ownership
 *
 */
void
PDM_part_mesh_bound_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 int                       i_group,
 PDM_bound_type_t          bound_type,
 int                      *pn_bound,
 int                     **pbound,
 PDM_g_num_t             **pbound_ln_to_gn,
 PDM_ownership_t           ownership
);


/**
 *
 * \brief Get the connection graph between partition for the requested entity type
 * \param [in]  pmesh                 Pointer to \ref PDM_part_mesh_t instance
 * \param [in]  i_part                Id of part
 * \param [in]  bound_type            Type of mesh entity
 * \param [out] ppart_bound_proc_idx  Partitioning boundary entities index from process (size = n_proc + 1)
 * \param [out] ppart_bound_part_idx  Partitioning boundary entities index from partition (size = n_total_part + 1)
 * \param [out] ppart_bound           Partitioning boundary entities (size = 4 * n_entity_part_bound)
 * \param [in]  ownership             Choice of ownership of the resulting arrays \ref PDM_ownership_t
 */
void
PDM_part_mesh_part_graph_comm_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_bound_type_t          bound_type,
 int                      *ppart_bound_proc_idx,
 int                      *ppart_bound_part_idx,
 int                      *ppart_bound,
 PDM_ownership_t           ownership
);

/**
 *
 * \brief Get the connection graph between partition for the requested entity type
 * \param [in]  pmesh                 Pointer to \ref PDM_part_mesh_t instance
 * \param [in]  i_part                Id of part
 * \param [in]  bound_type            Type of mesh entity
 * \param [in]  ppart_bound_proc_idx  Partitioning boundary entities index from process (size = n_proc + 1)
 * \param [in]  ppart_bound_part_idx  Partitioning boundary entities index from partition (size = n_total_part + 1)
 * \param [in]  ppart_bound           Partitioning boundary entities (size = 4 * n_entity_part_bound)
 * \param [in]  ownership             Choice of ownership of the resulting arrays \ref PDM_ownership_t
 */
void
PDM_part_mesh_part_graph_comm_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_bound_type_t          bound_type,
 int                     **ppart_bound_proc_idx,
 int                     **ppart_bound_part_idx,
 int                     **ppart_bound,
 PDM_ownership_t           ownership
);


/**
 *
 * \brief Set part_comm_graph onto part_mesh struct
 *
 * \param [in]  pmesh       Pointer to \ref PDM_part_mesh_t instance
 * \param [in]  pcg         Pointer to \ref PDM_part_comm_graph_t instance
 * \param [in]  entity_type Entity type (vertex, edge, face or cell)
 * \param [in]  ownership   part_mesh ownership on given part_comm_graph
 *
 */
void
PDM_part_mesh_part_comm_graph_set
(
  PDM_part_mesh_t       *pmesh,
  PDM_part_comm_graph_t *pcg,
  PDM_mesh_entities_t    entity_type,
  PDM_ownership_t        ownership
);


/**
 *
 * \brief Get part_mesh's part_comm_graph
 *
 * \param [in]   pmesh       Pointer to \ref PDM_part_mesh_t instance
 * \param [in]   entity_type Entity type (vertex, edge, face or cell)
 * \param [out]  pcg         Pointer to \ref PDM_part_comm_graph_t instance
 * \param [in]   ownership   part_mesh ownership on returned part_comm_graph
 *
 */
void
PDM_part_mesh_part_comm_graph_get
(
  PDM_part_mesh_t        *pmesh,
  PDM_mesh_entities_t     entity_type,
  PDM_part_comm_graph_t **pcg,
  PDM_ownership_t         ownership
);

/**
 * \brief Export a partitioned mesh in Ensight format
 *
 * \param [in] pmesh          Pointer to \ref PDM_part_mesh_t object
 * \param [in] directory      Output directory
 * \param [in] name           Output name
 * \param [in] export_bounds  Option to export bounds
 *
 */
void
PDM_part_mesh_dump_ensight
(
 PDM_part_mesh_t *pmesh,
 const char      *directory,
 const char      *name,
 PDM_bool_t       export_bounds
);

/**
 * \brief Compute the concatenate group information to fall back with old API
 *
 * \param [in] pmesh          Pointer to \ref PDM_part_mesh_t object
 * \param [in] i_part         Partition identifier
 * \param [in]  bound_type            Type of mesh entity
 *
 */
void
PDM_part_mesh_bound_concat_compute
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_bound_type_t          bound_type
);

/**
 *
 * \brief Set partition group
 *
 * \param [in]   pmesh                 PDM_part_mesh_t
 * \param [in]   i_part                part identifier
 * \param [in]   i_group               group identifier
 * \param [in]   bound_type            Kind of group
 * \param [in]   n_bound               Number of group
 * \param [in]   pbound_idx            Index of pbound for all groups (size = n_group+1)
 * \param [in]   pbound                List of entity in group (size = pbound_idx[n_group_entity])
 * \param [in]   pbound_ln_to_gn       List of global identifier in group (size = pbound_idx[n_group_entity])
 * \param [in]   ownership             Ownership
 *
 */
void
PDM_part_mesh_bound_concat_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_bound_type_t          bound_type,
 int                       n_bound,
 int                      *pbound_idx,
 int                      *pbound,
 PDM_g_num_t              *pbound_ln_to_gn,
 PDM_ownership_t           ownership
);

/**
 *
 * \brief Get partition group
 *
 * \param [in]   pmesh                 PDM_part_mesh_t
 * \param [in]   i_part                part identifier
 * \param [in]   i_group               group identifier
 * \param [in]   bound_type            Kind of group
 * \param [out]  pbound_idx            Index of pbound for all groups (size = n_group+1)
 * \param [out]  pbound                List of entity in group (size = pbound_idx[n_group_entity])
 * \param [out]  pbound_ln_to_gn       List of global identifier in group (size = pbound_idx[n_group_entity])
 * \param [in]   ownership             Ownership
 *
 */
void
PDM_part_mesh_bound_concat_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_bound_type_t          bound_type,
 int                     **pbound_idx,
 int                     **pbound,
 PDM_g_num_t             **pbound_ln_to_gn,
 PDM_ownership_t           ownership
);


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_MESH_H__ */
