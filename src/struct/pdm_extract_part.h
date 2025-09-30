/*
 * \file
 */

#ifndef __PDM_EXTRACT_PART_H__
#define __PDM_EXTRACT_PART_H__

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
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_to_part.h"

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

typedef struct _pdm_extract_part_t PDM_extract_part_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief Build an Extract Part instance
 *
 * \param [in]   dim                 Mesh dimension
 * \param [in]   n_part_in           Number of input  partitions
 * \param [in]   n_part_out          Number of output partitions
 * \param [in]   extract_kind        Extraction kind (local/reequilibrate/from target)
 * \param [in]   split_dual_method   Repartitioning method (used only in \ref PDM_EXTRACT_PART_KIND_REEQUILIBRATE mode)
 * \param [in]   compute_child_gnum  Enable generation of new global IDs for extraction
 * \param [in]   ownership           Ownership of extraction
 * \param [in]   comm                MPI communicator
 *
 * \return   Initialized \ref PDM_extract_part_t instance
 */
PDM_extract_part_t*
PDM_extract_part_create
(
  const int                     dim,
  const int                     n_part_in,
  const int                     n_part_out,
        PDM_extract_part_kind_t extract_kind,
        PDM_split_dual_t        split_dual_method,
        PDM_bool_t              compute_child_gnum,
        PDM_ownership_t         ownership,
        PDM_MPI_Comm            comm
);


/**
 * \brief Compute extraction
 *
 * \param [in]   extrp      \ref PDM_extract_part_t instance
 */
void
PDM_extract_part_compute
(
  PDM_extract_part_t        *extrp
);


/**
 * \brief Select local entities to extract.
 * \note Use only in \ref PDM_EXTRACT_PART_KIND_LOCAL or \ref PDM_EXTRACT_PART_KIND_REEQUILIBRATE mode
 *
 * \param [in]   extrp         \ref PDM_extract_part_t instance
 * \param [in]   i_part        Partition identifier
 * \param [in]   n_extract     Number of entities to extract
 * \param [in]   extract_lnum  Local IDs of entities to extract (1-based, size = \p n_extract)
 * \param [in]   ownership     Ownership
 */
void
PDM_extract_part_selected_lnum_set
(
  PDM_extract_part_t       *extrp,
  int                       i_part,
  int                       n_extract,
  int                      *extract_lnum,
  PDM_ownership_t           ownership
);


/**
 * \brief Set the target entities.
 * \note Use only in \ref PDM_EXTRACT_PART_KIND_FROM_TARGET mode.
 *
 * \param [in]   extrp             \ref PDM_extract_part_t instance
 * \param [in]   i_part            Partition identifier
 * \param [in]   n_target          Number of target entities
 * \param [in]   target_gnum       Global IDs of target entities (size = \p n_target)
 * \param [in]   target_location   Initial location of target entities (size = 3 * \p n_target or NULL)
 * \param [in]   ownership         Ownership
 */
void
PDM_extract_part_target_set
(
  PDM_extract_part_t       *extrp,
  int                       i_part,
  int                       n_target,
  PDM_g_num_t              *target_gnum,
  int                      *target_location,
  PDM_ownership_t           ownership
);


/**
 * \brief Set partition
 *
 * \param [in]   extrp             \ref PDM_extract_part_t instance
 * \param [in]   i_part            Partition identifer
 * \param [in]   n_cell            Number of cells
 * \param [in]   n_face            Number of faces
 * \param [in]   n_edge            Number of edges
 * \param [in]   n_vtx             Number of vertices
 * \param [in]   cell_face_idx     Index for cell→face connectivity (size = n_cell+1)
 * \param [in]   cell_face         Cell→face connectivity (size = cell_face_idx[n_cell])
 * \param [in]   face_edge_idx     Index for face→edge connectivity (size = n_face+1)
 * \param [in]   face_edge         Face→edge connectivity (size = face_edge_idx[n_face])
 * \param [in]   edge_vtx          Edge→vtx connectivity (size = 2*n_edge)
 * \param [in]   face_vtx_idx      Index for face→vtx connectivity (size = n_face+1)
 * \param [in]   face_vtx          Face→vtx connectivity (size = face_vtx_idx[n_face])
 * \param [in]   cell_ln_to_gn     Cell global IDs (size = n_cell)
 * \param [in]   face_ln_to_gn     Face global IDs (size = n_face)
 * \param [in]   edge_ln_to_gn     Edge global IDs (size = n_edge)
 * \param [in]   vtx_ln_to_gn      Vertex global IDs (size = n_vtx)
 * \param [in]   vtx_coord         Vertex coordinates (size = 3*n_vtx)
 */
void
PDM_extract_part_part_set
(
  PDM_extract_part_t        *extrp,
  int                       i_part,
  int                       n_cell,
  int                       n_face,
  int                       n_edge,
  int                       n_vtx,
  int                      *cell_face_idx,
  int                      *cell_face,
  int                      *face_edge_idx,
  int                      *face_edge,
  int                      *edge_vtx,
  int                      *face_vtx_idx,
  int                      *face_vtx,
  PDM_g_num_t              *cell_ln_to_gn,
  PDM_g_num_t              *face_ln_to_gn,
  PDM_g_num_t              *edge_ln_to_gn,
  PDM_g_num_t              *vtx_ln_to_gn,
  double                   *vtx_coord
);


/**
 * \brief Set number of groups
 *
 * \param [in]   extrp       \ref PDM_extract_part_t instance
 * \param [in]   bound_type  Kind of group
 * \param [in]   n_group     Number of groups
 */
void
PDM_extract_part_n_group_set
(
  PDM_extract_part_t        *extrp,
  PDM_bound_type_t           bound_type,
  int                        n_group
);


/**
 * \brief Set partition group
 *
 * \param [in]   extrp                  \ref PDM_extract_part_t instance
 * \param [in]   i_part                 Partition identifier
 * \param [in]   i_group                Group identifier
 * \param [in]   bound_type             Kind of group
 * \param [in]   n_group_entity         Number of entities in current group
 * \param [in]   group_entity           Local IDs of entities in group (size = n_group_entity)
 * \param [in]   group_entity_ln_to_gn  Group-specific global IDs of entities in group (size = n_group_entity)
 */
void
PDM_extract_part_part_group_set
(
  PDM_extract_part_t        *extrp,
  int                       i_part,
  int                       i_group,
  PDM_bound_type_t          bound_type,
  int                       n_group_entity,
  int                      *group_entity,
  PDM_g_num_t              *group_entity_ln_to_gn
);


/**
 * \brief Set PDM_part_mesh_nodal_t
 *
 * \param [in]   extrp            \ref PDM_extract_part_t instance
 * \param [in]   pmn              \ref PDM_part_mesh_nodal_t instance
 */
void
PDM_extract_part_part_nodal_set
(
  PDM_extract_part_t    *extrp,
  PDM_part_mesh_nodal_t *pmn
);


/**
 * \brief Set entity centers (useful for Hilbert numbering in \ref PDM_EXTRACT_PART_KIND_REEQUILIBRATE mode)
 *
 * \param [in]   extrp            \ref PDM_extract_part_t instance
 * \param [in]   i_part           Partition identifier
 * \param [in]   entity_center    Coordinates of entity centers
 */
void
PDM_extract_part_entity_center_set
(
  PDM_extract_part_t          *extrp,
  int                          i_part,
  double                      *entity_center
);


/**
 * \brief Get the number of entities of a given type in extraction
 *
 * \param [in]   extrp        \ref PDM_extract_part_t instance
 * \param [in]   i_part_out   Partition identifier
 * \param [in]   entity_type  Type of entity
 *
 * \return Number of entities
 */
int
PDM_extract_part_n_entity_get
(
  PDM_extract_part_t       *extrp,
  int                       i_part_out,
  PDM_mesh_entities_t       entity_type
);


/**
 * \brief Get connectivity in extraction
 *
 * \param [in]  extrp               \ref PDM_extract_part_t instance
 * \param [in]  i_part              Partition identifier
 * \param [in]  connectivity_type   Type of connectivity
 * \param [in]  connect             Connectivity array (size = connect_idx[n_entity])
 * \param [in]  connect_idx         Connectivity index (size = n_entity+1)
 * \param [in]  ownership           Ownership
 *
 * \return Number of leading entities
 */
int
PDM_extract_part_connectivity_get
(
  PDM_extract_part_t        *extrp,
  int                        i_part_out,
  PDM_connectivity_type_t    connectivity_type,
  int                      **connect,
  int                      **connect_idx,
  PDM_ownership_t            ownership
);


/**
 * \brief Get global IDs of entities in extraction
 *
 * \param [in]  extrp             \ref PDM_extract_part_t instance
 * \param [in]  i_part            Partition identifier
 * \param [in]  entity_type       Type of entity
 * \param [out] entity_ln_to_gn   Global IDs
 * \param [in]  ownership         Ownership
 *
 * \return Number of entities
 */
int
PDM_extract_part_ln_to_gn_get
(
  PDM_extract_part_t        *extrp,
  int                        i_part_out,
  PDM_mesh_entities_t        entity_type,
  PDM_g_num_t              **entity_ln_to_gn,
  PDM_ownership_t            ownership
);


/**
 * \brief Get color of entities in extraction (if a renumbering method was used)
 *
 * \param [in]  extrp            \ref PDM_extract_part_t instance
 * \param [in]  i_part           Partition identifier
 * \param [in]  entity_type      Type of entity
 * \param [out] entity_color     Entity color
 * \param [in]  ownership        Ownership
 *
 * \return Number of entities
 */
int
PDM_extract_part_color_get
(
  PDM_extract_part_t   *extrp,
  int                   i_part_out,
  PDM_mesh_entities_t   entity_type,
  int                 **entity_color,
  PDM_ownership_t       ownership
);


/**
 * \brief Get parent global IDs of entities in extraction
 *
 * \param [in]  extrp                   \ref PDM_extract_part_t instance
 * \param [in]  i_part                  Partition identifier
 * \param [in]  entity_type             Type of entity
 * \param [out] parent_entity_ln_to_gn  Parent global IDs
 * \param [in]  ownership               Ownership
 *
 * \return Number of entities
 */
int
PDM_extract_part_parent_ln_to_gn_get
(
  PDM_extract_part_t   *extrp,
  int                   i_part_out,
  PDM_mesh_entities_t   entity_type,
  PDM_g_num_t         **parent_entity_ln_to_gn,
  PDM_ownership_t       ownership
);


/**
 * \brief Get local IDs of parent entities.
 * \note Use only in \ref PDM_EXTRACT_PART_KIND_LOCAL mode
 *
 * \param [in]  extrp               \ref PDM_extract_part_t instance
 * \param [in]  i_part              Partition identifier
 * \param [in]  entity_type         Type of entity
 * \param [out] parent_entity_lnum  Parent local IDs
 * \param [in]  ownership           Ownership
 *
 * \return Number of entities
 */
int
PDM_extract_part_parent_lnum_get
(
  PDM_extract_part_t   *extrp,
  int                   i_part_out,
  PDM_mesh_entities_t   entity_type,
  int                 **parent_entity_lnum,
  PDM_ownership_t       ownership
);


/**
 * \brief Get the initial location of extracted entities
 *
 * \param [in]  extrp          \ref PDM_extract_part_t instance
 * \param [in]  i_part         Partition identifier
 * \param [in]  entity_type    Type of entity
 * \param [out] init_location  Initial location (rank, part, local ID) of extracted entities (size = n_entity * 3)
 * \param [in]  ownership      Ownership
 *
 * \return Number of entities
 */
int
PDM_extract_part_init_location_get
(
 PDM_extract_part_t   *extrp,
 int                   i_part_out,
 PDM_mesh_entities_t   entity_type,
 int                 **init_location,
 PDM_ownership_t       ownership
);


/**
 * \brief Get vertex coordinates in extraction
 *
 * \param [in]   extrp      \ref PDM_extract_part_t instance
 * \param [in]   i_part     Partition identifier
 * \param [out]  vtx_coord  Vertex coordinates (size = 3 * n_vtx)
 * \param [in]   ownership  Ownership
 *
 * \return Number of vertices in extraction
 */
int
PDM_extract_part_vtx_coord_get
(
 PDM_extract_part_t  *extrp,
 int                  i_part_out,
 double             **vtx_coord,
 PDM_ownership_t      ownership
);


/**
 * \brief Retrieve the extracted \ref PDM_part_mesh_nodal_t
 *
 * \param [in]  extrp        \ref PDM_extract_part_t instance
 * \param [out] extract_pmn  \ref PDM_part_mesh_nodal_t instance
 * \param [in]  ownership    Ownership
 */
void
PDM_extract_part_part_mesh_nodal_get
(
  PDM_extract_part_t     *extrp,
  PDM_part_mesh_nodal_t **extract_pmn,
  PDM_ownership_t         ownership
);


/**
 * \brief Free the structure
 *
 * \param [in]   extrp  \ref PDM_extract_part_t instance
 */
void
PDM_extract_part_free
(
  PDM_extract_part_t *extrp
);

/**
 * \brief Free all resulting arrays if not owner
 *
 * \note It is not necessary to call this function before calling \ref PDM_extract_part_free.
 *
 * \param [in]   extrp  \ref PDM_extract_part_t instance
 */
void
PDM_extract_part_partial_free
(
  PDM_extract_part_t  *extrp
);


/**
 * \brief Get the \ref PDM_part_to_part_t instance for a given entity type
 *
 * \note *Direct* exchanges go from extraction to input
 *       and *reverse* exchanges go from input to extraction.
 *
 * \param [in]   extrp        \ref PDM_extract_part_t instance
 * \param [in]   entity_type  Type of entity
 * \param [out]  ptp          \ref PDM_part_to_part_t instance
 * \param [in]   ownership    Ownership
 */
void
PDM_extract_part_part_to_part_get
(
        PDM_extract_part_t   *extrp,
  const PDM_mesh_entities_t   entity_type,
        PDM_part_to_part_t  **ptp,
        PDM_ownership_t       ownership
);


/**
 * \brief Get the Part-to-Part instance for a given group
 *
 * \note *Direct* exchanges go from extraction to input
 *       and *reverse* exchanges go from input to extraction.
 *
 * \param [in]   extrp        \ref PDM_extract_part_t instance
 * \param [in]   bound_type   Type of group
 * \param [in]   i_group      Group identifier
 * \param [out]  ptp          \ref PDM_part_to_part_t instance
 * \param [in]   ownership    Ownership
 */
void
PDM_extract_part_part_to_part_group_get
(
        PDM_extract_part_t   *extrp,
  const PDM_bound_type_t      bound_type,
        int                   i_group,
        PDM_part_to_part_t  **ptp,
        PDM_ownership_t       ownership
);


/**
 * \brief Get partition group
 *
 * \param [in]   extrp                                  \ref PDM_extract_part_t instance
 * \param [in]   bound_type                             Kind of group
 * \param [in]   i_part                                 Partition identifier
 * \param [in]   i_group                                Group identifier
 * \param [out]  pn_extract_group_entity                Number of entities in current group
 * \param [out]  pextract_group_entity                  Local IDs of entities in group (size = \p pn_extract_group_entity)
 * \param [out]  pextract_group_entity_ln_to_gn         Group-specific global IDs (in extraction) of entities in group (size = \p pn_extract_group_entity)
 * \param [out]  pextract_group_entity_parent_ln_to_gn  Group-specific global IDs of entities in group (size = \p pn_extract_group_entity)
 * \param [in]   ownership                              Ownership
 */
void
PDM_extract_part_group_get
(
        PDM_extract_part_t   *extrp,
  const PDM_bound_type_t      bound_type,
        int                   i_part,
        int                   i_group,
        int                  *pn_extract_group_entity,
        int                 **pextract_group_entity,
        PDM_g_num_t         **pextract_group_entity_ln_to_gn,
        PDM_g_num_t         **pextract_group_entity_parent_ln_to_gn,
        PDM_ownership_t       ownership
);


/**
 * \brief Set the reordering method to be used after partitioning
 *
 * \param [in]   extrp                    \ref PDM_extract_part_t instance
 * \param [in]   mesh_entity              Type of entity
 * \param [in]   renum_entity_method      Renumbering method
 * \param [in]   renum_entity_properties  Renumbering parameters for chosen method (can be *NULL*, check out \verbatim embed:rst:inline :ref:`this page<renumbering>` \endverbatim for more details)
 */
void
PDM_extract_part_renum_method_set
(
  PDM_extract_part_t  *extrp,
  PDM_mesh_entities_t  mesh_entity,
  const char          *renum_entity_method,
  const int           *renum_entity_properties
);


/**
 * \brief Get the extracted mesh as a \ref PDM_part_mesh_t instance
 *
 * \param [in]   extrp                   \ref PDM_extract_part_t instance
 * \param [out]  pmesh                   \ref PDM_part_mesh_t object
 * \param [in]   pmesh_takes_ownership   Wether ownerhip is transferred to \p pmesh
 */
void
PDM_extract_part_part_mesh_get
(
  PDM_extract_part_t  *extrp,
  PDM_part_mesh_t    **pmesh,
  PDM_bool_t           pmesh_takes_ownership
);


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_EXTRACT_PART_H__ */
