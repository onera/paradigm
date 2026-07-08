#ifndef __PDM_PART_MESH_NODAL_ALGORITHM_H__
#define __PDM_PART_MESH_NODAL_ALGORITHM_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_mpi.h"
#include "pdm_part_comm_graph.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*============================================================================
 * Types definition
 *============================================================================*/

/**
 *
 * \brief Compute internal part_comm_graph from vertices entity global ids.
 *
 * \param [in] pmn Pointer to \ref PDM_part_mesh_nodal_t instance
 *
 */
void
PDM_part_mesh_nodal_part_comm_graph_vtx_compute_from_gnum
(
  PDM_part_mesh_nodal_t *pmn
);

/**
 *
 * \brief Compute internal part_comm_graph from part_mesh_nodal entity global ids.
 *
 * \param [in]   pmn         Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]   geom_kind   Geometry kind (see \ref PDM_geometry_kind_t )
 *
 */
void
PDM_part_mesh_nodal_part_comm_graph_compute_from_gnum
(
  PDM_part_mesh_nodal_t *pmn,
  PDM_geometry_kind_t    geom_kind
);

/**
 *
 * \brief Compute vertex global ids from vtx part_comm_graph.
 *
 * \param [in] pmn Pointer to \ref PDM_part_mesh_nodal_t instance
 *
 */
void
PDM_part_mesh_nodal_gnum_vtx_compute_from_part_comm_graph
(
  PDM_part_mesh_nodal_t *pmn
);

/**
 *
 * \brief Compute global ids for ``geom_kind`` entities from related part_comm_graph.
 *
 * \param [in] pmn       Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in] geom_kind Geometry kind (see \ref PDM_geometry_kind_t )
 *
 */
void
PDM_part_mesh_nodal_gnum_compute_from_part_comm_graph
(
  PDM_part_mesh_nodal_t *pmn,
  PDM_geometry_kind_t    geom_kind
);

/**
 *
 * \brief Compute part_comm_graph for geom_kind using only connectivity and vertices part_comm_graph
 *
 * \param [in]  pmn         Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  geom_kind   Geometry kind (ridge or surface)
 *
 */
void
PDM_part_mesh_nodal_part_comm_graph_deduce_from_vtx
(
  PDM_part_mesh_nodal_t *pmn,
  PDM_geometry_kind_t    geom_kind
);

/**
 *
 * \brief Complete all internal \ref PDM_part_comm_graph_t inside a \ref PDM_part_mesh_nodal_t structure.
 * If vtx_ln_to_gn is only specified, compute the part_comm_graph for vertices, then deduces all other possible connectivity.
 * If only \ref PDM_part_comm_graph_t vertices is inside the \ref PDM_part_mesh_nodal we directly compute all other possible connectivity.
 * If not global ids or communication graph is provided, we raise an error.
 * This call can be is collective
 *
 * \param [in]   pmn         Pointer to \ref PDM_part_mesh_nodal_t instance
 *
 */
void
PDM_part_mesh_nodal_complete_part_comm_graph
(
  PDM_part_mesh_nodal_t *pmn
);

/**
 *
 * \brief Compute entities straddling at least two different groups from geom_kind entities.
 * Generated entities are stored in pmn structure.
 *
 * \param [in]  pmn           Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  geom_kind     Geometry kind (ridge or surface)
 * \param [in]  geom_kind_tgt Geometry kind of wanted generated entities (ridge or corner)
 */
void
PDM_part_mesh_nodal_compute_straddling_entities
(
  PDM_part_mesh_nodal_t  *pmn,
  PDM_geometry_kind_t     geom_kind,
  PDM_geometry_kind_t     geom_kind_tgt
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_MESH_NODAL_ALGORITHM_H__ */
