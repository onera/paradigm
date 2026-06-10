#ifndef __PDM_PART_MESH_GEOM_H__
#define __PDM_PART_MESH_GEOM_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_part_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Compute entity center with downward conectivity
 *
 * \param [in]  n_entity      Number of entities
 * \param [in]  connect_idx   Connectivity index between entity and coordinates (size = n_entity+1)
 * \param [in]  connect       Connectivity array between entity and coordinates (size = connect_idx[n_entity])
 * \param [in]  coord         Vertex coordinates (size = 3 * n_vtx)
 * \param [out] entity_center Entity center computed by mean with vertices coordinnates
 *
 */
void
PDM_compute_entity_center
(
  int     n_entity,
  int    *connect_idx,
  int    *connect,
  double *coord,
  double *entity_center
);

/**
 *
 * \brief Compute dual volume associated at each vertices (but without syncrhonise around partition)
 *
 * \param [in]  n_part          Number of partitions
 * \param [in]  n_face          Number of faces for each partitions (size = n_part)
 * \param [in]  n_edge          Number of edges for each partitions (size = n_part)
 * \param [in]  n_vtx           Vertex coordinates (size = 3 * n_vtx)
 * \param [in]  face_edge_idx   Faces edges connectivity index (size = n_face+1)
 * \param [in]  face_edge       Faces edges connectivity array (size = face_edge_idx[n_face])
 * \param [in]  edge_vtx        Edges vertices connectivity array (size = 2 * n_edge)
 * \param [in]  vtx_coord       Vertex coordinates (size = 3 * n_vtx)
 * \param [out] out_vtx_volume  Entity center computed by mean with vertices coordinnates (size = n_vtx)
 *
 */
void
PDM_compute_dual_volume_ngon_2d
(
  int       n_part,
  int      *n_face,
  int      *n_edge,
  int      *n_vtx,
  int     **face_edge_idx,
  int     **face_edge,
  int     **edge_vtx,
  double  **vtx_coord,
  double ***out_vtx_volume
);

/**
 *
 * \brief Compute dual volume associated at each vertices (but without syncrhonise around partition)
 *
 * \param [in]  n_part          Number of partitions
 * \param [in]  n_cell          Number of cells for each partitions (size = n_part)
 * \param [in]  n_face          Number of faces for each partitions (size = n_part)
 * \param [in]  n_edge          Number of edges for each partitions (size = n_part)
 * \param [in]  n_vtx           Vertex coordinates (size = 3 * n_vtx)
 * \param [in]  cell_face_idx   Cells faces connectivity index (size = n_cell+1)
 * \param [in]  cell_face       Cells faces connectivity array (size = face_edge_idx[n_cell])
 * \param [in]  face_edge_idx   Faces edges connectivity index (size = n_face+1)
 * \param [in]  face_edge       Faces edges connectivity array (size = face_edge_idx[n_face])
 * \param [in]  edge_vtx        Edges vertices connectivity array (size = 2 * n_edge)
 * \param [in]  vtx_coord       Vertex coordinates (size = 3 * n_vtx)
 * \param [out] out_vtx_volume  Entity center computed by mean with vertices coordinnates (size = n_vtx)
 *
 */
void
PDM_compute_dual_volume_ngon_3d
(
  int       n_part,
  int      *n_cell,
  int      *n_face,
  int      *n_edge,
  int      *n_vtx,
  int     **cell_face_idx,
  int     **cell_face,
  int     **face_edge_idx,
  int     **face_edge,
  int     **edge_vtx,
  double  **vtx_coord,
  double ***out_vtx_volume
);

/**
 *
 * \brief Compute dual volumes
 *
 * \param [in]  pm           Pointer to \ref PDM_part_mesh_t instance
 * \param [in]  synchronize  Enable synchronization at partition boundaries
 * \param [in]  dual_vol     For each part, dual volume for each vertex
 *
 */
void
PDM_part_mesh_dual_volume_compute
(
  PDM_part_mesh_t   *pm,
  PDM_bool_t         synchronize,
  double          ***out_dual_vol
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_MESH_GEOM_H__ */
