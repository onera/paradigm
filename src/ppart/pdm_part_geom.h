/*
 * \file
 */

#ifndef __PDM_PART_GEOM_H__
#define __PDM_PART_GEOM_H__

/*============================================================================
 * Mesh partitioning with geometric methods (SFC)
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_part.h"
#include "pdm_mpi.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Macro and type definitions
 *============================================================================*/

/**
 * \enum PDM_split_geom_t
 * \brief Geometric plit method
 *
 */

typedef enum {
  PDM_PART_GEOM_HILBERT = 1,
  PDM_PART_GEOM_MORTON  = 2,
} PDM_part_geom_t;


/*============================================================================
 * Public function definitions
 *============================================================================*/

void
PDM_dcompute_cell_center
(
  const PDM_MPI_Comm  comm,
  const int           dn_cell,
  const int          *dcell_face_idx,
  const PDM_g_num_t  *dcell_face,
  const int          *dface_vtx_idx,
  const PDM_g_num_t  *dface_vtx,
  const PDM_g_num_t  *dface_proc,
  const double       *dvtx_coord,
  const PDM_g_num_t  *dvtx_proc,
  double             *cell_center
);

/**
 *
 * \brief Perform geometric partitioning
 *
 * \param [in]   method          Geometric method
 * \param [in]   n_part          Number of partition to build on this process
 * \param [in]   comm            Communicator
 * \param [in]   dn_entity       Number of distributed cells
 * \param [in]   dentity_coord   Distributed entity coordinates (size : 3*dn_entity)
 * \param [in]   dcell_weight    Entity weight (size : dn_entity) or NULL
 * \param [inout] dentity_part   Distributed entity partitioning (size = dn_entity)
 */
void
PDM_part_entity_geom
(
 PDM_part_geom_t     method,
 const int           n_part,
 const PDM_MPI_Comm  comm,
 const PDM_g_num_t   dn_entity,
 const double       *dentity_coord,
 const double       *dentity_weight,
       int          *dentity_part
);


void
PDM_part_geom_0d
(
 PDM_part_geom_t     method,
 const int           n_part,
 const PDM_MPI_Comm  comm,
 const int           dn_vtx,
 const double       *dvtx_coord,
 const double       *dvtx_weight,
       int          *dvtx_part
);


void
PDM_part_geom_1d
(
 PDM_part_geom_t     method,
 const int           n_part,
 const PDM_MPI_Comm  comm,
 const int           dn_edge,
 const int           dn_vtx,
 const PDM_g_num_t  *dedge_vtx,
 const double       *dvtx_coord,
 const double       *dedge_weight,
       int          *dedge_part
);


void
PDM_part_geom_2d
(
 PDM_part_geom_t     method,
 const int           n_part,
 const PDM_MPI_Comm  comm,
 const int           dn_face,
 const int           dn_edge,
 const int           dn_vtx,
 const int          *dface_vtx_idx,
 const PDM_g_num_t  *dface_vtx,
 const int          *dface_edge_idx,
 const PDM_g_num_t  *dface_edge,
 const PDM_g_num_t  *dedge_vtx,
 const double       *dvtx_coord,
 const double       *dface_weight,
       int          *dface_part
);

/**
 *
 * \brief Perform geometric partitioning
 *
 * \param [in]   method         Geometric method
 * \param [in]   n_part          Number of partition to build on this process
 * \param [in]   comm           Communicator
 * \param [in]   dn_cell         Number of distributed cells
 * \param [in]   dn_face         Number of distributed faces
 * \param [in]   dn_vtx          Number of distributed vertices
 * \param [in]   dcell_face_idx   Distributed cell face connectivity index or NULL
 *                              (size : dn_cell + 1, numbering : 0 to n-1)
 * \param [in]   dcell_face      Distributed cell face connectivity or NULL
 *                              (size : dface_vtx_idx[dn_cell], numbering : 1 to n)
 * \param [in]   dcell_weight    Cell weight (size : n_cell) or NULL
 * \param [in]   dface_vtx_idx    Distributed face to vertex connectivity index
 *                              (size : dn_face + 1, numbering : 0 to n-1)
 * \param [in]   dface_vtx       Distributed face to vertex connectivity
 *                              (size : dface_vtx_idx[dn_face], numbering : 1 to n)
 * \param [in]   dvtx_coord      Distributed vertex coordinates
 *                              (size : 3*dn_vtx)
 * \param [inout]   dcell_part      Distributed cell partitioning
 *                              (size = dn_cell)
 *
 */

void
PDM_part_geom
(
 PDM_part_geom_t     method,
 const int           n_part,
 const PDM_MPI_Comm  comm,
 const int           dn_cell,
 const int          *dcell_face_idx,
 const PDM_g_num_t  *dcell_face,
 const int          *dcell_weight,
 const int          *dface_vtx_idx,
 const PDM_g_num_t  *dface_vtx,
 const PDM_g_num_t  *dface_proc,
 const double       *dvtx_coord,
 const PDM_g_num_t  *dvtx_proc,
 int                *dcell_part
);


void
PDM_dreorder_from_coords
(
 PDM_part_geom_t  method,
 int              dim,
 PDM_g_num_t     *distrib_vtx,
 double          *dvtx_coord,
 PDM_g_num_t     *vtx_ln_to_gn,
 PDM_MPI_Comm     comm
);

void
PDM_dreorder_from_length
(
 int              dim,
 PDM_g_num_t     *distrib_in,
 double          *length,
 PDM_g_num_t     *ln_to_gn,
 PDM_MPI_Comm     comm
);


/**
 * \brief Compute edge centers
 *
 * \note If \p selected_lnum is NULL, all edges are selected
 *
 * \param [in]   n_part         Number of partitions
 * \param [in]   n_selected     Number of selected edges (size = \p n_part)
 * \param [in]   selected_lnum  Local IDs of selected edges (optional, size = \p n_part or NULL)
 * \param [in]   pedge_vtx      Edge->Vertex connectivity (size = \p n_part, for each part, size = 2 * \p n_selected)
 * \param [in]   pvtx_coord     Vertex coordinates (size = \p n_part)
 * \param [out]  edge_center    Edge centers (size = \p n_part, for each part, size = 3 * \p n_selected)
 */
void
PDM_part_geom_edge_center
(
  int       n_part,
  int      *n_selected,
  int     **selected_lnum,
  int     **pedge_vtx,
  double  **pvtx_coord,
  double ***edge_center
);

/**
 * \brief Compute face centers from face->edge and edge->vtx connectivities
 *
 * \note If \p selected_lnum is NULL, all faces are selected
 *
 * \param [in]   n_part          Number of partitions
 * \param [in]   n_selected      Number of selected faces (size = \p n_part)
 * \param [in]   selected_lnum   Local IDs of selected faces (optional, size = \p n_part or NULL)
 * \param [in]   pface_edge_idx  Index for Face->Edge connectivity (size = \p n_part, for each part, size = \p n_selected + 1)
 * \param [in]   pface_edge      Face->Edge connectivity (size = \p n_part, for each part, size = \p pface_edge_idx[\p n_selected])
 * \param [in]   pedge_vtx       Edge->Vertex connectivity (size = \p n_part)
 * \param [in]   pvtx_coord      Vertex coordinates (size = \p n_part)
 * \param [out]  face_center     Face centers (size = \p n_part, for each part, size = 3 * \p n_selected)
 */
void
PDM_part_geom_face_center_from_edge
(
  int       n_part,
  int      *n_selected,
  int     **selected_lnum,
  int     **pface_edge_idx,
  int     **pface_edge,
  int     **pedge_vtx,
  double  **pvtx_coord,
  double ***face_center
);

/**
 * \brief Compute face centers from face->vtx connectivity
 *
 * \note If \p selected_lnum is NULL, all faces are selected
 *
 * \param [in]   n_part          Number of partitions
 * \param [in]   n_selected      Number of selected faces (size = \p n_part)
 * \param [in]   selected_lnum   Local IDs of selected faces (optional, size = \p n_part or NULL)
 * \param [in]   pface_vtx_idx   Index for Face->Vertex connectivity (size = \p n_part, for each part, size = \p n_selected + 1)
 * \param [in]   pface_vtx       Face->Vertex connectivity (size = \p n_part, for each part, size = \p pface_vtx_idx[\p n_selected])
 * \param [in]   pvtx_coord      Vertex coordinates (size = \p n_part)
 * \param [out]  face_center     Face centers (size = \p n_part, for each part, size = 3 * \p n_selected)
 */
void
PDM_part_geom_face_center_from_vtx
(
  int       n_part,
  int      *n_selected,
  int     **selected_lnum,
  int     **pface_vtx_idx,
  int     **pface_vtx,
  double  **pvtx_coord,
  double ***face_center
);

/**
 * \brief Compute cell centers
 *
 * \note If \p selected_lnum is NULL, all cells are selected
 *
 * \warning Either {\p pface_edge_idx, \p pface_edge, \p pedge_vtx} or {\p pface_vtx_idx, \p pface_vtx} must be provided
 *
 * \param [in]   n_part          Number of partitions
 * \param [in]   n_selected      Number of selected cells (size = \p n_part)
 * \param [in]   selected_lnum   Local IDs of selected cells (optional, size = \p n_part or NULL)
 * \param [in]   pcell_face_idx  Index for Cell->Face connectivity (size = \p n_part, for each part, size = \p n_selected + 1)
 * \param [in]   pcell_face      Cell->Face connectivity (size = \p n_part, for each part, size = \p pcell_face_idx[\p n_selected])
 * \param [in]   pface_edge_idx  Index for Face->Edge connectivity (optional)
 * \param [in]   pface_edge      Face->Edge connectivity (optional)
 * \param [in]   pface_vtx_idx   Index for Face->Vertex connectivity (optional)
 * \param [in]   pface_vtx       Face->Vertex connectivity (optional)
 * \param [in]   pedge_vtx       Edge->Vertex connectivity (optional)
 * \param [in]   pvtx_coord      Vertex coordinates (size = \p n_part)
 * \param [out]  cell_center     Cell centers (size = \p n_part, for each part, size = 3 * \p n_selected)
 */
void
PDM_part_geom_cell_center
(
  int       n_part,
  int      *n_selected,
  int     **selected_lnum,
  int     **pcell_face_idx,
  int     **pcell_face,
  int     **pface_edge_idx,
  int     **pface_edge,
  int     **pface_vtx_idx,
  int     **pface_vtx,
  int     **pedge_vtx,
  double  **pvtx_coord,
  double ***cell_center
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_GEOM_H__ */
