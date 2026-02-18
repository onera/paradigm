/*
 * \file
 */

#ifndef PDM_PART_CONNECTIVITY_TRANSFORM_H_
#define PDM_PART_CONNECTIVITY_TRANSFORM_H_

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/

#ifdef  __cplusplus
extern "C" {
#endif

/*============================================================================
 * Types definition
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/**
 *
 * \brief Compress connectivity by removing duplicate entries
 *
 * \p entity1_entity2 can be realloc'd after calling this function
 *
 * \param [in]     n_entity1           Number of entity1
 * \param [inout]  entity1_entity2_idx Connectivity index between entity1 and entity2 (size = \p n_entity1 + 1)
 * \param [inout]  entity1_entity2     Connectivity between entity1 and entity2 (size = \p entity1_entity2_idx[\p n_entity1])
 *
 */
void
PDM_compress_connectivity
(
  int  n_entity1,
  int *entity1_entity2_idx,
  int *entity1_entity2
);

/**
 *
 * \brief Combine \p entity1_entity2 and \p entity2_entity3 connectivities to get \p entity1_entity3
 *
 * \param [in]  n_entity1           Number of entity1
 * \param [in]  entity1_entity2_idx Connectivity index between entity1 and entity2 (size = \p n_entity1 + 1)
 * \param [in]  entity1_entity2     Connectivity between entity1 and entity2 (size = \p entity1_entity2_idx[\p n_entity1])
 * \param [in]  entity2_entity3_idx Connectivity index between entity2 and entity3 (size = \p n_entity2 + 1)
 * \param [in]  entity2_entity3     Connectivity between entity2 and entity3 (size = \p entity2_entity3_idx[\p n_entity2] )
 * \param [out] entity1_entity3_idx Connectivity index between entity1 and entity3 (size = \p n_entity1 + 1)
 * \param [out] entity1_entity3     Connectivity between entity1 and entity3 (size = \p entity1_entity3_idx[\p n_entity1])
 *
 */
void
PDM_combine_connectivity
(
  int   n_entity1,
  int  *entity1_entity2_idx,
  int  *entity1_entity2,
  int  *entity2_entity3_idx,
  int  *entity2_entity3,
  int **entity1_entity3_idx,
  int **entity1_entity3
);


/**
 *
 * \brief Transpose connectivity \p entity1_entity2 to get \p entity2_entity1
 *
 * \param [in]  n_entity1           Number of entity1
 * \param [in]  n_entity1           Number of entity2
 * \param [in]  entity1_entity2_idx Connectivity index between entity1 and entity2 (size = \p n_entity1 + 1)
 * \param [in]  entity1_entity2     Connectivity between entity1 and entity2 (size = \p entity1_entity2_idx[\p n_entity1])
 * \param [out] entity2_entity1_idx Connectivity index between entity2 and entity1 (size = \p n_entity2 + 1)
 * \param [out] entity2_entity1     Connectivity between entity2 and entity1 (size = \p entity2_entity1_idx[\p n_entity2] )
 *
 */
void
PDM_connectivity_transpose
(
  const int   n_entity1,
  const int   n_entity2,
        int  *entity1_entity2_idx,
        int  *entity1_entity2,
        int **entity2_entity1_idx,
        int **entity2_entity1
);

/**
 *
 * \brief Transpose connectivity \p entity1_entity2 to get \p entity2_entity1 for multiple partitions
 *
 * \param [in]  n_part              Number of partitions in current process
 * \param [in]  n_entity1           Number of entity1
 * \param [in]  n_entity1           Number of entity2
 * \param [in]  entity1_entity2_idx Connectivity index between entity1 and entity2 (size = \p n_entity1 + 1)
 * \param [in]  entity1_entity2     Connectivity between entity1 and entity2 (size = \p entity1_entity2_idx[\p n_entity1])
 * \param [out] entity2_entity1_idx Connectivity index between entity2 and entity1 (size = \p n_entity2 + 1)
 * \param [out] entity2_entity1     Connectivity between entity2 and entity1 (size = \p entity2_entity1_idx[\p n_entity2])
 *
 */
void
PDM_part_connectivity_transpose
(
  const int    n_part,
  const int   *n_entity1,
  const int   *n_entity2,
        int  **entity1_entity2_idx,
        int  **entity1_entity2,
        int ***entity2_entity1_idx,
        int ***entity2_entity1
);

/**
 *
 * \brief Combine connectivity between \p entity1_entity2 and \p entity2_entity3 to get \p entity1_entity3 for multiple partitions
 *
 * \param [in]  n_part              Number of partitions in current process
 * \param [in]  n_entity1           Number of entity1
 * \param [in]  entity1_entity2_idx Connectivity index between entity1 and entity2 (size = \p n_entity1 + 1)
 * \param [in]  entity1_entity2     Connectivity between entity1 and entity2 (size = \p entity1_entity2_idx[\p n_entity1])
 * \param [in]  entity2_entity3_idx Connectivity index between entity2 and entity3 (size = \p n_entity2 + 1)
 * \param [in]  entity2_entity3     Connectivity between entity2 and entity3 (size = \p entity2_entity3_idx[\p n_entity2])
 * \param [out] entity1_entity3_idx Connectivity index between entity1 and entity3 (size = \p n_entity1 + 1)
 * \param [out] entity1_entity3     Connectivity between entity1 and entity3 (size = \p entity1_entity3_idx[\p n_entity1])
 *
 */
void
PDM_part_combine_connectivity
(
  const int    n_part,
  int         *n_entity1,
  int        **entity1_entity2_idx,
  int        **entity1_entity2,
  int        **entity2_entity3_idx,
  int        **entity2_entity3,
  int       ***entity1_entity3_idx,
  int       ***entity1_entity3
);


/**
 *
 * \brief Build dual graph \p entity1_entity1 (and \p entity2_entity1) from \p entity1_entity1 connectivity (multiple partitions)
 *
 * \param [in]  n_part              Number of partitions in current process
 * \param [in]  n_entity1           Number of entity1
 * \param [in]  n_entity1           Number of entity2
 * \param [in]  entity1_entity2_idx Connectivity index between entity1 and entity2 (size = \p n_entity1 + 1)
 * \param [in]  entity1_entity2     Connectivity between entity1 and entity2 (size = \p entity1_entity2_idx[\p n_entity1] )
 * \param [in]  entity2_entity1_idx Connectivity index between entity2 and entity1 (size = \p n_entity2 + 1)
 * \param [in]  entity2_entity1     Connectivity between entity2 and entity1 (size = \p entity1_entity2_idx[\p n_entity2] )
 * \param [out] entity1_entity1_idx Connectivity index between entity1 and entity1 (size = \p n_entity1 + 1)
 * \param [out] entity1_entity1     Connectivity between entity1 and entity1 (size = \p entity1_entity1_idx[\p n_entity1] )
 *
 */
void
PDM_part_graph_dual
(
  int    n_part,
  int   *n_entity1,
  int   *n_entity2,
  int  **entity1_entity2_idx,
  int  **entity1_entity2,
  int ***entity2_entity1_idx,
  int ***entity2_entity1,
  int ***entity1_entity1_idx,
  int ***entity1_entity1
);


/**
 *
 * \brief Convert implicit pair connectivity, to a connectivity with index. Useful for converting face_cell or edge_vtx.
 *
 * \param [in]  n_part              Number of partitions in current process
 * \param [in]  n_entity1           Number of entity1
 * \param [in]  entity1_entity2_in  Implicit connectivity (face_cell for example with right cell is boundary face[2*i+1] == 0)
 * \param [out] entity1_entity2_idx Connectivity index between entity1 and entity2 (size = \p n_entity1 + 1)
 * \param [out] entity1_entity2     Connectivity between entity1 and entity2 (size = \p entity1_entity2_idx[\p n_entity1])
 *
 */
void
PDM_part_connectivity_to_connectivity_idx
(
  const int    n_part,
  const int   *n_entity1,
        int  **entity1_entity2_in,
        int ***entity1_entity2_idx,
        int ***entity1_entity2
);

/**
 *
 * \brief Generate face->vtx connectivity from (*signed*) face->edge and edge->vtx connectivities
 *
 * \param [in]  n_face        Number of faces
 * \param [in]  face_edge_idx Index for face->edge connectivity (size = \p n_face + 1)
 * \param [in]  face_edge     Face->edge connectivity (1-based, signed) (size = \p face_edge_idx[\p n_face])
 * \param [in]  edge_vtx      Edge->vertex connectivity (size = 2 * n_edge)
 * \param [out] face_vtx      Face->vertex (size = \p face_edge_idx[\p n_face])
 *
 */
void
PDM_compute_face_vtx_from_face_and_edge
(
  int   n_face,
  int  *face_edge_idx,
  int  *face_edge,
  int  *edge_vtx,
  int **face_vtx
);

/**
 *
 * \brief Generate face->vtx connectivity from (*unsigned*) face->edge and edge->vtx connectivities
 *
 * \param [in]  n_face        Number of faces
 * \param [in]  face_edge_idx Index for face->edge connectivity (size = \p n_face + 1)
 * \param [in]  face_edge     Face->edge connectivity (1-based, unsigned) (size = \p face_edge_idx[\p n_face])
 * \param [in]  edge_vtx      Edge->vertex connectivity (size = 2 * n_edge)
 * \param [out] face_vtx      Face->vertex (size = \p face_edge_idx[\p n_face])
 *
 */
void
PDM_compute_face_vtx_from_face_and_edge_unsigned
(
  int   n_face,
  int  *face_edge_idx,
  int  *face_edge,
  int  *edge_vtx,
  int **face_vtx
);


/**
 *
 * \brief Generate face->vtx with face->edge and edge->vtx using block-distributions (**Block-distributed**)
 *
 * \param [in]  comm            PDM_MPI communicator
 * \param [in]  distrib_face    Distribution of faces among process (size = n_rank+1)
 * \param [in]  distrib_edge    Distribution of faces among process (size = n_rank+1)
 * \param [in]  dface_edge_idx  Connectivity index between face and edge (size = dn_face + 1)
 * \param [in]  dface_edge      Connectivity between face and edge (size = \p dface_edge_idx[dn_face])
 * \param [out] dface_vtx       Connectivity between face and vtx (size = \p dface_edge_idx[dn_face])
 *
 */
void
PDM_compute_dface_vtx_from_edges_distrib
(
  PDM_MPI_Comm   comm,
  PDM_g_num_t   *distrib_face,
  PDM_g_num_t   *distrib_edge,
  int           *dface_edge_idx,
  PDM_g_num_t   *dface_edge,
  PDM_g_num_t   *dedge_vtx,
  PDM_g_num_t  **dface_vtx
);

/**
 *
 * \brief Generate face->vtx connectivity from face->edge and edge->vtx connectivities (**Block-distributed**)
 *
 * \param [in]  comm            PDM_MPI communicator
 * \param [in]  dn_face         Number of faces
 * \param [in]  dn_edge         Number of edges
 * \param [in]  dface_edge_idx  Index of Face->edge connectivity (size = \p dn_face + 1)
 * \param [in]  dface_edge      Face->edge connectivity (size = \p dface_edge_idx[\p dn_face])
 * \param [out] dface_vtx       Face->vertex connectivity (size = \p dface_edge_idx[\p dn_face])
 *
 */
void
PDM_compute_dface_vtx_from_edges
(
  PDM_MPI_Comm   comm,
  int            dn_face,
  int            dn_edge,
  int           *dface_edge_idx,
  PDM_g_num_t   *dface_edge,
  PDM_g_num_t   *dedge_vtx,
  PDM_g_num_t  **dface_vtx
);

/**
 *
 * \brief Sort and unique on a graph (obtained for example with a combination of part_connectivity_transpose / combine) to get vtx_vtx graph.
 *        The diagonal part is removed.
 *        This method is useful to prepare a graph in order to split it with METIS or Scotch.
 *        We can reallocate \p graph after the call.
 *        \p graph is 1-based in input and 0-based in output.
 *
 * \param [in]     n_entity        Number of entity is the current graph
 * \param [inout]  graph_idx       Index of graph (size = \p n_entity + 1)
 * \param [inout]  graph           Graph
 *
 */
void
PDM_graph_compress
(
  int  n_entity,
  int *graph_idx,
  int *graph
);

/**
 *
 * \brief Filter a connectivity \p entity1_to_entity2 based on flags on both entity types
 *
 * This function creates a sub-connectivity by keeping only entity1 having \p entity1_flag set to 1,
 * and within those, keeping only connections to entity2 having \p entity2_flag set to 1.
 * Arrays \p sub_entity1_to_entity2_idx and \p sub_entity1_to_entity2 are allocated
 * within the function and must be freed by the user.
 *
 * \param [in]  n_entity1                  Number of initial entity1
 * \param [in]  entity1_to_entity2_idx     Initial connectivity index (size = \p n_entity1 + 1)
 * \param [in]  entity1_to_entity2         Initial connectivity (size = \p entity1_to_entity2_idx[\p n_entity1])
 * \param [in]  entity1_flag               Filter flag for entity1 (1 to keep, 0 to skip, size = \p n_entity1)
 * \param [in]  entity2_flag               Filter flag for entity2 (1 to keep, 0 to skip, size = max(\p entity1_to_entity2))
 * \param [out] sub_entity1_to_entity2_idx Filtered connectivity index (size = return_value + 1)
 * \param [out] sub_entity1_to_entity2     Filtered connectivity (size = (*\p sub_entity1_to_entity2_idx)[return_value])
 *
 * \return Number of entity1 kept in the filtered connectivity
 */
int
PDM_connectivity_filter
(
  int   n_entity1,
  int  *entity1_to_entity2_idx,
  int  *entity1_to_entity2,
  int  *entity1_flag,
  int  *entity2_flag,
  int **sub_entity1_to_entity2_idx,
  int **sub_entity1_to_entity2
);

#ifdef  __cplusplus
}
#endif

#endif  /* PDM_PART_CONNECTIVITY_TRANSFORM_H_ */
