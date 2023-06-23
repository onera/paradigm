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

void
PDM_part_connectivity_to_connectity_idx
(
const int    n_part,
const int   *n_entity1,
      int  **entity1_entity2_in,
      int ***entity1_entity2_idx,
      int ***entity1_entity2
);

void
PDM_compute_face_vtx_from_face_and_edge
(
 int   n_face,
 int  *face_edge_idx,
 int  *face_edge,
 int  *edge_vtx,
 int **face_vtx
);

void
PDM_compute_face_vtx_from_face_and_edge_unsigned
(
 int   n_face,
 int  *face_edge_idx,
 int  *face_edge,
 int  *edge_vtx,
 int **face_vtx
);


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

#ifdef  __cplusplus
}
#endif

#endif  /* PDM_PART_CONNECTIVITY_TRANSFORM_H_ */
