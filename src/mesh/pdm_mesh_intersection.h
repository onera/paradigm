#ifndef PDM_MESH_INTERSECTION_H
#define PDM_MESH_INTERSECTION_H

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_overlay.h"
#include "pdm_part_to_part.h"

/*----------------------------------------------------------------------------*/

#ifdef  __cplusplus
extern "C" {
#if 0
} /* Fake brace */
#endif
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _pdm_mesh_intersection_t PDM_mesh_intersection_t;


/*============================================================================
 * Public function definitions
 *============================================================================*/

PDM_mesh_intersection_t*
PDM_mesh_intersection_create
(
 const PDM_mesh_intersection_kind_t intersection_kind,
 const int                          dim_mesh_a,
 const int                          dim_mesh_b,
 const int                          n_part_mesh_a,
 const int                          n_part_mesh_b,
 const double                       project_coeff,
       PDM_MPI_Comm                 comm
);

void
PDM_mesh_intersection_compute
(
  PDM_mesh_intersection_t  *mi
);

void
PDM_mesh_intersection_part_set
(
  PDM_mesh_intersection_t  *mi,
  PDM_ol_mesh_t             i_mesh,
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

void
PDM_mesh_intersection_tetraisation_pt_set
(
 PDM_mesh_intersection_t* mi,
 int     tetraisation_pt_type,
 double *tetraisation_pt_coord
);

void
PDM_mesh_intersection_stat_get
(
 PDM_mesh_intersection_t* mi,
 double *local_vol_A_B,
 double *global_vol_A_B,
 double *global_vol_A

);

void
PDM_mesh_intersection_free
(
 PDM_mesh_intersection_t* mi
);


/**
 * \brief Get part_to_part object to exchange data between the intersected meshes
 *
 * \param [in ] mi         Pointer to \ref PDM_mesh_intersection_t object
 * \param [out] ptp        Pointer to \ref PDM_part_to_part_t object
 * \param [in ] ownership  Ownership for ptp
 *
 */

void
PDM_mesh_intersection_part_to_part_get
(
 PDM_mesh_intersection_t  *mi,
 PDM_part_to_part_t      **ptp,
 PDM_ownership_t           ownership
 );

#ifdef  __cplusplus
}
#endif

#endif  /* PDM_MESH_INTERSECTION_H */
