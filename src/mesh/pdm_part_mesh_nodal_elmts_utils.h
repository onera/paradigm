#ifndef __PDM_PART_MESH_NODAL_ELMTS_UTILS_H__
#define __PDM_PART_MESH_NODAL_ELMTS_UTILS_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_io.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

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
PDM_part_mesh_nodal_elmts_sections_decompose_faces
(
  PDM_part_mesh_nodal_elmts_t *pmne,
  int                         *elmt_face_vtx_idx,
  PDM_g_num_t                 *elmt_face_vtx,
  PDM_g_num_t                 *elmt_face_cell,
  int                         *elmt_cell_face_idx,
  PDM_g_num_t                 *elmt_cell_face,
  int                         *parent_elmt_position
);



void
PDM_part_mesh_nodal_elmts_decompose_faces_get_size
(
 PDM_part_mesh_nodal_elmts_t *pmne,
 int                         *n_face_elt_tot,
 int                         *n_sum_vtx_face_tot
);



#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_MESH_NODAL_ELMTS_UTILS_H__ */
