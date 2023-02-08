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
PDM_part_mesh_nodal_tetra_decomposes_faces
(
       int          n_elt,
       int          order,
       int         *parent_node,
       int         *n_elt_current,
       int         *n_face_current,
 const PDM_g_num_t *vtx_ln_to_gn,
 const int         *connectivity_elmt_vtx,
       int         *elmt_face_vtx_idx,
       PDM_g_num_t *elmt_face_vtx,
       int         *elmt_cell_face_idx,
       int         *parent_elmt_position
);


void
PDM_part_mesh_nodal_std_decomposes_faces
(
       PDM_Mesh_nodal_elt_t  t_elt,
       int                   n_elt,
       int                   order,
       int                  *parent_node,
       int                  *n_elt_current,
       int                  *n_face_current,
 const PDM_g_num_t          *vtx_ln_to_gn,
 const int                  *connectivity_elmt_vtx,
       int                  *elmt_face_vtx_idx,
       PDM_g_num_t          *elmt_face_vtx,
       int                  *elmt_cell_face_idx,
       int                  *parent_elmt_position
);

void
PDM_part_mesh_nodal_elmts_sections_decompose_faces
(
  PDM_part_mesh_nodal_elmts_t  *pmne,
  PDM_g_num_t                 **vtx_ln_to_gn,
  int                          *elmt_face_vtx_idx,
  PDM_g_num_t                  *elmt_face_vtx,
  int                          *elmt_cell_face_idx,
  int                          *parent_elmt_position
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
