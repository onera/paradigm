#ifndef __PDM_DMESH_NODAL_ELEMENTS_UTILS_H__
#define __PDM_DMESH_NODAL_ELEMENTS_UTILS_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


void
PDM_hexa_section_decompose_elemt_to_face
(
      PDM_g_num_t  n_elmt,
const PDM_g_num_t *elmt_vtx,
      int         *elmt_face_vtx_idx,
      PDM_g_num_t *elmt_face_vtx,
      PDM_g_num_t *elmt_face_cell,
      PDM_g_num_t *elmt_cell_face
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MESH_NODAL_H__ */
