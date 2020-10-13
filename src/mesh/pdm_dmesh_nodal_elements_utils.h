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

/**
*
* \brief PDM_section_size_elt_faces_get
*
* \param [in]     mesh               Current mesh
* \param [in]     id_section         Section identifier
* \param [inout]  elt_face_vtx_idx   Index of element faces connectivity (preallocated)
* \param [inout]  elt_face_vtx       Element faces connectivity (preallocated)
*
*/
int
PDM_section_size_elt_faces_get
(
  PDM_DMesh_nodal_t *mesh,
  int               *s_elt_face_vtx_idx,
  int               *s_elt_face_vtx,
  int               *s_elt_face_cell
);


/**
*
* \brief PDM_section_size_elt_edges_get
*
* \param [in]     mesh               Current mesh
* \param [in]     id_section         Section identifier
* \param [inout]  elt_edge_vtx_idx   Index of element faces connectivity (preallocated)
* \param [inout]  elt_edge_vtx       Element faces connectivity (preallocated)
*
*/
int
PDM_section_size_elt_edges_get
(
  PDM_DMesh_nodal_t *mesh,
  int               *s_elt_edge_vtx_idx,
  int               *s_elt_edge_vtx,
  int               *s_elt_edge_cell
);


/**
 * \brief Return for standard elements the number of face that build this element
 *
 */
int
PDM_n_face_elt_per_elmt
(
  PDM_Mesh_nodal_elt_t t_elt
);

/**
 * \brief Return for standard elements the number of edge that build this element
 *
 */
int
PDM_n_nedge_elt_per_elmt
(
  PDM_Mesh_nodal_elt_t t_elt
);

/**
 * \brief Return for standard elements the total number of face vtx connectivity that build this element
 *
 */
int
PDM_n_sum_vtx_face_per_elmt
(
  PDM_Mesh_nodal_elt_t t_elt
);


/**
 * \brief Return for standard elements the total number of edge vtx connectivity that build this element
 *
 */
int
PDM_n_sum_vtx_edge_per_elmt
(
  PDM_Mesh_nodal_elt_t t_elt
);


void
PDM_sections_decompose_faces
(
  PDM_DMesh_nodal_t *mesh,
  int               *elmt_face_vtx_idx,
  PDM_g_num_t       *elmt_face_vtx,
  PDM_g_num_t       *elmt_face_cell,
  int               *elmt_cell_face_idx,
  PDM_g_num_t       *elmt_cell_face
);

void
PDM_sections_decompose_edges
(
  PDM_DMesh_nodal_t *mesh,
  int               *elmt_edge_vtx_idx,
  PDM_g_num_t       *elmt_edge_vtx,
  PDM_g_num_t       *elmt_edge_cell,
  int               *elmt_cell_edge_int,
  PDM_g_num_t       *elmt_cell_edge
);

/**
 *
 * \brief Decompose tetra cell_vtx connectivity to a flatten view of faces
 */
void
PDM_tetra_decomposes_faces
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_face_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_face_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
       int         *elmt_face_vtx_idx,
       PDM_g_num_t *elmt_face_vtx,
       PDM_g_num_t *elmt_face_cell,
       int         *elmt_cell_face_idx,
       PDM_g_num_t *elmt_cell_face
);

/**
 *
 * \brief Decompose pyra cell_vtx connectivity to a flatten view of faces
 */
void
PDM_pyra_decomposes_faces
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_face_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_face_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
       int         *elmt_face_vtx_idx,
       PDM_g_num_t *elmt_face_vtx,
       PDM_g_num_t *elmt_face_cell,
       int         *elmt_cell_face_idx,
       PDM_g_num_t *elmt_cell_face
);

/**
 *
 * \brief Decompose prism cell_vtx connectivity to a flatten view of faces
 */
void
PDM_prism_decomposes_faces
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_face_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_face_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
       int         *elmt_face_vtx_idx,
       PDM_g_num_t *elmt_face_vtx,
       PDM_g_num_t *elmt_face_cell,
       int         *elmt_cell_face_idx,
       PDM_g_num_t *elmt_cell_face
);

/**
 *
 * \brief Decompose hexa cell_vtx connectivity to a flatten view of faces
 */
void
PDM_hexa_decomposes_faces
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_face_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_face_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
       int         *elmt_face_vtx_idx,
       PDM_g_num_t *elmt_face_vtx,
       PDM_g_num_t *elmt_face_cell,
       int         *elmt_cell_face_idx,
       PDM_g_num_t *elmt_cell_face
);


/**
*
* \brief Decompose quad cell_vtx connectivity to a flatten view of edges
*/
void
PDM_tri_decomposes_edges
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_edge_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_face_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
       int         *elmt_edge_vtx_idx,
       PDM_g_num_t *elmt_edge_vtx,
       PDM_g_num_t *elmt_edge_cell,
       int         *elmt_cell_edge_idx,
       PDM_g_num_t *elmt_cell_edge
);

/**
*
* \brief Decompose quad cell_vtx connectivity to a flatten view of edges
*/
void
PDM_quad_decomposes_edges
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_edge_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_face_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
       int         *elmt_edge_vtx_idx,
       PDM_g_num_t *elmt_edge_vtx,
       PDM_g_num_t *elmt_edge_cell,
       int         *elmt_cell_edge_idx,
       PDM_g_num_t *elmt_cell_edge
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MESH_NODAL_H__ */
