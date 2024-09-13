/*
 * \file
 */

#ifndef __PDM_DMESH_NODAL_ELEMENTS_UTILS_H__
#define __PDM_DMESH_NODAL_ELEMENTS_UTILS_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include "pdm_dmesh_nodal.h"
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
 * \brief Return for standard elements the number of face that build this element
 *
 */
int
PDM_n_face_elt_per_elmt
(
  PDM_Mesh_nodal_elt_t t_elt
);

/**
 * \brief Return face->vtx connectivity of a standard elements decomposed into faces
 *
 * \param [in]  t_elt         Standard element type
 * \param [out] face_vtx_idx  Index for face->vertex connectivity
 * \param [out] face_vtx_idx  Face->vertex connectivity
 *
 * \return Number of faces per element
 */
int
PDM_face_vtx_per_elmt
(
  PDM_Mesh_nodal_elt_t   t_elt,
  const int            **face_vtx_idx,
  const int            **face_vtx
);

/**
 * \brief Return for standard elements the number of edge that build this element
 *
 */
int
PDM_n_edge_elt_per_elmt
(
  PDM_Mesh_nodal_elt_t t_elt
);

/**
 * \brief Return edge->vtx connectivity of a standard elements decomposed into edges
 *
 * \param [in]  t_elt         Standard element type
 * \param [out] edge_vtx_idx  Edge->vertex connectivity
 *
 * \return Number of edges per element
 */
int
PDM_edge_vtx_per_elmt
(
  PDM_Mesh_nodal_elt_t   t_elt,
  const int            **edge_vtx
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
PDM_poly2d_decomposes_faces
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_face_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_face_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
 const int         *connectivity_elmt_vtx_idx,
       int         *elmt_face_vtx_idx,
       PDM_g_num_t *elmt_face_vtx,
       PDM_g_num_t *elmt_face_cell,
       int         *elmt_cell_face_idx,
       PDM_g_num_t *elmt_cell_face,
       int         *parent_elmt_position
);

void
PDM_poly2d_decomposes_edges
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_edge_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_edge_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
 const int         *connectivity_elmt_vtx_idx,
       int         *elmt_edge_vtx_idx,
       PDM_g_num_t *elmt_edge_vtx,
       PDM_g_num_t *elmt_edge_cell,
       int         *elmt_cell_edge_idx,
       PDM_g_num_t *elmt_cell_edge,
       int         *parent_elmt_position
);

void
PDM_poly3d_decomposes_faces
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_face_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_face_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
 const PDM_g_num_t *connectivity_elmt_vtx_idx,
       int         *elmt_face_vtx_idx,
       PDM_g_num_t *elmt_face_vtx,
       PDM_g_num_t *elmt_face_cell,
       int         *elmt_cell_face_idx,
       PDM_g_num_t *elmt_cell_face,
       int         *parent_elmt_position
);

/**
*
* \brief Decompose quad cell_vtx connectivity to a flatten view of edges
*/
void
PDM_poly3d_decomposes_edges
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_edge_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_edge_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
 const PDM_g_num_t *connectivity_elmt_vtx_idx,
       int         *elmt_edge_vtx_idx,
       PDM_g_num_t *elmt_edge_vtx,
       PDM_g_num_t *elmt_edge_cell,
       int         *elmt_cell_edge_idx,
       PDM_g_num_t *elmt_cell_edge,
       int         *parent_elmt_position
);


void
PDM_sections_decompose_faces
(
  PDM_dmesh_nodal_elmts_t *dmn_elts,
  int                     *elmt_face_vtx_idx,
  PDM_g_num_t             *elmt_face_vtx,
  PDM_g_num_t             *elmt_face_cell,
  int                     *elmt_cell_face_idx,
  PDM_g_num_t             *elmt_cell_face,
  int                     *parent_elmt_position
);

void
PDM_sections_decompose_edges
(
  PDM_dmesh_nodal_elmts_t *dmn_elts,
  int                     *elmt_edge_vtx_idx,
  PDM_g_num_t             *elmt_edge_vtx,
  PDM_g_num_t             *elmt_edge_cell,
  int                     *elmt_cell_edge_idx,
  PDM_g_num_t             *elmt_cell_edge,
  int                     *parent_elmt_position
);

void
PDM_std_decomposes_faces
(
       PDM_Mesh_nodal_elt_t  t_elt,
       int                   order,
       int                  *parent_node,
       int                   n_elt,
       int                  *n_elt_current,
       int                  *n_dface_current,
       PDM_g_num_t           beg_gnum_elt_current,
       PDM_g_num_t           beg_gnum_face_current,
 const PDM_g_num_t          *connectivity_elmt_vtx,
       int                  *elmt_face_vtx_idx,
       PDM_g_num_t          *elmt_face_vtx,
       PDM_g_num_t          *elmt_face_cell,
       int                  *elmt_cell_face_idx,
       PDM_g_num_t          *elmt_cell_face,
       int                  *parent_elmt_position
);

void
PDM_std_decomposes_edges
(
       PDM_Mesh_nodal_elt_t  t_elt,
       int                   order,
       int                  *parent_node,
       int                   n_elt,
       int                  *n_elt_current,
       int                  *n_dedge_current,
       PDM_g_num_t           beg_gnum_elt_current,
       PDM_g_num_t           beg_gnum_edge_current,
 const PDM_g_num_t          *connectivity_elmt_vtx,
       int                  *elmt_edge_vtx_idx,
       PDM_g_num_t          *elmt_edge_vtx,
       PDM_g_num_t          *elmt_edge_cell,
       int                  *elmt_cell_edge_idx,
       PDM_g_num_t          *elmt_cell_edge,
       int                  *parent_elmt_position
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DMESH_NODAL_ELEMENTS_UTILS_H__ */
