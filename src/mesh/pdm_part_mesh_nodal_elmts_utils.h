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

/**
 * \brief Decompose a PartMeshNodalElmts into faces.
 *        (All sections and partitions are concatenated.)
 *
 * \note Array sizes are computed by \ref PDM_part_mesh_nodal_elmts_decompose_faces_get_size
 *
 * \param [in]     pmne                  Pointer to \ref PDM_part_mesh_nodal_elmts_t instance
 * \param [in]     vtx_ln_to_gn          Vertex global IDs (size = n_part)
 * \param [inout]  elmt_face_vtx_idx     Index for ElementFace->Vertex connectivity
 * \param [inout]  elmt_face_vtx         ElementFace->Vertex connectivity (global IDs)
 * \param [inout]  elmt_cell_face_idx    Index for Cell->ElementFace connectivity
 * \param [inout]  elmt_face_cell        ElementFace->Cell connectivity (Parent Element) (global IDs)
 * \param [inout]  parent_elmt_position  Position in parent element for each face
 *
 */
void
PDM_part_mesh_nodal_elmts_sections_decompose_faces
(
  PDM_part_mesh_nodal_elmts_t  *pmne,
  PDM_g_num_t                 **vtx_ln_to_gn,
  int                          *elmt_face_vtx_idx,
  PDM_g_num_t                  *elmt_face_vtx,
  int                          *elmt_cell_face_idx,
  PDM_g_num_t                  *elmt_face_cell,
  int                          *parent_elmt_position
);


/**
 * \brief Compute size of a PartMeshNodalElmts decomposition into faces
 *        (All sections and partitions are concatenated.)
 *
 * \param [in] pmne                 Pointer to \ref PDM_part_mesh_nodal_elmts_t instance
 * \param [out] n_elt_tot           Total number of Elements (sum over all sections and all partitions)
 * \param [out] n_face_elt_tot      Total number of ElementFaces (sum over all sections and all partitions)
 * \param [out] n_sum_vtx_face_tot  Total number of ElementFaceVertices (sum over all sections and all partitions)
 * \param [out] elmt_face_vtx_idx   Index for ElementFace->Vertex connectivity (size = \p n_face_elt_tot + 1)
 * \param [out] elmt_cell_face_idx  Index for Cell->ElementFace connectivity (size = \p n_elt_tot + 1)
 */
void
PDM_part_mesh_nodal_elmts_decompose_faces_get_size
(
 PDM_part_mesh_nodal_elmts_t  *pmne,
 int                          *n_elt_tot,
 int                          *n_face_elt_tot,
 int                          *n_sum_vtx_face_tot,
 int                         **elmt_face_vtx_idx,
 int                         **elmt_cell_face_idx
);


/**
 * \brief Decompose a standard section into edges
 *        (All sections and partitions are concatenated.)
 *
 * \note Array sizes are computed by \ref PDM_part_mesh_nodal_elmts_decompose_edges_get_size
 *
 * \param [in]     pmne                  Pointer to \ref PDM_part_mesh_nodal_elmts_t instance
 * \param [in]     vtx_ln_to_gn          Vertex global IDs (size = n_part)
 * \param [inout]  elmt_edge_vtx_idx     Index for ElementEdge->Vertex connectivity (trivial)
 * \param [inout]  elmt_edge_vtx         ElementEdge->Vertex connectivity (global IDs)
 * \param [inout]  elmt_cell_edge_idx    Index for Cell->ElementEdge connectivity
 * \param [inout]  elmt_edge_cell        ElementEdge->Cell connectivity (Parent Element) (global IDs)
 * \param [inout]  parent_elmt_position  Position in parent element for each edge
 *
 */
void
PDM_part_mesh_nodal_elmts_sections_decompose_edges
(
  PDM_part_mesh_nodal_elmts_t  *pmne,
  PDM_g_num_t                 **vtx_ln_to_gn,
  int                          *elmt_edge_vtx_idx,
  PDM_g_num_t                  *elmt_edge_vtx,
  int                          *elmt_cell_edge_idx,
  PDM_g_num_t                  *elmt_edge_cell,
  int                          *parent_elmt_position
);


/**
 * \brief Compute size of a PartMeshNodalElmts decomposition into edges
 *        (All sections and partitions are concatenated.)
 *
 * \param [in] pmne                 Pointer to \ref PDM_part_mesh_nodal_elmts_t instance
 * \param [out] n_elt_tot           Total number of Elements (sum over all sections and all partitions)
 * \param [out] n_edge_elt_tot      Total number of ElementEdges (sum over all sections and all partitions)
 * \param [out] n_sum_vtx_edge_tot  Total number of ElementEdgeVertices (sum over all sections and all partitions)
 * \param [out] elmt_edge_vtx_idx   Index for ElementEdge->Vertex connectivity (size = \p n_edge_elt_tot + 1)
 * \param [out] elmt_cell_edge_idx  Index for Cell->ElementEdge connectivity (size = \p n_elt_tot + 1)
 */
void
PDM_part_mesh_nodal_elmts_decompose_edges_get_size
(
 PDM_part_mesh_nodal_elmts_t *pmne,
 int                         *n_elt_tot,
 int                         *n_edge_elt_tot,
 int                         *n_sum_vtx_edge_tot
);


/**
 * \brief Decompose a poly2d section into edges
 *
 * \param [in]    n_elt                      Number of elements
 * \param [inout] n_elt_current              Current position in concatenated sections
 * \param [inout] n_edge_current             Current position in concatenated edges
 * \param [in]    vtx_ln_to_gn               Vertex global IDs
 * \param [in]    connectivity_elmt_vtx      Element->Vertex connectivity (local IDs, size = \p connectivity_elmt_vtx_idx[\p n_elt])
 * \param [in]    connectivity_elmt_vtx_idx  Index for Element->Vertex connectivity (size = \p n_elt + 1)
 * \param [in]    elmt_ln_to_gn              Element global IDs (size = \p n_elt)
 * \param [inout] elmt_edge_vtx_idx          Index for ElementEdge->Vertex connectivity
 * \param [inout] elmt_edge_vtx              ElementEdge->Vertex connectivity (global IDs)
 * \param [inout] elmt_cell_edge_idx         Index for Cell->ElementEdge connectivity
 * \param [inout] elmt_edge_cell             ElementEdge->Cell connectivity (Parent Element) (global IDs)
 * \param [inout] parent_elmt_position       Position in parent element for each edge
 *
 */
void
PDM_part_mesh_nodal_poly2d_decomposes_edges
(
       int                   n_elt,
       int                  *n_elt_current,
       int                  *n_edge_current,
 const PDM_g_num_t          *vtx_ln_to_gn,
 const int                  *connectivity_elmt_vtx,
 const int                  *connectivity_elmt_vtx_idx,
 const PDM_g_num_t          *elmt_ln_to_gn,
       int                  *elmt_edge_vtx_idx,
       PDM_g_num_t          *elmt_edge_vtx,
       int                  *elmt_cell_edge_idx,
       PDM_g_num_t          *elmt_edge_cell,
       int                  *parent_elmt_position
);


/**
 * \brief Decompose locally a standard element section into edges (single partition).
 *
 * \param [in]    t_elt                  Element type
 * \param [in]    n_elt                  Number of elements
 * \param [in]    order                  Element order
 * \param [in]    parent_node            Permutation of principal nodes (for high-order elements)
 * \param [inout] n_elt_current          Current position in concatenated sections
 * \param [inout] n_edge_current         Current position in concatenated edges
 * \param [in]    connectivity_elmt_vtx  Element->Vertex connectivity (local IDs)
 * \param [in]    parent_num                 Element parent numbering (or NULL)
 * \param [inout] elmt_edge_vtx_idx      Index for ElementEdge->Vertex connectivity
 * \param [inout] elmt_edge_vtx          ElementEdge->Vertex connectivity (local IDs)
 * \param [inout] elmt_cell_edge_idx     Index for ElementEdge->Edge connectivity
 * \param [inout] parent_elmt            ElementEdge->Cell connectivity (Parent Element) (local IDs)
 * \param [inout] parent_elmt_position   Position in parent element for each edge
 *
 */
void
PDM_part_mesh_nodal_std_decompose_local_edges
(
       PDM_Mesh_nodal_elt_t  t_elt,
       int                   n_elt,
       int                   order,
       int                  *parent_node,
       int                  *n_elt_current,
       int                  *n_edge_current,
 const int                  *connectivity_elmt_vtx,
 const int                  *parent_num,
       int                  *elmt_edge_vtx_idx,
       int                  *elmt_edge_vtx,
       int                  *elmt_cell_edge_idx,
       int                  *parent_elmt,
       int                  *parent_elmt_position
);


/**
 * \brief Decompose locally a poly2d section into edges (single partition).
 *
 * \param [in]    n_elt                      Number of elements
 * \param [inout] n_elt_current              Current position in concatenated sections
 * \param [inout] n_edge_current             Current position in concatenated edges
 * \param [in]    connectivity_elmt_vtx      Element->Vertex connectivity (local IDs)
 * \param [in]    connectivity_elmt_vtx_idx  Index for Element->Vertex connectivity
 * \param [in]    parent_num                 Element parent numbering (or NULL)
 * \param [inout] elmt_edge_vtx_idx          Index for ElementEdge->Vertex connectivity
 * \param [inout] elmt_edge_vtx              ElementEdge->Vertex connectivity (local IDs)
 * \param [inout] elmt_cell_edge_idx         Index for ElementEdge->Edge connectivity
 * \param [inout] parent_elmt                ElementEdge->Cell connectivity (Parent Element) (local IDs)
 * \param [inout] parent_elmt_position       Position in parent element for each edge
 *
 */
void
PDM_part_mesh_nodal_poly2d_decompose_local_edges
(
       int                   n_elt,
       int                  *n_elt_current,
       int                  *n_edge_current,
 const int                  *connectivity_elmt_vtx,
 const int                  *connectivity_elmt_vtx_idx,
 const int                  *parent_num,
       int                  *elmt_edge_vtx_idx,
       int                  *elmt_edge_vtx,
       int                  *elmt_cell_edge_idx,
       int                  *parent_elmt,
       int                  *parent_elmt_position
);


/**
 * \brief Decompose locally a PartMeshNodalElements into edges.
 *
 * \param [in]  pmne                       Pointer to \ref PDM_part_mesh_nodal_elmts_t instance
 * \param [out] out_n_decompose_elmt_edge  Number of ElementEdges (size = n_part)
 * \param [out] out_elmt_edge_idx          Index for Cell->ElementEdge connectivity (size = n_part)
 * \param [out] out_elmt_edge_vtx_idx      Index for ElementEdge->Vertex connectivity (size = n_part)
 * \param [out] out_elmt_edge_vtx          ElementEdge->Vertex connectivity (local IDs, size = n_part)
 * \param [out] out_parent_elmt            ElementEdge->Cell connectivity (Parent Element) (local IDs, size = n_part)
 * \param [out] out_parent_elmt_position   Position in parent element for each edge (size = n_part)
 *
 */
void
PDM_part_mesh_nodal_elmts_sections_local_decompose_edges
(
  PDM_part_mesh_nodal_elmts_t   *pmne,
  int                          **out_n_decompose_elmt_edge,
  int                         ***out_elmt_edge_idx,
  int                         ***out_elmt_edge_vtx_idx,
  int                         ***out_elmt_edge_vtx,
  int                         ***out_parent_elmt,
  int                         ***out_parent_elmt_position
);


/**
 * \brief Decompose locally a standard element section into faces (single partition).
 *
 * \param [in]    t_elt                  Element type
 * \param [in]    n_elt                  Number of elements
 * \param [in]    order                  Element order
 * \param [in]    parent_node            Permutation of principal nodes (for high-order elements)
 * \param [inout] n_elt_current          Current position in concatenated sections
 * \param [inout] n_face_current         Current position in concatenated faces
 * \param [in]    connectivity_elmt_vtx  Element->Vertex connectivity (local IDs)
 * \param [in]    parent_num             Element parent numbering (or NULL)
 * \param [inout] elmt_face_vtx_idx      Index for ElementFace->Vertex connectivity
 * \param [inout] elmt_face_vtx          ElementFace->Vertex connectivity (local IDs)
 * \param [inout] elmt_cell_face_idx     Index for ElementFace->Face connectivity
 * \param [inout] parent_elmt            ElementFace->Cell connectivity (Parent Element) (local IDs)
 * \param [inout] parent_elmt_position   Position in parent element for each edge
 *
 */
void
PDM_part_mesh_nodal_std_decompose_local_faces
(
       PDM_Mesh_nodal_elt_t  t_elt,
       int                   n_elt,
       int                   order,
       int                  *parent_node,
       int                  *n_elt_current,
       int                  *n_face_current,
 const int                  *connectivity_elmt_vtx,
 const int                  *parent_num,
       int                  *elmt_face_vtx_idx,
       int                  *elmt_face_vtx,
       int                  *elmt_cell_face_idx,
       int                  *parent_elmt,
       int                  *parent_elmt_position
);


/**
 * \brief Decompose locally a poly2d section into faces (single partition).
 *
 * \param [in]    n_elt                      Number of elements
 * \param [inout] n_elt_current              Current position in concatenated sections
 * \param [inout] n_face_current             Current position in concatenated faces
 * \param [in]    connectivity_elmt_vtx      Element->Vertex connectivity (local IDs)
 * \param [in]    connectivity_elmt_vtx_idx  Index for Element->Vertex connectivity
 * \param [in]    parent_num                 Element parent numbering (or NULL)
 * \param [inout] elmt_face_vtx_idx          Index for Elementface->Vertex connectivity
 * \param [inout] elmt_face_vtx              ElementFace->Vertex connectivity (local IDs)
 * \param [inout] elmt_cell_face_idx         Index for ElementFace->Face connectivity
 * \param [inout] parent_elmt                ElementFace->Cell connectivity (Parent Element) (local IDs)
 * \param [inout] parent_elmt_position       Position in parent element for each face
 *
 */
void
PDM_part_mesh_nodal_poly2d_decompose_local_faces
(
       int  n_elt,
       int *n_elt_current,
       int *n_face_current,
 const int *elt_vtx_idx,
 const int *elt_vtx,
 const int *parent_num,
       int *elmt_face_vtx_idx,
       int *elmt_face_vtx,
       int *elmt_cell_face_idx,
       int *parent_elmt,
       int *parent_elmt_position
);


/**
 * \brief Decompose locally a PartMeshNodalElements into faces.
 *
 * \param [in]  pmne                       Pointer to \ref PDM_part_mesh_nodal_elmts_t instance
 * \param [out] out_n_decompose_elmt_face  Number of ElementFaces (size = n_part)
 * \param [out] out_elmt_face_idx          Index for Cell->ElementFace connectivity (size = n_part)
 * \param [out] out_elmt_face_vtx_idx      Index for ElementFace->Vertex connectivity (size = n_part)
 * \param [out] out_elmt_face_vtx          ElementFace->Vertex connectivity (local IDs, size = n_part)
 * \param [out] out_parent_elmt            ElementFace->Cell connectivity (Parent Element) (local IDs, size = n_part)
 * \param [out] out_parent_elmt_position   Position in parent element for each face (size = n_part)
 *
 */
void
PDM_part_mesh_nodal_elmts_sections_local_decompose_faces
(
  PDM_part_mesh_nodal_elmts_t   *pmne,
  int                          **out_n_decompose_elmt_face,
  int                         ***out_elmt_face_idx,
  int                         ***out_elmt_face_vtx_idx,
  int                         ***out_elmt_face_vtx,
  int                         ***out_parent_elmt,
  int                         ***out_parent_elmt_position
);


/**
 * \brief Compute link between a parent and child PartMeshNodalElmts
 *
 * \note Dimension of \p pmne_parent must be greater than that of \p pmne_child
 *
 * \param [in]  pmne_parent               Pointer to parent \ref PDM_part_mesh_nodal_elmts_t instance
 * \param [in]  pmne_child                Pointer to child \ref PDM_part_mesh_nodal_elmts_t instance
 * \param [in]  only_child_link           Build only ascending link (child->parent)
 * \param [out] out_child_to_parent_idx   Index for child->parent connectivity (size = n_child+1 (total number of elt in \p pmne_child))
 * \param [out] out_child_to_parent       Child->parent connectivity (SIGN???)
 * \param [out] out_n_entity              Number of child entities in pmne_parent (if \p only_child_link is PDM_FALSE, size = \p n_part)
 * \param [out] out_entity_to_vtx_idx     Index for ChildEntity->Vertex connectivity (if \p only_child_link is PDM_FALSE, size = \p n_part and for each part, size = \p out_n_entity + 1)
 * \param [out] out_entity_to_vtx         ChildEntity->Vertex connectivity (if \p only_child_link is PDM_FALSE, size = \p n_part)
 * \param [out] out_parent_to_entity_idx  Index for ParentEntity->ChildEntity connectivity (if \p only_child_link is PDM_FALSE, size = \p n_part)
 * \param [out] out_parent_to_entity      ParentEntity->ChildEntity connectivity (if \p only_child_link is PDM_FALSE, size = \p n_part)
 */
void
PDM_part_mesh_nodal_elmts_compute_child_entities
(
 PDM_part_mesh_nodal_elmts_t   *pmne_parent,
 PDM_part_mesh_nodal_elmts_t   *pmne_child,
 PDM_bool_t                     only_child_link,
 int                         ***out_child_to_parent_idx,
 int                         ***out_child_to_parent,
 int                          **out_n_entity,
 int                         ***out_entity_to_vtx_idx,
 int                         ***out_entity_to_vtx,
 int                         ***out_parent_to_entity_idx,
 int                         ***out_parent_to_entity
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_MESH_NODAL_ELMTS_UTILS_H__ */
