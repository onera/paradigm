/*
 * \file
 */

#ifndef __PDM_MESH_NODAL_H__
#define __PDM_MESH_NODAL_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
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

typedef enum {

  PDM_BLOCK_ID_BLOCK_STD    = 0,
  PDM_BLOCK_ID_BLOCK_POLY2D = 1000000,
  PDM_BLOCK_ID_BLOCK_POLY3D = 2000000

} PDM_block_id_block_t;

/*----------------------------------------------------------------------------
 * Geometric element type
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_MESH_NODAL_POINT,
  PDM_MESH_NODAL_BAR2,
  PDM_MESH_NODAL_TRIA3,
  PDM_MESH_NODAL_QUAD4,
  PDM_MESH_NODAL_POLY_2D,
  PDM_MESH_NODAL_TETRA4,
  PDM_MESH_NODAL_PYRAMID5,
  PDM_MESH_NODAL_PRISM6,
  PDM_MESH_NODAL_HEXA8,
  PDM_MESH_NODAL_POLY_3D,
  PDM_MESH_NODAL_BARHO,
  PDM_MESH_NODAL_TRIAHO,
  PDM_MESH_NODAL_QUADHO,
  PDM_MESH_NODAL_TETRAHO,
  PDM_MESH_NODAL_PYRAMIDHO,
  PDM_MESH_NODAL_PRISMHO,
  PDM_MESH_NODAL_HEXAHO,
  PDM_MESH_NODAL_BARHO_BEZIER, // temporary add to visualize Bezier curves
  PDM_MESH_NODAL_TRIAHO_BEZIER, // temporary add to visualize Bezier triangles
  PDM_MESH_NODAL_N_ELEMENT_TYPES
} PDM_Mesh_nodal_elt_t;


typedef struct PDM_Mesh_nodal_vtx_t PDM_Mesh_nodal_vtx_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/**
 *
 * \brief Get the dimension of an element
 *
 * \param[in]  type    Element type
 *
 * \return     Dimension of the element
 *
 */

int
PDM_Mesh_nodal_elt_dim_get
(
 PDM_Mesh_nodal_elt_t type
 );

int
PDM_Mesh_nodal_is_2D_element
(
  PDM_Mesh_nodal_elt_t type
);

int
PDM_Mesh_nodal_is_3D_element
(
  PDM_Mesh_nodal_elt_t type
);


/**
 * \brief Get the number of vertices of an element type
 *
 * \param [in]   type     Element type
 * \param [in]   comm     Element order
 *
 * \return       Number of vertices
 *
 */

int
PDM_Mesh_nodal_n_vtx_elt_get
(
  PDM_Mesh_nodal_elt_t type,
  const int            order
);

/**
 * \brief Get if the element is consider as HO or not
 *
 * \param [in]   type     Element type
 *
 * \return       0 if element is not ho else 1
 *
 */
int
PDM_Mesh_nodal_elmt_is_ho
(
  PDM_Mesh_nodal_elt_t type
);


/**
 * \brief Returns the number of vertices in an element
 *
 * \param [in]  elt_type   Element type
 * \param [in]  order      Element order
 *
 * \return      Number of vertices in element
 *
 */
// A REMPLACER PAR PDM_Mesh_nodal_n_vtx_elt_get?

int
PDM_Mesh_nodal_n_vertices_element
(
 const PDM_Mesh_nodal_elt_t elt_type,
 const int                  order
);


void
PDM_Mesh_nodal_ho_parent_node
(
 const PDM_Mesh_nodal_elt_t  t_elt,
 const int                   order,
 const char                 *ho_ordering,
       int                  *parent_node
);


void
PDM_Mesh_nodal_reorder_elt_vtx
(
 const PDM_Mesh_nodal_elt_t  t_elt,
 const int                   order,
 const char                 *ho_ordering_in,
 const char                 *ho_ordering_out,
 const int                   n_elt,
       int                  *elt_vtx_in,
       int                  *elt_vtx_out
);



PDM_geometry_kind_t
PDM_Mesh_nodal_geom_kind_from_elt_type
(
 PDM_Mesh_nodal_elt_t t_elt
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


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MESH_NODAL_H__ */
