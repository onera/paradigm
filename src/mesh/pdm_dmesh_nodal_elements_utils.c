
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_error.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_dmesh_nodal_elmts.h"
#include "pdm_dmesh_nodal_elmts_priv.h"
#include "pdm_dmesh_nodal_elements_utils.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
*
* \brief Decompose standard elements to a flatten view of faces
*/
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
)
{
  PDM_UNUSED(beg_gnum_face_current);
  PDM_UNUSED(elmt_cell_face_idx);
  PDM_UNUSED(elmt_cell_face);

  int parent_node_std[8] = {0, 1, 2, 3, 4, 5, 6, 7};

  int *_parent_node;
  if (parent_node == NULL) {
    _parent_node = parent_node_std;
  } else {
    _parent_node = parent_node;
  }

  const int *elt_face_vtx_idx = NULL;
  const int *elt_face_vtx     = NULL;
  int n_face_elt = PDM_face_vtx_per_elmt(t_elt,
                                         &elt_face_vtx_idx,
                                         &elt_face_vtx);

  if (n_face_elt <= 0 || elt_face_vtx == NULL) {
    return;
  }

  int n_sum_vtx_face = elt_face_vtx_idx[n_face_elt];
  int n_vtx_elt      = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);

  int          _n_face_current            = *n_dface_current;
  int         *_current_elmt_face_vtx_idx = elmt_face_vtx_idx    + _n_face_current;
  PDM_g_num_t *_current_elmt_face_vtx     = elmt_face_vtx        + elmt_face_vtx_idx[_n_face_current];
  PDM_g_num_t *_current_elmt_face_cell    = elmt_face_cell       + _n_face_current;
  int         *_parent_elmt_position      = parent_elmt_position + _n_face_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int i_elt = 0; i_elt < n_elt; i_elt++) {

    for (int i_face = 0; i_face < n_face_elt; i_face++) {

      _current_elmt_face_vtx_idx[i_elt * n_face_elt + i_face + 1] = _current_elmt_face_vtx_idx[i_elt * n_face_elt + i_face] + elt_face_vtx_idx[i_face+1] - elt_face_vtx_idx[i_face];
      _current_elmt_face_cell   [i_elt * n_face_elt + i_face    ] =  beg_gnum_elt_current + i_elt + 1;
      _parent_elmt_position     [i_elt * n_face_elt + i_face    ] =  i_face;

      for (int i = elt_face_vtx_idx[i_face]; i < elt_face_vtx_idx[i_face+1]; i++) {
        int i_vtx = elt_face_vtx[i];
        _current_elmt_face_vtx[n_sum_vtx_face*i_elt + i] = connectivity_elmt_vtx[n_vtx_elt * i_elt + _parent_node[i_vtx]];
      }

    }

  }

  *n_elt_current   += n_elt;
  *n_dface_current += n_elt * n_face_elt;
}

/**
*
* \brief Decompose standard elements to a flatten view of edges
*/
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
)
{
  PDM_UNUSED(beg_gnum_edge_current);
  PDM_UNUSED(elmt_cell_edge_idx);
  PDM_UNUSED(elmt_cell_edge);

  int parent_node_std[8] = {0, 1, 2, 3, 4, 5, 6, 7};

  int *_parent_node;
  if (parent_node == NULL) {
    _parent_node = parent_node_std;
  } else {
    _parent_node = parent_node;
  }

  const int *elt_edge_vtx = NULL;
  int n_edge_elt = PDM_edge_vtx_per_elmt(t_elt,
                                         &elt_edge_vtx);

  if (n_edge_elt <= 0 || elt_edge_vtx == NULL) {
    return;
  }

  int n_sum_vtx_edge = n_edge_elt * 2;
  int n_vtx_elt      = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);

  int          _n_edge_current            = *n_dedge_current;
  int         *_current_elmt_edge_vtx_idx = elmt_edge_vtx_idx    + _n_edge_current;
  PDM_g_num_t *_current_elmt_edge_vtx     = elmt_edge_vtx        + elmt_edge_vtx_idx[_n_edge_current];
  PDM_g_num_t *_current_elmt_edge_cell    = elmt_edge_cell       + _n_edge_current;
  int         *_parent_elmt_position      = parent_elmt_position + _n_edge_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int i_elt = 0; i_elt < n_elt; i_elt++) {

    for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {

      _current_elmt_edge_vtx_idx[i_elt * n_edge_elt + i_edge + 1] = _current_elmt_edge_vtx_idx[i_elt * n_edge_elt + i_edge] + 2;
      _current_elmt_edge_cell   [i_elt * n_edge_elt + i_edge    ] =  beg_gnum_elt_current + i_elt + 1;
      _parent_elmt_position     [i_elt * n_edge_elt + i_edge    ] =  i_edge;

      for (int i = 2*i_edge; i < 2*(i_edge+1); i++) {
        int i_vtx = elt_edge_vtx[i];
        _current_elmt_edge_vtx[n_sum_vtx_edge*i_elt + i] = connectivity_elmt_vtx[n_vtx_elt * i_elt + _parent_node[i_vtx]];
      }

    }

  }

  *n_elt_current   += n_elt;
  *n_dedge_current += n_elt * n_edge_elt;
}


/**
*
* \brief Decompose polygons to a flatten view of edges
*/
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
)
{
  PDM_UNUSED(elmt_cell_face_idx);
  PDM_UNUSED(elmt_cell_face);
  PDM_UNUSED(beg_gnum_face_current);

  int          _n_face_current            = *n_face_current;
  int         *_current_elmt_face_vtx_idx = elmt_face_vtx_idx    + _n_face_current;
  PDM_g_num_t *_current_elmt_face_vtx     = elmt_face_vtx        + elmt_face_vtx_idx[_n_face_current];
  PDM_g_num_t *_current_elmt_face_cell    = elmt_face_cell       + _n_face_current;
  int         *_parent_elmt_position      = parent_elmt_position + _n_face_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  int idx = 0;
  for (int ielt = 0; ielt < n_elt; ielt++) {
    int beg = connectivity_elmt_vtx_idx[ielt];
    int n_vtx_on_face = connectivity_elmt_vtx_idx[ielt+1] - beg;
    _current_elmt_face_vtx_idx[idx + 1] = _current_elmt_face_vtx_idx[idx] + n_vtx_on_face;
    _current_elmt_face_cell   [idx    ] = beg_gnum_elt_current + ielt + 1;
    _parent_elmt_position     [idx    ] = 0;
    for(int ivtx = 0; ivtx < n_vtx_on_face; ++ivtx ) {
       _current_elmt_face_vtx[idx++] = connectivity_elmt_vtx[beg+ivtx];
    }
  }

  *n_elt_current  += n_elt;
  *n_face_current += n_elt;
}

/**
*
* \brief Decompose polygons to a flatten view of edges
*/
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
)
{
  PDM_UNUSED(elmt_cell_edge_idx);
  PDM_UNUSED(elmt_cell_edge);
  PDM_UNUSED(beg_gnum_edge_current);

  int          _n_edge_current            = *n_edge_current;
  int         *_current_elmt_edge_vtx_idx = elmt_edge_vtx_idx    + _n_edge_current;
  PDM_g_num_t *_current_elmt_edge_vtx     = elmt_edge_vtx        + elmt_edge_vtx_idx[_n_edge_current];
  PDM_g_num_t *_current_elmt_edge_cell    = elmt_edge_cell       + _n_edge_current;
  int         *_parent_elmt_position      = parent_elmt_position + _n_edge_current;
  /*
   * For each element we flaten all connectivities in one array
   */
  int idx = 0;
  for (int ielt = 0; ielt < n_elt; ielt++) {
    // Reminder for poly2d -> Number of vertex = Number of edge
    int n_edge_elt = connectivity_elmt_vtx_idx[ielt+1] - connectivity_elmt_vtx_idx[ielt];
    *n_edge_current += n_edge_elt;
    int idx2 = connectivity_elmt_vtx_idx[ielt];
    for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
      _current_elmt_edge_vtx_idx[idx + 1] = _current_elmt_edge_vtx_idx[idx] + 2;
      _current_elmt_edge_cell   [idx    ] = beg_gnum_elt_current + ielt + 1;
      _parent_elmt_position     [idx    ] = i_edge;

      int inext = (i_edge + 1) % n_edge_elt;
      _current_elmt_edge_vtx[2 * idx    ]  = connectivity_elmt_vtx[idx2 + i_edge];
      _current_elmt_edge_vtx[2 * idx + 1]  = connectivity_elmt_vtx[idx2 + inext ];

      idx += 1;
    }
  }

  *n_elt_current += n_elt;
}



/**
*
* \brief Decompose polyhedra to a flatten view of edges
*/
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
)
{
  PDM_UNUSED(n_elt);
  PDM_UNUSED(n_elt_current);
  PDM_UNUSED(n_face_current);
  PDM_UNUSED(beg_gnum_elt_current);
  PDM_UNUSED(beg_gnum_face_current);
  PDM_UNUSED(connectivity_elmt_vtx);
  PDM_UNUSED(connectivity_elmt_vtx_idx);
  PDM_UNUSED(elmt_face_vtx_idx);
  PDM_UNUSED(elmt_face_vtx);
  PDM_UNUSED(elmt_face_cell);
  PDM_UNUSED(elmt_cell_face_idx);
  PDM_UNUSED(elmt_cell_face);
  PDM_UNUSED(parent_elmt_position);

  PDM_error(__FILE__, __LINE__, 0, "PDM_poly3d_decomposes_faces : Not yet implemented\n");
}

/**
*
* \brief Decompose polyhedra to a flatten view of edges
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
)
{
  PDM_UNUSED(n_elt);
  PDM_UNUSED(n_elt_current);
  PDM_UNUSED(n_edge_current);
  PDM_UNUSED(beg_gnum_elt_current);
  PDM_UNUSED(beg_gnum_edge_current);
  PDM_UNUSED(connectivity_elmt_vtx);
  PDM_UNUSED(connectivity_elmt_vtx_idx);
  PDM_UNUSED(elmt_edge_vtx_idx);
  PDM_UNUSED(elmt_edge_vtx);
  PDM_UNUSED(elmt_edge_cell);
  PDM_UNUSED(elmt_cell_edge_idx);
  PDM_UNUSED(elmt_cell_edge);
  PDM_UNUSED(parent_elmt_position);
  PDM_error(__FILE__, __LINE__, 0, "PDM_poly3d_decomposes_edges : Not yet implemented\n");
}



/**
*
* \brief PDM_sections_decompose_faces
*
* \param [in]     mesh               Current mesh
* \param [inout]  elt_face_vtx_idx   Index of element faces connectivity (preallocated)
* \param [inout]  elt_face_vtx       Element faces connectivity (preallocated)
* \param [inout]  elmt_face_cell     Element faces connectivity (preallocated or NULL )
* \param [inout]  elmt_cell_face     Element faces connectivity (preallocated or NULL )
*
*/
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
)
{

  // A faire : local_num_in_parent_element

  int n_elt_current   = 0;
  int n_dface_current = 0;

  /* We need to loop over all sections in the good order */
  int n_section = dmn_elts->n_section;


  int parent_node[8];

  for (int i_section = 0; i_section < n_section; i_section++) {

    int id_section = dmn_elts->sections_id[i_section];

    const PDM_g_num_t* distrib = PDM_DMesh_nodal_elmts_distrib_section_get(dmn_elts, id_section);

    if(0 == 1) {
      printf("i_section = %i [%i] \n", i_section, n_section);
      printf("id_section = %i \n", id_section);
      printf("distrib[%i] = "PDM_FMT_G_NUM" \n", dmn_elts->i_rank, distrib[dmn_elts->i_rank]);
      printf("dmn_elts->section_distribution[%i] = "PDM_FMT_G_NUM" \n", i_section, dmn_elts->section_distribution[i_section]);
    }

    PDM_g_num_t beg_elmt_gnum = distrib[dmn_elts->i_rank] + dmn_elts->section_distribution[i_section];
    PDM_g_num_t beg_face_gnum = 0; // Useless in this context

    PDM_Mesh_nodal_elt_t t_elt = PDM_DMesh_nodal_elmts_section_type_get(dmn_elts, id_section);
    int n_elt                  = PDM_DMesh_nodal_elmts_section_n_elt_get(dmn_elts, id_section);


    if (t_elt == PDM_MESH_NODAL_POLY_2D) {
      /* Polygons */
      int         *connec_idx = NULL;
      PDM_g_num_t *connec     = NULL;
      PDM_DMesh_nodal_elmts_section_poly2d_get(dmn_elts,
                                               id_section,
                                               &connec_idx,
                                               &connec);

      PDM_poly2d_decomposes_faces(n_elt,
                                  &n_elt_current,
                                  &n_dface_current,
                                  beg_elmt_gnum,
                                  beg_face_gnum,
                                  connec,
                                  connec_idx,
                                  elmt_face_vtx_idx,
                                  elmt_face_vtx,
                                  elmt_face_cell,
                                  elmt_cell_face_idx,
                                  elmt_cell_face,
                                  parent_elmt_position);
    }

    else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
      /* Polyhedra */
      PDM_error(__FILE__, __LINE__, 0, "Error PDM_sections_decompose_faces : POLY_3D not yet supported\n");
    }

    else {
      /* Standard elements */
      int order = dmn_elts->sections_std[i_section]->order;
      int *_parent_node = NULL;

      PDM_g_num_t *connec = PDM_DMesh_nodal_elmts_section_std_get(dmn_elts, id_section);

      if (PDM_Mesh_nodal_elmt_is_ho(t_elt)) {
        const char* ho_ordering = dmn_elts->sections_std[i_section]->ho_ordering;
        PDM_Mesh_nodal_ho_parent_node(t_elt,
                                      order,
                                      ho_ordering,
                                      parent_node);
        _parent_node = parent_node;
      }
      else {
        order = 1;
      }

      PDM_std_decomposes_faces(t_elt,
                               order,
                               _parent_node,
                               n_elt,
                               &n_elt_current,
                               &n_dface_current,
                               beg_elmt_gnum,
                               beg_face_gnum,
                               connec,
                               elmt_face_vtx_idx,
                               elmt_face_vtx,
                               elmt_face_cell,
                               elmt_cell_face_idx,
                               elmt_cell_face,
                               parent_elmt_position);

    }
  }
}

/**
*
* \brief PDM_sections_decompose_faces
*
* \param [in]     mesh               Current mesh
* \param [inout]  elmt_edge_vtx_idx  Index of element faces connectivity (preallocated)
* \param [inout]  elmt_edge_vtx      Element faces connectivity (preallocated)
* \param [inout]  elmt_edge_cell     Element faces connectivity (preallocated or NULL )
* \param [inout]  elmt_cell_edge     Element faces connectivity (preallocated or NULL )
*
*/
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
)
{
  int n_elt_current   = 0;
  int n_dedge_current = 0;

  if(dmn_elts == NULL) {
    return;
  }

  int parent_node[8];
  /* We need to loop over all sections in the good order */
  int n_section = dmn_elts->n_section;

  for (int i_section = 0; i_section < n_section; i_section++) {

    int id_section = dmn_elts->sections_id[i_section];

    const PDM_g_num_t* distrib = PDM_DMesh_nodal_elmts_distrib_section_get(dmn_elts, id_section);

    if(0 == 1) {
      printf("i_section = %i [%i] \n", i_section, n_section);
      printf("id_section = %i \n", id_section);
      printf("distrib[%i] = "PDM_FMT_G_NUM" \n", dmn_elts->i_rank, distrib[dmn_elts->i_rank]);
      printf("dmn_elts->section_distribution[%i] = "PDM_FMT_G_NUM" \n", i_section, dmn_elts->section_distribution[i_section]);
    }

    PDM_g_num_t beg_elmt_gnum = distrib[dmn_elts->i_rank] + dmn_elts->section_distribution[i_section];
    PDM_g_num_t beg_edge_gnum = 0; // Useless in this context

    PDM_Mesh_nodal_elt_t t_elt = PDM_DMesh_nodal_elmts_section_type_get(dmn_elts, id_section);
    int n_elt                  = PDM_DMesh_nodal_elmts_section_n_elt_get(dmn_elts, id_section);

    if (t_elt == PDM_MESH_NODAL_POLY_2D) {
      /* Polygons */
      int         *connec_idx = NULL;
      PDM_g_num_t *connec     = NULL;
      PDM_DMesh_nodal_elmts_section_poly2d_get(dmn_elts,
                                               id_section,
                                               &connec_idx,
                                               &connec);

      PDM_poly2d_decomposes_edges(n_elt,
                                  &n_elt_current,
                                  &n_dedge_current,
                                  beg_elmt_gnum,
                                  beg_edge_gnum,
                                  connec,
                                  connec_idx,
                                  elmt_edge_vtx_idx,
                                  elmt_edge_vtx,
                                  elmt_edge_cell,
                                  elmt_cell_edge_idx,
                                  elmt_cell_edge,
                                  parent_elmt_position);
    }

    else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
      /* Polyhedra */
      PDM_error(__FILE__, __LINE__, 0, "Error PDM_sections_decompose_edges : POLY_3D not yet supported\n");
    }

    else {
      /* Standard elements */
      int order = dmn_elts->sections_std[i_section]->order;
      int *_parent_node = NULL;

      PDM_g_num_t *connec = PDM_DMesh_nodal_elmts_section_std_get(dmn_elts, id_section);

      if (PDM_Mesh_nodal_elmt_is_ho(t_elt)) {
        const char* ho_ordering = dmn_elts->sections_std[i_section]->ho_ordering;
        PDM_Mesh_nodal_ho_parent_node(t_elt,
                                      order,
                                      ho_ordering,
                                      parent_node);
        _parent_node = parent_node;
      }
      else {
        order = 1;
      }

      PDM_std_decomposes_edges(t_elt,
                               order,
                               _parent_node,
                               n_elt,
                               &n_elt_current,
                               &n_dedge_current,
                               beg_elmt_gnum,
                               beg_edge_gnum,
                               connec,
                               elmt_edge_vtx_idx,
                               elmt_edge_vtx,
                               elmt_edge_cell,
                               elmt_cell_edge_idx,
                               elmt_cell_edge,
                               parent_elmt_position);

    }
  }
}

