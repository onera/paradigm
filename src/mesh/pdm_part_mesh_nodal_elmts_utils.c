/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_mesh_nodal_priv.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_part_mesh_nodal_elmts_utils.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_geom_elem.h"
#include "pdm_array.h"
#include "pdm_vtk.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Maximum number of blocks depending of block type
 *----------------------------------------------------------------------------*/

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
* \brief Decompose tetra cell_vtx connectivity to a flatten view of faces
*/
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
 const PDM_g_num_t *elmt_ln_to_gn,
       int         *elmt_face_vtx_idx,
       PDM_g_num_t *elmt_face_vtx,
       int         *elmt_cell_face_idx,
       PDM_g_num_t *elmt_face_cell,
       int         *parent_elmt_position
)
{
  const int n_face_elt        = 4;
  const int n_sum_vtx_face    = 12;
  int n_sum_vtx_elt           = 4;

  if(order > 1) {
    n_sum_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TETRAHO, order);
  }

  int __parent_node[4] = {0, 1, 2, 3};
  int *_parent_node;
  if (parent_node == NULL) {
    _parent_node = __parent_node;
  } else {
    _parent_node = parent_node;
  }

  int _n_face_current = *n_face_current;
  int _n_elt_current  = *n_elt_current;
  int         *_current_elmt_face_vtx_idx = elmt_face_vtx_idx    + _n_face_current;
  PDM_g_num_t *_current_elmt_face_vtx     = elmt_face_vtx        + elmt_face_vtx_idx[_n_face_current];
  int         *_parent_elmt_position      = parent_elmt_position + _n_face_current;
  int         *_elmt_cell_face_idx        = elmt_cell_face_idx   + _n_elt_current;
  PDM_g_num_t *_elmt_face_cell            = elmt_face_cell       + _n_face_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    for (int i_face = 0; i_face < n_face_elt; i_face++) {
      _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face + 1] = _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face] + 3;
      _parent_elmt_position     [ielt * n_face_elt + i_face    ] = i_face;
      _elmt_face_cell           [ielt * n_face_elt + i_face    ] = elmt_ln_to_gn[ielt];
    }

    _elmt_cell_face_idx[ielt+1] = _elmt_cell_face_idx[ielt] + n_face_elt;

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 0]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 1]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 2]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 3]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 4]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 5]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 6]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 7]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 8]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 9]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 10] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 11] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];

  }

  *n_elt_current  += n_elt;
  *n_face_current += n_elt * n_face_elt;

}

void
PDM_part_mesh_nodal_hexa_decomposes_faces
(
       int          n_elt,
       int          order,
       int         *parent_node,
       int         *n_elt_current,
       int         *n_face_current,
 const PDM_g_num_t *vtx_ln_to_gn,
 const int         *connectivity_elmt_vtx,
 const PDM_g_num_t *elmt_ln_to_gn,
       int         *elmt_face_vtx_idx,
       PDM_g_num_t *elmt_face_vtx,
       int         *elmt_cell_face_idx,
       int         *elmt_face_cell,
       int         *parent_elmt_position
)
{
  const int n_face_elt     = 6;
  const int n_sum_vtx_face = 24;
  int n_sum_vtx_elt        = 8;
  if(order > 1) {
    n_sum_vtx_elt     = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_HEXAHO, order);
  }

  int __parent_node[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  int *_parent_node;
  if (parent_node == NULL) {
    _parent_node = __parent_node;
  } else {
    _parent_node = parent_node;
  }

  int _n_face_current = *n_face_current;
  int _n_elt_current  = *n_elt_current;
  int         *_current_elmt_face_vtx_idx = elmt_face_vtx_idx    + _n_face_current;
  PDM_g_num_t *_current_elmt_face_vtx     = elmt_face_vtx        + elmt_face_vtx_idx[_n_face_current];
  int         *_parent_elmt_position      = parent_elmt_position + _n_face_current;
  int         *_elmt_cell_face_idx        = elmt_cell_face_idx   + _n_elt_current;
  PDM_g_num_t *_elmt_face_cell            = elmt_face_cell       + _n_face_current;

  /*
   * For each element we flaten all connectivities in one array
   */


  for (int ielt = 0; ielt < n_elt; ielt++) {

    /* Store the face_cell */
    for (int i_face = 0; i_face < n_face_elt; i_face++) {
      _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face + 1] = _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face] + 4;
      _parent_elmt_position     [ielt * n_face_elt + i_face    ] = i_face;
      _elmt_face_cell           [ielt * n_face_elt + i_face    ] = elmt_ln_to_gn[ielt];
    }

    _elmt_cell_face_idx[ielt+1] = _elmt_cell_face_idx[ielt] + n_face_elt;

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 0]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 1]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 2]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 3]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 4]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[6]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 5]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[7]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 6]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[4]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 7]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[5]]-1];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 8]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[4]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 9]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[7]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 10] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 11] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 12] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[7]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 13] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[6]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 14] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 15] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 16] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 17] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[6]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 18] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[5]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 19] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 20] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 21] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[5]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 22] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[4]]-1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 23] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];

  }


  *n_elt_current  += n_elt;
  *n_face_current += n_elt * n_face_elt;
}

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
 const PDM_g_num_t          *elmt_ln_to_gn,
       int                  *elmt_face_vtx_idx,
       PDM_g_num_t          *elmt_face_vtx,
       int                  *elmt_cell_face_idx,
       PDM_g_num_t          *elmt_face_cell,
       int                  *parent_elmt_position
)
{
  switch (t_elt) {
   case PDM_MESH_NODAL_POINT:
     abort();
     break;
   case PDM_MESH_NODAL_BAR2:
   case PDM_MESH_NODAL_BARHO:
   case PDM_MESH_NODAL_BARHO_BEZIER:
     abort();
     break;
   case PDM_MESH_NODAL_TRIA3:
   case PDM_MESH_NODAL_TRIAHO:
    case PDM_MESH_NODAL_TRIAHO_BEZIER:
     // PDM_tri_decomposes_faces(n_elt,
     //                          order,
     //                          parent_node,
     //                          n_elt_current,
     //                          n_dface_current,
     //                          beg_gnum_elt_current,
     //                          beg_gnum_face_current,
     //                          connectivity_elmt_vtx,
     //                          elmt_face_vtx_idx,
     //                          elmt_face_vtx,
     //                          elmt_face_cell,
     //                          elmt_cell_face_idx,
     //                          elmt_cell_face,
     //                          parent_elmt_position);
     break;
   case PDM_MESH_NODAL_QUAD4:
   case PDM_MESH_NODAL_QUADHO:
     // PDM_quad_decomposes_faces(n_elt,
     //                           order,
     //                           parent_node,
     //                           n_elt_current,
     //                           n_dface_current,
     //                           beg_gnum_elt_current,
     //                           beg_gnum_face_current,
     //                           connectivity_elmt_vtx,
     //                           elmt_face_vtx_idx,
     //                           elmt_face_vtx,
     //                           elmt_face_cell,
     //                           elmt_cell_face_idx,
     //                           elmt_cell_face,
     //                           parent_elmt_position);
     break;
   case PDM_MESH_NODAL_TETRA4:
   case PDM_MESH_NODAL_TETRAHO:
     PDM_part_mesh_nodal_tetra_decomposes_faces(n_elt,
                                                order,
                                                parent_node,
                                                n_elt_current,
                                                n_face_current,
                                                vtx_ln_to_gn,
                                                connectivity_elmt_vtx,
                                                elmt_ln_to_gn,
                                                elmt_face_vtx_idx,
                                                elmt_face_vtx,
                                                elmt_cell_face_idx,
                                                elmt_face_cell,
                                                parent_elmt_position);
     break;
   case PDM_MESH_NODAL_PYRAMID5:
   case PDM_MESH_NODAL_PYRAMIDHO:
     // PDM_pyra_decomposes_faces(n_elt,
     //                           order,
     //                           parent_node,
     //                           n_elt_current,
     //                           n_dface_current,
     //                           beg_gnum_elt_current,
     //                           beg_gnum_face_current,
     //                           connectivity_elmt_vtx,
     //                           elmt_face_vtx_idx,
     //                           elmt_face_vtx,
     //                           elmt_face_cell,
     //                           elmt_cell_face_idx,
     //                           elmt_cell_face,
     //                           parent_elmt_position);
     break;
   case PDM_MESH_NODAL_PRISM6:
   case PDM_MESH_NODAL_PRISMHO:
     // PDM_prism_decomposes_faces(n_elt,
     //                            order,
     //                            parent_node,
     //                            n_elt_current,
     //                            n_dface_current,
     //                            beg_gnum_elt_current,
     //                            beg_gnum_face_current,
     //                            connectivity_elmt_vtx,
     //                            elmt_face_vtx_idx,
     //                            elmt_face_vtx,
     //                            elmt_face_cell,
     //                            elmt_cell_face_idx,
     //                            elmt_cell_face,
     //                            parent_elmt_position);
     break;
   case PDM_MESH_NODAL_HEXA8:
   case PDM_MESH_NODAL_HEXAHO:
     PDM_part_mesh_nodal_hexa_decomposes_faces(n_elt,
                                                order,
                                                parent_node,
                                                n_elt_current,
                                                n_face_current,
                                                vtx_ln_to_gn,
                                                connectivity_elmt_vtx,
                                                elmt_ln_to_gn,
                                                elmt_face_vtx_idx,
                                                elmt_face_vtx,
                                                elmt_cell_face_idx,
                                                elmt_face_cell,
                                                parent_elmt_position);
     break;
   default:
     PDM_error(__FILE__, __LINE__, 0, "Error PDM_sections_decompose_faces : Element type is supported\n");
  }
}





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
)
{

  int  n_section  = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *section_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  int n_elt_current  = 0;
  int n_face_current = 0;

  // int parent_node[8];


  for(int i_part = 0; i_part < pmne->n_part; ++i_part) {

    for (int i_section = 0; i_section < n_section; i_section++) {
      int id_section = section_id[i_section];
      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne, id_section);

      switch (t_elt) {
        case PDM_MESH_NODAL_POINT:
        case PDM_MESH_NODAL_BAR2:
        case PDM_MESH_NODAL_TRIA3:
        case PDM_MESH_NODAL_QUAD4:
        case PDM_MESH_NODAL_TETRA4:
        case PDM_MESH_NODAL_PYRAMID5:
        case PDM_MESH_NODAL_PRISM6:
        case PDM_MESH_NODAL_HEXA8:
        {
          int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne, id_section, i_part);

          int         *connec              = NULL;
          PDM_g_num_t *numabs              = NULL;
          int         *parent_num          = NULL;
          PDM_g_num_t *parent_entity_g_num = NULL;
          PDM_part_mesh_nodal_elmts_section_std_get(pmne,
                                                    id_section,
                                                    i_part,
                                                    &connec,
                                                    &numabs,
                                                    &parent_num,
                                                    &parent_entity_g_num);

          PDM_part_mesh_nodal_std_decomposes_faces(t_elt,
                                                   n_elt,
                                                   1,
                                                   NULL,
                                                   &n_elt_current,
                                                   &n_face_current,
                                                   vtx_ln_to_gn[i_part],
                                                   connec,
                                                   numabs,
                                                   elmt_face_vtx_idx,
                                                   elmt_face_vtx,
                                                   elmt_cell_face_idx,
                                                   elmt_face_cell,
                                                   parent_elmt_position);

          break;
        }
        case PDM_MESH_NODAL_BARHO:
        case PDM_MESH_NODAL_BARHO_BEZIER:
        case PDM_MESH_NODAL_TRIAHO:
        case PDM_MESH_NODAL_TRIAHO_BEZIER:
        case PDM_MESH_NODAL_QUADHO:
        case PDM_MESH_NODAL_TETRAHO:
        case PDM_MESH_NODAL_PYRAMIDHO:
        case PDM_MESH_NODAL_PRISMHO:
        case PDM_MESH_NODAL_HEXAHO:
        {
        break;
        }
        case PDM_MESH_NODAL_POLY_2D:
        {
          PDM_error(__FILE__, __LINE__, 0, "Error PDM_part_mesh_nodal_elmts_sections_decompose_faces : Element type is supported\n");
          break;
        }

        case PDM_MESH_NODAL_POLY_3D:
        {
          PDM_error(__FILE__, __LINE__, 0, "Error PDM_part_mesh_nodal_elmts_sections_decompose_faces : Element type is supported\n");
          break;
        }

        default:
          PDM_error(__FILE__, __LINE__, 0, "Error PDM_part_mesh_nodal_elmts_sections_decompose_faces : Element type is supported\n");
      }
    }
  }
}


void
PDM_part_mesh_nodal_elmts_decompose_faces_get_size
(
 PDM_part_mesh_nodal_elmts_t *pmne,
 int                         *n_elt_tot,
 int                         *n_face_elt_tot,
 int                         *n_sum_vtx_face_tot
)
{
  /* Get current structure to treat */
  *n_face_elt_tot     = 0;
  *n_sum_vtx_face_tot = 0;

  for(int i_part = 0; i_part < pmne->n_part; ++i_part) {

    for (int i_section = 0; i_section < pmne->n_section_std; i_section++) {

      int n_face_elt     = PDM_n_face_elt_per_elmt    (pmne->sections_std[i_section]->t_elt);
      int n_sum_vtx_face = PDM_n_sum_vtx_face_per_elmt(pmne->sections_std[i_section]->t_elt);

      *n_elt_tot          += pmne->sections_std[i_section]->n_elt[i_part];
      *n_face_elt_tot     += pmne->sections_std[i_section]->n_elt[i_part] * n_face_elt;
      *n_sum_vtx_face_tot += pmne->sections_std[i_section]->n_elt[i_part] * n_sum_vtx_face;

    }

    for (int i_section = 0; i_section < pmne->n_section_poly3d; i_section++) {
      int _n_face      = pmne->sections_poly3d[i_section]->n_face[i_part];
      *n_face_elt_tot += _n_face;
      int n_face_vtx   = pmne->sections_poly3d[i_section]->_facvtx_idx[i_part][_n_face];
      *n_sum_vtx_face_tot += n_face_vtx;

      *n_elt_tot          += pmne->sections_poly3d[i_section]->n_elt[i_part];
    }

    for (int i_section = 0; i_section < pmne->n_section_poly2d; i_section++) {
      int _n_face      = pmne->sections_poly2d[i_section]->n_elt[i_part];
      *n_face_elt_tot += _n_face;
      int n_edge_vtx   = pmne->sections_poly2d[i_section]->_connec_idx[i_part][_n_face];
      *n_sum_vtx_face_tot += n_edge_vtx;
      *n_elt_tot          += pmne->sections_poly2d[i_section]->n_elt[i_part];
    }

    assert(pmne->n_section_poly3d == 0); // Not implemented
    assert(pmne->n_section_poly2d == 0); // Not implemented
  }

  // printf("n_face_elt_tot     ::%i\n", *n_face_elt_tot   );
  // printf("n_sum_vtx_face_tot::%i\n" , *n_sum_vtx_face_tot);
}


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
