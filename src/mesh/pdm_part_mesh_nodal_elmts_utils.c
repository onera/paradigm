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
)
{

  for(int i_part = 0; i_part < pmne->n_part; ++i_part) {



  }
}


void
PDM_part_mesh_nodal_elmts_decompose_faces_get_size
(
 PDM_part_mesh_nodal_elmts_t *pmne,
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

      *n_face_elt_tot     += pmne->sections_std[i_section]->n_elt[i_part] * n_face_elt;
      *n_sum_vtx_face_tot += pmne->sections_std[i_section]->n_elt[i_part] * n_sum_vtx_face;

    }

    for (int i_section = 0; i_section < pmne->n_section_poly3d; i_section++) {
      int _n_face      = pmne->sections_poly3d[i_section]->n_face[i_part];
      *n_face_elt_tot += _n_face;
      int n_face_vtx   = pmne->sections_poly3d[i_section]->_facvtx_idx[i_part][_n_face];
      *n_sum_vtx_face_tot += n_face_vtx;
    }

    for (int i_section = 0; i_section < pmne->n_section_poly2d; i_section++) {
      int _n_face      = pmne->sections_poly2d[i_section]->n_elt[i_part];
      *n_face_elt_tot += _n_face;
      int n_edge_vtx   = pmne->sections_poly2d[i_section]->_connec_idx[i_part][_n_face];
      *n_sum_vtx_face_tot += n_edge_vtx;
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
