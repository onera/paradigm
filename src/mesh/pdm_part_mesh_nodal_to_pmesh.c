
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
#include "pdm_error.h"
#include "pdm_binary_search.h"
#include "pdm_dmesh_nodal_to_dmesh_priv.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_elmts_utils.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_logging.h"
#include "pdm_dmesh.h"
#include "pdm_gnum.h"
#include "pdm_distrib.h"
#include "pdm_quick_sort.h"
#include "pdm_para_graph_dual.h"
#include "pdm_array.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_part_mesh_nodal_to_pmesh.h"
#include "pdm_sort.h"

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

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/


PDM_part_mesh_t*
PDM_part_mesh_nodal_to_part_mesh
(
        PDM_part_mesh_nodal_t                      *pmn,
  const PDM_dmesh_nodal_to_dmesh_transform_t        transform_kind,
  const PDM_dmesh_nodal_to_dmesh_translate_group_t  transform_group_kind
)
{
  // PDM_part_mesh_nodal_elmts_t* pmne_vol = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn, PDM_GEOMETRY_KIND_VOLUMIC);
  PDM_UNUSED(transform_kind);
  PDM_UNUSED(transform_group_kind);

  int n_elmt_vol_tot         = 0;
  int n_face_elt_vol_tot     = 0;
  int n_sum_vtx_vol_face_tot = 0;

  int n_elmt_surf_tot         = 0;
  int n_face_elt_surf_tot     = 0;
  int n_sum_vtx_surf_face_tot = 0;

  PDM_part_mesh_nodal_elmts_decompose_faces_get_size(pmn->volumic,
                                                     &n_elmt_vol_tot,
                                                     &n_face_elt_vol_tot,
                                                     &n_sum_vtx_vol_face_tot );

  PDM_part_mesh_nodal_elmts_decompose_faces_get_size(pmn->surfacic,
                                                     &n_elmt_surf_tot,
                                                     &n_face_elt_surf_tot,
                                                     &n_sum_vtx_surf_face_tot);

  PDM_g_num_t **vtx_ln_to_gn = malloc(sizeof(PDM_g_num_t *));
  for(int i_part = 0; i_part < pmn->n_part; ++i_part) {
    vtx_ln_to_gn[i_part] = PDM_part_mesh_nodal_vtx_g_num_get(pmn, i_part);
  }

  int         *elmt_face_vtx_idx      = malloc((n_face_elt_vol_tot+1) * sizeof(int        ));
  PDM_g_num_t *elmt_face_vtx          = malloc(n_sum_vtx_vol_face_tot * sizeof(PDM_g_num_t));
  int         *elmt_cell_face_vtx_idx = malloc((n_elmt_vol_tot+1)     * sizeof(int        ));
  int         *parent_elmt_position   = malloc(n_face_elt_vol_tot     * sizeof(int        ));

  printf("n_face_elt_vol_tot     : %i\n", n_face_elt_vol_tot    );
  printf("n_sum_vtx_vol_face_tot : %i\n", n_sum_vtx_vol_face_tot);
  printf("n_elmt_vol_tot         : %i\n", n_elmt_vol_tot        );
  printf("n_face_elt_vol_tot     : %i\n", n_face_elt_vol_tot    );

  elmt_face_vtx_idx     [0] = 0;
  elmt_cell_face_vtx_idx[0] = 0;
  PDM_part_mesh_nodal_elmts_sections_decompose_faces(pmn->volumic,
                                                     vtx_ln_to_gn,
                                                     elmt_face_vtx_idx,
                                                     elmt_face_vtx,
                                                     elmt_cell_face_vtx_idx,
                                                     parent_elmt_position);

  // PDM_log_trace_array_int(elmt_face_vtx_idx, n_face_elt_vol_tot, "elmt_face_vtx_idx ::");


  free(vtx_ln_to_gn);
  free(elmt_face_vtx_idx     );
  free(elmt_face_vtx         );
  free(elmt_cell_face_vtx_idx);
  free(parent_elmt_position  );

  return NULL;
}




#ifdef __cplusplus
}
#endif /* __cplusplus */
