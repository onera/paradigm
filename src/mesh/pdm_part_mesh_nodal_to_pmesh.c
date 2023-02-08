
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

  int n_face_elt_vol_tot     = 0;
  int n_sum_vtx_vol_face_tot = 0;

  int n_face_elt_surf_tot     = 0;
  int n_sum_vtx_surf_face_tot = 0;

  PDM_part_mesh_nodal_elmts_decompose_faces_get_size(pmn->volumic,
                                                     &n_face_elt_vol_tot,
                                                     &n_sum_vtx_vol_face_tot );

  PDM_part_mesh_nodal_elmts_decompose_faces_get_size(pmn->surfacic,
                                                     &n_face_elt_surf_tot,
                                                     &n_sum_vtx_surf_face_tot);

  return NULL;
}




#ifdef __cplusplus
}
#endif /* __cplusplus */
