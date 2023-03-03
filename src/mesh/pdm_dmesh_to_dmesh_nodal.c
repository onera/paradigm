
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
#include "pdm_unique.h"
#include "pdm_quick_sort.h"
#include "pdm_para_graph_dual.h"
#include "pdm_array.h"
#include "pdm_order.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_part_mesh_nodal_to_pmesh.h"
#include "pdm_dmesh_to_dmesh_nodal.h"
#include "pdm_dmesh_to_dmesh_nodal_priv.h"
#include "pdm_sort.h"
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

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

static
PDM_dmesh_nodal_t*
_dmesh_to_dmesh_nodal
(
 PDM_MPI_Comm    comm,
 PDM_g_num_t    *distrib_cell,
 PDM_g_num_t    *distrib_face,
 PDM_g_num_t    *distrib_edge,
 PDM_g_num_t    *distrib_vtx,
 PDM_g_num_t    *dcell_face,
 int            *dcell_face_idx,
 PDM_g_num_t    *dface_edge,
 int            *dface_edge_idx,
 PDM_g_num_t    *dface_vtx,
 int            *dface_vtx_idx,
 PDM_g_num_t    *dedge_vtx,
 int            *n_bound,
 int           **dbound_idx,
 PDM_g_num_t   **dbound
)
{




}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_dmesh_to_dmesh_nodal_t*
PDM_dmesh_to_dmesh_nodal_create
(
 const int             n_mesh,
 const PDM_MPI_Comm    comm
)
{

  PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn = (PDM_dmesh_to_dmesh_nodal_t *) malloc(sizeof(PDM_dmesh_to_dmesh_nodal_t));

  dm_to_dmn->comm              = comm;
  dm_to_dmn->results_is_getted = PDM_FALSE;
  dm_to_dmn->n_mesh            = n_mesh;

  dm_to_dmn->dcell_face     = malloc(n_mesh * sizeof(PDM_g_num_t *));
  dm_to_dmn->dcell_face_idx = malloc(n_mesh * sizeof(int         *));
  dm_to_dmn->dface_edge     = malloc(n_mesh * sizeof(PDM_g_num_t *));
  dm_to_dmn->dface_edge_idx = malloc(n_mesh * sizeof(int         *));
  dm_to_dmn->dface_vtx      = malloc(n_mesh * sizeof(PDM_g_num_t *));
  dm_to_dmn->dface_vtx_idx  = malloc(n_mesh * sizeof(int         *));
  dm_to_dmn->dedge_vtx      = malloc(n_mesh * sizeof(PDM_g_num_t *));

  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {
    dm_to_dmn->dcell_face    [i_mesh] = NULL;
    dm_to_dmn->dcell_face_idx[i_mesh] = NULL;
    dm_to_dmn->dface_edge    [i_mesh] = NULL;
    dm_to_dmn->dface_edge_idx[i_mesh] = NULL;
    dm_to_dmn->dface_vtx     [i_mesh] = NULL;
    dm_to_dmn->dface_vtx_idx [i_mesh] = NULL;
    dm_to_dmn->dedge_vtx     [i_mesh] = NULL;
  }

  dm_to_dmn->distrib_cell = malloc(n_mesh * sizeof(PDM_g_num_t *));
  dm_to_dmn->distrib_face = malloc(n_mesh * sizeof(PDM_g_num_t *));
  dm_to_dmn->distrib_edge = malloc(n_mesh * sizeof(PDM_g_num_t *));
  dm_to_dmn->distrib_vtx  = malloc(n_mesh * sizeof(PDM_g_num_t *));

  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {
    dm_to_dmn->distrib_cell[i_mesh] = NULL;
    dm_to_dmn->distrib_face[i_mesh] = NULL;
    dm_to_dmn->distrib_edge[i_mesh] = NULL;
    dm_to_dmn->distrib_vtx [i_mesh] = NULL;
  }

  dm_to_dmn->n_bound         = malloc( n_mesh * sizeof(int          *) );
  dm_to_dmn->dbound          = malloc( n_mesh * sizeof(PDM_g_num_t **) );
  dm_to_dmn->dbound_idx      = malloc( n_mesh * sizeof(int         **) );

  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {

    dm_to_dmn->n_bound   [i_mesh] = malloc( PDM_BOUND_TYPE_MAX * sizeof(int          ) );
    dm_to_dmn->dbound    [i_mesh] = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_g_num_t *) );
    dm_to_dmn->dbound_idx[i_mesh] = malloc( PDM_BOUND_TYPE_MAX * sizeof(int         *) );

    for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i ) {
      dm_to_dmn->n_bound   [i_mesh][i] = 0;
      dm_to_dmn->dbound    [i_mesh][i] = NULL;
      dm_to_dmn->dbound_idx[i_mesh][i] = NULL;
    }
  }

  dm_to_dmn->dmn = malloc( n_mesh * sizeof(PDM_dmesh_nodal_t *) );
  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {
    dm_to_dmn->dmn[i_mesh] = NULL;
  }

  return dm_to_dmn;
}

void
PDM_dmesh_to_dmesh_nodal_compute
(
 PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn
)
{
  for(int i_mesh = 0; i_mesh < dm_to_dmn->n_mesh; ++i_mesh) {
    dm_to_dmn->dmn[i_mesh] = _dmesh_to_dmesh_nodal(dm_to_dmn->comm,
                                                   dm_to_dmn->distrib_cell  [i_mesh],
                                                   dm_to_dmn->distrib_face  [i_mesh],
                                                   dm_to_dmn->distrib_edge  [i_mesh],
                                                   dm_to_dmn->distrib_vtx   [i_mesh],
                                                   dm_to_dmn->dcell_face    [i_mesh],
                                                   dm_to_dmn->dcell_face_idx[i_mesh],
                                                   dm_to_dmn->dface_edge    [i_mesh],
                                                   dm_to_dmn->dface_edge_idx[i_mesh],
                                                   dm_to_dmn->dface_vtx     [i_mesh],
                                                   dm_to_dmn->dface_vtx_idx [i_mesh],
                                                   dm_to_dmn->dedge_vtx     [i_mesh],
                                                   dm_to_dmn->n_bound       [i_mesh],
                                                   dm_to_dmn->dbound_idx    [i_mesh],
                                                   dm_to_dmn->dbound        [i_mesh]);
  }


}

void
PDM_dmesh_to_dmesh_nodal_set_dmesh
(
        PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn,
  const int                         i_mesh,
        PDM_dmesh_t                *dm
)
{
  PDM_UNUSED(dm_to_dmn);
  PDM_UNUSED(i_mesh);
  PDM_UNUSED(dm);
  abort();
}


void
PDM_dmesh_to_dmesh_nodal_distribution_set
(
        PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn,
  const int                         i_mesh,
        PDM_g_num_t                *distrib_cell,
        PDM_g_num_t                *distrib_face,
        PDM_g_num_t                *distrib_edge,
        PDM_g_num_t                *distrib_vtx
)
{
  dm_to_dmn->distrib_cell[i_mesh] = distrib_cell;
  dm_to_dmn->distrib_face[i_mesh] = distrib_face;
  dm_to_dmn->distrib_edge[i_mesh] = distrib_edge;
  dm_to_dmn->distrib_vtx [i_mesh] = distrib_vtx;
}

void
PDM_dmesh_to_dmesh_nodal_connectivity_set
(
        PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn,
  const int                         i_mesh,
        int                        *dcell_face_idx,
        PDM_g_num_t                *dcell_face,
        PDM_g_num_t                *dface_edge,
        int                        *dface_edge_idx,
        PDM_g_num_t                *dedge_vtx,
        int                        *dface_vtx_idx,
        PDM_g_num_t                *dface_vtx
)
{
  dm_to_dmn->dcell_face_idx[i_mesh] = dcell_face_idx;
  dm_to_dmn->dcell_face    [i_mesh] = dcell_face;
  dm_to_dmn->dface_edge    [i_mesh] = dface_edge;
  dm_to_dmn->dface_edge_idx[i_mesh] = dface_edge_idx;
  dm_to_dmn->dedge_vtx     [i_mesh] = dedge_vtx;
  dm_to_dmn->dface_vtx_idx [i_mesh] = dface_vtx_idx;
  dm_to_dmn->dface_vtx     [i_mesh] = dface_vtx;
}

void
PDM_dmesh_to_dmesh_nodal_group_set
(
        PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn,
  const int                         i_mesh,
        PDM_bound_type_t            bound_type,
        int                         n_group,
        int                        *dbound_idx,
        PDM_g_num_t                *dbound
)
{
  dm_to_dmn->n_bound   [i_mesh][bound_type] = n_group;
  dm_to_dmn->dbound_idx[i_mesh][bound_type] = dbound_idx;
  dm_to_dmn->dbound    [i_mesh][bound_type] = dbound;
}


void
PDM_dmesh_to_dmesh_nodal_free
(
 PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn
)
{
  if(dm_to_dmn == NULL) {
    return;
  }

  free(dm_to_dmn->dcell_face    );
  free(dm_to_dmn->dcell_face_idx);
  free(dm_to_dmn->dface_edge    );
  free(dm_to_dmn->dface_edge_idx);
  free(dm_to_dmn->dface_vtx     );
  free(dm_to_dmn->dface_vtx_idx );
  free(dm_to_dmn->dedge_vtx     );

  free(dm_to_dmn->distrib_cell);
  free(dm_to_dmn->distrib_face);
  free(dm_to_dmn->distrib_edge);
  free(dm_to_dmn->distrib_vtx );

  for(int i = 0; i < dm_to_dmn->n_mesh; ++i ) {
    free(dm_to_dmn->n_bound   [i]);
    free(dm_to_dmn->dbound    [i]);
    free(dm_to_dmn->dbound_idx[i]);
  }

  free(dm_to_dmn->n_bound      );
  free(dm_to_dmn->dbound       );
  free(dm_to_dmn->dbound_idx   );

  if(dm_to_dmn->results_is_getted == PDM_FALSE) {
    for(int i_mesh = 0; i_mesh < dm_to_dmn->n_mesh; ++i_mesh) {
      PDM_DMesh_nodal_free(dm_to_dmn->dmn[i_mesh]);
    }
  }

  free(dm_to_dmn->dmn);

  free(dm_to_dmn);
  dm_to_dmn = NULL;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
