
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
  PDM_g_num_t _max_vtx_gnum = -1;
  for(int i_part = 0; i_part < pmn->n_part; ++i_part) {
    vtx_ln_to_gn[i_part] = PDM_part_mesh_nodal_vtx_g_num_get(pmn, i_part);
    int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
    for(int i_vtx = 0; i_vtx < n_vtx; ++i_vtx) {
      _max_vtx_gnum = PDM_MAX(_max_vtx_gnum, vtx_ln_to_gn[i_part][i_vtx]);
    }
  }
  PDM_g_num_t max_vtx_gnum = 0;
  PDM_MPI_Allreduce(&_max_vtx_gnum, &max_vtx_gnum, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, pmn->comm);


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

  PDM_log_trace_array_int(elmt_face_vtx_idx, n_face_elt_vol_tot, "elmt_face_vtx_idx ::");


  /*
   * Create hash table
   */
  PDM_g_num_t *key_ln_to_gn = malloc(n_face_elt_vol_tot * sizeof(PDM_g_num_t));
  double      *key_weight   = malloc(n_face_elt_vol_tot * sizeof(double     ));
  PDM_g_num_t key_mod = 4 * max_vtx_gnum;
  for(int i_face = 0; i_face < n_face_elt_vol_tot; ++i_face) {
    PDM_g_num_t key = 0;
    for(int idx = elmt_face_vtx_idx[i_face]; idx < elmt_face_vtx_idx[i_face+1]; ++idx) {
      key += elmt_face_vtx[idx];
    }
    // min_vtx =
    key_ln_to_gn[i_face] = key % key_mod + 1;
    key_weight  [i_face] = elmt_face_vtx_idx[i_face+1]-elmt_face_vtx_idx[i_face];
  }


  int n_rank;
  PDM_MPI_Comm_size (pmn->comm, &n_rank);
  int sampling_factor = 2;
  int n_iter_max      = 5;
  double tol = 0.10;
  PDM_g_num_t *distrib_key = NULL;
  PDM_distrib_weight(      sampling_factor,
                           n_rank,
                           1,
                           &n_face_elt_vol_tot,
    (const PDM_g_num_t **) &key_ln_to_gn,
    (const double      **) &key_weight,
                           n_iter_max,
                           tol,
                           pmn->comm,
                           &distrib_key);

  PDM_log_trace_array_int(distrib_key, n_rank+1, "distrib_key ::");

  /*
   * Prepare send
   */
  int *send_n   = malloc( n_rank    * sizeof(int));
  int *send_idx = malloc((n_rank+1) * sizeof(int));
  int *recv_n   = malloc( n_rank    * sizeof(int));
  int *recv_idx = malloc((n_rank+1) * sizeof(int));
  for(int i = 0; i < n_rank; ++i) {
    send_n[i] = 0;
  }

  for(int i_face = 0; i_face < n_face_elt_vol_tot; ++i_face) {
    PDM_g_num_t key = key_ln_to_gn[i_face];
    int t_rank = PDM_binary_search_gap_long(key-1, distrib_key, n_rank+1);
    send_n[t_rank]++;
  }

  PDM_MPI_Alltoall(send_n, 1, PDM_MPI_INT,
                   recv_n, 1, PDM_MPI_INT, pmn->comm);

  PDM_log_trace_array_int(send_n  , n_rank  , "send_n   ::");

  send_idx[0] = 0;
  recv_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    send_idx[i+1] = send_idx[i] + send_n[i];
    recv_idx[i+1] = recv_idx[i] + recv_n[i];
    send_n  [i] = 0;
  }

  PDM_log_trace_array_int(send_idx, n_rank+1, "send_idx ::");
  PDM_log_trace_array_int(recv_n  , n_rank  , "recv_n   ::");
  PDM_log_trace_array_int(recv_idx, n_rank+1, "recv_idx ::");

  /*
   * Prepare send
   */
  int *send_face_vtx_n = malloc(n_face_elt_vol_tot * sizeof(int));
  for(int i_face = 0; i_face < n_face_elt_vol_tot; ++i_face) {
    PDM_g_num_t key = key_ln_to_gn[i_face];
    int t_rank = PDM_binary_search_gap_long(key-1, distrib_key, n_rank+1);
    int idx_write = send_idx[t_rank] + send_n[t_rank]++;
    send_face_vtx_n[idx_write] = elmt_face_vtx_idx[i_face+1]-elmt_face_vtx_idx[i_face];
  }

  PDM_log_trace_array_int(send_face_vtx_n, n_face_elt_vol_tot, "send_face_vtx_n ::");


  free(send_face_vtx_n);
  free(send_n);
  free(send_idx);
  free(recv_n);
  free(recv_idx);

  free(distrib_key);
  free(key_ln_to_gn);
  free(key_weight);
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
