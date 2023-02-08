
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
#include "pdm_order.h"
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
  int *send_n    = malloc( n_rank              * sizeof(int));
  int *send_idx  = malloc((n_rank+1)           * sizeof(int));
  int *recv_n    = malloc( n_rank              * sizeof(int));
  int *recv_idx  = malloc((n_rank+1)           * sizeof(int));
  int *dest_rank = malloc(n_face_elt_vol_tot   * sizeof(int));
  for(int i = 0; i < n_rank; ++i) {
    send_n[i] = 0;
  }

  for(int i_face = 0; i_face < n_face_elt_vol_tot; ++i_face) {
    PDM_g_num_t key = key_ln_to_gn[i_face];
    int t_rank = PDM_binary_search_gap_long(key-1, distrib_key, n_rank+1);
    dest_rank[i_face] = t_rank;
    send_n[t_rank]++;
  }

  PDM_MPI_Alltoall(send_n, 1, PDM_MPI_INT,
                   recv_n, 1, PDM_MPI_INT, pmn->comm);

  PDM_log_trace_array_int(send_n  , n_rank  , "send_n   ::");

  int *send_s_face_vtx_n   = malloc( n_rank * sizeof(int));
  int *recv_s_face_vtx_n   = malloc( n_rank * sizeof(int));
  int *send_s_face_vtx_idx = malloc((n_rank+1) * sizeof(int));
  int *recv_s_face_vtx_idx = malloc((n_rank+1) * sizeof(int));
  send_idx[0] = 0;
  recv_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    send_idx         [i+1] = send_idx[i] + send_n[i];
    recv_idx         [i+1] = recv_idx[i] + recv_n[i];
    send_n           [i] = 0;
    send_s_face_vtx_n[i] = 0;
  }

  PDM_log_trace_array_int(send_idx, n_rank+1, "send_idx ::");
  PDM_log_trace_array_int(recv_n  , n_rank  , "recv_n   ::");
  PDM_log_trace_array_int(recv_idx, n_rank+1, "recv_idx ::");

  /*
   * Prepare send
   */
  int         *send_face_vtx_n   = malloc( n_face_elt_vol_tot * sizeof(int        ));
  PDM_g_num_t *send_face_key     = malloc( n_face_elt_vol_tot * sizeof(PDM_g_num_t));
  for(int i_face = 0; i_face < n_face_elt_vol_tot; ++i_face) {
    int t_rank = dest_rank[i_face];
    int idx_write = send_idx[t_rank] + send_n[t_rank]++;

    send_face_vtx_n[idx_write] = elmt_face_vtx_idx[i_face+1]-elmt_face_vtx_idx[i_face];
    send_face_key  [idx_write] = key_ln_to_gn[i_face];
    send_s_face_vtx_n[t_rank] += elmt_face_vtx_idx[i_face+1]-elmt_face_vtx_idx[i_face];

  }

  int         *recv_face_vtx_n = malloc(recv_idx[n_rank] * sizeof(int        ));
  PDM_g_num_t *recv_face_key   = malloc(recv_idx[n_rank] * sizeof(PDM_g_num_t));

  PDM_log_trace_array_int(send_face_vtx_n, n_face_elt_vol_tot, "send_face_vtx_n ::");
  PDM_MPI_Alltoallv(send_face_vtx_n,
                    send_n,
                    send_idx,
                    PDM_MPI_INT,
                    recv_face_vtx_n,
                    recv_n,
                    recv_idx,
                    PDM_MPI_INT,
                    pmn->comm);
  free(send_face_vtx_n);

  PDM_MPI_Alltoallv(send_face_key,
                    send_n,
                    send_idx,
                    PDM_MPI_INT,
                    recv_face_key,
                    recv_n,
                    recv_idx,
                    PDM_MPI_INT,
                    pmn->comm);
  free(send_face_key);

  /*
   * Exchange size of connectivity
   */
  PDM_MPI_Alltoall(send_s_face_vtx_n, 1, PDM_MPI_INT,
                   recv_s_face_vtx_n, 1, PDM_MPI_INT, pmn->comm);

  send_s_face_vtx_idx[0] = 0;
  recv_s_face_vtx_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    send_s_face_vtx_idx[i+1] = send_s_face_vtx_idx[i] + send_s_face_vtx_n[i];
    recv_s_face_vtx_idx[i+1] = recv_s_face_vtx_idx[i] + recv_s_face_vtx_n[i];
    send_s_face_vtx_n  [i] = 0;
  }

  /*
   * Exchange of connectivity
   */
  PDM_g_num_t* send_face_vtx = malloc(send_s_face_vtx_idx[n_rank] * sizeof(PDM_g_num_t));
  PDM_g_num_t* recv_face_vtx = malloc(recv_s_face_vtx_idx[n_rank] * sizeof(PDM_g_num_t));
  for(int i_face = 0; i_face < n_face_elt_vol_tot; ++i_face) {
    int t_rank = dest_rank[i_face];
    int idx_write = send_s_face_vtx_idx[t_rank];

    for(int k = elmt_face_vtx_idx[i_face]; k < elmt_face_vtx_idx[i_face+1]; ++k) {
      send_face_vtx[idx_write+send_s_face_vtx_n[t_rank]++] = elmt_face_vtx[k];
    }
  }


  PDM_MPI_Alltoallv(send_face_vtx,
                    send_s_face_vtx_n,
                    send_s_face_vtx_idx,
                    PDM_MPI_INT,
                    recv_face_vtx,
                    recv_s_face_vtx_n,
                    recv_s_face_vtx_idx,
                    PDM_MPI_INT,
                    pmn->comm);

  PDM_log_trace_array_int (recv_face_vtx_n, recv_idx[n_rank]           , "recv_face_vtx_n ::");
  PDM_log_trace_array_long(recv_face_key  , recv_idx[n_rank]           , "recv_face_key   ::");
  PDM_log_trace_array_long(recv_face_vtx  , recv_s_face_vtx_idx[n_rank], "recv_face_vtx   ::");

  int *recv_face_vtx_idx = malloc((recv_idx[n_rank]+1) * sizeof(int));
  recv_face_vtx_idx[0] = 0;
  int n_max_vtx = 0;
  for(int i = 0; i < recv_idx[n_rank]; ++i) {
    recv_face_vtx_idx[i+1] = recv_face_vtx_idx[i] + recv_face_vtx_n[i];
    n_max_vtx = PDM_MAX(n_max_vtx, recv_face_vtx_n[i]);
  }

  /*
   * All data are exchange we need to order the key and resolve conflit in hash table
   */
  int n_recv_key = recv_idx[n_rank];
  int *order = malloc(n_recv_key * sizeof(int));
  PDM_order_gnum_s(recv_face_key, 1, order, n_recv_key);

  int n_conflit_to_solve = 0;
  PDM_g_num_t last_gnum = -1;

  int *key_conflict_idx = malloc((n_recv_key+1) * sizeof(int));
  key_conflict_idx[0] = 0;
  for(int i = 0; i < n_recv_key; ++i) {
    if(recv_face_key[order[i]] != last_gnum){
      key_conflict_idx[n_conflit_to_solve+1] = key_conflict_idx[n_conflit_to_solve]+1;
      n_conflit_to_solve++;
      last_gnum = recv_face_key[order[i]];
    } else {
      key_conflict_idx[n_conflit_to_solve]++;
    }
  }

  PDM_log_trace_array_int(key_conflict_idx, n_conflit_to_solve, "key_conflict_idx ::  ");


  int n_max_entity_per_key = 0;
  for(int i = 0; i < n_conflit_to_solve; ++i) {
    n_max_entity_per_key = PDM_MAX(n_max_entity_per_key, key_conflict_idx[i+1]-key_conflict_idx[i]);
  }


  /*
   * Solve conflict
   */
  if(1 == 1) {
    for(int i = 0; i < n_conflit_to_solve; ++i) {
      log_trace(" ------ i = %i \n", i);
      for(int i_key = key_conflict_idx[i]; i_key < key_conflict_idx[i+1]; ++i_key) {
        int i_conflict = order[i_key];
        int beg = recv_face_vtx_idx[i_conflict];
        int n_vtx_in_face = recv_face_vtx_idx[i_conflict+1] - beg;
        log_trace(" \t i_key = %i | beg = %i / n = %i \n", recv_face_key[i_conflict], beg, n_vtx_in_face);
      }
    }
  }


  PDM_g_num_t* loc_entity_vtx_1 = (PDM_g_num_t *) malloc(  n_max_vtx               * sizeof(PDM_g_num_t) );
  PDM_g_num_t* loc_entity_vtx_2 = (PDM_g_num_t *) malloc(  n_max_vtx               * sizeof(PDM_g_num_t) );
  int*         already_treat    = (int         *) malloc(  n_max_entity_per_key    * sizeof(int        ) );
  int*         same_entity_idx  = (int         *) malloc( (n_max_entity_per_key+1) * sizeof(int        ) );
  int*         sens_entity      = (int         *) malloc(  n_max_entity_per_key    * sizeof(int        ) );
  PDM_g_num_t *tmp_parent       = (PDM_g_num_t *) malloc(n_max_entity_per_key      * sizeof(PDM_g_num_t) );
  int         *order_parent     = (int         *) malloc(n_max_entity_per_key      * sizeof(int        ) );

  for(int i = 0; i < n_conflit_to_solve; ++i) {

    int n_conflict_entitys = key_conflict_idx[i+1] - key_conflict_idx[i];

    for(int j = 0; j < n_conflict_entitys; ++j ) {
      already_treat[j] = -1;
      order_parent [j] = key_conflict_idx[i]+j;
    }

    for(int idx_entity = key_conflict_idx[i]; idx_entity < key_conflict_idx[i+1]; ++idx_entity) {
      // int i_entity = order[idx_entity];
      int i_entity = order_parent[idx_entity-key_conflict_idx[i]];

      int beg_1 = recv_face_vtx_idx[i_entity];
      int n_vtx_in_entity1 = recv_face_vtx_idx[i_entity+1] - beg_1;

      int idx_next_same_entity = 0;
      sens_entity    [idx_next_same_entity] = 1;
      same_entity_idx[idx_next_same_entity++] = i_entity;

      /* Only if not treated we search correspondance with other */
      if(already_treat[i_entity] != 1) {

        PDM_g_num_t key_1 = 0;
        int idx_min_1 = -1;
        PDM_g_num_t min_1 = max_vtx_gnum+1;
        for(int j = 0; j < n_vtx_in_entity1; ++j) {
          loc_entity_vtx_1[j] = recv_face_vtx[beg_1+j];
          key_1 += loc_entity_vtx_1[j];
          if(loc_entity_vtx_1[j] < min_1) {
            min_1 = loc_entity_vtx_1[j];
            idx_min_1 = j;
          };
        }
        PDM_quick_sort_long(loc_entity_vtx_1, 0, n_vtx_in_entity1-1);

        for(int idx_entity2 = key_conflict_idx[i]; idx_entity2 < key_conflict_idx[i+1]; ++idx_entity2) {
          int i_entity_next = order_parent[idx_entity2-key_conflict_idx[i]];

          if (i_entity_next == i_entity) {
            continue;
          }

          printf("conflict : i_entity = %d, i_entity_next = %d...\n", i_entity, i_entity_next);
          if(already_treat[i_entity_next] == 1) {
            continue;
          }

          int beg_2 = recv_face_vtx_idx[i_entity_next];
          int n_vtx_in_entity2 = recv_face_vtx_idx[i_entity_next+1] - beg_2;
          if(n_vtx_in_entity1 == n_vtx_in_entity2 ) {

            PDM_g_num_t key_2 = 0;
            int idx_min_2 = -1;
            PDM_g_num_t min_2 = max_vtx_gnum+1;
            for(int j = 0; j < n_vtx_in_entity1; ++j) {
              loc_entity_vtx_2[j] = recv_face_vtx[beg_2+j];
              key_2 += loc_entity_vtx_2[j];
              if(loc_entity_vtx_2[j] < min_2) {
                min_2 = loc_entity_vtx_2[j];
                idx_min_2 = j;
              };
            }
            PDM_quick_sort_long(loc_entity_vtx_2, 0, n_vtx_in_entity2-1);

            printf("key_1 : %d, key_2 = %d...\n", key_1, key_2);
            assert(key_1 == key_2);

            int is_same_entity = 1;
            for(int i_vtx = 0; i_vtx < n_vtx_in_entity1; ++i_vtx) {
              if(loc_entity_vtx_1[i_vtx] != loc_entity_vtx_2[i_vtx]) {
                is_same_entity = -1;
                break;
              }
            }
          }
        }

      } /* End if already_treat[i_entity] */
    } /* End conflict */
  }

  PDM_log_trace_array_int (order, n_recv_key , "order ::");

  /*
   * Reconstruction du face_ln_to_gn + cell_face local
   */




  /*
   * Renvoie de la connectivité dface_vtx via block_to_part
   */




  free(loc_entity_vtx_1);
  free(loc_entity_vtx_2);
  free(already_treat   );
  free(same_entity_idx );
  free(sens_entity     );
  free(tmp_parent      );
  free(order_parent    );

  free(key_conflict_idx);
  free(order);

  free(send_n);
  free(send_idx);
  free(recv_n);
  free(recv_idx);
  free(dest_rank);
  free(recv_face_vtx_n);
  free(recv_face_vtx_idx);

  free(send_face_vtx);
  free(recv_face_vtx);
  free(recv_face_key);

  free(send_s_face_vtx_n  );
  free(recv_s_face_vtx_n  );
  free(send_s_face_vtx_idx);
  free(recv_s_face_vtx_idx);

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
