#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_part.h"
#include "pdm_dcube_gen.h"
#include "pdm_logging.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_priv.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_sort.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_dmesh_nodal_to_dmesh.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

static
void
_exchange_point_list
(
 int            n_group_join,
 int           *group_join_to_join_opp,
 int           *dface_join_idx,
 PDM_g_num_t   *dface_join,
 PDM_g_num_t  **dface_join_opp,
 PDM_MPI_Comm   comm
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /*
   * We have for all extraction zone, we need to exchange id
   */
  PDM_g_num_t **distrib_join   = malloc(n_group_join * sizeof(PDM_g_num_t *));
  // PDM_g_num_t **dface_join_opp = malloc(n_group_join * sizeof(PDM_g_num_t *));
  for(int i_group_join = 0; i_group_join < n_group_join; ++i_group_join) {
    int dn_face_join = dface_join_idx[i_group_join+1] - dface_join_idx[i_group_join];
    distrib_join[i_group_join] = PDM_compute_entity_distribution(comm, dn_face_join);
  }

  *dface_join_opp = malloc(dface_join_idx[n_group_join] * sizeof(PDM_g_num_t));
  PDM_g_num_t* _dface_join_opp = *dface_join_opp;

  for(int i_group_join = 0; i_group_join < n_group_join; ++i_group_join) {

    int i_group_join_opp = group_join_to_join_opp[i_group_join];
    int dn_face_join = dface_join_idx[i_group_join+1] - dface_join_idx[i_group_join];
    PDM_g_num_t* distrib_join_cur = distrib_join[i_group_join    ];
    PDM_g_num_t* distrib_join_opp = distrib_join[i_group_join_opp];

    PDM_g_num_t *join_ln_to_gn = malloc(dn_face_join * sizeof(PDM_g_num_t));
    for(int i = 0; i < dn_face_join; ++i) {
      join_ln_to_gn[i] = distrib_join_cur[i_rank] + i + 1;
    }

    /*
     * Exchange
     */
    PDM_g_num_t *blk_dface_join_cur = &dface_join[dface_join_idx[i_group_join    ]];
    PDM_g_num_t *blk_dface_join_opp = &dface_join[dface_join_idx[i_group_join_opp]];

    PDM_block_to_part_t *btp = PDM_block_to_part_create(distrib_join_opp,
                                 (const PDM_g_num_t **) &join_ln_to_gn,
                                                        &dn_face_join,
                                                        1,
                                                        comm);

    int cst_stride = 1;
    PDM_g_num_t* sub_dface_join_opp = &_dface_join_opp[dface_join_idx[i_group_join]];
    PDM_block_to_part_exch(btp,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST,
                           &cst_stride,
                  (void *) blk_dface_join_opp,
                           NULL,
               (void ** ) &sub_dface_join_opp);

    if(1 == 1) {
      PDM_log_trace_array_long(blk_dface_join_cur, dn_face_join, "dface_join_cur :: ");
      PDM_log_trace_array_long(sub_dface_join_opp, dn_face_join, "dface_join_opp :: ");
    }

    PDM_block_to_part_free(btp);
    free(join_ln_to_gn);
  }

  for(int i_group_join = 0; i_group_join < n_group_join; ++i_group_join) {
    free(distrib_join[i_group_join]);
  }
  free(distrib_join);
}






static
void
_deduce_descending_join
(
 int            n_zone,
 int            n_group_join,
 int           *group_join_to_zone_cur,
 int           *group_join_to_zone_opp,
 int           *group_join_to_join_opp,
 int           *dface_join_idx,
 PDM_g_num_t   *dface_join,
 PDM_g_num_t   *dface_join_opp,
 int          **dface_vtx_idx,
 PDM_g_num_t  **dface_vtx,
 PDM_g_num_t  **extract_face_distribution,
 PDM_g_num_t  **extract_vtx_distribution,
 int          **dextract_face_vtx_idx,
 PDM_g_num_t  **dextract_face_vtx,
 int          **dextract_face_group_id,
 int          **dextract_face_group_sens,
 PDM_g_num_t  **dextract_face_join,
 PDM_g_num_t  **dextract_face_join_opp,
 PDM_g_num_t  **dparent_face_g_num,
 PDM_g_num_t  **dparent_vtx_g_num,
 PDM_g_num_t  **pextract_old_to_new,
 double       **dextract_vtx_coord,
 PDM_MPI_Comm   comm
)
{
  PDM_UNUSED(n_zone);
  PDM_UNUSED(dface_join_idx);
  PDM_UNUSED(dface_join);
  PDM_UNUSED(dface_vtx_idx);
  PDM_UNUSED(dface_vtx);
  PDM_UNUSED(extract_face_distribution);
  PDM_UNUSED(extract_vtx_distribution);
  PDM_UNUSED(dextract_face_vtx_idx);
  PDM_UNUSED(dextract_face_vtx);
  PDM_UNUSED(dparent_face_g_num);
  PDM_UNUSED(dparent_vtx_g_num);
  PDM_UNUSED(pextract_old_to_new);
  PDM_UNUSED(dextract_vtx_coord);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  PDM_g_num_t **dedge_distrib  = malloc(n_zone * sizeof(PDM_g_num_t *));
  int         **dedge_vtx_idx  = malloc(n_zone * sizeof(int         *));
  PDM_g_num_t **dedge_vtx      = malloc(n_zone * sizeof(PDM_g_num_t *));
  int         **dedge_face_idx = malloc(n_zone * sizeof(int         *));
  PDM_g_num_t **dedge_face     = malloc(n_zone * sizeof(PDM_g_num_t *));

  PDM_g_num_t **key_ln_to_gn    = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t*));

  /* Connectivity */
  int         **data_send_connect_n  = (int         **) malloc( n_zone * sizeof(int        *));
  int         **data_send_face_n     = (int         **) malloc( n_zone * sizeof(int        *));
  PDM_g_num_t **data_send_connect    = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t*));
  PDM_g_num_t **data_send_edge_g_num = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t*));
  PDM_g_num_t **data_send_group      = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t*));
  PDM_g_num_t **data_send_sens       = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t*));
  PDM_g_num_t **data_send_face_g_num = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t*));

  int         **stride_one      = (int         **) malloc( n_zone * sizeof(int        *));
  int         **zone_id         = (int         **) malloc( n_zone * sizeof(int        *));

  double      **weight        = (double      **) malloc( n_zone * sizeof(double     *));
  int          *dn_edge       = (int          *) malloc( n_zone * sizeof(int         ));

  // PDM_g_num_t key_mod = ; // Comment definir un bon key_mod ?
  for(int i_zone = 0; i_zone < n_zone; ++i_zone) {
    /*
     * Generate edge numbering
     */
    int dn_face = extract_face_distribution[i_zone][i_rank+1] - extract_face_distribution[i_zone][i_rank];
    int n_edge_elt_tot = dextract_face_vtx_idx[i_zone][dn_face];

    PDM_g_num_t* tmp_dface_edge         = (PDM_g_num_t *) malloc(     n_edge_elt_tot    * sizeof(PDM_g_num_t) );
    int*         tmp_parent_elmt_pos    = (int         *) malloc(     n_edge_elt_tot    * sizeof(int        ) );
    int*         tmp_dface_edge_vtx_idx = (int         *) malloc( ( n_edge_elt_tot + 1) * sizeof(int        ) );
    PDM_g_num_t* tmp_dface_edge_vtx     = (PDM_g_num_t *) malloc( 2 * n_edge_elt_tot    * sizeof(PDM_g_num_t) );

    int n_elmt_current = 0;
    int n_edge_current = 0;
    tmp_dface_edge_vtx_idx[0] = 0;
    PDM_poly2d_decomposes_edges(dn_face,
                                &n_elmt_current,
                                &n_edge_current,
                                extract_face_distribution[i_zone][i_rank],
                                -1,
                                dextract_face_vtx[i_zone],
                                dextract_face_vtx_idx[i_zone],
                                tmp_dface_edge_vtx_idx,
                                tmp_dface_edge_vtx,
                                tmp_dface_edge,
                                NULL,
                                NULL,
                                tmp_parent_elmt_pos);
    assert(n_edge_current == n_edge_elt_tot);

    /*
     *  Compute edges connectivity
     */
    PDM_generate_entitiy_connectivity_raw(comm,
                                          extract_vtx_distribution[i_zone][n_rank],
                                          n_edge_elt_tot,
                                          tmp_dface_edge,
                                          tmp_dface_edge_vtx_idx,
                                          tmp_dface_edge_vtx,
                                          &dn_edge[i_zone],
                                          &dedge_distrib [i_zone],
                                          &dedge_vtx_idx [i_zone],
                                          &dedge_vtx     [i_zone],
                                          &dedge_face_idx[i_zone],
                                          &dedge_face    [i_zone]);
    free(tmp_parent_elmt_pos    );

    if(1 == 1) {
      PDM_log_trace_array_long(dedge_vtx   [i_zone], dedge_vtx_idx [i_zone][dn_edge[i_zone]], "dedge_vtx :: ");
      PDM_log_trace_array_long(dedge_face  [i_zone], dedge_face_idx[i_zone][dn_edge[i_zone]], "dedge_face :: ");
      // PDM_log_trace_array_long(key_ln_to_gn[i_zone], dn_face, "key_ln_to_gn :: ");
    }

    /*
     * Echange des données pour construire la clé
     *   - exchange throw dedge_face
     */
    PDM_g_num_t *dedge_face_abs = (PDM_g_num_t *) malloc(dedge_face_idx[i_zone][dn_edge[i_zone]] * sizeof(PDM_g_num_t));
    int         *dedge_face_sgn = (int         *) malloc(dedge_face_idx[i_zone][dn_edge[i_zone]] * sizeof(int        ));
    for(int i = 0; i < dedge_face_idx[i_zone][dn_edge[i_zone]]; ++i) {
      dedge_face_abs[i] = PDM_ABS (dedge_face[i_zone][i]);
      dedge_face_sgn[i] = PDM_SIGN(dedge_face[i_zone][i]);
    }

    PDM_block_to_part_t *btp = PDM_block_to_part_create(extract_face_distribution[i_zone],
                                 (const PDM_g_num_t **) &dedge_face_abs,
                                                        &dedge_face_idx[i_zone][dn_edge[i_zone]],
                                                        1,
                                                        comm);

    /*
     * Exchange
     */
    int* dedge_face_group_id   = malloc(dedge_face_idx[i_zone][dn_edge[i_zone]] * sizeof(int        ));
    int* dedge_face_group_sens = malloc(dedge_face_idx[i_zone][dn_edge[i_zone]] * sizeof(int        ));
    int* dedge_face_join       = malloc(dedge_face_idx[i_zone][dn_edge[i_zone]] * sizeof(PDM_g_num_t));
    int* dedge_face_join_opp   = malloc(dedge_face_idx[i_zone][dn_edge[i_zone]] * sizeof(PDM_g_num_t));

    int cst_stride = 1;
    PDM_block_to_part_exch(btp,
                           sizeof(int),
                           PDM_STRIDE_CST,
                           &cst_stride,
                  (void *) dextract_face_group_id[i_zone],
                           NULL,
               (void ** ) &dedge_face_group_id);

    PDM_block_to_part_exch(btp,
                           sizeof(int),
                           PDM_STRIDE_CST,
                           &cst_stride,
                  (void *) dextract_face_group_sens[i_zone],
                           NULL,
               (void ** ) &dedge_face_group_sens);

    PDM_block_to_part_exch(btp,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST,
                           &cst_stride,
                  (void *) dextract_face_join[i_zone],
                           NULL,
               (void ** ) &dedge_face_join);

    PDM_block_to_part_exch(btp,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST,
                           &cst_stride,
                  (void *) dextract_face_join_opp[i_zone],
                           NULL,
               (void ** ) &dedge_face_join_opp);

    PDM_block_to_part_free(btp);

    free(dedge_face_abs);
    free(dedge_face_sgn);

    /*
     * Remplissage des clés + préparation buffer d'envoi
     *    - Une clé par edge pour unifier !!
     *    - Maybe besoin de signé le dgroup_face pour identifié les vtx aprés !!!
     *  Attention les l'ordre des faces dans le dfaces_group n'est pas le même que dans l'extraction !!!!!
     *  Il faut appliqué le pextract_old_to_new
     *
     *   -> Le mieux est peut-être de prétraité pour ne plus avoir le pb aprés ...
     */
    key_ln_to_gn[i_zone] = malloc(dn_edge[i_zone] * sizeof(PDM_g_num_t)); // Toutes les edges ont une clé car tout vient de l'extraction
    weight      [i_zone] = malloc(dn_edge[i_zone] * sizeof(double     )); // Toutes les edges ont une clé car tout vient de l'extraction

    stride_one          [i_zone] = malloc( (    dn_edge[i_zone]                        ) * sizeof(int        ));
    zone_id             [i_zone] = malloc( (    dn_edge[i_zone]                        ) * sizeof(int        ));
    data_send_connect_n [i_zone] = malloc( (    dn_edge[i_zone]                        ) * sizeof(int        ));
    data_send_face_n    [i_zone] = malloc( (    dn_edge[i_zone]                        ) * sizeof(int        ));
    data_send_connect   [i_zone] = malloc( (2 * dedge_face_idx[i_zone][dn_edge[i_zone]]) * sizeof(PDM_g_num_t));
    data_send_edge_g_num[i_zone] = malloc( (    dn_edge[i_zone]                        ) * sizeof(PDM_g_num_t));
    data_send_group     [i_zone] = malloc( (2 * dedge_face_idx[i_zone][dn_edge[i_zone]]) * sizeof(PDM_g_num_t));
    data_send_sens      [i_zone] = malloc( (    dedge_face_idx[i_zone][dn_edge[i_zone]]) * sizeof(PDM_g_num_t));
    data_send_face_g_num[i_zone] = malloc( (    dedge_face_idx[i_zone][dn_edge[i_zone]]) * sizeof(PDM_g_num_t));

    /*
     *  Let's build key
     */
    // int idx_write     = 0;
    int idx_write_connect    = 0;
    int idx_write_group      = 0;
    int idx_write_sens       = 0;
    int idx_write_face_g_num = 0;
    for(int i_edge = 0; i_edge < dn_edge[i_zone]; ++i_edge) {

      stride_one          [i_zone][i_edge] = 1;
      zone_id             [i_zone][i_edge] = i_zone;
      data_send_edge_g_num[i_zone][i_edge] = dedge_distrib [i_zone][i_rank] + i_edge + 1;
      data_send_connect_n [i_zone][i_edge] = 2 * (dedge_face_idx[i_zone][i_edge+1] - dedge_face_idx[i_zone][i_edge]);
      data_send_face_n    [i_zone][i_edge] = (dedge_face_idx[i_zone][i_edge+1] - dedge_face_idx[i_zone][i_edge]); // For sens and face_g_num

      PDM_g_num_t key = 0;
      for(int j = dedge_face_idx[i_zone][i_edge]; j < dedge_face_idx[i_zone][i_edge+1]; ++j) {
        key += (dedge_face_join[j] + dedge_face_join_opp[j]);
        data_send_connect[i_zone][idx_write_connect++] = dedge_face_join    [j];
        data_send_connect[i_zone][idx_write_connect++] = dedge_face_join_opp[j];
      }

      /* Send group id */
      for(int j = dedge_face_idx[i_zone][i_edge]; j < dedge_face_idx[i_zone][i_edge+1]; ++j) {
        data_send_group[i_zone][idx_write_group++] = dedge_face_group_id[j];
        int group_join_opp = group_join_to_join_opp[dedge_face_group_id[j]];
        data_send_group[i_zone][idx_write_group++] = group_join_opp;
      }

      /* Send face_g_num */
      for(int j = dedge_face_idx[i_zone][i_edge]; j < dedge_face_idx[i_zone][i_edge+1]; ++j) {
        data_send_face_g_num[i_zone][idx_write_face_g_num++] = dedge_face[i_zone][j];
      }

      /* Send group sens - Caution / 2 */
      for(int j = dedge_face_idx[i_zone][i_edge]; j < dedge_face_idx[i_zone][i_edge+1]; ++j) {
        data_send_sens[i_zone][idx_write_sens++] = dedge_face_group_sens[j];
      }

      // key_ln_to_gn[i_zone][i_edge] = key % key_mod + 1; // TODO
      key_ln_to_gn[i_zone][i_edge] = key;

      weight[i_zone][i_edge] = 2 * (dedge_face_idx[i_zone][i_edge+1] - dedge_face_idx[i_zone][i_edge]);

    }

    free(dedge_face_group_id);
    free(dedge_face_join);
    free(dedge_face_join_opp);
    free(dedge_face_group_sens);

    if(1 == 1) {
      PDM_log_trace_array_long(key_ln_to_gn[i_zone], dn_edge[i_zone], "key_ln_to_gn :: ");
    }
  }


  /*
   * Distributed hash table
   */
  PDM_part_to_block_t* ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      key_ln_to_gn,
                                                      weight,
                                                      dn_edge,
                                                      n_zone,
                                                      comm);


  int *blk_n_entity_per_key = NULL;
  int *blk_entity_n         = NULL;
  int blk_entity_n_size = PDM_part_to_block_exch(ptb,
                                                 sizeof(int),
                                                 PDM_STRIDE_VAR,
                                                 -1,
                                                 stride_one,
                                       (void **) data_send_connect_n,
                                                 &blk_n_entity_per_key,
                                       (void **) &blk_entity_n);

  int *tmp_blk_data_face_n = NULL;
  int *blk_data_face_n     = NULL;
  int blk_data_face_n_size = PDM_part_to_block_exch(ptb,
                                                    sizeof(int),
                                                    PDM_STRIDE_VAR,
                                                    -1,
                                                    stride_one,
                                          (void **) data_send_face_n,
                                                    &tmp_blk_data_face_n,
                                          (void **) &blk_data_face_n);
  free(tmp_blk_data_face_n);

  int *blk_zone_id_n = NULL;
  int *blk_zone_id         = NULL;
  PDM_part_to_block_exch(ptb,
                         sizeof(int),
                         PDM_STRIDE_VAR,
                         -1,
                         stride_one,
               (void **) zone_id,
                         &blk_zone_id_n,
               (void **) &blk_zone_id);
  free(blk_zone_id_n); // Same as blk_n_entity_per_key

  int*         tmp_blk_n_entity_per_key = NULL;
  PDM_g_num_t* blk_edge_g_num           = NULL;
  int blk_entity_n_edge_g_num = PDM_part_to_block_exch(ptb,
                                                       sizeof(PDM_g_num_t),
                                                       PDM_STRIDE_VAR,
                                                       -1,
                                                       stride_one,
                                             (void **) data_send_edge_g_num,
                                                       &tmp_blk_n_entity_per_key,
                                             (void **) &blk_edge_g_num);
  free(tmp_blk_n_entity_per_key); // Same as blk_n_entity_per_key

  int*         blk_data_connect_n = NULL;
  PDM_g_num_t* blk_data_connect   = NULL;
  int blk_data_connect_size = PDM_part_to_block_exch(ptb,
                                                     sizeof(PDM_g_num_t),
                                                     PDM_STRIDE_VAR,
                                                     -1,
                                                     data_send_connect_n,
                                           (void **) data_send_connect,
                                                     &blk_data_connect_n,
                                           (void **) &blk_data_connect);

  int*         blk_data_group_n = NULL;
  PDM_g_num_t* blk_data_group   = NULL;
  int blk_data_group_size = PDM_part_to_block_exch(ptb,
                                                   sizeof(int),
                                                   PDM_STRIDE_VAR,
                                                   -1,
                                                   data_send_connect_n,
                                         (void **) data_send_group,
                                                   &blk_data_group_n,
                                         (void **) &blk_data_group);

  int*         blk_data_sens_n = NULL;
  PDM_g_num_t* blk_data_sens   = NULL;
  int blk_data_sens_size = PDM_part_to_block_exch(ptb,
                                                  sizeof(int),
                                                  PDM_STRIDE_VAR,
                                                  -1,
                                                  data_send_face_n,
                                        (void **) data_send_sens,
                                                  &blk_data_sens_n,
                                        (void **) &blk_data_sens);

  int*         blk_data_face_g_num_n = NULL;
  PDM_g_num_t* blk_data_face_g_num   = NULL;
  PDM_part_to_block_exch(ptb,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR,
                         -1,
                         data_send_face_n,
               (void **) data_send_face_g_num,
                         &blk_data_face_g_num_n,
               (void **) &blk_data_face_g_num);
  free(blk_data_face_g_num_n); // Same as blk_data_sens_n


  /*
   * Free
   */
  for(int i_zone = 0; i_zone < n_zone; ++i_zone) {
    free(key_ln_to_gn        [i_zone]);
    free(data_send_connect_n [i_zone]);
    free(data_send_connect   [i_zone]);
    free(data_send_group     [i_zone]);
    free(data_send_edge_g_num[i_zone]);
    free(data_send_sens      [i_zone]);
    free(stride_one          [i_zone]);
    free(data_send_face_g_num[i_zone]);
    free(zone_id             [i_zone]);
    free(weight              [i_zone]);
  }
  free(key_ln_to_gn        );
  free(data_send_connect_n );
  free(data_send_connect   );
  free(data_send_group     );
  free(data_send_edge_g_num);
  free(data_send_sens      );
  free(data_send_face_g_num);
  free(stride_one          );
  free(zone_id             );
  free(weight              );


  int blk_size = PDM_part_to_block_n_elt_block_get(ptb);
  if(1 == 1) {
    PDM_log_trace_array_int(blk_n_entity_per_key, blk_size             , "blk_n_entity_per_key:: ");
    PDM_log_trace_array_int(blk_entity_n        , blk_entity_n_size    , "blk_entity_n        :: ");

    PDM_log_trace_array_int(blk_zone_id         , blk_entity_n_size    , "blk_zone_id         :: ");

    PDM_log_trace_array_long(blk_edge_g_num     , blk_entity_n_size    , "blk_edge_g_num      :: ");

    PDM_log_trace_array_int (blk_data_connect_n  , blk_size             , "blk_data_connect_n  :: ");
    PDM_log_trace_array_long(blk_data_connect    , blk_data_connect_size, "blk_data_connect    :: ");

    PDM_log_trace_array_long(blk_data_group_n    , blk_size             , "blk_data_group_n    :: ");
    PDM_log_trace_array_long(blk_data_group      , blk_data_group_size  , "blk_data_group      :: ");

    PDM_log_trace_array_long(blk_data_sens_n     , blk_size             , "blk_data_sens_n     :: ");

    PDM_log_trace_array_long(blk_data_face_n     , blk_data_face_n_size , "blk_data_face_n     :: ");

    PDM_log_trace_array_long(blk_data_sens       , blk_data_sens_size   , "blk_data_sens       :: ");
    PDM_log_trace_array_long(blk_data_face_g_num , blk_data_sens_size   , "blk_data_face_g_num :: ");
  }

  /*
   * Post-treatment
   */
  int n_max_entity_per_key = 0;
  int n_max_connec         = 0;
  for(int i = 0; i < blk_size; ++i) {
    n_max_entity_per_key = PDM_MAX(blk_n_entity_per_key[i], n_max_entity_per_key);
    n_max_connec         = PDM_MAX(blk_data_connect_n  [i], n_max_connec);
  }

  int* zone_id_n = (int * ) malloc(n_zone * sizeof(int));
  for(int i = 0; i < n_zone; ++i){
    zone_id_n[i] = 0;
  }

  int* blk_entity_idx    = malloc((blk_entity_n_size    + 1) * sizeof(int));
  int* blk_data_face_idx = malloc((blk_data_face_n_size + 1) * sizeof(int));
  blk_entity_idx[0] = 0;
  for(int i = 0; i < blk_entity_n_size; ++i) {
    blk_entity_idx[i+1] = blk_entity_idx[i] + blk_entity_n[i];
    zone_id_n[blk_zone_id[i]]++;
  }

  blk_data_face_idx[0] = 0;
  for(int i = 0; i < blk_data_face_n_size; ++i) {
    blk_data_face_idx[i+1] = blk_data_face_idx[i] + blk_data_face_n[i];
  }

  PDM_log_trace_array_int(zone_id_n, n_zone, "zone_id_n :: ");
  PDM_log_trace_array_int(blk_entity_idx, blk_entity_n_size+1, "blk_entity_idx :: ");
  PDM_log_trace_array_int(blk_data_face_idx, blk_data_face_n_size+1, "blk_data_face_idx :: ");

  printf("n_max_entity_per_key = %i \n", n_max_entity_per_key);
  printf("n_max_connec         = %i \n", n_max_connec);

  // PDM_g_num_t** loc_connec       = (PDM_g_num_t *) malloc(  n_max_entity_per_key * sizeof(PDM_g_num_t *) );
  int*          already_treat    = (int         *) malloc(  n_max_entity_per_key * sizeof(int          ) );
  int*          same_entity_idx  = (int         *) malloc(  n_max_entity_per_key * sizeof(int          ) );
  int*          sens_entity      = (int         *) malloc(  n_max_entity_per_key * sizeof(int          ) );

  PDM_g_num_t **ln_to_gn_res     = malloc(n_zone * sizeof(PDM_g_num_t *));
  PDM_g_num_t **results_edge     = malloc(n_zone * sizeof(PDM_g_num_t *));
  PDM_g_num_t **results_edge_opp = malloc(n_zone * sizeof(PDM_g_num_t *));
  int         **results_zone_opp = malloc(n_zone * sizeof(int         *));
  for(int i = 0; i < n_zone; ++i) {
    ln_to_gn_res    [i] = malloc(40 * zone_id_n[i] * sizeof(PDM_g_num_t));
    results_edge    [i] = malloc(40 * zone_id_n[i] * sizeof(PDM_g_num_t));
    results_edge_opp[i] = malloc(40 * zone_id_n[i] * sizeof(PDM_g_num_t));
    results_zone_opp[i] = malloc(40 * zone_id_n[i] * sizeof(int        ));
  }

  /* Reset to fill */
  for(int i = 0; i < n_zone; ++i){
    zone_id_n[i] = 0;
  }

  int idx  = 0;
  int idx2 = 0;
  for(int i_key = 0; i_key < blk_size; ++i_key) {

    int n_conflict_entitys = blk_n_entity_per_key[i_key];

    log_trace(" i_key = %i | n_conflict_entitys = %i \n", i_key, n_conflict_entitys);

    /* Reset */
    for(int j = 0; j < n_conflict_entitys; ++j) {
      already_treat[j] = -1;
    }

    /* Loop over all entitys in conflict and sort all */
    for(int i_entity = 0; i_entity < n_conflict_entitys; ++i_entity) {

      int beg    = blk_entity_idx[idx+i_entity  ];
      int end    = blk_entity_idx[idx+i_entity+1];
      int l_size = end - beg;

      // Caution inplace sort !!!!
      PDM_sort_long(&blk_data_connect[beg], NULL, l_size);
      PDM_sort_int (&blk_data_group  [beg], NULL, l_size);

      if(0 == 1) {
        PDM_log_trace_array_long(&blk_data_connect[beg], l_size, "blk_data_connect (sort) :: ");
        PDM_log_trace_array_int (&blk_data_group  [beg], l_size, "blk_data_group   (sort) :: ");
      }
    }

    /*
     *  Identify pair or invilid other
     */
    for(int i_entity = 0; i_entity < n_conflict_entitys; ++i_entity) {
      int beg1    = blk_entity_idx[idx+i_entity  ];
      int end1    = blk_entity_idx[idx+i_entity+1];
      int l_size1 = end1 - beg1;

      int n_same         = 0;
      int i_entity2_same = -1;

      if(already_treat[i_entity] == 1) {
        continue;
      }

      for(int i_entity2 = i_entity+1; i_entity2 < n_conflict_entitys; ++i_entity2) {

        if(already_treat[i_entity2] == -1) {
          int beg2    = blk_entity_idx[idx+i_entity2  ];
          int end2    = blk_entity_idx[idx+i_entity2+1];
          int l_size2 = end2 - beg2;

          if(l_size1 != l_size2) {
            continue;
          }

          // Compare
          int have_same_group_pair = 1;
          for(int k = 0; k < l_size1; ++k) {
            if(blk_data_group[beg1+k] != blk_data_group[beg2+k]){
              have_same_group_pair = 0;
              break;
            }
          }

          if(!have_same_group_pair) {
            continue;
          }

          int have_same_entity_list = 1;
          for(int k = 0; k < l_size1; ++k) {
            if(blk_data_connect[beg1+k] != blk_data_connect[beg2+k]){
              have_same_entity_list = 0;
              break;
            }
          }

          if(!have_same_entity_list) {
            continue;
          }

          already_treat[i_entity2] = 1;
          same_entity_idx[n_same++] = i_entity2;
          i_entity2_same = i_entity2;

          if(n_same > 1) {
            i_entity2_same = -1;
          }

        }
      } /* End for i_entity2 */

      log_trace("i_entity = %i | i_entity2_same = %i | n_same = %i \n", i_entity, i_entity2_same, n_same);

      assert(n_same >= 1);
      log_trace("idx2 = %i \n", idx2);
      /*
       * Renvoie des resultats :
       *    - Par edge --> Il faut donc également le numero de zones
       */
      int beg_face_1 = blk_data_face_idx[idx2+i_entity  ];
      int end_face_1 = blk_data_face_idx[idx2+i_entity+1];
      int n_face_1   = end_face_1 - beg_face_1;

      for(int k = 0; k < n_face_1; ++k) {
        PDM_g_num_t face_cur   = blk_data_face_g_num[beg_face_1+k];
        PDM_g_num_t edge_cur   = blk_edge_g_num     [idx+i_entity];
        int         i_zone_cur = blk_zone_id        [idx+i_entity];

        ln_to_gn_res    [i_zone_cur][zone_id_n[i_zone_cur]  ] = face_cur;
        results_edge    [i_zone_cur][zone_id_n[i_zone_cur]  ] = edge_cur;
        results_edge_opp[i_zone_cur][zone_id_n[i_zone_cur]  ] = edge_cur;
        results_zone_opp[i_zone_cur][zone_id_n[i_zone_cur]++] = i_zone_cur; // In order to manage corner we need this

        for(int i_same = 0; i_same < n_same; ++i_same) {

          int i_entity2  = same_entity_idx[i_same];
          int beg_face_2 = blk_data_face_idx[idx2+i_entity2  ];
          int end_face_2 = blk_data_face_idx[idx2+i_entity2+1];
          int n_face_2   = end_face_2 - beg_face_2;

          for(int l = 0; l < n_face_2; ++l) {
            PDM_g_num_t face_opp   = blk_data_face_g_num[beg_face_2+l];
            PDM_g_num_t edge_opp   = blk_edge_g_num     [idx+i_entity2];
            int         i_zone_opp = blk_zone_id        [idx+i_entity2];

            // Symetric fill
            log_trace("face_cur = %i \n", face_cur);
            ln_to_gn_res    [i_zone_cur][zone_id_n[i_zone_cur]  ] = face_cur;
            results_edge    [i_zone_cur][zone_id_n[i_zone_cur]  ] = edge_cur;
            results_edge_opp[i_zone_cur][zone_id_n[i_zone_cur]  ] = edge_opp;
            results_zone_opp[i_zone_cur][zone_id_n[i_zone_cur]++] = i_zone_opp; // In order to manage corner we need this

            // Symetric fill
            ln_to_gn_res    [i_zone_opp][zone_id_n[i_zone_opp]  ] = face_opp;
            results_edge    [i_zone_opp][zone_id_n[i_zone_opp]  ] = edge_opp;
            results_edge_opp[i_zone_opp][zone_id_n[i_zone_opp]  ] = edge_cur;
            results_zone_opp[i_zone_opp][zone_id_n[i_zone_opp]++] = i_zone_cur; // In order to manage corner we need this
          }
        }
      }

      // Renvoi de tout les edges candidats à travers les faces ???
      already_treat[i_entity] = 1;
    }

    idx  += n_conflict_entitys;
    idx2 += n_conflict_entitys;
    // idx2 += (blk_data_face_idx[idx2+n_conflict_entitys] - blk_data_face_idx[idx2]);
    // idx2 = blk_data_face_idx[idx2+n_conflict_entitys];

  }

  /*
   *  Part to block avec face_extraction_dsitrib
   */
  for(int i = 0; i < n_zone; ++i) {
    PDM_log_trace_array_long(ln_to_gn_res    [i], zone_id_n[i], "ln_to_gn_res     :: ");
    PDM_log_trace_array_long(results_edge    [i], zone_id_n[i], "results_edge     :: ");
    PDM_log_trace_array_long(results_edge_opp[i], zone_id_n[i], "results_edge_opp :: ");
    PDM_log_trace_array_int (results_zone_opp[i], zone_id_n[i], "results_zone_opp :: ");
  }


  /*
   * A verifier que les edges qui sont indeterminé pointe sur une seule face !!!
   *   Quand on decompose il faut également envoyer le numero de la face (de l'extraction)
   */
  free(zone_id_n);

  for(int i = 0; i < n_zone; ++i) {
    free(ln_to_gn_res    [i]);
    free(results_edge    [i]);
    free(results_edge_opp[i]);
    free(results_zone_opp[i]);
  }
  free(ln_to_gn_res    );
  free(results_edge    );
  free(results_edge_opp);
  free(results_zone_opp);

  // free(loc_entity_1    );
  // free(loc_entity_2    );
  free(already_treat   );
  free(same_entity_idx );
  free(sens_entity     );

  free(blk_n_entity_per_key);
  free(blk_entity_n        );

  free(blk_data_connect_n);
  free(blk_data_connect  );
  free(blk_entity_idx);
  free(blk_data_face_idx);
  free(blk_data_face_n);

  free(blk_data_group_n);
  free(blk_data_group  );
  free(blk_data_sens_n);
  free(blk_data_sens);
  free(blk_edge_g_num);
  free(blk_zone_id);
  free(blk_data_face_g_num);


  PDM_part_to_block_free(ptb);

  /*
   * Renvoie des resultats
   *   btp avec les key_ln_to_gn
   *   Attention aux block partiels !!!
   */

  /*
   * Pour les edges recalitrant
   *   --> On extract les faces dont les edges_opp ne sont pas referencé
   *   --> Normalement il y au moins un edge par face qui est connecté
   *   --> On échange le dface_edge + l'information du sens devrait permmettre de tout reorienté
   *       ==> On trouve le edge commun puis on tourne dans le sens de l'orientation respective
   */


  /*
   * Pour les coins coins il faut probablement un echange en plus entre les graph 2 à 2 puis unifier
   * Maybe utiliser le PDM_part_dentity_group_to_pentity_group
   * A la fin de l'algo un vtx est potentiellement lié à plusieurs groupes.
   * On peut merger avec un part_to_block tout les tags d'un vertex en appliquant sur les PL et PLD
   * --> Ca gère le coin
   */
  for(int i_zone = 0; i_zone < n_zone; ++i_zone) {
    free(dedge_distrib  [i_zone]);
    free(dedge_vtx_idx  [i_zone]);
    free(dedge_vtx      [i_zone]);
    free(dedge_face_idx [i_zone]);
    free(dedge_face     [i_zone]);
  }

  free(dedge_distrib);
  free(dedge_vtx_idx);
  free(dedge_vtx);
  free(dedge_face_idx);
  free(dedge_face);
  free(dn_edge);
}

/**
 *
 * \brief  Usage
 *
 */

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t  *n_vtx_seg,
           double        *length)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  PDM_g_num_t        n_vtx_seg = 10;
  double             length    = 1.;
  int                n_zone    = 2;


  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length);

  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm     comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /* Alloc for distributed mesh */
  int *dn_cell       = (int *) malloc(n_zone * sizeof(int));
  int *dn_face       = (int *) malloc(n_zone * sizeof(int));
  int *dn_vtx        = (int *) malloc(n_zone * sizeof(int));
  int *n_face_group  = (int *) malloc(n_zone * sizeof(int));
  int *dface_vtx_s   = (int *) malloc(n_zone * sizeof(int));
  int *dface_group_s = (int *) malloc(n_zone * sizeof(int));

  PDM_g_num_t  **dface_cell      = (PDM_g_num_t **) malloc(n_zone * sizeof(PDM_g_num_t *));
  int          **dface_vtx_idx   = (int         **) malloc(n_zone * sizeof(int         *));
  PDM_g_num_t  **dface_vtx       = (PDM_g_num_t **) malloc(n_zone * sizeof(PDM_g_num_t *));
  double       **dvtx_coord      = (double      **) malloc(n_zone * sizeof(double      *));
  int          **dface_group_idx = (int         **) malloc(n_zone * sizeof(int         *));
  PDM_g_num_t  **dface_group     = (PDM_g_num_t **) malloc(n_zone * sizeof(PDM_g_num_t *));
  int          **dface_bnd_idx   = (int         **) malloc(n_zone * sizeof(int         *));
  PDM_g_num_t  **dface_bnd       = (PDM_g_num_t **) malloc(n_zone * sizeof(PDM_g_num_t *));

  PDM_g_num_t  **extract_face_distribution = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t *));
  PDM_g_num_t  **extract_vtx_distribution  = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t *));
  int          **dextract_face_vtx_idx     = (int         **) malloc( n_zone * sizeof(int         *));
  PDM_g_num_t  **dextract_face_vtx         = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t *));
  PDM_g_num_t  **dparent_face_g_num        = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t *));
  PDM_g_num_t  **dparent_vtx_g_num         = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t *));
  PDM_g_num_t  **pextract_old_to_new       = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t *));
  double       **dextract_vtx_coord        = (double      **) malloc( n_zone * sizeof(double      *));

  int n_group_join = 2*(n_zone-1);
  int *group_join_to_zone_cur = (int *) malloc( n_group_join * sizeof(int));
  int *group_join_to_zone_opp = (int *) malloc( n_group_join * sizeof(int));
  int *group_join_to_join_opp = (int *) malloc( n_group_join * sizeof(int));

  int          *dface_join_idx  = (int         *) malloc((n_group_join + 1) * sizeof(int        ));
  PDM_g_num_t  *dface_join      = NULL; // A allouer propremet

  int tmp_i_zone = 0;
  for (int i_join = 0; i_join < n_group_join; i_join++) {
    if (i_join % 2 == 0) {
      group_join_to_join_opp[i_join] = i_join + 1;
      group_join_to_zone_opp[i_join] = tmp_i_zone + 1;
    } else {
      group_join_to_join_opp[i_join] = i_join - 1;
      group_join_to_zone_opp[i_join] = tmp_i_zone++;
    }
  }

  PDM_log_trace_array_int(group_join_to_join_opp, n_group_join, "group_join_to_join_opp :: ");
  PDM_log_trace_array_int(group_join_to_zone_opp, n_group_join, "group_join_to_zone_opp :: ");

  for(int i_group_join = 0; i_group_join < n_group_join+1; ++i_group_join) {
    dface_join_idx[i_group_join] = 0;
  }

  PDM_dcube_t **dcube = (PDM_dcube_t **) malloc(n_zone * sizeof(PDM_dcube_t *));
  int tmp_i_group_join = 0;
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {

    dcube[i_zone] = PDM_dcube_gen_init(comm, n_vtx_seg, length, i_zone, 0., 0., PDM_OWNERSHIP_KEEP);
    PDM_dcube_gen_dim_get(dcube         [i_zone],
                          &n_face_group [i_zone],
                          &dn_cell      [i_zone],
                          &dn_face      [i_zone],
                          &dn_vtx       [i_zone],
                          &dface_vtx_s  [i_zone],
                          &dface_group_s[i_zone]);

    PDM_dcube_gen_data_get(dcube          [i_zone],
                          &dface_cell     [i_zone],
                          &dface_vtx_idx  [i_zone],
                          &dface_vtx      [i_zone],
                          &dvtx_coord     [i_zone],
                          &dface_group_idx[i_zone],
                          &dface_group    [i_zone]);

    /*
     * Les faces groups du dcube sont : zmin, zmax, xmin, xmax, ymin, ymax
     * Il faut les séparer en faces de bords et faces raccord, sachant que
     * les zones sont alignées selon X
     */
    int n_bnd = 4;
    int n_jn  = 2;
    if (i_zone == 0){
      n_bnd++;
      n_jn-- ;
    }
    if (i_zone == n_zone-1){
      n_bnd++;
      n_jn-- ;
    }

    // Join numbering (left to right, increasing i_zone)
    printf("n_jn = %i \n", n_jn);
    printf("tmp_i_group_join = %i \n", tmp_i_group_join);

    dface_bnd_idx [i_zone] = (int *) malloc((n_bnd        + 1) * sizeof(int));

    // First pass to count and allocate
    int i_bnd = 1;
    int i_jn  = tmp_i_group_join+1;
    dface_bnd_idx[i_zone][0]  = 0;
    // dface_join_idx[0] = 0;
    for (int igroup = 0; igroup < n_face_group[i_zone]; igroup++) {
      int copy_to_bnd = (igroup != 2 || (igroup == 2 && i_zone == 0)) && (igroup != 3 || (igroup == 3 && i_zone == n_zone-1));
      int group_size = dface_group_idx[i_zone][igroup+1] - dface_group_idx[i_zone][igroup];
      if (copy_to_bnd) { //Its a boundary
        dface_bnd_idx[i_zone][i_bnd++] = group_size;
      } else { //Its a join
        group_join_to_zone_cur[i_jn-1] = i_zone;
        dface_join_idx[i_jn++] = group_size;
      }
    }
    for (int i = 0; i < n_bnd; i++) {
      dface_bnd_idx[i_zone][i+1] = dface_bnd_idx[i_zone][i+1] + dface_bnd_idx[i_zone][i];
    }

    printf("i_jn = %i \n", i_jn);

    for (int i = tmp_i_group_join; i < i_jn-1; i++) {
      dface_join_idx[i+1] = dface_join_idx[i+1] + dface_join_idx[i];
    }

    PDM_log_trace_array_int(dface_join_idx, n_group_join+1, "dface_join_idx :: ");

    /* A bit stupid but complicated to made it in ohter way for a clear test */
    dface_join = realloc(dface_join, dface_join_idx[i_jn-1] * sizeof(PDM_g_num_t));

    // Second pass to copy
    dface_bnd [i_zone] = (PDM_g_num_t *) malloc(dface_bnd_idx [i_zone][n_bnd        ] * sizeof(PDM_g_num_t));
    i_bnd = 0;
    i_jn  = tmp_i_group_join;
    for (int igroup = 0; igroup < n_face_group[i_zone]; igroup++) {
      int copy_to_bnd = (igroup != 2 || (igroup == 2 && i_zone == 0)) && (igroup != 3 || (igroup == 3 && i_zone == n_zone-1));
      if (copy_to_bnd){ //Its a boundary
        for (int i = dface_group_idx[i_zone][igroup]; i < dface_group_idx[i_zone][igroup+1]; i++) {
          dface_bnd[i_zone][i_bnd++] = dface_group[i_zone][i];
        }
      } else { //Its a join
        int k = 0;
        for (int i = dface_group_idx[i_zone][igroup]; i < dface_group_idx[i_zone][igroup+1]; i++) {
          dface_join[dface_join_idx[i_jn]+k++] = dface_group[i_zone][i];
        }
        i_jn++;
      }
    }

    /*
     *  Go to nexts join
     */
    tmp_i_group_join += n_jn;
  }

  PDM_log_trace_array_int(dface_join_idx, n_group_join+1, "dface_join_idx (END) :: ");
  PDM_log_trace_array_int(group_join_to_zone_cur, n_group_join, "group_join_to_zone_cur :: ");

  /*
   * Setup dface_join_opp + group_id in current layout
   */
  PDM_g_num_t *dface_join_opp = NULL;
  _exchange_point_list(n_group_join,
                       group_join_to_join_opp,
                       dface_join_idx,
                       dface_join,
                       &dface_join_opp,
                       comm);



  /*
   * Algorithm begin - Extract faces
   */
  int         **dextract_face_group_id   = malloc(n_zone * sizeof(int         *));
  int         **dextract_face_group_sens = malloc(n_zone * sizeof(int         *));
  PDM_g_num_t **dextract_face_join       = malloc(n_zone * sizeof(PDM_g_num_t *));
  PDM_g_num_t **dextract_face_join_opp   = malloc(n_zone * sizeof(PDM_g_num_t *));

  tmp_i_group_join = 0;
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    int n_jn  = 2;
    if (i_zone == 0       ){n_jn--;}
    if (i_zone == n_zone-1){n_jn--;}

    /*
     *  Now we have all joins create we need to extract them
     */
    PDM_g_num_t* face_distribution = PDM_compute_entity_distribution(comm, dn_face[i_zone]);
    PDM_g_num_t* vtx_distribution  = PDM_compute_entity_distribution(comm, dn_vtx [i_zone]);

    /*
     *  By construction of actual join, each join of same zone are concatenate
     */
    int dn_l_face_join = dface_join_idx[tmp_i_group_join+n_jn] - dface_join_idx[tmp_i_group_join];
    PDM_g_num_t* l_dface_join = &dface_join[dface_join_idx[tmp_i_group_join]];

    PDM_dconnectivity_to_extract_dconnectivity(comm,
                                               dn_l_face_join,
                                               l_dface_join,
                                               face_distribution,
                                               dface_vtx_idx[i_zone],
                                               dface_vtx[i_zone],
                                               &extract_face_distribution[i_zone],
                                               &extract_vtx_distribution[i_zone],
                                               &dextract_face_vtx_idx[i_zone],
                                               &dextract_face_vtx[i_zone],
                                               &dparent_face_g_num[i_zone],
                                               &dparent_vtx_g_num[i_zone],
                                               &pextract_old_to_new[i_zone]);

    /*
     * Equilibrate the pextract_old_to_new
     */
    PDM_g_num_t* dface_group_init_distrib = PDM_compute_entity_distribution(comm, dn_l_face_join);
    int dn_extract_face = extract_face_distribution[i_zone][i_rank+1] - extract_face_distribution[i_zone][i_rank];
    PDM_g_num_t* extract_face_ln_to_gn = malloc(dn_extract_face * sizeof(PDM_g_num_t));
    for(int i = 0; i < dn_extract_face; ++i) {
      extract_face_ln_to_gn[i] = extract_face_distribution[i_zone][i_rank] + i + 1;
    }

    PDM_g_num_t** tmp_dextract_old_to_new = NULL;
    PDM_part_dfield_to_pfield(comm,
                              1,
                              sizeof(PDM_g_num_t),
                              dface_group_init_distrib,
      (unsigned char    *)    pextract_old_to_new[i_zone],
                              &dn_extract_face,
      (const PDM_g_num_t **)  &extract_face_ln_to_gn,
      (unsigned char ***)     &tmp_dextract_old_to_new);
    PDM_g_num_t *dextract_old_to_new = tmp_dextract_old_to_new[0];
    free(tmp_dextract_old_to_new);

    int dn_extract_vtx  = extract_vtx_distribution[i_zone][i_rank+1] - extract_vtx_distribution[i_zone][i_rank];

    double** tmp_dextract_vtx_coord = NULL;
    PDM_part_dcoordinates_to_pcoordinates(comm,
                                          1,
                                          vtx_distribution,
                                          dvtx_coord[i_zone],
                                          &dn_extract_vtx,
                   (const PDM_g_num_t **) &dparent_vtx_g_num[i_zone],
                                          &tmp_dextract_vtx_coord);

    /*
     *  Pour chaque face extraite :
     *     - Echange du numero de face/face_opposé (Attention on est obligé de le faire car PDM_dconnectivity_to_extract_dconnectivity reordonne tout)
     *     - Exchange i_group also
     */
    int* dface_group_tag  = malloc(dn_l_face_join * sizeof(int));
    int* dface_group_sens = malloc(dn_l_face_join * sizeof(int)); // 1 = Exterior / -1 = Interior
    int k = 0;
    for(int i_group = tmp_i_group_join; i_group < tmp_i_group_join+n_jn; ++i_group) {
      for(int i_face = dface_join_idx[i_group]; i_face < dface_join_idx[i_group+1]; ++i_face) {
        dface_group_tag [k  ] = i_group;
        dface_group_sens[k++] = 1;
      }
    }

    if(0 == 1) {
      log_trace("dn_extract_face = %i \n", dn_extract_face);
      log_trace("dn_l_face_join  = %i \n", dn_l_face_join);
    }

    /*
     * Exchange tag
     */
    int** tmp_dextract_face_group_id = NULL;
    PDM_part_dfield_to_pfield(comm,
                              1,
                              sizeof(int),
                              dface_group_init_distrib,
      (unsigned char    *)    dface_group_tag,
                              &dn_extract_face,
      (const PDM_g_num_t **)  &dextract_old_to_new,
      (unsigned char ***)     &tmp_dextract_face_group_id);
    dextract_face_group_id[i_zone] = tmp_dextract_face_group_id[0];
    free(tmp_dextract_face_group_id);

    /*
     * Exchange sens
     */
    int** tmp_dextract_face_group_sens = NULL;
    PDM_part_dfield_to_pfield(comm,
                              1,
                              sizeof(int),
                              dface_group_init_distrib,
      (unsigned char    *)    dface_group_sens,
                              &dn_extract_face,
      (const PDM_g_num_t **)  &dextract_old_to_new,
      (unsigned char ***)     &tmp_dextract_face_group_sens);
    dextract_face_group_sens[i_zone] = tmp_dextract_face_group_sens[0];
    free(tmp_dextract_face_group_sens);

    /*
     *  Exchange dface_join
     */
    PDM_g_num_t* sub_dface_join = &dface_join[dface_join_idx[tmp_i_group_join]];
    int** tmp_dextract_face_join = NULL;
    PDM_part_dfield_to_pfield(comm,
                              1,
                              sizeof(PDM_g_num_t),
                              dface_group_init_distrib,
      (unsigned char    *)    sub_dface_join,
                              &dn_extract_face,
      (const PDM_g_num_t **)  &dextract_old_to_new,
      (unsigned char ***)     &tmp_dextract_face_join);
    dextract_face_join[i_zone] = tmp_dextract_face_join[0];
    free(tmp_dextract_face_join);

    /*
     *  Exchange dface_join_opp
     */
    PDM_g_num_t* sub_dface_join_opp = &dface_join_opp[dface_join_idx[tmp_i_group_join]];
    int** tmp_dextract_face_join_opp = NULL;
    PDM_part_dfield_to_pfield(comm,
                              1,
                              sizeof(PDM_g_num_t),
                              dface_group_init_distrib,
      (unsigned char    *)    sub_dface_join_opp,
                              &dn_extract_face,
      (const PDM_g_num_t **)  &dextract_old_to_new,
      (unsigned char ***)     &tmp_dextract_face_join_opp);
    dextract_face_join_opp[i_zone] = tmp_dextract_face_join_opp[0];
    free(tmp_dextract_face_join_opp);


    free(extract_face_ln_to_gn);

    dextract_vtx_coord[i_zone] = tmp_dextract_vtx_coord[0];
    free(tmp_dextract_vtx_coord);

    free(face_distribution);
    free(vtx_distribution);
    free(dface_group_init_distrib);
    free(dface_group_tag);
    free(dface_group_sens);
    free(dextract_old_to_new); // Maybe to keep ??

    /*
     * To keep
     */
    if(1 == 1) {
      PDM_log_trace_array_long(dextract_face_group_id  [i_zone], dn_extract_face, "dextract_face_group_id :: ");
      PDM_log_trace_array_long(dextract_face_group_sens[i_zone], dn_extract_face, "dextract_face_group_id :: ");
      PDM_log_trace_array_long(dextract_face_join      [i_zone], dn_extract_face, "dextract_face_join     :: ");
      PDM_log_trace_array_long(dextract_face_join_opp  [i_zone], dn_extract_face, "dextract_face_join_opp :: ");
    }

    /*
     *  Go to nexts join
     */
    tmp_i_group_join += n_jn;

  }

  /*
   * Real algorithm
   */
  _deduce_descending_join(n_zone,
                          n_group_join,
                          group_join_to_zone_cur,
                          group_join_to_zone_opp,
                          group_join_to_join_opp,
                          dface_join_idx,
                          dface_join,
                          dface_join_opp,
                          dface_vtx_idx,
                          dface_vtx,
                          extract_face_distribution,
                          extract_vtx_distribution,
                          dextract_face_vtx_idx,
                          dextract_face_vtx,
                          dextract_face_group_id,
                          dextract_face_group_sens,
                          dextract_face_join,
                          dextract_face_join_opp,
                          dparent_face_g_num,
                          dparent_vtx_g_num,
                          pextract_old_to_new,
                          dextract_vtx_coord,
                          comm);

  /* Free memory */
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    free(dface_bnd_idx [i_zone]);
    free(dface_bnd     [i_zone]);
    free(extract_face_distribution[i_zone]);
    free(extract_vtx_distribution [i_zone]);
    free(dextract_face_vtx_idx    [i_zone]);
    free(dextract_face_vtx        [i_zone]);
    free(dparent_face_g_num       [i_zone]);
    free(dparent_vtx_g_num        [i_zone]);
    free(pextract_old_to_new      [i_zone]);
    free(dextract_vtx_coord       [i_zone]);
    free(dextract_face_group_id   [i_zone]);
    free(dextract_face_group_sens [i_zone]);
    free(dextract_face_join       [i_zone]);
    free(dextract_face_join_opp   [i_zone]);
    PDM_dcube_gen_free(dcube[i_zone]);
  }
  free(dcube);
  free(dn_cell);
  free(dn_face);
  free(dn_vtx);
  free(n_face_group);
  free(dface_group_s);
  free(dface_vtx_s);
  free(dface_cell);
  free(dface_vtx_idx);
  free(dface_vtx);
  free(dvtx_coord);
  free(dface_group_idx);
  free(dface_group);
  free(dface_bnd_idx);
  free(dface_bnd);
  free(dface_join_idx);
  free(dface_join);
  free(group_join_to_zone_cur);
  free(group_join_to_zone_opp);
  free(group_join_to_join_opp);
  free(dextract_face_group_id);
  free(dextract_face_group_sens);
  free(dextract_face_join    );
  free(dextract_face_join_opp);
  free(dface_join_opp);

  free(extract_face_distribution);
  free(extract_vtx_distribution );
  free(dextract_face_vtx_idx    );
  free(dextract_face_vtx        );
  free(dparent_face_g_num       );
  free(dparent_vtx_g_num        );
  free(pextract_old_to_new      );
  free(dextract_vtx_coord       );

  PDM_MPI_Finalize();

  return 0;
}
