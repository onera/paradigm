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
#include "pdm_multi_block_to_part.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_sort.h"
#include "pdm_array.h"
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

static PDM_g_num_t* _per_block_offset(int n_block, PDM_g_num_t *sizes, PDM_MPI_Comm comm) {
  PDM_g_num_t *per_block_offset = (PDM_g_num_t *) malloc((n_block+1) * sizeof(PDM_g_num_t));
  per_block_offset[0] = 0;
  PDM_MPI_Allreduce(sizes, &per_block_offset[1], n_block, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
  PDM_array_accumulate_gnum(per_block_offset, n_block+1);
  return per_block_offset;
}


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

    if(0 == 1) {
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


static void _split_paired_connectivity(int n_item, int *in_array_idx, PDM_g_num_t *in_array,
    int **out_array_idx, PDM_g_num_t **out_array1, PDM_g_num_t **out_array2)
{

  assert (n_item % 2 == 0);
  assert (in_array_idx[n_item] % 2 == 0);

  *out_array_idx = malloc(((n_item / 2 )+1)               * sizeof(int));
  *out_array1    = malloc((in_array_idx[n_item] / 2) * sizeof(PDM_g_num_t));
  *out_array2    = malloc((in_array_idx[n_item] / 2) * sizeof(PDM_g_num_t));

  int         *_out_array_idx = *out_array_idx;
  PDM_g_num_t *_out_array1    = *out_array1;
  PDM_g_num_t *_out_array2    = *out_array2;

  _out_array_idx[0] = 0;
  int recv_data_idx = 0;
  for (int i_item = 0; i_item < n_item/2; i_item++) {
    int n_elt = in_array_idx[2*i_item+1] - in_array_idx[2*i_item];
    _out_array_idx[i_item+1] = _out_array_idx[i_item] + n_elt;
    memcpy(&_out_array1[_out_array_idx[i_item]], &in_array[recv_data_idx]      , n_elt * sizeof(PDM_g_num_t));
    memcpy(&_out_array2[_out_array_idx[i_item]], &in_array[recv_data_idx+n_elt], n_elt * sizeof(PDM_g_num_t));
    recv_data_idx += 2*n_elt;
  }
  assert (recv_data_idx == in_array_idx[n_item]);
}

static int _extract_and_shift_jn_faces
(
 int           n_zone,
 PDM_g_num_t  *dn_face,
 PDM_g_num_t  *dn_vtx,
 int           n_group_join,
 int          *group_join_to_join_opp,
 int          *group_join_to_zone_cur,
 int          *group_join_to_zone_opp,
 int          *dface_join_idx,
 PDM_g_num_t  *dface_join,
 PDM_g_num_t  *dface_join_opp,
 int         **dface_vtx_idx,
 PDM_g_num_t **dface_vtx,
 int         **face_vtx_both_idx,
 PDM_g_num_t **face_vtx_both,
 PDM_MPI_Comm  comm
)
{
  PDM_g_num_t *face_per_block_offset = _per_block_offset(n_zone, dn_face, comm);
  PDM_g_num_t *vtx_per_block_offset  = _per_block_offset(n_zone, dn_vtx,  comm);
  
  int n_face_join = 0; // Put faces only once
  for (int i_join = 0; i_join < n_group_join; i_join++) {
    if (i_join <= group_join_to_join_opp[i_join])
      n_face_join += 2*(dface_join_idx[i_join+1] - dface_join_idx[i_join]);
  }

  PDM_g_num_t *multi_gnum = malloc(n_face_join * sizeof(PDM_g_num_t));

  int idx = 0;
  for (int i_join = 0; i_join < n_group_join; i_join++) {
    int i_join_opp = group_join_to_join_opp[i_join];
    int i_zone_cur = group_join_to_zone_cur[i_join];
    int i_zone_opp = group_join_to_zone_opp[i_join];

    log_trace("Treat jn %i\n", i_join);
    if (i_join <= i_join_opp) {
      for (int i_face_jn = dface_join_idx[i_join]; i_face_jn < dface_join_idx[i_join+1]; i_face_jn++) {
        multi_gnum[idx++] = dface_join[i_face_jn] + face_per_block_offset[i_zone_cur];
        multi_gnum[idx++] = dface_join_opp[i_face_jn] + face_per_block_offset[i_zone_opp];
      }
      /*for (int i_face_jn = dface_join_idx[i_join]; i_face_jn < dface_join_idx[i_join+1]; i_face_jn++) {*/
        /*multi_gnum[idx++] = dface_join_opp[i_face_jn] + face_per_block_offset[i_zone_opp];*/
      /*}*/
    }
  }
  assert (idx == n_face_join);

  PDM_log_trace_array_long(face_per_block_offset, n_zone+1, "face_per_block_offset :: ");
  PDM_log_trace_array_long(vtx_per_block_offset,  n_zone+1, "vtx_per_block_offset :: ");
  PDM_log_trace_array_long(multi_gnum, n_face_join, "multi_gnum :: ");

  PDM_g_num_t **all_face_distribution = malloc(n_zone * sizeof(PDM_g_num_t*));
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    all_face_distribution[i_zone] = PDM_compute_entity_distribution(comm, dn_face[i_zone]);
  }

  PDM_multi_block_to_part_t *mptb = PDM_multi_block_to_part_create(face_per_block_offset,
                                                                   n_zone,
                                            (const PDM_g_num_t **) all_face_distribution,
                                            (const PDM_g_num_t **)&multi_gnum,
                                                                  &n_face_join,
                                                                   1,
                                                                   comm);
  //Prepare data to send : face -> vtx connectivity 
  int **face_vtx_n       = malloc(n_zone * sizeof(int*));
  PDM_g_num_t **face_vtx_shifted = malloc(n_zone * sizeof(int*));
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    face_vtx_n[i_zone]       = malloc(dn_face[i_zone] * sizeof(int));
    face_vtx_shifted[i_zone] = malloc(dface_vtx_idx[i_zone][dn_face[i_zone]] * sizeof(PDM_g_num_t));

    for (int i_face = 0; i_face < dn_face[i_zone]; i_face++) {
      face_vtx_n[i_zone][i_face] = dface_vtx_idx[i_zone][i_face+1] - dface_vtx_idx[i_zone][i_face];
      for (int j = dface_vtx_idx[i_zone][i_face]; j < dface_vtx_idx[i_zone][i_face+1]; j++) {
        face_vtx_shifted[i_zone][j] = dface_vtx[i_zone][j] + vtx_per_block_offset[i_zone];
      }
    }
    //PDM_log_trace_array_long(face_vtx_shifted[i_zone], dface_vtx_idx[i_zone][dn_face[i_zone]], "face_vtx_shifted :: ");
  }
  //int **face_vtx = malloc(n_zone * sizeof(int*));



  int         **part_stride = NULL;
  PDM_g_num_t **part_data   = NULL;
  PDM_multi_block_to_part_exch2(mptb,
                                sizeof(PDM_g_num_t),
                                PDM_STRIDE_VAR,
                                face_vtx_n,
                 (void **)      face_vtx_shifted,
                               &part_stride,
                 (void ***)    &part_data);

  *face_vtx_both_idx = PDM_array_new_idx_from_sizes_int(part_stride[0], n_face_join);
  *face_vtx_both     = part_data[0];

  free(part_data);
  free(part_stride);

  PDM_multi_block_to_part_free(mptb);

  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    free(all_face_distribution[i_zone]);
    free(face_vtx_n[i_zone]);
    free(face_vtx_shifted[i_zone]);
  }
  free(all_face_distribution);
  free(face_vtx_n);
  free(face_vtx_shifted);
  free(multi_gnum);

  free(face_per_block_offset);
  free(vtx_per_block_offset);

  return n_face_join;
}

int _generate_edge_face
(
 int            n_face,
 int           *face_vtx_idx,
 PDM_g_num_t   *face_vtx,
 PDM_g_num_t  **dedge_distrib,
 int          **dedge_vtx_idx,
 PDM_g_num_t  **dedge_vtx,
 int          **dedge_face_idx,
 PDM_g_num_t  **dedge_face,
 PDM_MPI_Comm  comm
)
{

  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  
  // 1. Get the number of unique vertex
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                     &face_vtx,
                                                      NULL,
                                                     &face_vtx_idx[n_face],
                                                      1,
                                                      comm);
  int n_vtx;
  int n_vtx_loc = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_MPI_Allreduce(&n_vtx_loc, &n_vtx, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
  PDM_part_to_block_free(ptb);

  // 2. Get (unmerged) edges
  int n_edge_unmerged = face_vtx_idx[n_face];
  int n_elmt_current  = 0;
  int n_edge_current  = 0;

  PDM_g_num_t *face_distri         = PDM_compute_entity_distribution(comm, n_face);
  PDM_g_num_t *tmp_edge_face       = (PDM_g_num_t *) malloc(n_edge_unmerged       * sizeof(PDM_g_num_t));
  int         *tmp_parent_elmt_pos = (int         *) malloc(n_edge_unmerged       * sizeof(int        ));
  int         *tmp_edge_vtx_idx    = (int         *) malloc((n_edge_unmerged + 1) * sizeof(int        ));
  PDM_g_num_t *tmp_edge_vtx        = (PDM_g_num_t *) malloc(2*n_edge_unmerged     * sizeof(PDM_g_num_t));

  tmp_edge_vtx_idx[0] = 0;
  PDM_poly2d_decomposes_edges(n_face,
                             &n_elmt_current,
                             &n_edge_current,
                              face_distri[i_rank],
                              -1,
                              face_vtx,
                              face_vtx_idx,
                              tmp_edge_vtx_idx, //Numéro de sommet des edges
                              tmp_edge_vtx,
                              tmp_edge_face,         // Numéro des edges pour une face
                              NULL,
                              NULL,
                              tmp_parent_elmt_pos);
  assert(n_edge_current == n_edge_unmerged);

  /*PDM_log_trace_array_long(tmp_dface_edge, n_edge_unmerged, "dface_edge :: ");*/
  /*PDM_log_trace_connectivity_long(tmp_dface_edge_vtx_idx, tmp_dface_edge_vtx, n_edge_unmerged, "dface_edge :: ");*/

  // 3. Merge shared edges
  PDM_g_num_t dn_edge;
  PDM_generate_entitiy_connectivity_raw(comm,
                                        n_vtx,
                                        n_edge_unmerged,
                                        tmp_edge_face,
                                        tmp_edge_vtx_idx,
                                        tmp_edge_vtx,
                                       &dn_edge,
                                        dedge_distrib,
                                        dedge_vtx_idx,
                                        dedge_vtx,
                                        dedge_face_idx,
                                        dedge_face);
  free(tmp_parent_elmt_pos);
  free(face_distri);
  return dn_edge;
}

static int _match_internal_edges
(
 int            dn_edge,
 int           *dedge_distrib,
 int           *dedge_face_idx,
 PDM_g_num_t   *dedge_face,
 int           *dedge_face_group_id,
 int           *dedge_face_group_id_opp,
 PDM_g_num_t   *dedge_face_join,
 PDM_g_num_t   *dedge_face_join_opp,
 PDM_g_num_t  **dedge_gnum,
 PDM_g_num_t  **dedge_gnum_opp,
 PDM_MPI_Comm  comm
)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  
  // 0. Count the number of internal edges
  int dn_internal_edge = 0;
  for(int i_edge = 0; i_edge < dn_edge; ++i_edge) {
    int n_face_this_edge = dedge_face_idx[i_edge+1] - dedge_face_idx[i_edge];
    dn_internal_edge += (int) (n_face_this_edge > 1);
  }
  log_trace("dn internal edges is %i \n", dn_internal_edge);
  log_trace("dn external edges is %i \n", dn_edge - dn_internal_edge);

  // 1. Build hash keys
  PDM_g_num_t *key_ln_to_gn = malloc(dn_internal_edge * sizeof(PDM_g_num_t)); 
  int         *stride_one   = malloc(dn_internal_edge * sizeof(int        ));
  int         *stride_two   = malloc(dn_internal_edge * sizeof(int        ));
  int         *stride_four  = malloc(dn_internal_edge * sizeof(int        ));

  PDM_g_num_t *data_send_connect    = malloc(4*dn_internal_edge * sizeof(PDM_g_num_t));
  PDM_g_num_t *data_send_edge_g_num = malloc(  dn_internal_edge * sizeof(PDM_g_num_t));
  PDM_g_num_t *data_send_group      = malloc(4*dn_internal_edge * sizeof(PDM_g_num_t));
  PDM_g_num_t *data_send_sens       = malloc(2*dn_internal_edge * sizeof(PDM_g_num_t));
  PDM_g_num_t *data_send_face_g_num = malloc(2*dn_internal_edge * sizeof(PDM_g_num_t));

  int i_int_edge = 0;
  int idx_write2 = 0;
  int idx_write4 = 0;
  for(int i_edge = 0; i_edge < dn_edge; ++i_edge) {

    int n_face_this_edge = dedge_face_idx[i_edge+1] - dedge_face_idx[i_edge];
    if (n_face_this_edge == 1) {
      continue;
    }
    assert (n_face_this_edge == 2);

    stride_one [i_int_edge] = 1;
    stride_two [i_int_edge] = 2;
    stride_four[i_int_edge] = 4;
    //Retrive zone id using group id of any of two faces
    data_send_edge_g_num[i_int_edge] = dedge_distrib[i_rank] + i_edge + 1;

    int key = 0;
    for(int j = dedge_face_idx[i_edge]; j < dedge_face_idx[i_edge+1]; ++j) { //Do it for the two faces data
      key += (dedge_face_join[j] + dedge_face_join_opp[j]);

      //data_send_face_g_num[idx_write2]   = dedge_face[j];
      //data_send_sens      [idx_write2++] = dedge_face_group_sens[j];
      idx_write2++;

      data_send_connect[idx_write4] = dedge_face_join    [j];
      data_send_group[idx_write4++] = dedge_face_group_id[j];

      data_send_connect[idx_write4] = dedge_face_join_opp[j];
      data_send_group[idx_write4++] = dedge_face_group_id_opp[j];
    }
    key_ln_to_gn[i_int_edge] = key;

    i_int_edge++;
  }
  assert(idx_write2 == 2*dn_internal_edge);
  assert(idx_write4 == 4*dn_internal_edge);

  PDM_log_trace_array_long(key_ln_to_gn, dn_internal_edge, "key_ln_to_gn :: ");

  // 2. Exchange data over hash key
  //Attention, pb d'équilibrage car les clés sont réparties vers la fin ... Un proc risque
  // de se retrouver avec tt les clés
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                    (PDM_g_num_t **) &key_ln_to_gn,
                                                      NULL,
                                                     &dn_internal_edge,
                                                      1,
                                                      comm);
  // Get protocol data
  int blk_size = PDM_part_to_block_n_elt_block_get(ptb);
  // This one will contain 2 for non conflicting hash key, more otherwise
  int *gnum_n_occurences = PDM_part_to_block_block_gnum_count_get(ptb);

  int gnum_n_occurences_tot = 0;
  for (int k = 0; k < blk_size; k++) {
    gnum_n_occurences_tot += gnum_n_occurences[k];
  }
  
  int         *unused_recv_stride = NULL;
  PDM_g_num_t *blk_edge_g_num     = NULL;
  int exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR,
                                     -1,
                           (int **)  &stride_one,
                           (void **) &data_send_edge_g_num,
                                     &unused_recv_stride,
                           (void **) &blk_edge_g_num);
  free(unused_recv_stride); // Same as gnum_n_occurences 
  assert (exch_size == gnum_n_occurences_tot);

  /*
  PDM_g_num_t* blk_data_sens   = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(int),
                                     PDM_STRIDE_VAR,
                                     -1,
                           (int **)  &stride_two,
                           (void **) &data_send_sens,
                                     &unused_recv_stride,
                           (void **) &blk_data_sens);
  free(unused_recv_stride); // Same as 2*gnum_n_occurences 
  assert (exch_size == 2*gnum_n_occurences_tot);
  */

  /*
  PDM_g_num_t* blk_data_face_g_num   = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR,
                                     -1,
                           (int **)  &stride_two,
                           (void **) &data_send_face_g_num,
                                     &unused_recv_stride,
                           (void **) &blk_data_face_g_num);
  free(unused_recv_stride); // Same as 2*gnum_n_occurences
  assert (exch_size == 2*gnum_n_occurences_tot);
  */

  PDM_g_num_t *blk_data_connect = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR,
                                     -1,
                          (int **)  &stride_four,
                          (void **) &data_send_connect,
                                    &unused_recv_stride,
                          (void **) &blk_data_connect);
  free(unused_recv_stride); // Same as 4*gnum_n_occurences 
  assert (exch_size == 4*gnum_n_occurences_tot);

  PDM_g_num_t *blk_data_group = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(int),
                                     PDM_STRIDE_VAR,
                                     -1,
                           (int **)  &stride_four,
                           (void **) &data_send_group,
                                     &unused_recv_stride,
                           (void **) &blk_data_group);
  free(unused_recv_stride); // Same as 4*gnum_n_occurences 
  assert (exch_size == 4*gnum_n_occurences_tot);


  PDM_log_trace_array_int(gnum_n_occurences   , blk_size               , "gnum_n_occurences   :: ");
  PDM_log_trace_array_long(blk_edge_g_num     , gnum_n_occurences_tot  , "blk_edge_g_num      :: ");
  //PDM_log_trace_array_long(blk_data_face_g_num, 2*gnum_n_occurences_tot, "blk_data_face_g_num :: ");
  //PDM_log_trace_array_long(blk_data_sens      , 2*gnum_n_occurences_tot, "blk_data_sens       :: ");
  PDM_log_trace_array_long(blk_data_connect   , 4*gnum_n_occurences_tot, "blk_data_connect    :: ");
  PDM_log_trace_array_long(blk_data_group     , 4*gnum_n_occurences_tot, "blk_data_group      :: ");


  free(key_ln_to_gn        );
  free(data_send_connect   );
  free(data_send_group     );
  free(data_send_sens      );
  free(data_send_face_g_num);
  free(data_send_edge_g_num);
  free(stride_one          );
  free(stride_two          );
  free(stride_four         );



  // 3. Post treatemement : resolve conflicting keys
  PDM_g_num_t *results_edge     = malloc(gnum_n_occurences_tot * sizeof(PDM_g_num_t));
  PDM_g_num_t *results_edge_opp = malloc(gnum_n_occurences_tot * sizeof(PDM_g_num_t));

  int n_max_entity_per_key = 0;
  for(int i = 0; i < blk_size; ++i) {
    n_max_entity_per_key = PDM_MAX(gnum_n_occurences   [i], n_max_entity_per_key);
  }
  int n_max_connec = 4*n_max_entity_per_key; //todo clean me


  log_trace("n_max_entity_per_key = %i \n", n_max_entity_per_key);
  log_trace("n_max_connec         = %i \n", n_max_connec);
  int *already_treat   = (int *) malloc(n_max_entity_per_key * sizeof(int));
  int *same_entity_idx = (int *) malloc(n_max_entity_per_key * sizeof(int));
  // int *sens_entity     = (int *) malloc(n_max_entity_per_key * sizeof(int));


  int idx  = 0;
  int idx_w = 0;
  for(int i_key = 0; i_key < blk_size; ++i_key) {

    int n_matching_edge = gnum_n_occurences[i_key];

    log_trace(" i_key = %i | n_matching_edge = %i \n", i_key, n_matching_edge);

    /* Reset */
    PDM_array_reset_int(already_treat, n_matching_edge, -1);

    /* Loop over all entitys in conflict and sort all */
    for(int i_entity = 0; i_entity < n_matching_edge; ++i_entity) {

      //Each internal edge comes with 2 faces -> 4 data
      int beg    = 4*(idx+i_entity);

      // Caution inplace sort !!!!
      PDM_sort_long(&blk_data_connect[beg], NULL, 4);
      PDM_sort_int (&blk_data_group  [beg], NULL, 4);

      if(0 == 1) {
        PDM_log_trace_array_long(&blk_data_connect[beg], 4, "blk_data_connect (sort) :: ");
        PDM_log_trace_array_int (&blk_data_group  [beg], 4, "blk_data_group   (sort) :: ");
      }
    }

    /*
     *  Identify pair or invalid other //todo : shortcut only 2 occurences
     */
    for(int i_entity = 0; i_entity < n_matching_edge; ++i_entity) {
      int beg1    = 4*(idx+i_entity);

      int n_same         = 0;
      int i_entity2_same = -1;

      if(already_treat[i_entity] == 1) {
        continue;
      }

      for(int i_entity2 = i_entity+1; i_entity2 < n_matching_edge; ++i_entity2) {

        if(already_treat[i_entity2] == -1) {
          int beg2    = 4*(idx+i_entity2);

          if(!PDM_array_are_equal_int(&blk_data_group[beg1], &blk_data_group[beg2], 4)) {
            continue;
          }

          if(!PDM_array_are_equal_int(&blk_data_connect[beg1], &blk_data_connect[beg2], 4)) {
            continue;
          }

          already_treat[i_entity2] = 1;
          same_entity_idx[n_same++] = i_entity2;
          i_entity2_same = i_entity2;

        }
      } /* End for i_entity2 */
      assert(n_same == 1);


      int edge_idx     = (idx+i_entity);
      int edge_idx_opp = (idx+i_entity2_same);

      // Set data for edge
      results_edge    [idx_w] = blk_edge_g_num[edge_idx];
      results_edge_opp[idx_w++] = blk_edge_g_num[edge_idx_opp];

      // Set data for opposite edge
      results_edge    [idx_w] = blk_edge_g_num[edge_idx_opp];
      results_edge_opp[idx_w++] = blk_edge_g_num[edge_idx];

      already_treat[i_entity] = 1;
    }
    idx  += n_matching_edge;
  }
  assert (idx_w == gnum_n_occurences_tot);
  free(already_treat);
  free(same_entity_idx);

  free(blk_edge_g_num);
  free(blk_data_connect);
  free(blk_data_group);
  PDM_part_to_block_free(ptb); // Needed until here for gnum_n_occurences

  log_trace("Conflict resolved, gnum are\n");
  PDM_log_trace_array_long(results_edge, gnum_n_occurences_tot, "edge gnum ::");
  PDM_log_trace_array_long(results_edge_opp, gnum_n_occurences_tot, "edge gnum opp::");



  // 4. Send back result on edge distribution
                       ptb = PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_NOTHING,
                                                       1.,
                                                      &results_edge,
                                                       dedge_distrib,
                                                      &gnum_n_occurences_tot,
                                                       1,
                                                       comm);

  assert (PDM_part_to_block_n_elt_block_get(ptb) == dn_internal_edge);
  *dedge_gnum = malloc(dn_internal_edge * sizeof(PDM_g_num_t));

  PDM_g_num_t *dedge_gnum_tmp = PDM_part_to_block_block_gnum_get(ptb);
  memcpy(*dedge_gnum, dedge_gnum_tmp, dn_internal_edge*sizeof(PDM_g_num_t));

  PDM_part_to_block_exch(ptb,
                        sizeof(PDM_g_num_t),
                        PDM_STRIDE_CST,
                        1,
                        NULL,
             (void **) &(results_edge_opp),
                        NULL,
              (void **) dedge_gnum_opp);
  PDM_part_to_block_free(ptb);

  free(results_edge);
  free(results_edge_opp);

  return dn_internal_edge;
}

void _match_all_edges_from_faces
(
  int          dn_face,
  int         *face_edge_idx,
  PDM_g_num_t *face_edge,
  PDM_g_num_t *face_edge_wopp,
  PDM_g_num_t *pedge_vtx,
  PDM_g_num_t *p_all_vtx,
  PDM_g_num_t *p_all_vtx_opp
)
{
  // Avoid multiple malloc using max size
  int max_face_len = 0;
  for (int i_face = 0; i_face < dn_face/2; i_face++) {
    max_face_len = PDM_MAX(max_face_len, face_edge_idx[2*i_face+1] - face_edge_idx[2*i_face]);
  }
  PDM_g_num_t *ordered_edge     = malloc(max_face_len * sizeof(PDM_g_num_t));
  PDM_g_num_t *ordered_edge_opp = malloc(max_face_len * sizeof(PDM_g_num_t));
  PDM_g_num_t *ordered_vtx      = malloc(max_face_len * sizeof(PDM_g_num_t));
  PDM_g_num_t *ordered_vtx_opp  = malloc(max_face_len * sizeof(PDM_g_num_t));

  // Avec la construction des faces de bord, on a des paires faces / face opp
  assert (dn_face%2 == 0);
  assert (face_edge_idx[dn_face]%2 == 0);
  int glob_idx     = 0;
  int glob_idx_opp = face_edge_idx[dn_face]/2;
  for (int i_face = 0; i_face < dn_face/2; i_face++) {
    int start_idx     = face_edge_idx[2*i_face];
    int start_idx_opp = face_edge_idx[2*i_face+1];
    int face_len = face_edge_idx[2*i_face+1] - face_edge_idx[2*i_face];
    if (0 == 1) {
      log_trace("iface %d :\n", i_face);
      log_trace("face data\n");
      // PDM_log_trace_array_long(&face_vtx[face_edge_idx[i_face]], face_len, "face_vtx");
      PDM_log_trace_array_long(&face_edge[start_idx], face_len, "face_edge");
      PDM_log_trace_array_long(&pedge_vtx[2*start_idx], 2*face_len, "edge vertices");
      PDM_log_trace_array_long(&face_edge_wopp[start_idx], face_len, "face_edge_wopp");
      log_trace("face opp data\n");
      // PDM_log_trace_array_long(&face_vtx_opp[face_edge_idx[i_face]], face_len, "face_vtx");
      PDM_log_trace_array_long(&face_edge[start_idx_opp], face_len, "face_opp_edge_opp");
      PDM_log_trace_array_long(&pedge_vtx[2*start_idx_opp], 2*face_len, "edge vertices opp");
    }

    // Search any received edge (we should have at least one)
    PDM_g_num_t opp_edge_key = 0;
    int idx = -1;
    while (opp_edge_key == 0)
      opp_edge_key = face_edge_wopp[start_idx + 1 + idx++];
    assert (idx < face_len);

    //Search idx of opposite edge in opposite 
    PDM_g_num_t candidate = 0;
    int opp_idx = -1;
    while (candidate != opp_edge_key) {
      candidate = face_edge[start_idx_opp + 1 + opp_idx++];
      candidate = PDM_ABS(candidate);
    }
    assert (opp_idx < face_len);

    // Find starting point using the two edge indices and signs
    int sign     = PDM_SIGN(face_edge[start_idx + idx]);
    int opp_sign = PDM_SIGN(face_edge[start_idx_opp + opp_idx]);

    int next_vtx, next_vtx_opp, cur_vtx, cur_vtx_opp;
    if (sign == 1) {
      cur_vtx  = pedge_vtx[2*(start_idx + idx)];
      next_vtx = pedge_vtx[2*(start_idx + idx)+1];
    }
    else {
      cur_vtx  = pedge_vtx[2*(start_idx + idx)+1];
      next_vtx = pedge_vtx[2*(start_idx + idx)];
    }
    if (opp_sign == 1) {//Invert looping order for opposite face
      cur_vtx_opp  = pedge_vtx[2*(start_idx_opp + opp_idx)+1];
      next_vtx_opp = pedge_vtx[2*(start_idx_opp + opp_idx)];
    }
    else {
      cur_vtx_opp  = pedge_vtx[2*(start_idx_opp + opp_idx)];
      next_vtx_opp = pedge_vtx[2*(start_idx_opp + opp_idx)+1];
    }

    for (int i = 0; i < face_len; i++) {
      //Fill
      // log_trace("Cur vtx %d and opp %d \n", cur_vtx, cur_vtx_opp);
      ordered_edge[i]     = PDM_ABS(face_edge[start_idx+idx]);
      ordered_edge_opp[i] = PDM_ABS(face_edge[start_idx_opp+opp_idx]);
      ordered_vtx[i]      = cur_vtx;
      ordered_vtx_opp[i]  = cur_vtx_opp;
      //This is for face
      for (int j = 0; j < face_len; j++) {
        if (j != idx) {
          int vtx1 = pedge_vtx[2*(start_idx + j)];
          int vtx2 = pedge_vtx[2*(start_idx + j)+1];
          if (vtx1 == next_vtx) {
            idx = j;
            cur_vtx  = vtx1;
            next_vtx = vtx2;
            break;
          }
          else if (vtx2 == next_vtx) {
            idx = j;
            cur_vtx  = vtx2;
            next_vtx = vtx1;
            break;
          }
        }
      }
      //This is for opposite face
      for (int j = 0; j < face_len; j++) {
        if (j != opp_idx) {
          int vtx1 = pedge_vtx[2*(start_idx_opp+j)];
          int vtx2 = pedge_vtx[2*(start_idx_opp+j)+1];
          if (vtx1 == next_vtx_opp) {
            opp_idx = j;
            cur_vtx_opp  = vtx1;
            next_vtx_opp = vtx2;
            break;
          }
          else if (vtx2 == next_vtx_opp) {
            opp_idx = j;
            cur_vtx_opp  = vtx2;
            next_vtx_opp = vtx1;
            break;
          }
        }
      }
    }
    if (0 == 1) {
      PDM_log_trace_array_long(ordered_edge,     face_len, "ordered edges");
      PDM_log_trace_array_long(ordered_edge_opp, face_len, "ordered edges_opp");
      PDM_log_trace_array_long(ordered_vtx,     face_len, "ordered vtx");
      PDM_log_trace_array_long(ordered_vtx_opp, face_len, "ordered vtx_opp");
    }
    //Copy results for this face
    memcpy(&p_all_vtx[glob_idx],      ordered_vtx,      face_len*sizeof(PDM_g_num_t));
    memcpy(&p_all_vtx_opp[glob_idx],  ordered_vtx_opp,  face_len*sizeof(PDM_g_num_t));
    //memcpy(&p_all_edge[glob_idx],     ordered_edge      face_len*sizeof(PDM_g_num_t));
    //memcpy(&p_all_edge_opp[glob_idx], ordered_edge_opp, face_len*sizeof(PDM_g_num_t));
    //Copy results for opposite face
    memcpy(&p_all_vtx[glob_idx_opp],       ordered_vtx_opp,  face_len*sizeof(PDM_g_num_t));
    memcpy(&p_all_vtx_opp[glob_idx_opp],   ordered_vtx,      face_len*sizeof(PDM_g_num_t));
    //memcpy(&p_all_edge[glob_idx_opp],      ordered_edge_opp, face_len*sizeof(PDM_g_num_t));
    //memcpy(&p_all_edge_opp[glob_idx_opp],  ordered_edge_,    face_len*sizeof(PDM_g_num_t));
    glob_idx     += face_len;
    glob_idx_opp += face_len;
  }
  free(ordered_edge);
  free(ordered_edge_opp);
  free(ordered_vtx);
  free(ordered_vtx_opp);
}

static void
_create_vtx_join
(
int            n_group_join,
int            p_all_vtx_n,
int           *p_all_vtx,
int           *p_all_vtx_opp,
int           *p_all_vtx_group,
int          **dvtx_group_idx,
PDM_g_num_t  **dvtx_group,
PDM_g_num_t  **dvtx_group_opp,
PDM_MPI_Comm   comm
)
{
  //Todo : we could shift back to position of vtx in extraction to have a better
  //balance of edge distribution
  int *stride_one = PDM_array_const_int(p_all_vtx_n, 1);
  //Merge vertices using gnum -- exchange opp_gnum & group id
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                     &p_all_vtx,
                                                      NULL,
                                                     &p_all_vtx_n,
                                                      1,
                                                      comm);

  int blk_size = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_g_num_t *dall_vtx   = PDM_part_to_block_block_gnum_get(ptb);

  int         *recv_stride    = NULL;
  PDM_g_num_t *dall_vtx_opp   = NULL;
  PDM_g_num_t *dall_vtx_group = NULL;
  int exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR,
                                     -1,
                           (int  **) &stride_one,
                           (void **) &p_all_vtx_opp,
                                     &recv_stride,
                           (void **) &dall_vtx_opp);

  int *unused_recv_stride = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR,
                                     -1,
                           (int  **) &stride_one,
                           (void **) &p_all_vtx_group,
                                     &unused_recv_stride, //Same  than recv stride
                           (void **) &dall_vtx_group);
  free(unused_recv_stride);

  PDM_log_trace_array_long(dall_vtx,       blk_size,  "dall_vtx            :");
  PDM_log_trace_array_long(recv_stride,    blk_size,  "recv stride         :");
  PDM_log_trace_array_long(dall_vtx_opp,   exch_size, "recv dall_vtx_opp   :");
  PDM_log_trace_array_long(dall_vtx_group, exch_size, "recv dall_vtx_group :");
  free(stride_one);

  //Post treat vertex data
  int *all_order    = malloc(exch_size*sizeof(int));
  int *group_id_tmp = malloc(exch_size * sizeof(int));
  memcpy(group_id_tmp, dall_vtx_group, exch_size*sizeof(int));

  int start_vtx = 0;
  for (int i_vtx =  0; i_vtx < blk_size; i_vtx++) {
    int n_recv = recv_stride[i_vtx];
    //order array must be initialized
    for (int j = 0; j < n_recv; j++)
      all_order[start_vtx + j] = j;
    PDM_sort_int(&group_id_tmp[start_vtx], &all_order[start_vtx], n_recv);
    start_vtx += n_recv;
  }
  free(group_id_tmp);
  PDM_log_trace_array_int(all_order, exch_size, "all order ");
  
  //First pass to count
  int *_dvtx_group_n   = PDM_array_zeros_int(n_group_join);
  start_vtx = 0;
  for (int i_vtx = 0; i_vtx < blk_size; i_vtx++) {
    int n_recv = recv_stride[i_vtx];
    int *_dall_vtx_group = &dall_vtx_group[start_vtx];
    int *loc_order       = &all_order     [start_vtx];

    _dvtx_group_n[_dall_vtx_group[loc_order[0]]]++; //Vtx should be related to at least one face, so its ok
    for (int i = 1; i < n_recv; i++) {
      if (_dall_vtx_group[loc_order[i]] != _dall_vtx_group[loc_order[i-1]])
        _dvtx_group_n[_dall_vtx_group[loc_order[i]]]++;
    }
    start_vtx += n_recv;
  }

  //Second pass to fill
  int         *_dvtx_group_idx = PDM_array_new_idx_from_sizes_int(_dvtx_group_n, n_group_join);
  PDM_g_num_t *_dvtx_group     = malloc(_dvtx_group_idx[n_group_join]*sizeof(PDM_g_num_t));
  PDM_g_num_t *_dvtx_group_opp = malloc(_dvtx_group_idx[n_group_join]*sizeof(PDM_g_num_t));
  PDM_array_reset_int(_dvtx_group_n, n_group_join, 0); // reset to count

  start_vtx = 0;
  for (int i_vtx =  0; i_vtx < blk_size; i_vtx++) {
    int n_recv = recv_stride[i_vtx];
    int *_dall_vtx_group = &dall_vtx_group[start_vtx];
    int *loc_order       = &all_order     [start_vtx];

    PDM_g_num_t vtx_gnum     = dall_vtx[i_vtx];  //These one will be the same even if vertex appears more than once
    PDM_g_num_t opp_vtx_gnum = dall_vtx_opp[start_vtx];

    int i_group = _dall_vtx_group[loc_order[0]];
    _dvtx_group    [_dvtx_group_idx[i_group] + _dvtx_group_n[i_group]] = vtx_gnum;
    _dvtx_group_opp[_dvtx_group_idx[i_group] + _dvtx_group_n[i_group]] = opp_vtx_gnum;
    _dvtx_group_n[i_group]++;
    for (int i = 1; i < n_recv; i++) {
      if (_dall_vtx_group[loc_order[i]] != _dall_vtx_group[loc_order[i-1]]) {
        i_group = _dall_vtx_group[loc_order[i]];
        _dvtx_group    [_dvtx_group_idx[i_group] + _dvtx_group_n[i_group]] = vtx_gnum;
        _dvtx_group_opp[_dvtx_group_idx[i_group] + _dvtx_group_n[i_group]] = opp_vtx_gnum;
        _dvtx_group_n[i_group]++;
      }
    }
    start_vtx += n_recv;
  }

  free(_dvtx_group_n);
  free(all_order);
  free(recv_stride);
  free(dall_vtx_opp);
  free(dall_vtx_group);
  PDM_part_to_block_free(ptb);

  *dvtx_group_idx = _dvtx_group_idx;
  *dvtx_group     = _dvtx_group;
  *dvtx_group_opp = _dvtx_group_opp;
}



static void
merge_zone
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
 PDM_g_num_t  **dface_cell,
 PDM_MPI_Comm   comm
)
{

  /*
   * On a plusieurs connectivités distribué  :
   *   -> On veut unifier ce qui est commun
   *   -> Et actualiser les connectivités
   *
   *  PDM_multi_block_to_part
   *    -> Avec face_ln_to_gn = implicite face_ln_to_gn - face à supprimer
   *          --> PDM_redistribute
   *
   */


  /*
   * Par connectivité ascendante ou descendane ?
   *    --> face_cell
   *    --> Par ascendance --> Trop dure je pense
   *   Dans pdm_mesh_adpation c'est en ascendant mais on part des vertex .
   *   Subtile on reconstruit dans ce sens --> vtx -> edge -> face -> cell
   *   Mais on construit les connectivités descendantes --> edge_vtx -> face_edge...
   */


  /*
   * Par descendance :
   *   --> Calcul du dcell_face --> Attention au signe !!!!!!!
   *   --> Quand on transpose la connectvity il faut d'aborder transformé le dface_cell en orienté
   *   --> This one : PDM_setup_connectivity_idx
   *
   *   PMD_multiblock_to_part avec dcell_face + cell_ln_to_gn concateante
   *     --> Donc un dcell_face re-repartie mais faux en face
   *     --> Si on est des ouf --> Reorder_block_hilbert
   *   --> Pour update le dcell_face
   *         --> mbtp avec ln_to_gn = face_to_keep
   *   --> On utilise le dface_jon + dface_join_opp
   *   Attention au desequilibre de charge possibleeeeeeeee
   *   Pour l'update des numero de face --> block_to_part avec un block variable de numero de faces old_to_new
   *     --> Si la face est supprimé strid = 0
   *
   */

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

  int n_group_join = 2*(n_zone-1);
  int          *dface_join_idx  = (int         *) malloc((n_group_join + 1) * sizeof(int        ));
  PDM_g_num_t  *dface_join      = NULL; // A allouer propremet

  //int n_group_join = 2*2*(n_zone-1); // Try 2*2 jns
  int *group_join_to_zone_cur = (int *) malloc( n_group_join * sizeof(int));
  int *group_join_to_zone_opp = (int *) malloc( n_group_join * sizeof(int));
  int *group_join_to_join_opp = (int *) malloc( n_group_join * sizeof(int));


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
  // Test 2*2 jns
  /*group_join_to_join_opp[0] = 2;*/
  /*group_join_to_join_opp[1] = 3;*/
  /*group_join_to_join_opp[2] = 0;*/
  /*group_join_to_join_opp[3] = 1;*/
  /*group_join_to_zone_opp[0] = 1;*/
  /*group_join_to_zone_opp[1] = 1;*/
  /*group_join_to_zone_opp[2] = 0;*/
  /*group_join_to_zone_opp[3] = 0;*/


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
    //Try different setup with 2*2 jns
    /*
    n_jn = 2;
    n_face_group[i_zone] += 1;
    int *dface_group_idx_new = malloc((n_face_group[i_zone]+2) * sizeof(int));
    //For zone0, jn is 3 ; for zone1, jn is 2
    int jn_group = -1;
    if (i_zone == 0) jn_group = 3;
    else if (i_zone == 1) jn_group = 2;
    for (int i = 0; i < n_face_group[i_zone]+1; i++) {
      if (i <= jn_group) {
        dface_group_idx_new[i] = dface_group_idx[i_zone][i];
      }
      else if (i==jn_group+1) {
        int nface_this_group = dface_group_idx[i_zone][i] - dface_group_idx[i_zone][i-1];
        dface_group_idx_new[i] = dface_group_idx[i_zone][i-1] + (nface_this_group / 2);
        dface_group_idx_new[i+1] = dface_group_idx[i_zone][i-1] + nface_this_group;
      }
      else {
        dface_group_idx_new[i] = dface_group_idx[i_zone][i-1];
      }
    }
    PDM_log_trace_array_int(dface_group_idx[i_zone], 6+1, "dfacegroupidx :");
    PDM_log_trace_array_int(dface_group_idx_new, 7+1, "dfacegroupidxnew :");
    PDM_log_trace_array_long(dface_group[i_zone], dface_group_idx[i_zone][6], "dfacegroup :");
    free(dface_group_idx[i_zone]);
    dface_group_idx[i_zone] = dface_group_idx_new;
    */
    //

    // Join numbering (left to right, increasing i_zone)

    dface_bnd_idx [i_zone] = (int *) malloc((n_bnd        + 1) * sizeof(int));

    // First pass to count and allocate
    int i_bnd = 1;
    int i_jn  = tmp_i_group_join+1;
    dface_bnd_idx[i_zone][0]  = 0;
    // dface_join_idx[0] = 0;
    for (int igroup = 0; igroup < n_face_group[i_zone]; igroup++) {
      int copy_to_bnd = (igroup != 2 || (igroup == 2 && i_zone == 0)) && (igroup != 3 || (igroup == 3 && i_zone == n_zone-1));
      /*int copy_to_bnd = -1; //Try 2*2 jns*/
      /*if (i_zone == 0) copy_to_bnd = (igroup != 3) && (igroup !=4);*/
      /*if (i_zone == 1) copy_to_bnd = (igroup != 2) && (igroup !=3);*/
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


    for (int i = tmp_i_group_join; i < i_jn-1; i++) {
      dface_join_idx[i+1] = dface_join_idx[i+1] + dface_join_idx[i];
    }

    /* A bit stupid but complicated to made it in ohter way for a clear test */
    dface_join = realloc(dface_join, dface_join_idx[i_jn-1] * sizeof(PDM_g_num_t));

    // Second pass to copy
    dface_bnd [i_zone] = (PDM_g_num_t *) malloc(dface_bnd_idx [i_zone][n_bnd        ] * sizeof(PDM_g_num_t));
    i_bnd = 0;
    i_jn  = tmp_i_group_join;
    for (int igroup = 0; igroup < n_face_group[i_zone]; igroup++) {
      int copy_to_bnd = (igroup != 2 || (igroup == 2 && i_zone == 0)) && (igroup != 3 || (igroup == 3 && i_zone == n_zone-1));
      /*int copy_to_bnd = -1; //Try 2*2 jns*/
      /*if (i_zone == 0) copy_to_bnd = (igroup != 3) && (igroup !=4);*/
      /*if (i_zone == 1) copy_to_bnd = (igroup != 2) && (igroup !=3);*/
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

  log_trace("Global join data (%d)\n", n_group_join);
  PDM_log_trace_array_int(group_join_to_zone_cur, n_group_join, "group_join_to_zone_cur :: ");
  PDM_log_trace_array_int(group_join_to_join_opp, n_group_join, "group_join_to_join_opp :: ");
  PDM_log_trace_array_int(group_join_to_zone_opp, n_group_join, "group_join_to_zone_opp :: ");

  PDM_log_trace_array_int(dface_join_idx, n_group_join+1, "dface_join_idx :: ");
  PDM_log_trace_array_long(dface_join    , dface_join_idx[n_group_join], "dface_join     :: ");
  PDM_log_trace_array_long(dface_join_opp, dface_join_idx[n_group_join], "dface_join_opp :: ");

  // New version begins

  
  PDM_g_num_t *face_per_block_offset = _per_block_offset(n_zone, dn_face, comm);
  PDM_g_num_t *vtx_per_block_offset  = _per_block_offset(n_zone, dn_vtx,  comm);

  int         *face_vtx_both_idx = NULL;
  PDM_g_num_t *face_vtx_both     = NULL;
  int n_face_join = _extract_and_shift_jn_faces(n_zone,
                              dn_face,
                              dn_vtx,
                              n_group_join,
                              group_join_to_join_opp,
                              group_join_to_zone_cur,
                              group_join_to_zone_opp,
                              dface_join_idx,
                              dface_join,
                              dface_join_opp,
                              dface_vtx_idx,
                              dface_vtx,
                             &face_vtx_both_idx,
                             &face_vtx_both,
                              comm);

  PDM_g_num_t *extracted_face_distri = PDM_compute_entity_distribution(comm, n_face_join);
  
  // Transport some data to the extracted faces
  int         *dextract_face_group_id = malloc(n_face_join*sizeof(int));
  PDM_g_num_t *dextract_face_join     = malloc(n_face_join*sizeof(PDM_g_num_t));
  PDM_g_num_t *dextract_face_join_opp = malloc(n_face_join*sizeof(PDM_g_num_t));
  int idx = 0;
  for (int i_join = 0; i_join < n_group_join; i_join++) {
    int i_join_opp = group_join_to_join_opp[i_join];
    int i_zone_cur = group_join_to_zone_cur[i_join];
    int i_zone_opp = group_join_to_zone_opp[i_join];

    if (i_join <= i_join_opp) {
      for (int i_face_jn = dface_join_idx[i_join]; i_face_jn < dface_join_idx[i_join+1]; i_face_jn++) {
        //Face data
        dextract_face_join    [idx] = dface_join[i_face_jn]; //Here unshifted
        dextract_face_join_opp[idx] = dface_join_opp[i_face_jn];
        dextract_face_group_id[idx++] = i_join;
        //Opp face data
        dextract_face_join    [idx] = dface_join_opp[i_face_jn];
        dextract_face_join_opp[idx] = dface_join[i_face_jn];
        dextract_face_group_id[idx++] = i_join_opp;
      }
    }
  }
  assert (idx == n_face_join);

  log_trace("Face vtx received after MBTP\n");
  PDM_log_trace_array_int (face_vtx_both_idx, n_face_join+1, "face_vtx_idx :: ");
  PDM_log_trace_array_long(face_vtx_both, face_vtx_both_idx[n_face_join], "face_vtx :: ");


  //This is probably useless
  int         *face_vtx_idx = NULL;
  PDM_g_num_t *face_vtx     = NULL;
  PDM_g_num_t *face_vtx_opp = NULL;
  _split_paired_connectivity(n_face_join, face_vtx_both_idx, face_vtx_both, &face_vtx_idx, &face_vtx, &face_vtx_opp);

  log_trace("Split face_vtx & face_vtx donor \n");
  PDM_log_trace_array_int(face_vtx_idx, (n_face_join / 2)+1,            "face_vtx_idx :: ");
  PDM_log_trace_array_long(face_vtx,     face_vtx_idx[n_face_join / 2], "face_vtx :: ");
  PDM_log_trace_array_long(face_vtx_opp, face_vtx_idx[n_face_join / 2], "face_vtx_opp :: ");

  //Now we have some pairs of faces (each pair appears only one) + face_vtx for this pairs

  //Generate edge numbering

  PDM_g_num_t *dedge_distrib  = NULL;
  int         *dedge_vtx_idx  = NULL;
  PDM_g_num_t *dedge_vtx      = NULL;
  int         *dedge_face_idx = NULL;
  PDM_g_num_t *dedge_face     = NULL;
  int dn_edge = _generate_edge_face(n_face_join,
                                    face_vtx_both_idx,
                                    face_vtx_both,
                                   &dedge_distrib,
                                   &dedge_vtx_idx,
                                   &dedge_vtx,
                                   &dedge_face_idx,
                                   &dedge_face,
                                    comm);
  
  int ext_dn_face = n_face_join;
  log_trace("Edges rebuild\n");
  PDM_log_trace_array_long(dedge_distrib, n_rank+1, "dedge_distri ::");
  //PDM_log_trace_array_long(dedge_vtx, 2*dn_edge, "dedge_vtx ::");
  //PDM_log_trace_connectivity_long(dedge_face_idx, dedge_face, dn_edge, "dedge_face ::");

  // Transport face data to edges

  //Prepare numbering
  PDM_g_num_t *dedge_face_abs = (PDM_g_num_t *) malloc(dedge_face_idx[dn_edge] * sizeof(PDM_g_num_t));
  int         *dedge_face_sgn = (int         *) malloc(dedge_face_idx[dn_edge] * sizeof(int        ));
  for(int i = 0; i < dedge_face_idx[dn_edge]; ++i) {
    dedge_face_abs[i] = PDM_ABS (dedge_face[i]);
    dedge_face_sgn[i] = PDM_SIGN(dedge_face[i]);
  }


  PDM_g_num_t *dedge_face_join         = malloc(dedge_face_idx[dn_edge] * sizeof(PDM_g_num_t));
  PDM_g_num_t *dedge_face_join_opp     = malloc(dedge_face_idx[dn_edge] * sizeof(PDM_g_num_t));
  int         *dedge_face_group_id     = malloc(dedge_face_idx[dn_edge] * sizeof(int        ));
  int         *dedge_face_group_id_opp = malloc(dedge_face_idx[dn_edge] * sizeof(int        ));
  PDM_block_to_part_t *btp = PDM_block_to_part_create(extracted_face_distri,
                               (const PDM_g_num_t **) &dedge_face_abs,
                                                      &dedge_face_idx[dn_edge],
                                                      1,
                                                      comm);
  int cst_stride = 1;
  PDM_block_to_part_exch(btp,
                         sizeof(int),
                         PDM_STRIDE_CST,
                         &cst_stride,
                (void *) dextract_face_group_id,
                         NULL,
             (void ** ) &dedge_face_group_id);

  //Group id can be recomputed instead of exchanged
  for (int i = 0; i < dedge_face_idx[dn_edge]; i++) {
    dedge_face_group_id_opp[i] = group_join_to_join_opp[dedge_face_group_id[i]];
  }
  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST,
                         &cst_stride,
                (void *) dextract_face_join,
                         NULL,
             (void ** ) &dedge_face_join);

  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST,
                         &cst_stride,
                (void *) dextract_face_join_opp,
                         NULL,
             (void ** ) &dedge_face_join_opp);

  PDM_block_to_part_free(btp);
  free(dedge_face_abs);
  free(dedge_face_sgn);


  log_trace("Transport data on edges\n");
  PDM_log_trace_array_int (dedge_face_idx,      dn_edge+1,               "dedge_face_idx      ::");
  PDM_log_trace_array_int (dedge_face_group_id, dedge_face_idx[dn_edge], "dedge_face_group_id ::");
  PDM_log_trace_array_long(dedge_face_join    , dedge_face_idx[dn_edge], "dedge_face_join     ::");
  PDM_log_trace_array_long(dedge_face_join_opp, dedge_face_idx[dn_edge], "dedge_face_join_opp ::");

  // Match internal edges
  PDM_g_num_t *dedge_gnum     = NULL;
  PDM_g_num_t *dedge_gnum_opp = NULL;
  int dn_internal_edge = _match_internal_edges(dn_edge,
                                               dedge_distrib,
                                               dedge_face_idx,
                                               dedge_face,
                                               dedge_face_group_id,
                                               dedge_face_group_id_opp,
                                               dedge_face_join,
                                               dedge_face_join_opp,
                                              &dedge_gnum,
                                              &dedge_gnum_opp,
                                               comm);

  log_trace("Internal edge matches after conflict resolution \n");
  PDM_log_trace_array_long(dedge_gnum, dn_internal_edge, "dedge gnum :: ");
  PDM_log_trace_array_long(dedge_gnum_opp, dn_internal_edge, "dedge gnum_opp :: ");

  // To match external edges, we will work on face distribution. We need face->edge connectivity

  log_trace("Generate dface->edge from dedge->face\n");
  int         *dface_edge_idx = NULL;
  PDM_g_num_t *dface_edge     = NULL;
  PDM_dconnectivity_transpose(comm,
                              dedge_distrib,
                              extracted_face_distri,
                              dedge_face_idx,
                              dedge_face,
                              1,
                             &dface_edge_idx,
                             &dface_edge);

  PDM_log_trace_connectivity_long(dface_edge_idx, dface_edge, ext_dn_face, "dface_edge :: ");

  //Transfert some data from the edges to the edges know by the faces

  //Prepare gnum and stride
  PDM_g_num_t *dface_edge_abs = (PDM_g_num_t *) malloc(dface_edge_idx[ext_dn_face] * sizeof(PDM_g_num_t));
  for (int i = 0; i < dface_edge_idx[ext_dn_face]; i++) {
    dface_edge_abs[i] = PDM_ABS(dface_edge[i]);
  }
  idx = 0;
  int *dedge_gnum_n = malloc(dn_edge*sizeof(int));
  for (int i = 0; i < dn_edge; i++) {
    if (i + dedge_distrib[i_rank] + 1 == dedge_gnum[idx]) {
      dedge_gnum_n[i] = 1;
      idx++;
    }
    else {
      dedge_gnum_n[i] = 0;
    }
  }
  //PDM_log_trace_array_int(dedge_gnum_n, dn_edge, "dedge_gnum_n");

                       btp = PDM_block_to_part_create(dedge_distrib,
                               (const PDM_g_num_t **) &dface_edge_abs,
                                                      &dface_edge_idx[ext_dn_face],
                                                      1,
                                                      comm);


  int         **recv_stride_tmp = NULL;
  PDM_g_num_t **recv_data_tmp   = NULL;
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          dedge_gnum_n,
                          dedge_gnum_opp,
                         &recv_stride_tmp,
              (void ***) &recv_data_tmp);
  int         *pedge_gnum_n   = recv_stride_tmp[0];
  PDM_g_num_t *pedge_gnum_opp = recv_data_tmp[0];
  free(recv_stride_tmp);
  free(recv_data_tmp);

  int stride2 = 2;
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST,
                         &stride2,
                          dedge_vtx,
                          NULL,
              (void ***) &recv_data_tmp);
  PDM_g_num_t *pedge_vtx = recv_data_tmp[0];
  free(recv_data_tmp);
  PDM_log_trace_array_long(pedge_vtx, 2*dface_edge_idx[ext_dn_face], "pedge_vtx");



  PDM_g_num_t *face_edge_wopp = malloc(dface_edge_idx[ext_dn_face]*sizeof(PDM_g_num_t));
  idx = 0;
  for (int i = 0; i < dface_edge_idx[ext_dn_face]; i++) {
    if (pedge_gnum_n[i] == 1)
      face_edge_wopp[i] = pedge_gnum_opp[idx++];
    else
      face_edge_wopp[i] = 0;
  }

  PDM_log_trace_connectivity_long(dface_edge_idx, face_edge_wopp, ext_dn_face, "dface_edge :: ");

  free(dface_edge_abs);
  free(dedge_gnum_n);
  free(pedge_gnum_n);
  free(pedge_gnum_opp);
  PDM_block_to_part_free(btp);

  //Match external edges
  PDM_g_num_t *p_all_vtx      = malloc(dface_edge_idx[ext_dn_face] * sizeof(PDM_g_num_t));
  PDM_g_num_t *p_all_vtx_opp  = malloc(dface_edge_idx[ext_dn_face] * sizeof(PDM_g_num_t));
  int *p_all_vtx_group  = malloc(dface_edge_idx[ext_dn_face] * sizeof(int));
  /*PDM_g_num_t *p_all_edge_gnum     = malloc(dface_edge_idx[ext_dn_face] * sizeof(PDM_g_num_t));*/
  /*PDM_g_num_t *p_all_edge_gnum_opp = malloc(dface_edge_idx[ext_dn_face] * sizeof(PDM_g_num_t));*/
  _match_all_edges_from_faces(ext_dn_face,
                              dface_edge_idx,
                              dface_edge,
                              face_edge_wopp,
                              pedge_vtx,
                              p_all_vtx,
                              p_all_vtx_opp);

  //Copy group from face to vertices
  int glob_idx = 0;
  for (int i_face = 0; i_face < ext_dn_face/2; i_face++) {
    int face_len = dface_edge_idx[2*i_face+1] - dface_edge_idx[2*i_face];
    for (int i = 0; i < face_len; i++) {
      p_all_vtx_group[glob_idx] = dextract_face_group_id[2*i_face];
      p_all_vtx_group[glob_idx + (dface_edge_idx[ext_dn_face]/2)] = dextract_face_group_id[2*i_face+1];
      glob_idx ++;
    }
  }


  free(face_edge_wopp);
  log_trace("edge matching on face distribution\n");
  /*PDM_log_trace_array_long(p_all_edge_gnum,     dface_edge_idx[ext_dn_face], "p_all_edge_gnum     ::");*/
  /*PDM_log_trace_array_long(p_all_edge_gnum_opp, dface_edge_idx[ext_dn_face], "p_all_edge_gnum_opp ::");*/
  PDM_log_trace_array_long(p_all_vtx,     dface_edge_idx[ext_dn_face], "p_all_vtx     ::");
  PDM_log_trace_array_long(p_all_vtx_opp, dface_edge_idx[ext_dn_face], "p_all_vtx_opp ::");

  int p_all_vtx_n = dface_edge_idx[ext_dn_face];
  int         *dvtx_group_idx = NULL;
  PDM_g_num_t *dvtx_group     = NULL;
  PDM_g_num_t *dvtx_group_opp = NULL;
  _create_vtx_join(n_group_join,
                   p_all_vtx_n,
                   p_all_vtx,
                   p_all_vtx_opp,
                   p_all_vtx_group,
                  &dvtx_group_idx,
                  &dvtx_group,
                  &dvtx_group_opp,
                   comm);


  //Ultimate step : go back to original vtx numbering. All we have to do is retrieve zone
  // and substract zone offset  for (int i_vtx = 0; i_vtx < dvtx_group_idx[n_group_join]; i_vtx++) {
  for (int i_join = 0; i_join < n_group_join; i_join++) {
    int zone_cur_offset = vtx_per_block_offset[group_join_to_zone_cur[i_join]];
    int zone_opp_offset = vtx_per_block_offset[group_join_to_zone_opp[i_join]];
    for (int i_vtx = dvtx_group_idx[i_join]; i_vtx < dvtx_group_idx[i_join+1]; i_vtx++) {
      dvtx_group[i_vtx]     -= zone_cur_offset;
      dvtx_group_opp[i_vtx] -= zone_opp_offset;
    }
  }
  log_trace("Vtx & vtx opp per group in vtx global numbering\n");
  PDM_log_trace_array_int(dvtx_group_idx, n_group_join+1, "dvtx_group_idx ::");
  PDM_log_trace_array_long(dvtx_group, dvtx_group_idx[n_group_join], "dvtx_group     ::");
  PDM_log_trace_array_long(dvtx_group_opp, dvtx_group_idx[n_group_join], "dvtx_group_opp ::");



  free(dface_edge_idx);
  free(dface_edge);

  free(dextract_face_join);
  free(dextract_face_join_opp);
  free(dextract_face_group_id);

  free(dedge_face_join);
  free(dedge_face_join_opp);
  free(dedge_face_group_id);
  free(dedge_face_group_id_opp);


  free(extracted_face_distri);

  free(face_vtx_both_idx);
  free(face_vtx_both);



  free(face_vtx_idx);
  free(face_vtx);
  free(face_vtx_opp);

  free(dedge_distrib);
  free(dedge_vtx_idx);
  free(dedge_vtx);
  free(dedge_face_idx);
  free(dedge_face);

  free(dedge_gnum);
  free(dedge_gnum_opp);

  free(p_all_vtx);
  free(p_all_vtx_opp);
  free(p_all_vtx_group);
  free(face_per_block_offset);
  free(vtx_per_block_offset);

  /* Free memory */
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    free(dface_bnd_idx [i_zone]);
    free(dface_bnd     [i_zone]);
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
  free(dface_join_opp);

  free(group_join_to_zone_cur);
  free(group_join_to_zone_opp);
  free(group_join_to_join_opp);




  PDM_MPI_Finalize();

  return 0;
}
