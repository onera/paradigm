/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/
#include <sys/resource.h>

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
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_doctree_priv.h"
#include "pdm_doctree.h"
#include "pdm_timer.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_sort.h"
#include "pdm_box.h"
#include "pdm_box_tree.h"
#include "pdm_box_priv.h"
#include "pdm_point_tree_seq_priv.h"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_binary_search.h"
#include "pdm_point_tree_seq.h"
#include "pdm_vtk.h"


#ifdef __cplusplus
extern "C"
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

static
void
_redistribute_pts_geom
(
 PDM_doctree_t        *doct,
 PDM_part_to_block_t **ptb_out,
 double              **dpts_coords_out
)
{
  /*
   * Redistribute all pts and impose hilbert ordering
   */
  int **stride_one = malloc(doct->n_part_cloud * sizeof(int *));
  int **weight     = malloc(doct->n_part_cloud * sizeof(int *));
  for(int i_part = 0; i_part < doct->n_part_cloud; ++i_part) {
    weight    [i_part] = malloc(doct->n_point_cloud[i_part] * sizeof(int));
    stride_one[i_part] = malloc(doct->n_point_cloud[i_part] * sizeof(int));
    for(int i = 0; i < doct->n_point_cloud[i_part]; ++i) {
      weight    [i_part][i] = 1;
      stride_one[i_part][i] = 1;
    }
  }

  PDM_part_to_block_t* ptb = PDM_part_to_block_geom_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_CLEANUP, // A voir avec merge mais attention au init_location
                                                           1.,
                                                           PDM_PART_GEOM_HILBERT,
                                                           doct->pts_coords,
                                                           doct->pts_g_num,
                                                           weight,
                                                           doct->n_point_cloud,
                                                           doct->n_part_cloud,
                                                           doct->comm);

  for(int i_part = 0; i_part < doct->n_part_cloud; ++i_part) {
    free(weight[i_part]);
  }
  free(weight);

  int dn_pts = PDM_part_to_block_n_elt_block_get(ptb);

  /* Il faut le faire en MERGE mais enlever les doublons de coords */
  int    *blk_coord_n   = NULL;
  double *tmp_blk_pts_coord = NULL;
  PDM_part_to_block_exch(ptb,
                         3 * sizeof(double),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         stride_one,
               (void **) doct->pts_coords,
                         &blk_coord_n,
               (void **) &tmp_blk_pts_coord);


  // Copy
  double *blk_pts_coord = malloc(3 * dn_pts * sizeof(double));
  int idx_read  = 0;
  for(int i = 0; i < dn_pts; ++i) {

    blk_pts_coord[3*i  ] = tmp_blk_pts_coord[3*idx_read  ];
    blk_pts_coord[3*i+1] = tmp_blk_pts_coord[3*idx_read+1];
    blk_pts_coord[3*i+2] = tmp_blk_pts_coord[3*idx_read+2];

    idx_read += blk_coord_n[i];
  }
  free(blk_coord_n);
  free(tmp_blk_pts_coord);

  /* Transport init_location - Attention au merge du ptb à faire */
  int have_init_location = 1;
  for(int i = 0; i < doct->n_part_cloud; ++i) {
    if(doct->pts_init_location[i] == NULL) {
      have_init_location = 0;
    }
  }

  if(have_init_location == 1) {
    int *blk_init_location_pts_n = NULL;
    int *blk_init_location_pts   = NULL;

    PDM_part_to_block_exch(ptb,
                           3 * sizeof(int),
                           PDM_STRIDE_VAR_INTERLACED,
                           1,
                           stride_one,
                 (void **) doct->pts_init_location,
                           &blk_init_location_pts_n,
                 (void **) &blk_init_location_pts);


    free(blk_init_location_pts_n);
    free(blk_init_location_pts);

  }

  for(int i_part = 0; i_part < doct->n_part_cloud; ++i_part) {
    free(stride_one[i_part]);
  }
  free(stride_one);

  *ptb_out         = ptb;
  *dpts_coords_out = blk_pts_coord;

}

// void
// _exchange_global_tree
// (
//  PDM_doctree_t        *doct,
//  int                   n_coarse_box,
//  double               *coarse_box_extents,
//  int                  *box_n_pts
// )
// {

// }


/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_doctree_t*
PDM_doctree_create
(
 PDM_MPI_Comm              comm,
 int                       dim,
 int                       n_part_cloud,
 double                   *global_extents,
 PDM_doctree_local_tree_t  local_tree_kind
)
{
  PDM_doctree_t* doct = (PDM_doctree_t *) malloc(sizeof(PDM_doctree_t));

  doct->comm = comm;
  doct->dim  = dim;

  doct->local_tree_kind = local_tree_kind;

  if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_OCTREE) {
    doct->coarse_depth_max          = 5;
    doct->coarse_points_in_leaf_max = 60;

    doct->local_depth_max          = 31;
    doct->local_points_in_leaf_max = 30;
    doct->local_tolerance          = 1e-6;

  } else if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_KDTREE){
    doct->coarse_depth_max          = 6;
    doct->coarse_points_in_leaf_max = 60;

    doct->local_depth_max          = 31;
    doct->local_points_in_leaf_max = 30;
    doct->local_tolerance          = 1e-6;
  } else {
    abort();
  }

  doct->coarse_tree   = NULL;
  doct->local_tree    = NULL;
  doct->shmem_tree    = NULL;

  doct->comm_dist_graph = PDM_MPI_COMM_NULL;
  doct->n_degree_in     = 0;
  doct->neighbor_in     = NULL;

  doct->comm_shared     = PDM_MPI_COMM_NULL;

  PDM_UNUSED(global_extents);

  doct->n_part_cloud      = n_part_cloud;
  doct->n_point_cloud     = malloc(n_part_cloud * sizeof(int          ));
  doct->pts_g_num         = malloc(n_part_cloud * sizeof(PDM_g_num_t *));
  doct->pts_coords        = malloc(n_part_cloud * sizeof(double      *));
  doct->pts_init_location = malloc(n_part_cloud * sizeof(int         *));

  for(int i = 0; i < n_part_cloud; ++i) {
    doct->n_point_cloud    [i] = 0;
    doct->pts_g_num        [i] = NULL;
    doct->pts_coords       [i] = NULL;
    doct->pts_init_location[i] = NULL;
  }

  // doct->solicitation_kind    = solicitation_kind;
  // doct->n_part               = n_part;
  // doct->n_entity             = malloc(n_part * sizeof(int          ));
  // doct->entity_gnum          = malloc(n_part * sizeof(PDM_g_num_t *));
  // doct->entity_coords        = malloc(n_part * sizeof(double      *));
  // doct->init_location_entity = malloc(n_part * sizeof(int         *));

  // for(int i = 0; i < n_part; ++i) {
  //   doct->n_entity            [i] = 0;
  //   doct->entity_gnum         [i] = NULL;
  //   doct->entity_coords       [i] = NULL;
  //   doct->init_location_entity[i] = NULL;
  // }

  return doct;
}


void
PDM_doctree_build
(
 PDM_doctree_t     *doct
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank (doct->comm, &i_rank);
  PDM_MPI_Comm_size (doct->comm, &n_rank);

  /*
   * Prepare graphe comm for hybrid MPI-MPI
   */
  PDM_MPI_setup_hybrid_dist_comm_graph(doct->comm,
                                       &doct->comm_shared,
                                       &doct->comm_dist_graph,
                                       &doct->n_degree_in,
                                       &doct->neighbor_in);

  int n_rank_in_shm, i_rank_in_shm;
  PDM_MPI_Comm_rank (doct->comm_shared, &i_rank_in_shm);
  PDM_MPI_Comm_size (doct->comm_shared, &n_rank_in_shm);

  /*
   * Redistribute all pts
   */
  double              *blk_pts_coord = NULL;
  PDM_part_to_block_t *ptb           = NULL;
  _redistribute_pts_geom(doct, &ptb, &blk_pts_coord);

  PDM_g_num_t* distrib_pts = PDM_part_to_block_distrib_index_get(ptb);

  /*
   * Step 2 : Build local coarse tree
   */
  int dn_pts = distrib_pts[i_rank+1] - distrib_pts[i_rank];
  if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_OCTREE ||
     doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_KDTREE) {

    assert(doct->coarse_tree == NULL);
    doct->coarse_tree = PDM_point_tree_seq_create(doct->local_tree_kind,
                                                  doct->coarse_depth_max,
                                                  doct->coarse_points_in_leaf_max,
                                                  doct->local_tolerance);

    PDM_point_tree_seq_point_cloud_set(doct->coarse_tree,
                                       dn_pts,
                                       blk_pts_coord);

    PDM_point_tree_seq_build(doct->coarse_tree);
    if(1 == 1) {
      char filename[999];
      sprintf(filename, "out_coarse_tree_%i.vtk", i_rank);
      PDM_point_tree_seq_write_nodes(doct->coarse_tree, filename);
    }

  } else {
    abort();
  }

  /*
   * Extract extents on all local_tree
   */
  int n_coarse_box = 0;
  int    *coarse_box_id      = NULL;
  double *coarse_box_extents = NULL;
  int    *coarse_box_n_pts   = NULL; // Number of point in boxes
  int     n_depth_per_proc   = 0;
  if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_OCTREE) {
    n_depth_per_proc = 2;
  } else if (doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_KDTREE){
    n_depth_per_proc = 6; // 2^(depth)
  } else {
    abort();
  }

  PDM_point_tree_seq_extract_nodes(doct->coarse_tree,
                                   0,
                                   n_depth_per_proc,
                                   &n_coarse_box,
                                   &coarse_box_id,
                                   &coarse_box_extents,
                                   &coarse_box_n_pts);

  /*
   * Equilibrate among nodes/numa - To reduce memory footprint we set up data in shared memory
   */
  int *lrecv_count = malloc(doct->n_degree_in * sizeof(int));
  PDM_MPI_Neighbor_allgather(&n_coarse_box , 1, PDM_MPI_INT,
                             lrecv_count   , 1, PDM_MPI_INT, doct->comm_dist_graph);

  PDM_mpi_win_shared_t* wshared_local_nodes_n   = PDM_mpi_win_shared_create(n_rank  , sizeof(int), doct->comm_shared);
  PDM_mpi_win_shared_t* wshared_local_nodes_idx = PDM_mpi_win_shared_create(n_rank+1, sizeof(int), doct->comm_shared);
  int *shared_local_nodes_n   = PDM_mpi_win_shared_get(wshared_local_nodes_n);
  int *shared_local_nodes_idx = PDM_mpi_win_shared_get(wshared_local_nodes_idx);
  PDM_mpi_win_shared_lock_all (0, wshared_local_nodes_n  );
  PDM_mpi_win_shared_lock_all (0, wshared_local_nodes_idx);

  for(int i = 0; i < doct->n_degree_in; ++i) {
    shared_local_nodes_n[doct->neighbor_in[i]] = lrecv_count[i];
  }
  PDM_MPI_Barrier(doct->comm_shared);

  if(i_rank_in_shm == 0) {
    shared_local_nodes_idx[0] = 0;
    for(int i = 0; i < n_rank; ++i) {
      shared_local_nodes_idx[i+1] = shared_local_nodes_idx[i] + shared_local_nodes_n[i];
    }
  }
  PDM_MPI_Barrier(doct->comm_shared);

  if(1 == 1) {
    PDM_log_trace_array_int(shared_local_nodes_n  , n_rank , "shared_local_nodes_n   ::");
    PDM_log_trace_array_int(shared_local_nodes_idx, n_rank+1 , "shared_local_nodes_idx ::");
    PDM_log_trace_array_int(doct->neighbor_in, doct->n_degree_in , "doct->neighbor_in ::");
  }

  // Hook local recv_shift
  int *recv_shift = malloc(doct->n_degree_in * sizeof(int));
  for(int i = 0; i < doct->n_degree_in; ++i) {
    recv_shift[i] = shared_local_nodes_idx[doct->neighbor_in[i]];
  }

  /*
   * Exchange extents
   */
  PDM_mpi_win_shared_t* wshared_coarse_box_n_pts   = PDM_mpi_win_shared_create(    shared_local_nodes_idx[n_rank], sizeof(int)   , doct->comm_shared);
  PDM_mpi_win_shared_t* wshared_coarse_box_extents = PDM_mpi_win_shared_create(6 * shared_local_nodes_idx[n_rank], sizeof(double), doct->comm_shared);
  int    *shared_coarse_box_n_pts   = PDM_mpi_win_shared_get(wshared_coarse_box_n_pts  );
  double *shared_coarse_box_extents = PDM_mpi_win_shared_get(wshared_coarse_box_extents);
  PDM_mpi_win_shared_lock_all (0, wshared_coarse_box_n_pts  );
  PDM_mpi_win_shared_lock_all (0, wshared_coarse_box_extents);

  PDM_MPI_Neighbor_allgatherv(coarse_box_n_pts       , n_coarse_box, PDM_MPI_INT,
                              shared_coarse_box_n_pts, lrecv_count  , recv_shift, PDM_MPI_INT, doct->comm_dist_graph);
  PDM_MPI_Barrier(doct->comm_shared);


  /* Update */
  for(int i = 0; i < doct->n_degree_in; ++i) {
    recv_shift [i]  = shared_local_nodes_idx[doct->neighbor_in[i]]*6;
    lrecv_count[i] *= 6;
  }

  PDM_MPI_Neighbor_allgatherv(coarse_box_extents, 6 * n_coarse_box, PDM_MPI_DOUBLE,
                              shared_coarse_box_extents, lrecv_count        , recv_shift, PDM_MPI_DOUBLE, doct->comm_dist_graph);
  PDM_MPI_Barrier(doct->comm_shared);

  /*
   * Solicitate
   */
  int n_shared_boxes = shared_local_nodes_idx[n_rank];
  PDM_g_num_t* distrib_shared_boxes = PDM_compute_uniform_entity_distribution(doct->comm_shared, n_shared_boxes);

  PDM_box_set_t  *box_set   = NULL;
  PDM_box_tree_t *bt_shared = NULL;
  int   max_boxes_leaf_shared = 10; // Max number of boxes in a leaf for coarse shared BBTree
  int   max_tree_depth_shared = 6; // Max tree depth for coarse shared BBTree
  float max_box_ratio_shared  = 5; // Max ratio for local BBTree (nConnectedBoxe < ratio * nBoxes)

  PDM_MPI_Comm comm_alone;
  PDM_MPI_Comm_split(doct->comm, i_rank, 0, &(comm_alone));
  const int n_info_location = 3;
  int *init_location_proc = PDM_array_zeros_int (n_info_location * n_shared_boxes);

  PDM_mpi_win_shared_t* wshared_coarse_boxes_gnum = PDM_mpi_win_shared_create(    n_shared_boxes, sizeof(PDM_g_num_t), doct->comm_shared);
  PDM_mpi_win_shared_t* wshared_coarse_box_center        = PDM_mpi_win_shared_create(3 * n_shared_boxes, sizeof(double     ), doct->comm_shared);
  PDM_g_num_t    *shared_coarse_boxes_gnum   = PDM_mpi_win_shared_get(wshared_coarse_boxes_gnum  );
  double         *shared_box_center          = PDM_mpi_win_shared_get(wshared_coarse_box_center  );

  PDM_mpi_win_shared_lock_all (0, wshared_coarse_boxes_gnum  );
  PDM_mpi_win_shared_lock_all (0, wshared_coarse_box_center  );
  for (int i = distrib_shared_boxes[i_rank_in_shm]; i < distrib_shared_boxes[i_rank_in_shm+1]; i++) {
    shared_coarse_boxes_gnum[i] = i + 1;

    shared_box_center[3*i  ] = 0.5 * (shared_coarse_box_extents[6*i  ] + shared_coarse_box_extents[6*i+3]);
    shared_box_center[3*i+1] = 0.5 * (shared_coarse_box_extents[6*i+1] + shared_coarse_box_extents[6*i+4]);
    shared_box_center[3*i+2] = 0.5 * (shared_coarse_box_extents[6*i+2] + shared_coarse_box_extents[6*i+5]);

  }
  PDM_MPI_Barrier(doct->comm_shared);


  box_set = PDM_box_set_create(3,
                               0,  // No normalization to preserve initial extents
                               0,  // No projection to preserve initial extents
                               n_shared_boxes,
                               shared_coarse_boxes_gnum,
                               shared_coarse_box_extents,
                               1,
                               &n_shared_boxes,
                               init_location_proc,
                               comm_alone);

  bt_shared = PDM_box_tree_create (max_tree_depth_shared,
                                   max_boxes_leaf_shared,
                                   max_box_ratio_shared);

  PDM_box_tree_set_boxes (bt_shared,
                          box_set,
                          PDM_BOX_TREE_ASYNC_LEVEL);

  free(init_location_proc);

  int* coarse_tree_box_to_box_idx = NULL;
  int* coarse_tree_box_to_box     = NULL;
  if(doct->solicitation_kind == PDM_TREE_SOLICITATION_BOXES_POINTS) {
    assert(doct->n_part == 1);
    PDM_box_tree_intersect_boxes_boxes2(bt_shared,
                                        -1,
                                        doct->n_entity[0],
                                        doct->entity_coords[0],
                                        &coarse_tree_box_to_box_idx,
                                        &coarse_tree_box_to_box);

    if(0 == 1) {
      PDM_log_trace_connectivity_int(coarse_tree_box_to_box_idx,
                                     coarse_tree_box_to_box,
                                     n_shared_boxes,
                                     "coarse_tree_box_to_box : ");
    }
  } else {
    abort();
  }

  /*
   * Pour chaque shared box on connait le poids de la solitation
   *    -> part_to_block sur les shared
   *  Optim = can shared coarse_boxes_gnum -> distribution par noeuds
   */
  int    *weight = malloc(n_shared_boxes * sizeof(int));
  for(int i = 0; i < n_shared_boxes; ++i) {
    weight[i] = coarse_tree_box_to_box_idx[i+1] - coarse_tree_box_to_box_idx[i];
  }

  /*
   * Equilibrate boxes leaf with solicitation
   */
  PDM_part_to_block_t* ptb_equi_box = PDM_part_to_block_geom_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                    PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                                    1.,
                                                                    PDM_PART_GEOM_HILBERT,
                                                                    &shared_box_center,
                                                                    &shared_coarse_boxes_gnum,
                                                                    &weight,
                                                                    &n_shared_boxes,
                                                                    1,
                                                                    doct->comm);

  int          dn_equi_tree     = PDM_part_to_block_n_elt_block_get  (ptb_equi_box);
  PDM_g_num_t* parent_tree_gnum = PDM_part_to_block_block_gnum_get   (ptb_equi_box);
  PDM_g_num_t* distrib_tree     = PDM_part_to_block_distrib_index_get(ptb_equi_box);

  if(0 == 1) {
    PDM_log_trace_array_long(parent_tree_gnum, dn_equi_tree , "parent_tree_gnum :: ");
    PDM_log_trace_array_long(distrib_tree    , n_rank+1, "distrib_tree : ");
  }

  free(weight);

  PDM_mpi_win_shared_unlock_all(wshared_coarse_boxes_gnum);
  PDM_mpi_win_shared_unlock_all(wshared_coarse_box_center);
  PDM_mpi_win_shared_free (wshared_coarse_boxes_gnum);
  PDM_mpi_win_shared_free (wshared_coarse_box_center);

  /*
   * Update permutation of pre-solicitation
   */
  PDM_MPI_Neighbor_allgather(&dn_equi_tree, 1, PDM_MPI_INT,
                             lrecv_count  , 1, PDM_MPI_INT, doct->comm_dist_graph);

  for(int i = 0; i < doct->n_degree_in; ++i) {
    shared_local_nodes_n[doct->neighbor_in[i]] = lrecv_count[i];
  }
  PDM_MPI_Barrier(doct->comm_shared);

  if(i_rank_in_shm == 0) {
    shared_local_nodes_idx[0] = 0;
    for(int i = 0; i < n_rank; ++i) {
      shared_local_nodes_idx[i+1] = shared_local_nodes_idx[i] + shared_local_nodes_n[i];
    }
  }
  PDM_MPI_Barrier(doct->comm_shared);

  // Hook local recv_shift
  for(int i = 0; i < doct->n_degree_in; ++i) {
    recv_shift[i] = shared_local_nodes_idx[doct->neighbor_in[i]];
  }

  /*
   * Exchange dparent_gnum
   */
  PDM_mpi_win_shared_t* wshared_old_to_new_box_rank = PDM_mpi_win_shared_create(shared_local_nodes_idx[n_rank], sizeof(int)   , doct->comm_shared);
  int *shared_old_to_new_box_rank                   = PDM_mpi_win_shared_get(wshared_old_to_new_box_rank);
  PDM_mpi_win_shared_lock_all (0, wshared_old_to_new_box_rank);

  int *proc_id = malloc(dn_equi_tree * sizeof(int));
  for(int i = 0; i < dn_equi_tree; ++i) {
    proc_id[i] = i_rank;
  }

  // PDM_MPI_Neighbor_allgatherv(parent_tree_gnum       , dn_equi_tree, PDM__PDM_MPI_G_NUM,
  PDM_MPI_Neighbor_allgatherv(proc_id                , dn_equi_tree, PDM_MPI_INT,
                              shared_old_to_new_box_rank, lrecv_count , recv_shift, PDM_MPI_INT, doct->comm_dist_graph);
  PDM_MPI_Barrier(doct->comm_shared);

  if(1 == 1) {
    PDM_log_trace_array_int(shared_old_to_new_box_rank, n_shared_boxes, "shared_old_to_new_box_rank :");
  }
  free(proc_id);


  /*
   * Setup partitioning
   */
  PDM_g_num_t* impli_distrib_tree = PDM_compute_entity_distribution(doct->comm, n_coarse_box);

  if(1 == 1) {
    PDM_log_trace_array_long(impli_distrib_tree    , n_rank+1, "impli_distrib_tree : ");
  }
  PDM_block_to_part_t* btp = PDM_block_to_part_create(impli_distrib_tree,
                               (const PDM_g_num_t **) &parent_tree_gnum,
                                                      &dn_equi_tree,
                                                      1,
                                                      doct->comm);

  /*
   * Prepare buffer
   */
  int n_pts_tot = 0;
  int point_range[2];
  for(int i = 0; i < n_coarse_box; ++i ) {
    int node_id = coarse_box_id[i];
    if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_OCTREE ||
       doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_KDTREE) {
      n_pts_tot += PDM_point_tree_seq_point_range_get(doct->coarse_tree, node_id, point_range);
    } else {
      abort();
    }
  }

  int dn_blk = impli_distrib_tree[i_rank+1] - impli_distrib_tree[i_rank];

  PDM_g_num_t* blk_impli_pts_gnum = PDM_part_to_block_block_gnum_get(ptb);
  assert(dn_blk == n_coarse_box);

  double      *reorder_blk_coord_send = malloc(3 * n_pts_tot * sizeof(double     ));
  PDM_g_num_t *reorder_blk_pts_gnum   = malloc(    n_pts_tot * sizeof(PDM_g_num_t));
  int idx_write = 0;
  for(int i = 0; i < n_coarse_box; ++i ) {
    int node_id = coarse_box_id[i];

    double *sorted_tree_coord = NULL;
    int    *new_to_old_pts    = NULL;

    if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_OCTREE ||
       doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_KDTREE) {
      PDM_point_tree_seq_point_range_get(doct->coarse_tree, node_id, point_range);
      PDM_point_tree_seq_sorted_points_get(doct->coarse_tree, &sorted_tree_coord);
      PDM_point_tree_seq_point_new_to_old_get(doct->coarse_tree, &new_to_old_pts);
    } else {
      abort();
    }

    for(int i_pt = point_range[0]; i_pt < point_range[1]; ++i_pt) {
      reorder_blk_coord_send[3*idx_write  ] = sorted_tree_coord[3*i_pt  ];
      reorder_blk_coord_send[3*idx_write+1] = sorted_tree_coord[3*i_pt+1];
      reorder_blk_coord_send[3*idx_write+2] = sorted_tree_coord[3*i_pt+2];

      reorder_blk_pts_gnum[idx_write] = blk_impli_pts_gnum[new_to_old_pts[i_pt]];

      idx_write++;
    }
  }

  /*
   * Coordinates
   */
  double **tmp_equi_pts_coords = NULL;
  int    **tmp_equi_n_pts      = NULL;
  PDM_block_to_part_exch(btp,
                         3 * sizeof(double),
                         PDM_STRIDE_VAR_INTERLACED,
                         coarse_box_n_pts,
                         reorder_blk_coord_send,
                         &tmp_equi_n_pts,
              (void ***) &tmp_equi_pts_coords);
  double *equi_pts_coords = tmp_equi_pts_coords[0];
  int    *equi_n_pts      = tmp_equi_n_pts[0];
  free(tmp_equi_pts_coords);
  free(tmp_equi_n_pts);


  int equi_n_pts_tot = 0;
  for(int i = 0; i < dn_equi_tree; ++i) {
    equi_n_pts_tot += equi_n_pts[i];
  }

  /*
   *  Compute size for directly have shared_memory pts_gnum
   */
  int *shm_equi_pts_tot_idx = malloc( (n_rank_in_shm + 1) * sizeof(int));
  shm_equi_pts_tot_idx[0] = 0;
  PDM_MPI_Allgather(&equi_n_pts_tot         , 1, PDM_MPI_INT,
                    &shm_equi_pts_tot_idx[1], 1, PDM_MPI_INT, doct->comm_shared);

  for(int i = 0; i < n_rank_in_shm; ++i) {
    shm_equi_pts_tot_idx[i+1] += shm_equi_pts_tot_idx[i];
  }
  int n_equi_pts_shared_tot = shm_equi_pts_tot_idx[n_rank_in_shm];

  PDM_mpi_win_shared_t* wequi_pts_gnum = PDM_mpi_win_shared_create(n_equi_pts_shared_tot, sizeof(PDM_g_num_t), doct->comm_shared);
  PDM_g_num_t *equi_pts_gnum           = PDM_mpi_win_shared_get(wequi_pts_gnum);

  PDM_mpi_win_shared_lock_all (0, wequi_pts_gnum);
  PDM_g_num_t *lequi_pts_gnum = &equi_pts_gnum[shm_equi_pts_tot_idx[i_rank_in_shm]];

  /*
   * g_num
   */
  // PDM_g_num_t** tmp_equi_pts_gnum = NULL;
  PDM_block_to_part_exch_in_place(btp,
                                  sizeof(PDM_g_num_t),
                                  PDM_STRIDE_VAR_INTERLACED,
                                  coarse_box_n_pts,
                                  reorder_blk_pts_gnum,
                                  &equi_n_pts,
                        (void **) &lequi_pts_gnum);
  // PDM_g_num_t *equi_pts_gnum = tmp_equi_pts_gnum[0];
  // free(tmp_equi_pts_gnum);
  // free(tmp_equi_n_pts[0]);
  // free(tmp_equi_n_pts);

  /*
   * Init location -> Attention si gnum de points dupliqué -> variable
   */
  PDM_g_num_t* equi_pts_init_location = NULL;

  if(0 == 1) {
    char filename[999];
    sprintf(filename, "out_equi_pts_%i.vtk", i_rank);
    PDM_vtk_write_point_cloud(filename, equi_n_pts_tot, equi_pts_coords, equi_pts_gnum, NULL);
  }


  free(reorder_blk_coord_send);
  free(reorder_blk_pts_gnum);
  PDM_block_to_part_free(btp);
  PDM_part_to_block_free(ptb_equi_box);

  free(distrib_shared_boxes);
  free(impli_distrib_tree);

  PDM_MPI_Comm_free(&comm_alone);

  PDM_box_set_destroy (&box_set);
  PDM_box_tree_destroy (&bt_shared);

  PDM_mpi_win_shared_unlock_all(wshared_coarse_box_n_pts);
  PDM_mpi_win_shared_unlock_all(wshared_coarse_box_extents);


  PDM_mpi_win_shared_free(wshared_coarse_box_n_pts);
  PDM_mpi_win_shared_free(wshared_coarse_box_extents);

  free(recv_shift);
  free(lrecv_count);

  free(coarse_box_id);
  free(coarse_box_n_pts);
  free(coarse_box_extents);

  free(blk_pts_coord);

  PDM_part_to_block_free(ptb);

  if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_OCTREE ||
     doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_KDTREE) {
    if(doct->coarse_tree != NULL) {
      PDM_point_tree_seq_free(doct->coarse_tree);
    }
  } else {
    abort();
  }

  // Preparation of send count and box_rank/box_rank_idx
  //
  int *send_entity_n = PDM_array_zeros_int(n_rank);
  int *recv_entity_n = malloc (sizeof(int) * n_rank);
  for(int i = 0; i < n_rank; ++i) {
    for(int j = shared_local_nodes_idx[i]; j < shared_local_nodes_idx[i+1]; ++j) {
      send_entity_n[shared_old_to_new_box_rank[j]] += coarse_tree_box_to_box_idx[j+1] - coarse_tree_box_to_box_idx[j];
    }
  }
  PDM_MPI_Alltoall (send_entity_n, 1, PDM_MPI_INT,
                    recv_entity_n, 1, PDM_MPI_INT,
                    doct->comm);

  /*
   * Setup shared
   */
  int *send_entity_idx = malloc ( ( n_rank + 1) * sizeof(int));
  int *recv_entity_idx = malloc ( ( n_rank + 1) * sizeof(int));
  send_entity_idx[0] = 0;
  recv_entity_idx[0] = 0;

  int* shared_recv_count  = malloc(n_rank_in_shm * sizeof(int));

  int n_tot_recv = 0;
  send_entity_idx[0] = 0;
  recv_entity_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    n_tot_recv += recv_entity_n[i];
    send_entity_idx[i+1] = send_entity_idx[i] + send_entity_n[i];
    recv_entity_idx[i+1] = recv_entity_idx[i] + recv_entity_n[i];
  }

  PDM_MPI_Allgather(&n_tot_recv,       1, PDM_MPI_INT,
                    shared_recv_count, 1, PDM_MPI_INT,
                    doct->comm_shared);


  int *shared_recv_idx = malloc((n_rank_in_shm+1) * sizeof(int));
  shared_recv_idx[0] = 0;
  for(int i = 0; i < n_rank_in_shm; ++i) {
    shared_recv_idx[i+1] = shared_recv_idx[i] + shared_recv_count[i];
  }

  if(0 == 1) {
    PDM_log_trace_array_int(send_entity_n, n_rank, "send_entity_n :: ");
    PDM_log_trace_array_int(recv_entity_n, n_rank, "recv_entity_n :: ");
  }

  int n_tot_recv_shared = shared_recv_idx[n_rank_in_shm];
  PDM_mpi_win_shared_t* wshared_entity_coord         = NULL;
  PDM_mpi_win_shared_t* wshared_entity_gnum          = NULL;
  PDM_mpi_win_shared_t* wshared_entity_init_location = NULL;

  PDM_MPI_Request req_entity_gnum          = -1;
  PDM_MPI_Request req_entity_coord         = -1;
  PDM_MPI_Request req_entity_init_location = -1;

  PDM_g_num_t *send_g_num         = malloc (    send_entity_idx[n_rank] * sizeof(PDM_g_num_t));
  double      *send_extents       = malloc (6 * send_entity_idx[n_rank] * sizeof(double     ));
  int         *send_init_location = malloc (3 * send_entity_idx[n_rank] * sizeof(double     ));

  PDM_MPI_Datatype mpi_entity_type;
  PDM_MPI_Datatype mpi_init_location_type;

  PDM_MPI_Type_create_contiguous(3, PDM_MPI_INT, &mpi_init_location_type);
  PDM_MPI_Type_commit(&mpi_init_location_type);

  if(doct->solicitation_kind == PDM_TREE_SOLICITATION_BOXES_POINTS) {

    PDM_MPI_Type_create_contiguous(6, PDM_MPI_DOUBLE, &mpi_entity_type);
    PDM_MPI_Type_commit(&mpi_entity_type);

    wshared_entity_coord         = PDM_mpi_win_shared_create(6 * n_tot_recv_shared, sizeof(double     ), doct->comm_shared);
    wshared_entity_gnum          = PDM_mpi_win_shared_create(    n_tot_recv_shared, sizeof(PDM_g_num_t), doct->comm_shared);
    wshared_entity_init_location = PDM_mpi_win_shared_create(3 * n_tot_recv_shared, sizeof(int        ), doct->comm_shared);

    PDM_mpi_win_shared_lock_all (0, wshared_entity_coord);
    PDM_mpi_win_shared_lock_all (0, wshared_entity_gnum);
    PDM_mpi_win_shared_lock_all (0, wshared_entity_init_location);

    PDM_g_num_t *shared_entity_gnum          = PDM_mpi_win_shared_get(wshared_entity_gnum );
    double      *shared_entity_coord         = PDM_mpi_win_shared_get(wshared_entity_coord);
    int         *shared_entity_init_location = PDM_mpi_win_shared_get(wshared_entity_init_location);

    // Prepare send
    for(int i = 0; i < n_rank; ++i) {
      send_entity_n[i] = 0;
    }

    assert(doct->n_part == 1);
    PDM_g_num_t *box_g_num         = doct->entity_gnum         [0];
    double      *box_extents       = doct->entity_coords       [0];
    int         *box_init_location = doct->init_location_entity[0];

    for(int i = 0; i < n_rank; ++i) {
      for(int j = shared_local_nodes_idx[i]; j < shared_local_nodes_idx[i+1]; ++j) {
        int t_rank = shared_old_to_new_box_rank[j];
        for(int k = coarse_tree_box_to_box_idx[j]; k < coarse_tree_box_to_box_idx[j+1]; ++k) {
          int lnum = coarse_tree_box_to_box[k];
          int idx_lwrite = send_entity_idx[t_rank] + send_entity_n[t_rank]++;
          send_g_num[idx_lwrite] = box_g_num[lnum];
          for (int l = 0; l < 6; l++) {
            send_extents[6*idx_lwrite + l] = box_extents[6*lnum + l];
          }

          send_init_location[3*idx_lwrite  ] = box_init_location[3*lnum  ];
          send_init_location[3*idx_lwrite+1] = box_init_location[3*lnum+1];
          send_init_location[3*idx_lwrite+2] = box_init_location[3*lnum+2];

        }
      }
    }

    PDM_g_num_t *lrecv_gnum          = &shared_entity_gnum         [    shared_recv_idx[i_rank_in_shm]];
    double      *lrecv_extents       = &shared_entity_coord        [6 * shared_recv_idx[i_rank_in_shm]];
    int         *lrecv_init_location = &shared_entity_init_location[3 * shared_recv_idx[i_rank_in_shm]];
    PDM_MPI_Ialltoallv(send_g_num, send_entity_n, send_entity_idx, PDM__PDM_MPI_G_NUM,
                       lrecv_gnum, recv_entity_n, recv_entity_idx, PDM__PDM_MPI_G_NUM,
                       doct->comm,
                       &req_entity_gnum);

    PDM_MPI_Ialltoallv(send_extents , send_entity_n, send_entity_idx, mpi_entity_type,
                       lrecv_extents, recv_entity_n, recv_entity_idx, mpi_entity_type,
                       doct->comm,
                       &req_entity_coord);

    PDM_MPI_Ialltoallv(send_init_location , send_entity_n, send_entity_idx, mpi_init_location_type,
                       lrecv_init_location, recv_entity_n, recv_entity_idx, mpi_init_location_type,
                       doct->comm,
                       &req_entity_init_location);

  } else {
    abort();
  }


  /*
   * Step 3 : All pts are redistribute to equilibrate solicitation
   *   We can now rebuild a finer tree to finalize solicitation
   */
  if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_OCTREE ||
     doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_KDTREE) {

    assert(doct->local_tree == NULL);
    doct->local_tree = PDM_point_tree_seq_create(doct->local_tree_kind,
                                                 doct->local_depth_max,
                                                 doct->local_points_in_leaf_max,
                                                 doct->local_tolerance);

    PDM_point_tree_seq_point_cloud_set(doct->local_tree,
                                       equi_n_pts_tot,
                                       equi_pts_coords);

    PDM_point_tree_seq_build(doct->local_tree);

    if(1 == 1) {
      char filename[999];
      sprintf(filename, "out_local_tree_%i.vtk", i_rank);
      PDM_point_tree_seq_write_nodes(doct->local_tree, filename);
    }

  } else {
    abort();
  }
  // free(equi_pts_gnum);
  free(equi_pts_coords);
  free(equi_n_pts);

  /*
   * Setup shared
   */
  doct->shmem_tree = PDM_point_tree_make_shared(doct->local_tree,
                                                doct->comm_shared);

  PDM_point_tree_seq_free(doct->local_tree);
  doct->local_tree = NULL;

  /*
   * Make shared the equi pts_gnum
   */


  /*
   * Wait message and free all useless buffer
   */
  PDM_MPI_Wait(&req_entity_gnum);
  PDM_MPI_Wait(&req_entity_coord);
  PDM_MPI_Wait(&req_entity_init_location);
  PDM_MPI_Type_free(&mpi_entity_type);
  PDM_MPI_Type_free(&mpi_init_location_type);
  free(send_g_num);
  free(send_extents);
  free(send_init_location);
  free(send_entity_n);
  free(recv_entity_n);
  free(send_entity_idx);
  free(recv_entity_idx);

  PDM_g_num_t *shared_entity_gnum          = PDM_mpi_win_shared_get(wshared_entity_gnum );
  double      *shared_entity_coord         = PDM_mpi_win_shared_get(wshared_entity_coord);
  int         *shared_entity_init_location = PDM_mpi_win_shared_get(wshared_entity_init_location);

  PDM_g_num_t* distrib_search = PDM_compute_uniform_entity_distribution(doct->comm_shared, n_tot_recv_shared);
  int  dn_shared_box = distrib_search[i_rank_in_shm+1] - distrib_search[i_rank_in_shm];

  if(0 == 1) {
    char filename[999];
    sprintf(filename, "equi_boxes_for_solicitate_%i.vtk", i_rank);

    int beg    = distrib_search[i_rank_in_shm  ];

    double      *ptr_shared_box_extents = &shared_entity_coord[6*beg];
    PDM_g_num_t *ptr_shared_box_gnum    = &shared_entity_gnum [  beg];

    PDM_vtk_write_boxes(filename,
                        dn_shared_box,
                        ptr_shared_box_extents,
                        ptr_shared_box_gnum);
  }


  free(shared_recv_count);
  PDM_mpi_win_shared_unlock_all (wshared_local_nodes_n  );
  PDM_mpi_win_shared_unlock_all (wshared_local_nodes_idx);

  PDM_mpi_win_shared_free (wshared_local_nodes_n  );
  PDM_mpi_win_shared_free (wshared_local_nodes_idx);

  free(coarse_tree_box_to_box_idx);
  free(coarse_tree_box_to_box);
  PDM_mpi_win_shared_unlock_all (wshared_old_to_new_box_rank  );
  PDM_mpi_win_shared_free (wshared_old_to_new_box_rank  );

  /*
   * Finalize solicitation
   */
  int *distrib_search_by_rank_idx = malloc((n_rank_in_shm+1) * sizeof(int));

  for(int i = 0; i < n_rank_in_shm+1; ++i) {
    distrib_search_by_rank_idx[i] = 0;
  }

  for(int i = distrib_search[i_rank_in_shm]; i < distrib_search[i_rank_in_shm+1]; ++i) {
    int t_rank = PDM_binary_search_gap_int(i, shared_recv_idx, n_rank_in_shm+1);
    distrib_search_by_rank_idx[t_rank+1]++;
  }

  for(int i = 0; i < n_rank_in_shm; ++i) {
    distrib_search_by_rank_idx[i+1] += distrib_search_by_rank_idx[i];
  }
  free(distrib_search);
  free(shared_recv_idx);

  int n_part_out = n_rank_in_shm;
  int          *part_n_box         = malloc (sizeof(int          ) * n_part_out);
  int         **box_pts_idx        = malloc (sizeof(int         *) * n_part_out);
  int         **box_pts_l_num      = malloc (sizeof(int         *) * n_part_out);
  PDM_g_num_t **res_box_g_num      = malloc (sizeof(PDM_g_num_t *) * n_part_out);
  int         **res_box_strid      = malloc (sizeof(int         *) * n_part_out);
  double      **res_box_weight     = malloc (sizeof(double      *) * n_part_out);
  double      **res_box_pts_coords = malloc (sizeof(double      *) * n_part_out);
  PDM_g_num_t **res_box_pts_gnum   = malloc (sizeof(PDM_g_num_t *) * n_part_out);

  if(doct->solicitation_kind == PDM_TREE_SOLICITATION_BOXES_POINTS) {

    for(int i_shm = 0; i_shm < n_rank_in_shm; ++i_shm) {

      int beg    = distrib_search_by_rank_idx[i_shm  ];
      int n_lbox = distrib_search_by_rank_idx[i_shm+1] - beg;

      part_n_box[i_shm] = n_lbox;
      // PDM_g_num_t *lbox_gnum    = &shared_entity_gnum [  beg];
      double      *lbox_extents = &shared_entity_coord[6*beg];

      res_box_g_num[i_shm] = &shared_entity_gnum[beg];

      if(1 == 1) {
        char filename[999];
        sprintf(filename, "equi_boxes_for_solicitate_%i_%i.vtk", i_shm, i_rank);
        PDM_vtk_write_boxes(filename,
                            part_n_box[i_shm],
                            lbox_extents,
                            res_box_g_num[i_shm]);
      }

      PDM_point_tree_seq_points_inside_boxes_shared(doct->shmem_tree,
                                                    i_shm,
                                                    part_n_box[i_shm],
                                                    lbox_extents,
                                                    // lbox_gnum,
                                                    &(box_pts_idx[i_shm]),
                                                    &(box_pts_l_num[i_shm]));

      res_box_weight[i_shm] = malloc(n_lbox * sizeof(double));
      res_box_strid [i_shm] = malloc(n_lbox * sizeof(int   ));

      for(int i = 0; i < n_lbox; ++i ){
        res_box_strid [i_shm][i] = box_pts_idx[i_shm][i+1] - box_pts_idx[i_shm][i];
        res_box_weight[i_shm][i] = box_pts_idx[i_shm][i+1] - box_pts_idx[i_shm][i];
      }

      double *sorted_tree_coord = NULL;
      int    *new_to_old_pts    = NULL;
      PDM_point_tree_seq_shm_sorted_points_get   (doct->shmem_tree, i_shm, &sorted_tree_coord);
      PDM_point_tree_seq_shm_point_new_to_old_get(doct->shmem_tree, i_shm, &new_to_old_pts);

      /*
       * Extract point and gnum
       */
      int *_box_pts_idx   = box_pts_idx  [i_shm];
      int *_box_pts_l_num = box_pts_l_num[i_shm];
      res_box_pts_coords[i_shm] = malloc(3 * _box_pts_idx[n_lbox] * sizeof(double     ));
      res_box_pts_gnum  [i_shm] = malloc(    _box_pts_idx[n_lbox] * sizeof(PDM_g_num_t));

      PDM_g_num_t *shm_equi_pts_gnum = &equi_pts_gnum[shm_equi_pts_tot_idx[i_shm]];

      for(int i = 0;  i < _box_pts_idx[n_lbox]; ++i) {
        int l_num = _box_pts_l_num[i];
        res_box_pts_gnum  [i_shm][i] = shm_equi_pts_gnum[l_num];
        for (int k = 0; k < 3; k++) {
          res_box_pts_coords[i_shm][3*i + k] = sorted_tree_coord[3*l_num + k];
        }
      }
    }
  } else {
    abort();
  }

  free(distrib_search_by_rank_idx);
  free(shm_equi_pts_tot_idx);

  PDM_mpi_win_shared_unlock_all (wequi_pts_gnum);
  PDM_mpi_win_shared_free (wequi_pts_gnum);

  /*
   * Equilibrate unitary works
   */
  doct->ptb_unit_op_equi = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                     PDM_PART_TO_BLOCK_POST_MERGE,
                                                     1.,
                                    (PDM_g_num_t **) res_box_g_num,
                                                     res_box_weight,
                                                     part_n_box,
                                                     n_rank_in_shm,
                                                     doct->comm);

  for(int i_shm = 0; i_shm < n_rank_in_shm; ++i_shm) {
    free(res_box_weight[i_shm]);
    free(box_pts_idx   [i_shm]);
    free(box_pts_l_num [i_shm]);
  }
  free(box_pts_idx   );
  free(box_pts_l_num );

  /*
   * Exchange of gnum
   */
  int request_gnum = -1;
  int         *block_pts_in_box_n     = NULL;
  PDM_g_num_t *block_pts_in_box_g_num = NULL;
  PDM_part_to_block_iexch (doct->ptb_unit_op_equi,
                           PDM_MPI_COMM_KIND_COLLECTIVE,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_VAR_INTERLACED,
                           1,
                           res_box_strid,
                 (void **) res_box_pts_gnum,
                           &block_pts_in_box_n,
                 (void **) &block_pts_in_box_g_num,
                           &request_gnum);

  int request_coord = -1;
  int    *block_stride           = NULL;
  double *block_pts_in_box_coord = NULL;
  PDM_part_to_block_iexch (doct->ptb_unit_op_equi,
                           PDM_MPI_COMM_KIND_COLLECTIVE,
                           3 * sizeof(double),
                           PDM_STRIDE_VAR_INTERLACED,
                           1,
                           res_box_strid,
                 (void **) res_box_pts_coords,
                           &block_stride,
                 (void **) &block_pts_in_box_coord,
                           &request_coord);

  PDM_part_to_block_iexch_wait(doct->ptb_unit_op_equi, request_gnum);
  PDM_part_to_block_iexch_wait(doct->ptb_unit_op_equi, request_coord);

  free(block_stride);



  //-->>
  /* Remove doubles */
  int n_unit_op_equi_elt_block = PDM_part_to_block_n_elt_block_get (doct->ptb_unit_op_equi);
  if (1) {
    int max_n = 0;
    for (int i = 0; i < n_unit_op_equi_elt_block; i++) {
      max_n = PDM_MAX (max_n, block_pts_in_box_n[i]);
    }

    int *order = malloc (sizeof(int) * max_n);
    double *tmp_coord = malloc (sizeof(double) * max_n * 3);
    int idx1 = 0, idx2 = 0;
    for (int i = 0; i < n_unit_op_equi_elt_block; i++) {
      if (block_pts_in_box_n[i] == 0) continue;

      PDM_g_num_t *_g_num1 = block_pts_in_box_g_num + idx1;
      double      *_coord1 = block_pts_in_box_coord + idx1*3;
      PDM_g_num_t *_g_num2 = block_pts_in_box_g_num + idx2;
      double      *_coord2 = block_pts_in_box_coord + idx2*3;

      memcpy (tmp_coord, _coord1, sizeof(double) * block_pts_in_box_n[i] * 3);

      for (int j = 0; j < block_pts_in_box_n[i]; j++) {
        order[j] = j;
      }
      PDM_sort_long (_g_num1, order, block_pts_in_box_n[i]);

      _g_num2[0] = _g_num1[0];
      for (int k = 0; k < 3; k++) {
        _coord2[k] = tmp_coord[3*order[0] + k];
      }
      int tmp_n = 1;
      for (int j = 1; j < block_pts_in_box_n[i]; j++) {
        if (_g_num1[j] != _g_num2[tmp_n-1]) {
          _g_num2[tmp_n] = _g_num1[j];
          for (int k = 0; k < 3; k++) {
            _coord2[3*tmp_n + k] = tmp_coord[3*order[j] + k];
          }
          tmp_n++;
        }
      }

      idx1 += block_pts_in_box_n[i];
      idx2 += tmp_n;
      block_pts_in_box_n[i] = tmp_n;
    }
    free (order);
    free (tmp_coord);
  }
  //<<--

  doct->block_pts_in_box_n     = block_pts_in_box_n;
  doct->block_pts_in_box_g_num = block_pts_in_box_g_num;
  doct->block_pts_in_box_coord = block_pts_in_box_coord;

  for(int i_shm = 0; i_shm < n_rank_in_shm; ++i_shm) {
    free(res_box_pts_coords[i_shm]);
    free(res_box_pts_gnum  [i_shm]);
    free(res_box_strid     [i_shm]);
  }
  free(part_n_box    );

  free(res_box_g_num );
  free(res_box_strid );
  free(res_box_weight);
  free(res_box_pts_coords);
  free(res_box_pts_gnum  );

  PDM_mpi_win_shared_unlock_all (wshared_entity_coord  );
  PDM_mpi_win_shared_unlock_all (wshared_entity_gnum);
  PDM_mpi_win_shared_unlock_all (wshared_entity_init_location);

  PDM_mpi_win_shared_free (wshared_entity_coord  );
  PDM_mpi_win_shared_free (wshared_entity_gnum);
  PDM_mpi_win_shared_free (wshared_entity_init_location);


  PDM_MPI_Comm_free(&doct->comm_dist_graph);
  PDM_MPI_Comm_free(&doct->comm_shared);
  free(doct->neighbor_in);
}

void
PDM_doctree_point_set
(
 PDM_doctree_t     *doct,
 const int          i_part_cloud,
 const int          n_points,
 const int         *pts_init_location,
 const PDM_g_num_t *pts_g_num,
 const double      *pts_coords
)
{
  assert(i_part_cloud < doct->n_part_cloud);

  doct->n_point_cloud    [i_part_cloud] = n_points;
  doct->pts_g_num        [i_part_cloud] = (PDM_g_num_t *) pts_g_num;
  doct->pts_coords       [i_part_cloud] = (double      *) pts_coords;
  doct->pts_init_location[i_part_cloud] = (int         *) pts_init_location;
}

void
PDM_doctree_solicitation_set
(
 PDM_doctree_t             *doct,
 PDM_tree_solicitation_t    solicitation_kind,
 int                        n_part,
 int                       *n_entity,
 int                      **init_location_entity,
 PDM_g_num_t              **entity_gnum,
 double                   **entity_coords
)
{
  doct->solicitation_kind    = solicitation_kind;
  doct->n_part               = n_part;
  doct->n_entity             = n_entity;
  doct->init_location_entity = init_location_entity;
  doct->entity_gnum          = entity_gnum;
  doct->entity_coords        = entity_coords;
}

void
PDM_doctree_results_in_orig_frame_get
(
 PDM_doctree_t       *doct,
 int                  n_boxes,
 PDM_g_num_t         *box_g_num,
 int                **box_pts_idx,
 PDM_g_num_t        **box_pts,
 double             **pts_coord
)
{

  /*
   *  Block to part
   */
  PDM_g_num_t* blk_gnum = PDM_part_to_block_block_gnum_get(doct->ptb_unit_op_equi);
  int n_unit_op_equi_elt_block = PDM_part_to_block_n_elt_block_get (doct->ptb_unit_op_equi);

  PDM_block_to_part_t *btp = PDM_block_to_part_create_from_sparse_block(blk_gnum,
                                                                        n_unit_op_equi_elt_block,
                                                (const PDM_g_num_t **) &box_g_num,
                                                                       &n_boxes,
                                                                       1,
                                                                       doct->comm);

  int         **_tmp_pts_in_box_n     = NULL;
  PDM_g_num_t **_tmp_pts_in_box_g_num = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         doct->block_pts_in_box_n,
                         doct->block_pts_in_box_g_num,
                         &_tmp_pts_in_box_n,
              (void ***) &_tmp_pts_in_box_g_num);

  int *pts_in_box_n = _tmp_pts_in_box_n[0];
  free(_tmp_pts_in_box_n);

  *box_pts_idx = PDM_array_new_idx_from_sizes_int(pts_in_box_n, n_boxes);
  free(pts_in_box_n);

  PDM_g_num_t *pts_in_box_g_num = _tmp_pts_in_box_g_num[0];
  free(_tmp_pts_in_box_g_num);

  double **tmp_pts_in_box_coord = NULL;

  PDM_block_to_part_exch(btp,
                          3 * sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          doct->block_pts_in_box_n,
                 (void *) doct->block_pts_in_box_coord,
                          &_tmp_pts_in_box_n,
               (void ***) &tmp_pts_in_box_coord);
  free(_tmp_pts_in_box_n[0]);
  free(_tmp_pts_in_box_n);

  PDM_block_to_part_free(btp);

  *box_pts   = pts_in_box_g_num;
  *pts_coord = tmp_pts_in_box_coord[0];
  free(tmp_pts_in_box_coord);
}


void
PDM_doctree_free
(
  PDM_doctree_t   *doct
)
{
  free(doct->n_point_cloud    );
  free(doct->pts_g_num        );
  free(doct->pts_coords       );
  free(doct->pts_init_location);

  free(doct->block_pts_in_box_g_num);
  free(doct->block_pts_in_box_coord);
  free(doct->block_pts_in_box_n    );

  PDM_part_to_block_free(doct->ptb_unit_op_equi);

  if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_OCTREE ||
     doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_KDTREE) {
    if(doct->local_tree != NULL) {
      PDM_point_tree_seq_free(doct->local_tree);
    }
    if(doct->shmem_tree != NULL) {
      PDM_point_tree_seq_shm_free(doct->shmem_tree);
    }
  } else {
    abort();
  }

  free(doct);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
