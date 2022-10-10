#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>

#include "pdm.h"
#include "pdm_priv.h"

#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_point_cloud_gen.h"
#include "pdm_octree.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_doctree.h"
#include "pdm_box_gen.h"
#include "pdm_sphere_surf_gen.h"
#include "pdm_part_to_block.h"
#include "pdm_point_tree_seq.h"
#include "pdm_box.h"
#include "pdm_box_tree.h"
#include "pdm_block_to_part.h"
#include "pdm_box_priv.h"
#include "pdm_binary_search.h"
#include "pdm_point_tree_seq_priv.h"
#include "pdm_array.h"
#include "pdm_unique.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/





    // /*
    //  * Unique id of boxes
    //  */
    // int *shared_box_tag = malloc(g_extract_boxes_idx[n_rank] * sizeof(int)); // Suralloc
    // int *shared_box_id  = malloc(g_extract_boxes_idx[n_rank] * sizeof(int)); // Suralloc
    // for(int i = 0; i < g_extract_boxes_idx[n_rank]; ++i) {
    //   shared_box_tag[i] = -1;
    // }

    // int n_shared_box_in = 0;

    // int idx_read = 0;
    // for(int i = 0; i < dn_equi_box; ++i) {
    //   for(int j = 0; j < blk_box_to_coarse_box_pts_n[i]; ++j) {
    //     if(shared_box_tag[blk_box_to_coarse_box_pts[idx_read]-1] == -1) {
    //       shared_box_tag[blk_box_to_coarse_box_pts[idx_read]-1] = blk_box_to_coarse_box_pts[idx_read]-1;
    //       shared_box_id[n_shared_box_in++] = blk_box_to_coarse_box_pts[idx_read]-1;
    //     }
    //     idx_read++;
    //   }
    // }



    // PDM_log_trace_array_int(shared_box_id, n_shared_box_in, "shared_box_id :");

    // // Ok now we know with wich proc we communicate
    // int *neighbor_in  = malloc(n_rank * sizeof(int)); // A adpater avec le nombre de range voisin courant
    // // int *neighbor_tag = malloc(n_rank * sizeof(int)); // A adpater avec le nombre de range voisin courant
    // for(int i = 0; i < n_rank; ++i) {
    //   neighbor_tag[i] = -1;
    // }

    // int *shared_id_rank = malloc(n_shared_box_in * sizeof(int));

    // for(int i = 0; i < n_shared_box_in; ++i) {
    //   int t_rank = PDM_binary_search_gap_long((PDM_g_num_t) shared_box_id[i], g_extract_boxes_idx, n_rank+1);
    //   if(neighbor_tag[t_rank] == -1) {
    //     neighbor_tag[t_rank] = n_neighbor_in;
    //     neighbor_in [n_neighbor_in++] = t_rank;
    //   }
    //   shared_id_rank[i] = neighbor_tag[t_rank]; // Donc le numero de neighbor
    // }


    // PDM_log_trace_array_int(neighbor_in, n_neighbor_in, "neighbor_in :");
    // PDM_log_trace_array_int(shared_id_rank, n_shared_box_in, "shared_id_rank :");

    // PDM_MPI_Comm comm_dist_graph;
    // PDM_MPI_setup_dist_graph_from_neighbor_in(comm, n_neighbor_in, neighbor_in, &comm_dist_graph);

    // int n_sources      = 0;
    // int n_destinations = 0;
    // int is_weight      = 0;
    // PDM_MPI_Dist_graph_neighbors_count(comm_dist_graph, &n_sources, &n_destinations, &is_weight);

    // int *sources      = malloc(n_sources      * sizeof(int));
    // int *destinations = malloc(n_destinations * sizeof(int));

    // PDM_MPI_Dist_graph_neighbors(comm_dist_graph, n_sources, sources, n_destinations, destinations);

    // PDM_log_trace_array_int(sources     , n_sources     , "sources ::");
    // PDM_log_trace_array_int(destinations, n_destinations, "destinations ::");

    // // We need the reverse ones ...
    // PDM_MPI_Comm comm_dist_graph_reverse;
    // PDM_MPI_Dist_graph_create_adjacent(comm,
    //                                    n_destinations,
    //                                    destinations,
    //                                    n_sources,
    //                                    sources,
    //                                    0,
    //                                    &comm_dist_graph_reverse);

    // /*
    //  * On connait les rang qui ont la data qu'on souhaite
    //  */
    // int *send_n = malloc(   n_sources   * sizeof(int));
    // int *recv_n = malloc(n_destinations * sizeof(int));

    // for(int i = 0; i < n_sources; ++i) {
    //   send_n[i] = 0;
    // }

    // for(int i = 0; i < n_shared_box_in; ++i) {
    //   send_n[shared_id_rank[i]]++;
    // }

    // PDM_MPI_Neighbor_alltoall(send_n, 1, PDM_MPI_INT,
    //                           recv_n, 1, PDM_MPI_INT, comm_dist_graph_reverse);

    // PDM_log_trace_array_int(send_n, n_sources     , "send_n ::");
    // PDM_log_trace_array_int(recv_n, n_destinations, "recv_n ::");

    // /*
    //  * On connait le nombre de request
    //  */
    // int* send_idx = malloc((n_sources     +1) * sizeof(int));
    // int* recv_idx = malloc((n_destinations+1) * sizeof(int));
    // send_idx[0] = 0;
    // recv_idx[0] = 0;
    // for(int i = 0; i < n_sources; ++i) {
    //   send_idx[i+1] = send_idx[i] + send_n[i];
    //   send_n[i] = 0;
    // }
    // for(int i = 0; i < n_destinations; ++i) {
    //   recv_idx[i+1] = recv_idx[i] + recv_n[i];
    // }

    // // for(int i = 0; i < n_shared_box_in; ++i) {
    // //   int idx_write = send_idx[shared_id_rank[i]] + send_n[shared_id_rank[i]]++;
    // //   send_node_id[idx_write] = shared_box_id[i];
    // // }


    // free(send_n);
    // free(recv_n);

    // free(sources);
    // free(destinations);
    // free(neighbor_in);
    // free(neighbor_tag);

    // free(g_extract_boxes_idx);
    // free(n_g_coarse_pts_box);
    // free(shared_box_tag);
    // free(shared_box_id);

    // PDM_MPI_Comm_free(&comm_dist_graph);
/**
 *
 * \brief  Usage
 *
 */

static void
_usage
(
int exit_code
)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of points (default : 10).\n\n"
     "  -radius <level>  Radius of domain (default : 10).\n\n"
     "  -local           Number of points is local (default : global).\n\n"
     "  -rand            Random definition of point coordinates (default : false).\n\n"
     "  -h               This message.\n\n");

  exit (exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]    argc   Number of arguments
 * \param [in]    argv   Arguments
 * \param [inout] nPts   Number of points
 * \param [inout] ls     Low scalability
 * \param [inout] length Length of domains
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 PDM_g_num_t   *nPts,
 double        *radius,
 int           *local,
 int           *rand
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _nPts = atol(argv[i]);
        *nPts = (PDM_g_num_t) _nPts;
      }
    }

    else if (strcmp(argv[i], "-radius") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *radius = atof(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-local") == 0) {
      *local = 1;
    }

    else if (strcmp(argv[i], "-rand") == 0) {
      *rand = 1;
    }

    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}

static
void
_adaptative_tree2
(
  int           n_pts,
  double       *pts_coord,
  PDM_g_num_t  *pts_gnum,
  int           n_box,
  double       *box_extents,
  PDM_g_num_t  *box_gnum,
  PDM_MPI_Comm  comm
)
{
  int i_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);
  /*
   * Hilbert of points AND boxes
   */
  int    *weight_pts = malloc(    n_pts * sizeof(int   ));
  int    *weight_box = malloc(    n_box * sizeof(int   ));
  double *box_center = malloc(3 * n_box * sizeof(double));
  for(int i = 0; i < n_pts; ++i) {
    weight_pts[i] = 1;
  }
  for(int i = 0; i < n_box; ++i) {
    weight_box[i] = 1;
    box_center[3*i  ] = 0.5 * (box_extents[6*i  ] + box_extents[6*i+3]);
    box_center[3*i+1] = 0.5 * (box_extents[6*i+1] + box_extents[6*i+4]);
    box_center[3*i+2] = 0.5 * (box_extents[6*i+2] + box_extents[6*i+5]);
  }

  PDM_part_to_block_t* ptb_pts = PDM_part_to_block_geom_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                               PDM_PART_TO_BLOCK_POST_CLEANUP, // A voir avec merge mais attention au init_location
                                                               1.,
                                                               PDM_PART_GEOM_HILBERT,
                                                               &pts_coord,
                                                               &pts_gnum,
                                                               &weight_pts,
                                                               &n_pts,
                                                               1,
                                                               comm);
  free(weight_pts);

  PDM_part_to_block_t* ptb_box = PDM_part_to_block_geom_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                               PDM_PART_TO_BLOCK_POST_CLEANUP, // A voir avec merge mais attention au init_location
                                                               1.,
                                                               PDM_PART_GEOM_HILBERT,
                                                               &box_center,
                                                               &box_gnum,
                                                               &weight_box,
                                                               &n_box,
                                                               1,
                                                               comm);
  free(weight_box);
  free(box_center);


  double *blk_pts_coord = NULL;
  PDM_part_to_block_exch(ptb_pts,
                         3 * sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         1,
                         NULL,
               (void **) &pts_coord,
                         NULL,
               (void **) &blk_pts_coord);


  double *blk_box_extents = NULL;
  PDM_part_to_block_exch(ptb_box,
                         6 * sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         1,
                         NULL,
               (void **) &box_extents,
                         NULL,
               (void **) &blk_box_extents);

  int dn_pts = PDM_part_to_block_n_elt_block_get(ptb_pts);
  int dn_box = PDM_part_to_block_n_elt_block_get(ptb_box);

  // Plus besoin du gnum pour l'instant ....
  PDM_part_to_block_free(ptb_pts);
  PDM_part_to_block_free(ptb_box);

  PDM_MPI_Comm comm_alone;
  PDM_MPI_Comm_split(comm, i_rank, 0, &(comm_alone));

  PDM_MPI_Datatype mpi_extent_type;
  PDM_MPI_Type_create_contiguous(6, PDM_MPI_DOUBLE, &mpi_extent_type);
  PDM_MPI_Type_commit(&mpi_extent_type);

  // Build only octree
  PDM_point_tree_seq_t* coarse_tree_pts = PDM_point_tree_seq_create(PDM_DOCTREE_LOCAL_TREE_OCTREE,
                                                                    4, // depth_max
                                                                    1,
                                                                    1e-8);
  PDM_point_tree_seq_point_cloud_set(coarse_tree_pts, dn_pts, blk_pts_coord);
  PDM_point_tree_seq_build(coarse_tree_pts);

  if(1 == 1) {
    char filename[999];
    sprintf(filename, "out_coarse_tree_%i.vtk", i_rank);
    PDM_point_tree_seq_write_nodes(coarse_tree_pts, filename);
  }

  if(1 == 1) {
    char filename[999];
    sprintf(filename, "blk_pts_coord_%i.vtk", i_rank);
    PDM_vtk_write_point_cloud(filename,
                              dn_pts,
                              blk_pts_coord,
                              NULL,
                              NULL);
  }


  int dn_node = coarse_tree_pts->n_nodes;
  PDM_g_num_t* distrib_pts_box = PDM_compute_entity_distribution(comm, dn_node);

  // PDM_log_trace_array_long(distrib_pts_box, n_rank+1, "distrib_pts_box :");

  /*
   * Initialisation : Extract extents on all local_tree
   */
  int     n_coarse_pts_box       = 0;
  int    *coarse_pts_box_id      = NULL;
  double *coarse_pts_box_extents = NULL;
  int    *coarse_pts_box_n_pts   = NULL; // Number of point in boxes

  PDM_point_tree_seq_extract_nodes(coarse_tree_pts,
                                   0,
                                   1, // Depth
                                   &n_coarse_pts_box,
                                   &coarse_pts_box_id,
                                   &coarse_pts_box_extents,
                                   &coarse_pts_box_n_pts);

  int *n_g_coarse_pts_box = malloc(n_rank * sizeof(int));
  PDM_MPI_Allgather (&n_coarse_pts_box , 1, PDM_MPI_INT,
                     n_g_coarse_pts_box, 1, PDM_MPI_INT, comm);

  int *g_extract_boxes_idx = (int *) malloc (sizeof(int) * (n_rank+1));
  g_extract_boxes_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    g_extract_boxes_idx[i+1] = g_extract_boxes_idx[i] + n_g_coarse_pts_box[i];
  }
  double *g_coarse_pts_box_extents = malloc(6 * g_extract_boxes_idx[n_rank] * sizeof(double));

  PDM_MPI_Allgatherv(coarse_pts_box_extents  , n_coarse_pts_box, mpi_extent_type,
                     g_coarse_pts_box_extents, n_g_coarse_pts_box,
                     g_extract_boxes_idx,
                     mpi_extent_type, comm);

  int *g_coarse_pts_box_id = malloc( g_extract_boxes_idx[n_rank] * sizeof(int));
  PDM_MPI_Allgatherv(coarse_pts_box_id  , n_coarse_pts_box, PDM_MPI_INT,
                     g_coarse_pts_box_id, n_g_coarse_pts_box,
                     g_extract_boxes_idx,
                     PDM_MPI_INT, comm);

  for(int i = 0; i < n_rank; ++i) {
    for(int j = g_extract_boxes_idx[i]; j < g_extract_boxes_idx[i+1]; ++j) {
      g_coarse_pts_box_id[j] += distrib_pts_box[i];
    }
  }

  free(coarse_pts_box_id     );
  free(coarse_pts_box_extents);
  free(coarse_pts_box_n_pts  );
  free(n_g_coarse_pts_box);

  int n_neighbor_current = n_rank;

  int n_iter = 3;
  for(int i_iter = 0; i_iter < n_iter; ++i_iter) {

    /*
     * Build tree
     */
    PDM_g_num_t *coarse_pts_box_gnum         = malloc(    g_extract_boxes_idx[n_neighbor_current] * sizeof(PDM_g_num_t));
    int         *init_location_coase_pts_box = malloc(3 * g_extract_boxes_idx[n_neighbor_current] * sizeof(int        ));
    for(int i = 0; i < g_extract_boxes_idx[n_neighbor_current]; ++i) {
      coarse_pts_box_gnum[i] = g_coarse_pts_box_id[i]+1;
      init_location_coase_pts_box[3*i  ] = 0;
      init_location_coase_pts_box[3*i+1] = 0;
      init_location_coase_pts_box[3*i+2] = i;
    }
    free(g_coarse_pts_box_id);

    // PDM_log_trace_connectivity_long(g_extract_boxes_idx, coarse_pts_box_gnum, n_neighbor_current, "coarse_pts_box_gnum ::");

    PDM_box_set_t  *coarse_pts_box_set = PDM_box_set_create(3,
                                                            1,  // No normalization to preserve initial extents
                                                            0,  // No projection to preserve initial extents
                                                            g_extract_boxes_idx[n_neighbor_current],
                                                            coarse_pts_box_gnum,
                                                            g_coarse_pts_box_extents,
                                                            1,
                                                            &g_extract_boxes_idx[n_neighbor_current],
                                                            init_location_coase_pts_box,
                                                            comm_alone);

    int   max_boxes_leaf_coarse = 1;   // Max number of boxes in a leaf for coarse coarse BBTree
    int   max_tree_depth_coarse = 31;  // Max tree depth for coarse coarse BBTree
    float max_box_ratio_coarse  = 5;   // Max ratio for local BBTree (nConnectedBoxe < ratio * nBoxes)
    PDM_box_tree_t* coarse_pts_bt_shared = PDM_box_tree_create (max_tree_depth_coarse,
                                                                max_boxes_leaf_coarse,
                                                                max_box_ratio_coarse);

    PDM_box_tree_set_boxes (coarse_pts_bt_shared,
                            coarse_pts_box_set,
                            PDM_BOX_TREE_ASYNC_LEVEL);
    const PDM_g_num_t *bt_box_pts_gnum   = PDM_box_set_get_g_num (coarse_pts_box_set);
    // const int         *bt_box_pts_origin = PDM_box_set_origin_get(coarse_pts_box_set);

    free(coarse_pts_box_gnum);
    free(init_location_coase_pts_box);

    if(1 == 1) {
      char filename[999];
      sprintf(filename, "coarse_pts_box_set_%i.vtk", i_rank);
      PDM_box_tree_write_vtk(filename, coarse_pts_bt_shared, -1, 0);

      sprintf(filename, "coarse_pts_box_%i.vtk", i_rank);
      PDM_vtk_write_boxes(filename,
                          g_extract_boxes_idx[n_neighbor_current],
                          g_coarse_pts_box_extents,
                          NULL);
    }


    // Intersect naîvily
    int *box_to_coarse_box_pts_idx = NULL;
    int *box_to_coarse_box_pts     = NULL;
    PDM_box_tree_intersect_boxes_boxes(coarse_pts_bt_shared,
                                       -1,
                                       dn_box,
                                       blk_box_extents,
                                       &box_to_coarse_box_pts_idx,
                                       &box_to_coarse_box_pts);

    // PDM_log_trace_connectivity_int(box_to_coarse_box_pts_idx,
    //                                box_to_coarse_box_pts,
    //                                dn_box, "box_to_coarse_box_pts ::");


    free(g_coarse_pts_box_extents);


    // On fait block_to_part sur les box en refaisant une distrib implicit !!
    // Attention il faut envoyer le gnum original également
    // On compresse l'info des boites de pts intersecter
    //    --> On cherche
    int n_extract = 0;
    int n_extract_box_to_coarse_box_pts_tot = box_to_coarse_box_pts_idx[dn_box];
    for(int i = 0; i < dn_box; ++i) {
      if(box_to_coarse_box_pts_idx[i+1] - box_to_coarse_box_pts_idx[i] == 0) {
        continue;
      }
      n_extract++;
    }

    PDM_g_num_t *extract_box_gnum                = malloc(    n_extract    * sizeof(PDM_g_num_t));
    int          *weight                          = malloc(    n_extract    * sizeof(double     ));
    double      *extract_box_extents             = malloc(6 * n_extract    * sizeof(double     ));
    double      *extract_box_center              = malloc(3 * n_extract    * sizeof(double     ));
    int         *extract_box_to_coarse_box_pts_n = malloc(6 * n_extract    * sizeof(int        ));
    int         *extract_box_to_coarse_box_pts   = malloc(n_extract_box_to_coarse_box_pts_tot * sizeof(int        ));

    int idx_write = 0;
    n_extract = 0;
    for(int i = 0; i < dn_box; ++i) {

      if(box_to_coarse_box_pts_idx[i+1] - box_to_coarse_box_pts_idx[i] == 0) {
        continue;
      }

      extract_box_gnum[n_extract] = n_extract+1;
      weight          [n_extract] = box_to_coarse_box_pts_idx[i+1] - box_to_coarse_box_pts_idx[i];
      for(int k = 0; k < 6; ++k) {
        extract_box_extents[6*n_extract+k] = blk_box_extents[6*i+k];
      }
      extract_box_center[3*n_extract  ] = 0.5 * ( extract_box_extents[6*n_extract  ] + extract_box_extents[6*n_extract+3]);
      extract_box_center[3*n_extract+1] = 0.5 * ( extract_box_extents[6*n_extract+1] + extract_box_extents[6*n_extract+4]);
      extract_box_center[3*n_extract+2] = 0.5 * ( extract_box_extents[6*n_extract+2] + extract_box_extents[6*n_extract+5]);


      extract_box_to_coarse_box_pts_n[n_extract] = box_to_coarse_box_pts_idx[i+1] - box_to_coarse_box_pts_idx[i];
      for(int j = box_to_coarse_box_pts_idx[i]; j < box_to_coarse_box_pts_idx[i+1]; ++j) {
        extract_box_to_coarse_box_pts[idx_write++] = (int) bt_box_pts_gnum[box_to_coarse_box_pts[j]];
      }

      n_extract++;
    }

    /*
     * Chaque boite peut connecter a un ou plusieurs rangs
     *  Au final vu qu'on fait du brute force sur les box, il faut garder le lien nouveau_proc -> boites anciennement connecté
     *  En faisant manuelement un dest_proc on peut envoyer cette info
     *  Attention il faut transmetre le numero de la boites de l'arbre, pas juste le numero de rank
     *  Mais il faut trié avant !!!
     */
    PDM_g_num_t _n_extract = n_extract;
    PDM_g_num_t offset = 0;
    PDM_MPI_Exscan(&_n_extract, &offset, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);

    for(int i = 0; i < n_extract; ++i) {
      extract_box_gnum[i] = extract_box_gnum[i] + offset;
    }
    PDM_log_trace_array_long(extract_box_gnum, n_extract, "extract_box_gnum ::");

    free(box_to_coarse_box_pts_idx);
    free(box_to_coarse_box_pts);

    // PDM_part_to_block_t* ptb_equi_box = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
    //                                                              PDM_PART_TO_BLOCK_POST_CLEANUP,
    //                                                              1.,
    //                                                              &extract_box_gnum,
    //                                                              &weight,
    //                                                              &n_extract,
    //                                                              1,
    //                                                              comm);
    PDM_part_to_block_t* ptb_equi_box = PDM_part_to_block_geom_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                      PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                                      1.,
                                                                      PDM_PART_GEOM_HILBERT,
                                                                      &extract_box_center,
                                                                      &extract_box_gnum,
                                                                      &weight,
                                                                      &n_extract,
                                                                      1,
                                                                      comm);
    free(extract_box_center);
    free(weight);
    free(extract_box_gnum);

    int *blk_box_to_coarse_box_pts_n = NULL;
    int *blk_box_to_coarse_box_pts   = NULL;
    int blk_size = PDM_part_to_block_exch(ptb_equi_box,
                                          sizeof(int),
                                          PDM_STRIDE_VAR_INTERLACED,
                                          1,
                                (int  **) &extract_box_to_coarse_box_pts_n,
                                (void **) &extract_box_to_coarse_box_pts,
                                          &blk_box_to_coarse_box_pts_n,
                                (void **) &blk_box_to_coarse_box_pts);
    free(extract_box_to_coarse_box_pts_n);
    free(extract_box_to_coarse_box_pts);

    // A faire -> Echange des gnum

    double *tmp_blk_box_extents = NULL;
    int     request_box_extents = -1;
    PDM_part_to_block_iexch(ptb_equi_box,
                            PDM_MPI_COMM_KIND_COLLECTIVE,
                            6 * sizeof(double),
                            PDM_STRIDE_CST_INTERLACED,
                            1,
                            NULL,
                  (void **) &extract_box_extents,
                            NULL,
                  (void **) &tmp_blk_box_extents,
                            &request_box_extents);

    /*
     * Make sort
     */
    int n_unique_box_pts = PDM_inplace_unique(blk_box_to_coarse_box_pts, 0, blk_size-1);
    blk_box_to_coarse_box_pts = realloc(blk_box_to_coarse_box_pts, n_unique_box_pts * sizeof(int));

    if(0 == 1) {
      PDM_log_trace_array_int(blk_box_to_coarse_box_pts, n_unique_box_pts, "blk_box_to_coarse_box_pts :");
    }

    PDM_part_to_block_iexch_wait(ptb_equi_box, request_box_extents);
    free(extract_box_extents);

    /*
     * Setup boxes
     */
    int          dn_equi_box      = PDM_part_to_block_n_elt_block_get(ptb_equi_box);

    free(blk_box_extents);
    dn_box = dn_equi_box;
    blk_box_extents = tmp_blk_box_extents;
    PDM_part_to_block_free(ptb_equi_box);

    /*
     * Creation du nouveau graphe de comm
     */
    int n_neighbor_in = 0;
    int *neighbor_tag = malloc(n_rank * sizeof(int));
    int *neighbor_in  = malloc(n_rank * sizeof(int));
    for(int i = 0; i < n_rank; ++i) {
      neighbor_tag[i] = -1;
    }

    // C'est sort donc on peux faire une range
    int *send_request_pts_box_n   = PDM_array_zeros_int(n_rank+1);
    int *send_request_pts_box_idx = PDM_array_zeros_int(n_rank+1);
    for(int i = 0; i < n_unique_box_pts; ++i) {
      // Meme pas besoin du binary search gap long car deja trié
      int t_rank = PDM_binary_search_gap_long(blk_box_to_coarse_box_pts[i]-1, distrib_pts_box, n_rank+1);
      if(neighbor_tag[t_rank] == -1) {
        neighbor_in[n_neighbor_in++] = t_rank;
        neighbor_tag[t_rank] = n_neighbor_in;
      }
      send_request_pts_box_n[n_neighbor_in-1]++; // Car tout est trié
    }

    send_request_pts_box_idx[0] = 0;
    for(int i = 0; i < n_neighbor_in; ++i) {
      send_request_pts_box_idx[i+1] = send_request_pts_box_idx[i] + send_request_pts_box_n[i];
    }

    if(1 == 1) {
      PDM_log_trace_array_int(neighbor_in, n_neighbor_in  , "neighbor_in :");
      PDM_log_trace_array_int(send_request_pts_box_idx, n_neighbor_in+1, "send_request_pts_box_idx :");
      PDM_log_trace_array_int(send_request_pts_box_n, n_neighbor_in  , "send_request_pts_box_n :");
    }

    PDM_MPI_Comm comm_dist_graph;
    PDM_MPI_setup_dist_graph_from_neighbor_in(comm, n_neighbor_in, neighbor_in, &comm_dist_graph);

    int n_sources      = 0;
    int n_destinations = 0;
    int is_weight      = 0;
    PDM_MPI_Dist_graph_neighbors_count(comm_dist_graph, &n_sources, &n_destinations, &is_weight);

    int *sources      = malloc(n_sources      * sizeof(int));
    int *destinations = malloc(n_destinations * sizeof(int));

    PDM_MPI_Dist_graph_neighbors(comm_dist_graph, n_sources, sources, n_destinations, destinations);

    if(1 == 1) {
      PDM_log_trace_array_int(sources     , n_sources     , "sources ::");
      PDM_log_trace_array_int(destinations, n_destinations, "destinations ::");
    }

    PDM_MPI_Comm comm_dist_graph_reverse;
    PDM_MPI_Dist_graph_create_adjacent(comm,
                                       n_destinations,
                                       destinations,
                                       n_sources,
                                       sources,
                                       0,
                                       &comm_dist_graph_reverse);

    int *recv_request_pts_box_n = malloc(n_destinations * sizeof(int));

    PDM_MPI_Neighbor_alltoall(send_request_pts_box_n, 1, PDM_MPI_INT,
                              recv_request_pts_box_n, 1, PDM_MPI_INT, comm_dist_graph_reverse);


    int *recv_request_pts_box_idx = malloc((n_destinations+1) * sizeof(int));
    recv_request_pts_box_idx[0] = 0;
    for(int i = 0; i < n_destinations; ++i) {
      recv_request_pts_box_idx[i+1] = recv_request_pts_box_idx[i] + recv_request_pts_box_n[i];
    }


    if(1 == 1) {
      PDM_log_trace_array_int(recv_request_pts_box_idx, n_destinations+1, "recv_request_pts_box_idx :");
      PDM_log_trace_array_int(recv_request_pts_box_n, n_destinations  , "recv_request_pts_box_n :");
    }

    int *recv_request_pts_box = malloc(recv_request_pts_box_idx[n_destinations] * sizeof(int));

    PDM_MPI_Neighbor_alltoallv(blk_box_to_coarse_box_pts, send_request_pts_box_n, send_request_pts_box_idx, PDM_MPI_INT,
                               recv_request_pts_box     , recv_request_pts_box_n, recv_request_pts_box_idx, PDM_MPI_INT, comm_dist_graph_reverse);

    // Unshift
    int n_node_to_extract = recv_request_pts_box_idx[n_destinations];
    for(int i = 0; i < n_node_to_extract; ++i) {
      recv_request_pts_box[i] -= (distrib_pts_box[i_rank] + 1);
    }

    if(1 == 1) {
      // PDM_log_trace_connectivity_long(recv_request_pts_box_idx, recv_request_pts_box, n_destinations, "recv_request_pts_box ::");
    }

    /*
     * Extract extents on all local_tree
     */
    int     n_extract_child   = 0;
    int    *extract_child_id  = NULL;
    int    *extract_is_leaf   = NULL;
    double *extract_extents   = NULL;
    int    *node_to_child_idx = NULL;

    PDM_point_tree_seq_extract_extents_by_child_ids(coarse_tree_pts,
                                                    n_node_to_extract,
                                                    recv_request_pts_box,
                                                    &n_extract_child,
                                                    &node_to_child_idx,
                                                    &extract_child_id,
                                                    &extract_is_leaf,
                                                    &extract_extents);

    if(0 == 1) {
      PDM_log_trace_array_int(extract_child_id, n_extract_child, "extract_child_id ::");
      PDM_log_trace_connectivity_int(node_to_child_idx, extract_child_id, n_node_to_extract, "node_to_child ::");
    }

    for(int i = 0; i < n_extract_child; ++i) {
      extract_child_id[i] += distrib_pts_box[i_rank];
    }

    // Accumulate to prepare send
    int *send_child_extract_n = malloc(n_destinations * sizeof(int));
    for(int i = 0; i < n_destinations; ++i) {
      send_child_extract_n[i] = 0;
      for(int j = recv_request_pts_box_idx[i]; j < recv_request_pts_box_idx[i+1]; ++j) {
        send_child_extract_n[i] += node_to_child_idx[j+1] - node_to_child_idx[j];
      }
    }

    int *send_child_extract_idx = malloc((n_destinations+1) * sizeof(int));
    send_child_extract_idx[0] = 0;
    for(int i = 0; i < n_destinations; ++i) {
      send_child_extract_idx[i+1] = send_child_extract_idx[i] + send_child_extract_n[i];
    }

    int *recv_child_extract_n = malloc(n_sources * sizeof(int));
    PDM_MPI_Neighbor_alltoall (send_child_extract_n, 1, PDM_MPI_INT,
                               recv_child_extract_n, 1, PDM_MPI_INT, comm_dist_graph);

    int *recv_child_extract_idx = malloc((n_sources+1) * sizeof(int));
    recv_child_extract_idx[0] = 0;
    for(int i = 0; i < n_sources; ++i) {
      recv_child_extract_idx[i+1] = recv_child_extract_idx[i] + recv_child_extract_n[i];
    }

    if(1 == 1) {
      PDM_log_trace_array_int(send_child_extract_idx, n_destinations+1, "send_child_extract_idx ::");
      PDM_log_trace_array_int(recv_child_extract_idx, n_sources+1     , "recv_child_extract_idx ::");
    }

    int *recv_extract_id = malloc(recv_child_extract_idx[n_sources] * sizeof(int));
    PDM_MPI_Neighbor_alltoallv (extract_child_id, send_child_extract_n, send_child_extract_idx, PDM_MPI_INT,
                                recv_extract_id , recv_child_extract_n, recv_child_extract_idx, PDM_MPI_INT, comm_dist_graph);



    double *recv_extract_extents = malloc(6 * recv_child_extract_idx[n_sources] * sizeof(double));
    PDM_MPI_Neighbor_alltoallv (extract_extents     , send_child_extract_n, send_child_extract_idx, mpi_extent_type,
                                recv_extract_extents, recv_child_extract_n, recv_child_extract_idx, mpi_extent_type, comm_dist_graph);

    if(1 == 1) {
      char filename[999];
      sprintf(filename, "recv_extract_extents_%i_%i.vtk", i_iter, i_rank);
      PDM_vtk_write_boxes(filename,
                          recv_child_extract_idx[n_sources],
                          recv_extract_extents,
                          NULL);
    }

    free(send_child_extract_idx);
    free(recv_child_extract_n);
    free(send_child_extract_n);

    free(node_to_child_idx);
    free(extract_child_id);
    free(extract_is_leaf );
    free(extract_extents );
    free(blk_box_to_coarse_box_pts_n);
    free(blk_box_to_coarse_box_pts);

    free(send_request_pts_box_idx);
    free(send_request_pts_box_n);
    free(recv_request_pts_box_n);
    free(recv_request_pts_box_idx);
    free(recv_request_pts_box);

    PDM_MPI_Comm_free(&comm_dist_graph);
    PDM_MPI_Comm_free(&comm_dist_graph_reverse);

    // Il faut faire attention a ne pas trop descendre dans l'arboresence !!!
    if(1 == 1) {
      char filename[999];
      sprintf(filename, "boxes_iter_%i_%i.vtk", i_iter, i_rank);
      PDM_vtk_write_boxes(filename,
                          dn_box,
                          blk_box_extents,
                          NULL);
    }

    // Ok on prepare pour le coup d'après
    free(g_extract_boxes_idx);

    n_neighbor_current       = n_sources;
    g_extract_boxes_idx      = recv_child_extract_idx;
    g_coarse_pts_box_id      = recv_extract_id;
    g_coarse_pts_box_extents = recv_extract_extents;

    PDM_box_set_destroy (&coarse_pts_box_set);
    PDM_box_tree_destroy(&coarse_pts_bt_shared);

    free(neighbor_tag);
    free(neighbor_in);
    free(sources);
    free(destinations);
  }

  // Extraction final pour équilibrer les points
  // sum(box/points) = sur chaque points

  // block_to_part sur l'octree initiale
  PDM_block_to_part_t* btp_new_pts = PDM_block_to_part_create(distrib_pts_box,
                                      (const PDM_g_num_t  **) &g_coarse_pts_box_id,
                                                              &g_extract_boxes_idx[n_neighbor_current],
                                                              1,
                                                              comm);

  int *blk_strid = malloc(dn_node * sizeof(int));
  for(int i = 0; i < dn_node; ++i) {
    blk_strid[i] = coarse_tree_pts->nodes->range[2*i+1] - coarse_tree_pts->nodes->range[2*i];
  }
  double *dpts_coord = NULL;
  PDM_point_tree_seq_sorted_points_get(coarse_tree_pts, &dpts_coord);

  int    **tmp_blk_strid = NULL;
  double **tmp_pts_coord = NULL;
  PDM_block_to_part_exch(btp_new_pts,
                         3 * sizeof(double),
                         PDM_STRIDE_VAR_INTERLACED,
                         blk_strid,
                         dpts_coord,
                         &tmp_blk_strid,
            (void ***)   &tmp_pts_coord);
  free(blk_strid);
  int    *new_pts_n      = tmp_blk_strid[0];
  double *new_pts_coords = tmp_pts_coord[0];
  free(tmp_pts_coord);
  free(tmp_blk_strid);

  int n_new_pts_tot = 0;
  for(int i = 0; i < g_extract_boxes_idx[n_neighbor_current]; ++i) {
    n_new_pts_tot += new_pts_n[i];
  }

  if(1 == 1) {
    char filename[999];
    sprintf(filename, "new_pts_coords_%i.vtk", i_rank);
    PDM_vtk_write_point_cloud(filename,
                              n_new_pts_tot,
                              new_pts_coords,
                              NULL,
                              NULL);
  }


  free(new_pts_coords);
  free(new_pts_n);


  PDM_block_to_part_free(btp_new_pts);

  free(g_extract_boxes_idx);
  free(g_coarse_pts_box_id);
  free(g_coarse_pts_box_extents);

  PDM_point_tree_seq_free(coarse_tree_pts);

  PDM_MPI_Comm_free(&comm_alone);
  PDM_MPI_Type_free(&mpi_extent_type);
  free(distrib_pts_box);

  free(blk_pts_coord);
  free(blk_box_extents);
}



/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief  Main
 *
 */

int
main
(
int argc,
char *argv[]
)
{

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  PDM_g_num_t nPts   = 10;
  double radius = 0.5;
  int local = 0;
  int rand = 0;

  _read_args(argc,
             argv,
             &nPts,
             &radius,
             &local,
             &rand);

  /* Initialize random */

  if (rand) {
    srand(time(NULL));
  }
  else {
    srand(i_rank);
  }

  /* Random point cloud */
  /* Generate src and tgt point clouds */
  int          n_src     = 0;
  double      *src_coord = NULL;

  const double x_center = -0.15;
  const double y_center = 0.15;
  const double z_center = 0.85;

  int         *dback_face_vtx_idx = NULL;
  PDM_g_num_t *dback_face_vtx     = NULL;
  PDM_g_num_t *back_distrib_vtx   = NULL;
  PDM_g_num_t *back_distrib_face  = NULL;
  PDM_sphere_surf_icosphere_gen(comm,
                                nPts,
                                x_center,
                                y_center,
                                z_center,
                                radius,
                                &src_coord,
                                &dback_face_vtx_idx,
                                &dback_face_vtx,
                                &back_distrib_vtx,
                                &back_distrib_face);

  n_src = back_distrib_vtx[i_rank+1] - back_distrib_vtx[i_rank];
  PDM_g_num_t *src_g_num = malloc(n_src * sizeof(PDM_g_num_t));
  for(int i = 0; i < n_src; ++i) {
    src_g_num[i] = back_distrib_vtx[i_rank] + i + 1;
  }

  free(dback_face_vtx_idx);
  free(dback_face_vtx    );
  free(back_distrib_vtx  );
  free(back_distrib_face );

  if(1 == 1) {
    char filename[999];
    sprintf(filename, "sphere_cloud_%i.vtk", i_rank);
    PDM_vtk_write_point_cloud(filename,
                              n_src,
                              src_coord,
                              src_g_num,
                              NULL);
  }




  // PDM_doctree_local_tree_t local_tree_kind = PDM_DOCTREE_LOCAL_TREE_OCTREE;
  // PDM_doctree_local_tree_t local_tree_kind = PDM_DOCTREE_LOCAL_TREE_KDTREE;
  // PDM_doctree_t *doct = PDM_doctree_create(comm,
  //                                          3,
  //                                          1,
  //                                          NULL, // global_extents
  //                                          local_tree_kind);

  int *init_location_pts = malloc(3 * n_src * sizeof(int));
  for(int i = 0; i < n_src; ++i) {
    init_location_pts[3*i  ] = i_rank;
    init_location_pts[3*i+1] = 0; // i_part
    init_location_pts[3*i+2] = i;
  }

  // PDM_doctree_point_set(doct,
  //                       0,
  //                       n_src,
  //                       init_location_pts,
  //                       src_g_num,
  //                       src_coord);
  int n_box   = 0;
  int n_vtx_x = 12;
  int n_vtx_y = 12;
  int n_vtx_z = 12;
  double      *box_extents = NULL;
  PDM_g_num_t *box_gnum    = NULL;
  PDM_box_gen_cartesian(comm,
                        n_vtx_x,
                        n_vtx_y,
                        n_vtx_z,
                        -0., -0., -0.,
                        2., 2., 2.,
                        &n_box,
                        &box_extents,
                        &box_gnum);

  if(1 == 1) {
    char filename[999];
    sprintf(filename, "boxes_%i.vtk", i_rank);
    PDM_vtk_write_boxes(filename,
                        n_box,
                        box_extents,
                        box_gnum);
  }

  int *init_location_box = malloc(3 * n_box * sizeof(int));
  for(int i = 0; i < n_box; ++i) {
    init_location_box[3*i  ] = i_rank;
    init_location_box[3*i+1] = 0; // i_part
    init_location_box[3*i+2] = i;
  }

  // PDM_doctree_solicitation_set(doct,
  //                              PDM_TREE_SOLICITATION_BOXES_POINTS,
  //                              1,
  //                              &n_box,
  //                              &init_location_box,
  //                              &box_gnum,
  //                              &box_extents);

  // PDM_doctree_build(doct);

  _adaptative_tree2(n_src,
                    src_coord,
                    src_g_num,
                    n_box,
                    box_extents,
                    box_gnum,
                    comm);

  // if(1 == 1) {

  //   int         *box_pts_idx = NULL;
  //   PDM_g_num_t *box_pts     = NULL;
  //   double      *pts_coord   = NULL;
  //   PDM_doctree_results_in_orig_frame_get(doct,
  //                                         n_box,
  //                                         box_gnum,
  //                                         &box_pts_idx,
  //                                         &box_pts,
  //                                         &pts_coord);

  //   free(box_pts_idx);
  //   free(box_pts    );
  //   free(pts_coord  );
  // }


  // PDM_doctree_free(doct);


  free(box_gnum);
  free(box_extents);
  free(init_location_box);
  free(init_location_pts);

  free (src_coord);
  free (src_g_num);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }


  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);

  PDM_MPI_Finalize ();

  return 0;

}



// static
// void
// _adaptative_tree
// (
//   int           n_pts,
//   double       *pts_coord,
//   PDM_g_num_t  *pts_gnum,
//   int           n_box,
//   double       *box_extents,
//   PDM_g_num_t  *box_gnum,
//   PDM_MPI_Comm  comm
// )
// {
//   int i_rank;
//   PDM_MPI_Comm_rank (comm, &i_rank);

//   int n_rank;
//   PDM_MPI_Comm_size (comm, &n_rank);
//   /*
//    * Hilbert of points AND boxes
//    */
//   int    *weight_pts = malloc(    n_pts * sizeof(int   ));
//   int    *weight_box = malloc(    n_box * sizeof(int   ));
//   double *box_center = malloc(3 * n_box * sizeof(double));
//   for(int i = 0; i < n_pts; ++i) {
//     weight_pts[i] = 1;
//   }
//   for(int i = 0; i < n_box; ++i) {
//     weight_box[i] = 1;
//     box_center[3*i  ] = 0.5 * (box_extents[6*i  ] + box_extents[6*i+3]);
//     box_center[3*i+1] = 0.5 * (box_extents[6*i+1] + box_extents[6*i+4]);
//     box_center[3*i+2] = 0.5 * (box_extents[6*i+2] + box_extents[6*i+5]);
//   }

//   PDM_part_to_block_t* ptb_pts = PDM_part_to_block_geom_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
//                                                                PDM_PART_TO_BLOCK_POST_CLEANUP, // A voir avec merge mais attention au init_location
//                                                                1.,
//                                                                PDM_PART_GEOM_HILBERT,
//                                                                &pts_coord,
//                                                                &pts_gnum,
//                                                                &weight_pts,
//                                                                &n_pts,
//                                                                1,
//                                                                comm);
//   free(weight_pts);

//   PDM_part_to_block_t* ptb_box = PDM_part_to_block_geom_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
//                                                                PDM_PART_TO_BLOCK_POST_CLEANUP, // A voir avec merge mais attention au init_location
//                                                                1.,
//                                                                PDM_PART_GEOM_HILBERT,
//                                                                &box_center,
//                                                                &box_gnum,
//                                                                &weight_box,
//                                                                &n_box,
//                                                                1,
//                                                                comm);
//   free(weight_box);
//   free(box_center);


//   double *blk_pts_coord = NULL;
//   PDM_part_to_block_exch(ptb_pts,
//                          3 * sizeof(double),
//                          PDM_STRIDE_CST_INTERLACED,
//                          1,
//                          NULL,
//                (void **) &pts_coord,
//                          NULL,
//                (void **) &blk_pts_coord);


//   double *blk_box_extents = NULL;
//   PDM_part_to_block_exch(ptb_box,
//                          6 * sizeof(double),
//                          PDM_STRIDE_CST_INTERLACED,
//                          1,
//                          NULL,
//                (void **) &box_extents,
//                          NULL,
//                (void **) &blk_box_extents);

//   int dn_pts = PDM_part_to_block_n_elt_block_get(ptb_pts);
//   int dn_box = PDM_part_to_block_n_elt_block_get(ptb_box);

//   // Plus besoin du gnum pour l'instant ....
//   PDM_part_to_block_free(ptb_pts);
//   PDM_part_to_block_free(ptb_box);

//   PDM_MPI_Comm comm_alone;
//   PDM_MPI_Comm_split(comm, i_rank, 0, &(comm_alone));

//   PDM_MPI_Datatype mpi_extent_type;
//   PDM_MPI_Type_create_contiguous(6, PDM_MPI_DOUBLE, &mpi_extent_type);
//   PDM_MPI_Type_commit(&mpi_extent_type);

//   /*
//    * Iterative algorithm
//    */
//   int n_iter = 2;
//   for(int i_iter = 0; i_iter < n_iter; ++i_iter) {

//     /* Global extents */
//     double g_global_extents[6];
//     double local_extents[6];

//     for (int i = 0; i < 3; i++) {
//       local_extents[i]     =  DBL_MAX;
//       local_extents[i + 3] = -DBL_MAX;
//     }

//     for (int i = 0; i < dn_box; i++) {
//       for (int j = 0; j < 3; j++) {
//         local_extents[j  ] = PDM_MIN(local_extents[j    ], blk_box_extents[i*3*2 + j    ]);
//         local_extents[j+3] = PDM_MAX(local_extents[j + 3], blk_box_extents[i*3*2 + j + 3]);
//       }
//     }

//     for (int i = 0; i < dn_pts; i++) {
//       for (int  j = 0; j < 3; j++) {
//         if (blk_pts_coord[i*3 + j] < local_extents[j]) {
//           local_extents[j] = blk_pts_coord[i*3 + j];
//         }
//         if (blk_pts_coord[i*3 + j] > local_extents[j + 3]) {
//           local_extents[j + 3] = blk_pts_coord[i*3 + j];
//         }
//       }
//     }

//     PDM_MPI_Allreduce(local_extents  , g_global_extents    , 3, PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);
//     PDM_MPI_Allreduce(local_extents+3, g_global_extents + 3, 3, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

//     double s[3];
//     double d[3];
//     for (int j = 0; j < 3; j++) {
//       s[j] = g_global_extents[j];
//       d[j] = g_global_extents[j+3] - g_global_extents[j];
//     }

//     if(1 == 1) {
//       char filename[999];
//       sprintf(filename, "g_global_extents.vtk");
//       PDM_g_num_t one = 1;
//       PDM_vtk_write_boxes(filename,
//                           1,
//                           g_global_extents,
//                           &one);
//     }

//     // _adaptative_tree_intersect(dn_pts,
//     //                            blk_pts_coord,
//     //                            dn_box,
//     //                            blk_box_extents,
//     //                            comm);

//     /*
//      * Create octree of point
//      */
//     PDM_point_tree_seq_t* coarse_tree_pts = PDM_point_tree_seq_create(PDM_DOCTREE_LOCAL_TREE_OCTREE,
//                                                                       1, // depth_max
//                                                                       1,
//                                                                       1e-8);
//     PDM_point_tree_seq_point_cloud_set(coarse_tree_pts, dn_pts, blk_pts_coord);
//     PDM_point_tree_seq_build(coarse_tree_pts);
//     if(1 == 1) {
//       char filename[999];
//       sprintf(filename, "out_coarse_tree_%i.vtk", i_rank);
//       PDM_point_tree_seq_write_nodes(coarse_tree_pts, filename);
//     }

//     /*
//      * Creation of box tree
//      */
//     PDM_g_num_t *blk_box_gnum      = malloc(dn_box * sizeof(PDM_g_num_t));
//     int         *init_location_box = malloc(3 * dn_box * sizeof(int));
//     for(int i = 0; i < dn_box; ++i) {
//       blk_box_gnum[i] = i + 1;
//       init_location_box[3*i  ] = 0;
//       init_location_box[3*i+1] = 0;
//       init_location_box[3*i+2] = i;
//     }

//     log_trace("dn_box = %i \n", dn_box);
//     PDM_box_set_t  *box_set = PDM_box_set_create(3,
//                                                  0,  // No normalization to preserve initial extents
//                                                  0,  // No projection to preserve initial extents
//                                                  dn_box,
//                                                  blk_box_gnum,
//                                                  blk_box_extents,
//                                                  1,
//                                                  &dn_box,
//                                                  init_location_box,
//                                                  comm_alone);
//     memcpy (box_set->d, d, sizeof(double) * 3);
//     memcpy (box_set->s, s, sizeof(double) * 3);

//     int   max_boxes_leaf_shared = 10; // Max number of boxes in a leaf for coarse shared BBTree
//     int   max_tree_depth_shared = 1;  // Max tree depth for coarse shared BBTree
//     float max_box_ratio_shared  = 5;  // Max ratio for local BBTree (nConnectedBoxe < ratio * nBoxes)
//     PDM_box_tree_t* bt_shared = PDM_box_tree_create (max_tree_depth_shared,
//                                                      max_boxes_leaf_shared,
//                                                      max_box_ratio_shared);

//     PDM_box_tree_set_boxes (bt_shared,
//                             box_set,
//                             PDM_BOX_TREE_ASYNC_LEVEL);

//     free(blk_box_gnum);
//     free(init_location_box);

//     if(1 == 1) {
//       char filename[999];
//       sprintf(filename, "bt_shared_%i.vtk", i_rank);
//       PDM_box_tree_write_vtk(filename, bt_shared, -1, 0);
//     }


//     /*
//      * On a tout les arbres locaux, maintenant on doit interoger les arbres en parallèle
//      */


//     /*
//      * Extract extents on all local_tree
//      */
//     int     n_coarse_pts_box       = 0;
//     int    *coarse_pts_box_id      = NULL;
//     double *coarse_pts_box_extents = NULL;
//     int    *coarse_pts_box_n_pts   = NULL; // Number of point in boxes

//     PDM_point_tree_seq_extract_nodes(coarse_tree_pts,
//                                      0,
//                                      1, // Depth
//                                      &n_coarse_pts_box,
//                                      &coarse_pts_box_id,
//                                      &coarse_pts_box_extents,
//                                      &coarse_pts_box_n_pts);

//     int *n_g_coarse_pts_box = malloc(n_rank * sizeof(int));
//     PDM_MPI_Allgather (&n_coarse_pts_box , 1, PDM_MPI_INT,
//                        n_g_coarse_pts_box, 1, PDM_MPI_INT, comm);

//     int *g_extract_boxes_idx = (int *) malloc (sizeof(int) * (n_rank+1));
//     g_extract_boxes_idx[0] = 0;
//     for(int i = 0; i < n_rank; ++i) {
//       g_extract_boxes_idx[i+1] = g_extract_boxes_idx[i] + n_g_coarse_pts_box[i];
//     }
//     double *g_coarse_pts_box_extents = malloc(6 * g_extract_boxes_idx[n_rank] * sizeof(double));

//     PDM_MPI_Allgatherv(coarse_pts_box_extents  , n_coarse_pts_box, mpi_extent_type,
//                        g_coarse_pts_box_extents, n_g_coarse_pts_box,
//                        g_extract_boxes_idx,
//                        mpi_extent_type, comm);

//     int *g_coarse_pts_box_id = malloc( g_extract_boxes_idx[n_rank] * sizeof(int));
//     PDM_MPI_Allgatherv(coarse_pts_box_id  , n_coarse_pts_box, PDM_MPI_INT,
//                        g_coarse_pts_box_id, n_g_coarse_pts_box,
//                        g_extract_boxes_idx,
//                        PDM_MPI_INT, comm);

//     for(int i = 0; i < n_rank; ++i) {
//       for(int j = g_extract_boxes_idx[i]; j < g_extract_boxes_idx[i+1]; ++j) {
//         g_coarse_pts_box_id[j] += g_extract_boxes_idx[i];
//       }
//     }


//     /*
//      * Build tree
//      */
//     PDM_g_num_t *coarse_pts_box_gnum         = malloc(    g_extract_boxes_idx[n_rank] * sizeof(PDM_g_num_t));
//     int         *init_location_coase_pts_box = malloc(3 * g_extract_boxes_idx[n_rank] * sizeof(int        ));
//     for(int i = 0; i < g_extract_boxes_idx[n_rank]; ++i) {
//       // coarse_pts_box_gnum[i] = g_coarse_pts_box_id[i] + 1;
//       coarse_pts_box_gnum[i] = g_coarse_pts_box_id[i]; // On suppose que root = 0 donc g_id lineraire a partir de 1
//       init_location_coase_pts_box[3*i  ] = 0;
//       init_location_coase_pts_box[3*i+1] = 0;
//       init_location_coase_pts_box[3*i+2] = i+1;
//     }
//     free(g_coarse_pts_box_id);

//     PDM_log_trace_connectivity_long(g_extract_boxes_idx, coarse_pts_box_gnum, n_rank, "coarse_pts_box_gnum ::");

//     PDM_box_set_t  *coarse_pts_box_set = PDM_box_set_create(3,
//                                                             0,  // No normalization to preserve initial extents
//                                                             0,  // No projection to preserve initial extents
//                                                             g_extract_boxes_idx[n_rank],
//                                                             coarse_pts_box_gnum,
//                                                             g_coarse_pts_box_extents,
//                                                             1,
//                                                             &g_extract_boxes_idx[n_rank],
//                                                             init_location_coase_pts_box,
//                                                             comm_alone);
//     memcpy (coarse_pts_box_set->d, d, sizeof(double) * 3);
//     memcpy (coarse_pts_box_set->s, s, sizeof(double) * 3);

//     int   max_boxes_leaf_coarse = 1;   // Max number of boxes in a leaf for coarse coarse BBTree
//     int   max_tree_depth_coarse = 31;  // Max tree depth for coarse coarse BBTree
//     float max_box_ratio_coarse  = 5;   // Max ratio for local BBTree (nConnectedBoxe < ratio * nBoxes)
//     PDM_box_tree_t* coarse_pts_bt_shared = PDM_box_tree_create (max_tree_depth_coarse,
//                                                                 max_boxes_leaf_coarse,
//                                                                 max_box_ratio_coarse);

//     PDM_box_tree_set_boxes (coarse_pts_bt_shared,
//                             coarse_pts_box_set,
//                             PDM_BOX_TREE_ASYNC_LEVEL);

//     free(coarse_pts_box_gnum);
//     free(init_location_coase_pts_box);

//     if(1 == 1) {
//       char filename[999];
//       sprintf(filename, "coarse_pts_box_set_%i.vtk", i_rank);
//       PDM_box_tree_write_vtk(filename, coarse_pts_bt_shared, -1, 0);

//       sprintf(filename, "coarse_pts_box_%i.vtk", i_rank);
//       PDM_vtk_write_boxes(filename,
//                           g_extract_boxes_idx[n_rank],
//                           g_coarse_pts_box_extents,
//                           NULL);
//     }


//     /*
//      * Extract coarse box for box_tree
//      */
//     int n_coarse_box_box = 0;
//     double *coarse_box_extents = NULL;
//     int  n_extract_child = 0;
//     int *extract_child_id = NULL;
//     PDM_box_tree_extract_extents(bt_shared,
//                                  0,
//                                  1,
//                                  &n_coarse_box_box,
//                                  &coarse_box_extents,
//                                  &n_extract_child,
//                                  &extract_child_id);
//     log_trace("n_coarse_box_box = %i \n", n_coarse_box_box);
//     log_trace("n_extract_child  = %i \n", n_extract_child);
//     PDM_log_trace_array_int(extract_child_id, n_extract_child, "extract_child_id ::");

//     if(1 == 1) {
//       char filename[999];
//       sprintf(filename, "coarse_box_box_%i.vtk", i_rank);
//       PDM_vtk_write_boxes(filename,
//                           n_coarse_box_box,
//                           coarse_box_extents,
//                           NULL);
//     }

//     /*
//      * Maintenant on veut interoger l'arbre global des points avec la solicitation des boites
//      *  On peut avoir la connectivité : coarse_pts_box_to_coarse_box
//      *  On devra dans tout les cas échangé des choses
//      */
//     int *coarse_box_to_coarse_box_pts_idx = NULL;
//     int *coarse_box_to_coarse_box_pts = NULL;
//     PDM_box_tree_intersect_boxes_boxes(coarse_pts_bt_shared,
//                                        -1,
//                                        n_coarse_box_box,
//                                        coarse_box_extents,
//                                        &coarse_box_to_coarse_box_pts_idx,
//                                        &coarse_box_to_coarse_box_pts);

//     free(coarse_box_extents);

//     // PDM_log_trace_connectivity_int(coarse_box_to_coarse_box_pts_idx, coarse_box_to_coarse_box_pts, n_coarse_box_box, "coarse_box_to_coarse_box_pts ::");

//     const int         *box_origin   = PDM_box_set_origin_get(coarse_pts_box_set);
//     const PDM_g_num_t *box_pts_gnum = PDM_box_set_get_g_num (coarse_pts_box_set);
//     PDM_log_trace_array_int(box_origin, 3 * g_extract_boxes_idx[n_rank], "box_origin") ;

//     /*
//      * Equilibrate coarse_box_pts -> map on box_gnum_id (of pts )
//      *  Puis chaque connexion correspond à une feuille de l'arbre de boites,
//      *  Qu'on envoie au futur boites de pts
//      *   Donc il faut créer une partition de boites à partir du box_tree -> On veut tout les boxes_id par feuilles en gros
//      *    On est malin -> Si aucune boites intersecte l'arbre de pts on met stride = 0
//      *  Il faut également update le coarse_box_to_coarse_box_pts (avec le nouveau numero de feuilles et l'envoyé)
//      *
//      *
//      *  A la reception des boites on connait donc deja leur lien avec l'octree de pts grossier
//      */
//     int n_connect_box_coarse_box_pts = coarse_box_to_coarse_box_pts_idx[n_coarse_box_box];
//     PDM_g_num_t  *coarse_box_to_coarse_box_pts_gnum = malloc(n_connect_box_coarse_box_pts * sizeof(PDM_g_num_t));
//     double       *weight                            = malloc(n_connect_box_coarse_box_pts * sizeof(double     ));

//     for(int i = 0; i < n_connect_box_coarse_box_pts; ++i) {
//       coarse_box_to_coarse_box_pts_gnum[i] = box_pts_gnum[coarse_box_to_coarse_box_pts[i]]; // + 1;
//       weight[i] = 1.;
//     }

//     PDM_log_trace_connectivity_long(coarse_box_to_coarse_box_pts_idx,
//                                     coarse_box_to_coarse_box_pts_gnum,
//                                     n_coarse_box_box,
//                                     "coarse_box_to_coarse_box_pts_gnum ::");


//     PDM_part_to_block_t* ptb_equi_pts_box = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
//                                                                      PDM_PART_TO_BLOCK_POST_MERGE,
//                                                                      1.,
//                                                                      &coarse_box_to_coarse_box_pts_gnum,
//                                                                      &weight,
//                                                                      &n_connect_box_coarse_box_pts,
//                                                                      1,
//                                                                      comm);

//     // Envoie de boites locals vers les nouvelles feuilles !
//     //
//     int *target_rank = PDM_part_to_block_destination_get(ptb_equi_pts_box);
//     int **boxes_ids = malloc(n_coarse_box_box * sizeof(int *));
//     int  *n_boxes_send = malloc(n_coarse_box_box * sizeof(int));

//     int *send_box_n   = malloc( n_rank    * sizeof(int));
//     int *recv_box_n   = malloc( n_rank    * sizeof(int));
//     int *send_box_idx = malloc((n_rank+1) * sizeof(int));
//     int *recv_box_idx = malloc((n_rank+1) * sizeof(int));
//     send_box_idx[0] = 0;
//     recv_box_idx[0] = 0;
//     for(int i = 0; i < n_rank; ++i) {
//       send_box_n[i] = 0;
//       recv_box_n[i] = 0;
//     }

//     int idx = 0;
//     for(int i_coarse_box = 0; i_coarse_box < n_coarse_box_box; ++i_coarse_box) {

//       if(coarse_box_to_coarse_box_pts_idx[i_coarse_box+1] - coarse_box_to_coarse_box_pts_idx[i_coarse_box] == 0) {
//         continue;
//       }

//       // Pas le choix on doit dededoubler
//       int n_boxes = PDM_box_tree_get_box_ids(bt_shared, extract_child_id[i_coarse_box], &boxes_ids[i_coarse_box]);
//       n_boxes_send[i_coarse_box] = n_boxes;

//       for(int j = coarse_box_to_coarse_box_pts_idx[i_coarse_box]; j < coarse_box_to_coarse_box_pts_idx[i_coarse_box+1]; ++j) {
//         int t_rank = target_rank[idx++];

//         // if(box intersect octree) -> Elimine les boites qui rentre pas dans la feuille de l'octree

//         send_box_n[t_rank] += n_boxes;
//       }
//     }

//     PDM_MPI_Alltoall(send_box_n, 1, PDM_MPI_INT,
//                      recv_box_n, 1, PDM_MPI_INT, comm);

//     for(int i = 0; i < n_rank; ++i) {
//       send_box_idx[i+1] = send_box_idx[i] + send_box_n[i];
//       recv_box_idx[i+1] = recv_box_idx[i] + recv_box_n[i];
//       send_box_n  [i  ] = 0;
//     }

//     double *send_box_extents = malloc(6 * send_box_idx[n_rank] * sizeof(double));
//     double *recv_box_extents = malloc(6 * recv_box_idx[n_rank] * sizeof(double));

//     // PDM_g_num_t *send_box_gnum = malloc(send_box_idx[n_rank] * sizeof(PDM_g_num_t));
//     // PDM_g_num_t *recv_box_gnum = malloc(recv_box_idx[n_rank] * sizeof(PDM_g_num_t));

//     idx = 0;
//     for(int i_coarse_box = 0; i_coarse_box < n_coarse_box_box; ++i_coarse_box) {

//       if(coarse_box_to_coarse_box_pts_idx[i_coarse_box+1] - coarse_box_to_coarse_box_pts_idx[i_coarse_box] == 0) {
//         continue;
//       }

//       for(int j = coarse_box_to_coarse_box_pts_idx[i_coarse_box]; j < coarse_box_to_coarse_box_pts_idx[i_coarse_box+1]; ++j) {
//         int t_rank = target_rank[idx++];

//         int idx_write = send_box_idx[t_rank] + send_box_n[t_rank];

//         for(int k = 0; k < n_boxes_send[i_coarse_box]; ++k) {
//           int box_id = boxes_ids[i_coarse_box][k];
//           for(int p = 0; p < 6; ++p) {
//             send_box_extents[6*idx_write+p] = box_set->local_boxes->extents[6*box_id + p];
//           }

//           // send_box_gnum[idx_write] = box_set->local_boxes->g_num[box_id];

//         }

//         send_box_n[t_rank] += n_boxes_send[i_coarse_box];
//       }
//       free(boxes_ids[i_coarse_box]);
//     }
//     free(boxes_ids);
//     free(n_boxes_send);

//     PDM_log_trace_array_int(send_box_n, n_rank, "send_box_n ::");
//     PDM_log_trace_array_int(recv_box_n, n_rank, "recv_box_n ::");
//     PDM_log_trace_array_int(send_box_idx, n_rank, "send_box_idx ::");
//     PDM_log_trace_array_int(recv_box_idx, n_rank, "recv_box_idx ::");

//     // PDM_MPI_Alltoallv(send_box_gnum, send_box_n, send_box_idx, PDM__PDM_MPI_G_NUM,
//     //                   recv_box_gnum, recv_box_n, recv_box_idx, PDM__PDM_MPI_G_NUM,
//     //                   comm);
//     log_trace("Avat \n");
//     PDM_MPI_Alltoallv(send_box_extents, send_box_n, send_box_idx, mpi_extent_type,
//                       recv_box_extents, recv_box_n, recv_box_idx, mpi_extent_type,
//                       comm);
//     log_trace("Apres \n");

//     /* Replace */
//     free(blk_box_extents);
//     // free(blk_box_gnum);
//     blk_box_extents = recv_box_extents;
//     // blk_box_gnum    = recv_box_gnum;
//     dn_box          = recv_box_idx[n_rank];

//     free(send_box_n);
//     free(send_box_idx);
//     free(recv_box_n);
//     free(recv_box_idx);
//     free(send_box_extents);
//     // free(recv_box_extents);
//     // free(send_box_gnum);
//     // free(recv_box_gnum);

//     free(extract_child_id);

//     free(coarse_box_to_coarse_box_pts_gnum);
//     free(weight);

//     int n_next_coarse_box_pts = PDM_part_to_block_n_elt_block_get(ptb_equi_pts_box);
//     PDM_g_num_t* next_gnum_box_pts = PDM_part_to_block_block_gnum_get(ptb_equi_pts_box);

//     PDM_log_trace_array_long(next_gnum_box_pts, n_next_coarse_box_pts, "next_gnum_box_pts :");

//     /*
//      * Il faut maintenant recupérer les pts des boites de pts !
//      *   block_to_part avec gnum = block_to_part_get_gnum()
//      */
//     PDM_g_num_t* distrib_coarse_box_pts = PDM_compute_entity_distribution(comm, n_coarse_pts_box);
//     PDM_block_to_part_t* btp = PDM_block_to_part_create(distrib_coarse_box_pts,
//                                 (const PDM_g_num_t **)  &next_gnum_box_pts,
//                                                         &n_next_coarse_box_pts,
//                                                         1,
//                                                         comm);

//     /*
//      * Extraction depuis l'octree et envoi
//      */
//     double* tree_coords = NULL;
//     PDM_point_tree_seq_sorted_points_get(coarse_tree_pts, &tree_coords);

//     int point_range[2];
//     int *n_pts_in_leaf = malloc(n_coarse_pts_box * sizeof(int));
//     PDM_log_trace_array_int(coarse_pts_box_id, n_coarse_pts_box, "coarse_pts_box_id ::");
//     assert(n_coarse_pts_box == 8);
//     for(int i = 0; i < n_coarse_pts_box; ++i) {
//       n_pts_in_leaf[coarse_pts_box_id[i]-1] = PDM_point_tree_seq_point_range_get(coarse_tree_pts, coarse_pts_box_id[i], point_range);
//       log_trace(" coarse_pts_box_id[%i] = %i | range = %i / %i \n", i, coarse_pts_box_id[i], point_range[0], point_range[1]);
//     }

//     if(1 == 1) {
//       char filename[999];
//       sprintf(filename, "tree_coords_%i.vtk", i_rank);
//       PDM_vtk_write_point_cloud(filename,
//                                 dn_pts,
//                                 tree_coords,
//                                 NULL,
//                                 NULL);
//     }

//     int    **tmp_next_pts_coords_n = NULL;
//     double **tmp_next_pts_coords   = NULL;
//     PDM_block_to_part_exch(btp,
//                            3 * sizeof(double),
//                            PDM_STRIDE_VAR_INTERLACED,
//                            n_pts_in_leaf,
//                            tree_coords,
//                            &tmp_next_pts_coords_n,
//                 (void ***) &tmp_next_pts_coords);
//     int    *next_pts_coords_n = tmp_next_pts_coords_n[0];
//     double *next_pts_coords   = tmp_next_pts_coords  [0];
//     free(tmp_next_pts_coords_n);
//     free(tmp_next_pts_coords  );

//     dn_pts = 0;
//     for(int i = 0; i < n_next_coarse_box_pts; ++i) {
//       dn_pts += next_pts_coords_n[i];
//     }

//     free(next_pts_coords_n);
//     free(n_pts_in_leaf);

//     free(blk_pts_coord);
//     blk_pts_coord = next_pts_coords;

//     if(1 == 1) {

//       char filename[999];
//       sprintf(filename, "blk_pts_coord_%i.vtk", i_rank);
//       PDM_vtk_write_point_cloud(filename,
//                                 dn_pts,
//                                 blk_pts_coord,
//                                 NULL,
//                                 NULL);
//     }



//     // exit(1);

//     PDM_point_tree_seq_free(coarse_tree_pts);

//     /*
//      * Envoie du lien implicite entre les boites et les pts --> Sous communicateur par feuilles ?
//      *  --> Construire un neighbor --> plus smart :) !!!
//      */
//     free(distrib_coarse_box_pts);
//     PDM_block_to_part_free(btp);

//     PDM_part_to_block_free(ptb_equi_pts_box);


//     PDM_box_set_destroy (&coarse_pts_box_set);
//     PDM_box_tree_destroy(&coarse_pts_bt_shared);

//     free(g_extract_boxes_idx);
//     free(coarse_pts_box_id     );
//     free(coarse_pts_box_extents);
//     free(coarse_pts_box_n_pts  );
//     free(n_g_coarse_pts_box);
//     free(g_coarse_pts_box_extents);
//     free(coarse_box_to_coarse_box_pts_idx);
//     free(coarse_box_to_coarse_box_pts);
//     PDM_box_set_destroy (&box_set);
//     PDM_box_tree_destroy(&bt_shared);

//     /*
//      * Pour l'octree adaptative, il faut faire converger par rang MPI le ratio box / pts
//      *  Donc on sait que il faut rajouter une profondeur sur une feuille en particulier
//      *  Dans l'échange de ptb_equi_pts_box on rajoute le nombre de points par boites de boites
//      *  Si c'est inférieur à une tolérence -> Ok la feuille est deja bien decoupé
//      *  Sinon ben c'est elle qu'on deglingue d'un niveau
//      *  Il faut enlever le truc sur les global_extents et laisser l'octree se faire separement sur les pts / boites
//      *  --> Sinon on contraint l'un avec l'autre -> C'est pas l'idée ici
//      *
//      */


//   }

//   PDM_MPI_Comm_free(&comm_alone);
//   PDM_MPI_Type_free(&mpi_extent_type);

//   free(blk_pts_coord);
//   free(blk_box_extents);

// }
