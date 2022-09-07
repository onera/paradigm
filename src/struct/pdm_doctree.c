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
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_binary_search.h"
#include "pdm_octree_seq.h"
#include "pdm_kdtree_seq.h"
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
  int **weight = malloc(doct->n_part_cloud * sizeof(int *));
  for(int i_part = 0; i_part < doct->n_part_cloud; ++i_part) {
    weight[i_part] = malloc(doct->n_point_cloud[i_part] * sizeof(int));
    for(int i = 0; i < doct->n_point_cloud[i_part]; ++i) {
      weight[i_part][i] = 1;
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

  double *blk_pts_coord = NULL;
  PDM_part_to_block_exch(ptb,
                         3 * sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         1,
                         NULL,
               (void **) doct->pts_coords,
                         NULL,
               (void **) &blk_pts_coord);

  /* Transport init_location - Attention au merge du ptb à faire */

  // int          n_parent    = PDM_part_to_block_n_elt_block_get  (ptb);
  // PDM_g_num_t* parent_gnum = PDM_part_to_block_block_gnum_get   (ptb);
  // PDM_g_num_t* distrib_pts = PDM_part_to_block_distrib_index_get(ptb);

  *ptb_out         = ptb;
  *dpts_coords_out = blk_pts_coord;

}



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

  doct->global_depth_max          = 5;
  doct->global_points_in_leaf_max = 60;

  doct->local_depth_max          = 5;
  doct->local_points_in_leaf_max = 30;
  doct->local_tolerance          = 1e-6;

  doct->local_tree_kind = local_tree_kind;
  doct->global_octree   = NULL;
  doct->local_octree    = NULL;
  doct->shmem_octree    = NULL;

  doct->comm_shared   = PDM_MPI_COMM_NULL;

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

  return doct;
}


void
PDM_doctree_build
(
 PDM_doctree_t     *doct
)
{
  int i_rank;
  PDM_MPI_Comm_rank (doct->comm, &i_rank);
  int n_rank;
  PDM_MPI_Comm_size (doct->comm, &n_rank);

  double              *blk_pts_coord = NULL;
  PDM_part_to_block_t *ptb           = NULL;
  _redistribute_pts_geom(doct, &ptb, &blk_pts_coord);

  PDM_g_num_t* distrib_pts = PDM_part_to_block_distrib_index_get(ptb);

  /*
   * Step 2 : Build local coarse tree
   */
  int dn_pts = distrib_pts[i_rank+1] - distrib_pts[i_rank];
  PDM_octree_seq_t *coarse_octree = NULL;
  PDM_kdtree_seq_t *coarse_kdtree = NULL;
  if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_OCTREE) {

    assert(coarse_octree == NULL);
    coarse_octree = PDM_octree_seq_create(1, // n_point_cloud
                                          doct->local_depth_max,
                                          doct->local_points_in_leaf_max,
                                          doct->local_tolerance);

    PDM_octree_seq_point_cloud_set(coarse_octree,
                                   0,
                                   dn_pts,
                                   blk_pts_coord);

    PDM_octree_seq_build(coarse_octree);


  } else if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_KDTREE){

    assert(coarse_kdtree == NULL);
    coarse_kdtree = PDM_kdtree_seq_create(1, // n_point_cloud
                                          doct->local_depth_max,
                                          doct->local_points_in_leaf_max,
                                          doct->local_tolerance);

    PDM_kdtree_seq_point_cloud_set(coarse_kdtree,
                                   0,
                                   dn_pts,
                                   blk_pts_coord);

    PDM_kdtree_seq_build(coarse_kdtree);
  }

  /*
   * Extract extents on all local_tree
   */
  int n_coarse_box = 0;
  double *coarse_box_extents = NULL;
  if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_OCTREE) {
    int n_depth_per_proc = 2;
    PDM_octree_seq_extract_extent(coarse_octree,
                                  0,
                                  n_depth_per_proc,
                                  &n_coarse_box,
                                  &coarse_box_extents);
  } else if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_KDTREE){
    int n_depth_per_proc = 16;
    PDM_kdtree_seq_extract_extent(coarse_kdtree,
                                  0,
                                  n_depth_per_proc,
                                  &n_coarse_box,
                                  &coarse_box_extents);
  }

  /*
   * Step 2 : Create coarse octree to equilibrate leaf
   *   --> Il faut la solicitation
   *   Creation bt_shared temporaire
   *     puis solicitation
   *   On connait le lien shared_to_box
   */






   // PDM_box_tree_intersect_boxes_boxes2(bt_shared,
   //                                      -1,
   //                                      n_boxes,
   //                                      box_extents,
   //                                      &shared_to_box_idx,
   //                                      &shared_to_box);
   // Preparation of send count and box_rank/box_rank_idx
   // for(int i = 0; i < n_rank; ++i) {
   //   for(int j = shared_all_rank_idx[i]; j < shared_all_rank_idx[i+1]; ++j) {
   //     send_count[i] += shared_to_box_idx[j+1] - shared_to_box_idx[j];
   //   }
   // }

   /*
    * Le nouveau tri donne le lien old_to_new_rank (car on permutera les bbox a peu de choses près)
    * Une fois qu'on connait le tri, on peut faire l'échange en asynchrone pdt que l'arbre se construit ?
    *
    */

  /*
   * Build a box_tree
   */
  int *g_coarse_box_n = (int *) malloc (n_rank * sizeof(int));
  PDM_MPI_Allgather (&n_coarse_box , 1, PDM_MPI_INT,
                     g_coarse_box_n, 1, PDM_MPI_INT,
                     doct->comm);

  int *g_coarse_box_idx = malloc((n_rank + 1 ) *sizeof(int));
  g_coarse_box_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    g_coarse_box_idx[i+1] =  g_coarse_box_idx[i] + g_coarse_box_n[i];
  }

  if(1 == 1) {
    PDM_log_trace_array_int(g_coarse_box_idx, n_rank+1, "g_coarse_box_idx ::");
  }

  double      *gcoarse_box_extents = malloc(g_coarse_box_idx[n_rank] * 6 * sizeof(double     ));
  PDM_g_num_t *g_coarse_box_color  = malloc(g_coarse_box_idx[n_rank] *     sizeof(PDM_g_num_t));

  for(int i = 0; i < n_rank; ++i) {
    for(int j = g_coarse_box_idx[i]; j < g_coarse_box_idx[i+1]; ++j) {
      g_coarse_box_color[j] = i;
    }
  }

  for(int i = 0; i < n_rank; ++i) {
    g_coarse_box_n  [i] *= 6;
    g_coarse_box_idx[i] *= 6;
  }
  PDM_MPI_Allgatherv (coarse_box_extents , 6 * n_coarse_box, PDM_MPI_DOUBLE,
                      gcoarse_box_extents, g_coarse_box_n  , g_coarse_box_idx, PDM_MPI_DOUBLE,
                      doct->comm);

  free(coarse_box_extents);

  if(1 == 1 && i_rank == 0) {
    char filename[999];
    sprintf(filename, "gcoarse_box_extents_%i.vtk", i_rank);
    PDM_vtk_write_boxes(filename,
                        g_coarse_box_idx[n_rank],
                        gcoarse_box_extents,
                        g_coarse_box_color);
  }


  /*
   *  On pourrait faire du hilbert sur des niveaux de feuilles avec gnum = node_id implicitement odered
   *  Avec du sampling
   *  On encode en hilbert le millieu des feuilles -> meme en kd-tree ca marchera
   *  On fait la pré-soliciation pour avoir des poids par boîtes
   *  On fait part_to_block sur des child_id implictement hilbert puis on echange le contenu des noeuds
   *   La stride = le nombre de points -> Ca fait des echanges mais osef !!!
   *   L'algo ressemble ENORMEMENT à PDM_part_assemble_partitions
   *   La pré-solicitation nous permet également d'avoir le lien grossier
   *   Si on le preserve on n'a plus a interoger l'octree !!!
   *   Pdt le transfert de la pré-solicitation --> On construit l'arbre fin
   *   On peut également utiliser l'info du g_child_id dans lequel on est solicité pour preconditionné la recherche local (on gagnera 3/4 niveaux)
   *   A affiner avec le double niveau node / numa
   *   Reprendre le Allgatherv du para_octree sur le comm circulaire
   */
  int n_shared_boxes = g_coarse_box_idx[n_rank];
  const int n_info_location = 3;
  int *init_location_proc = PDM_array_zeros_int (n_info_location * n_shared_boxes);

  int   max_boxes_leaf_shared = 10; // Max number of boxes in a leaf for coarse shared BBTree
  int   max_tree_depth_shared = 4; // Max tree depth for coarse shared BBTree
  float max_box_ratio_shared  = 5; // Max ratio for local BBTree (nConnectedBoxe < ratio * nBoxes)

  PDM_MPI_Comm bt_comm;
  PDM_MPI_Comm_split (doct->comm, i_rank, 0, &bt_comm);
  PDM_box_set_t  *box_set   = PDM_box_set_create(3,             // dim
                                                 1,             // normalize
                                                 0,             // allow_projection
                                                 n_shared_boxes,
                                                 g_coarse_box_color,
                                                 gcoarse_box_extents,
                                                 1,
                                                 &n_shared_boxes,
                                                 init_location_proc,
                                                 doct->comm);

  PDM_box_tree_t* bt_shared = PDM_box_tree_create (max_tree_depth_shared,
                                                   max_boxes_leaf_shared,
                                                   max_box_ratio_shared);

  PDM_box_tree_set_boxes (bt_shared,
                          box_set,
                          PDM_BOX_TREE_ASYNC_LEVEL);

  PDM_box_set_destroy (&box_set);
  PDM_box_tree_destroy(&bt_shared);
  PDM_MPI_Comm_free (&bt_comm);


  free(g_coarse_box_n);
  free(g_coarse_box_idx);
  free(gcoarse_box_extents);
  free(g_coarse_box_color);
  free(init_location_proc);

  free(blk_pts_coord);

  PDM_part_to_block_free(ptb);

  /*
   * Si l'octree ne change pas l'ordre des points -> On fait une allocation shared
   */

  /*
   *  L'octree seq fait un tri indirect sur les pts, c'est maybe possible de le faire shared inplace !
   *  Sinon on fait un chapeau pour le faire shared !!!!!
   */




  /*
   * A revoir quand on fera l'équilibrage
   */
  // PDM_octree_seq_free(coarse_octree);
  doct->local_octree = coarse_octree;
  doct->local_kdtree = coarse_kdtree;
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

}

// void
// PDM_doctree_results_get
// (
//  PDM_doctree_t      *doct,
//  int                *init_location_entity,
//  int                *dentity1_entity2_idx,
//  PDM_g_num_t        *dentity1_entity2,
//  int                *
//  double             *entity_coords
// )
// {

// }


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

  if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_OCTREE) {
    if(doct->local_octree != NULL) {
      PDM_octree_seq_free(doct->local_octree);
    }
  } else if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_KDTREE){
    if(doct->local_kdtree != NULL) {
      PDM_kdtree_seq_free(doct->local_kdtree);
    }
  }

  free(doct);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
