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

/*============================================================================
 * Type definitions
 *============================================================================*/



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

  /*
   * Redistribute all pts and impose hilbert ordering
   */
  int **weight = malloc(doct->n_part_cloud * sizeof(int *));
  for(int i_part = 0; i_part < doct->n_part_cloud; ++i_part) {
    weight[i_part] = malloc(doct->n_point_cloud[i_part] * sizeof(int));
    for(int i = 0; i < doct->n_point_cloud[i]; ++i) {
      weight[i_part][i] = 1;
    }
  }


  PDM_part_to_block_t* ptb = PDM_part_to_block_geom_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_MERGE,
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

  int          n_parent    = PDM_part_to_block_n_elt_block_get  (ptb);
  PDM_g_num_t* parent_gnum = PDM_part_to_block_block_gnum_get   (ptb);
  PDM_g_num_t* distrib_pts = PDM_part_to_block_distrib_index_get(ptb);

  free(blk_pts_coord);

  /*
   * Step 2 : Create coarse octree to equilibrate leaf
   *   --> Il faut la solicitation
   */



  /*
   * Step 3 : Build local octree
   */
  int dn_pts = distrib_pts[i_rank+1] - distrib_pts[i_rank];

  if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_OCTREE) {

    assert(doct->local_octree == NULL);
    doct->local_octree = PDM_octree_seq_create(1, // n_point_cloud
                                               doct->local_depth_max,
                                               doct->local_points_in_leaf_max,
                                               doct->local_tolerance);

    PDM_octree_seq_point_cloud_set(doct->local_octree,
                                   0,
                                   dn_pts,
                                   blk_pts_coord);

    PDM_octree_seq_build(doct->local_octree);


  } else {
    abort();
  }


  /*
   * Setup global tree to orien resarch in parallel
   *    We take the first two level of the tree
   */
  int n_depth_per_proc = 2;

  /*
   * Extract extents on all local_tree
   */
  int n_coarse_box = 0;
  double *coarse_box_extents = NULL;
  if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_OCTREE) {
    PDM_octree_seq_extract_extent(doct->local_octree,
                                  0,
                                  n_depth_per_proc,
                                  &n_coarse_box,
                                  &coarse_box_extents);
  } else {
    abort();
  }


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

  double *gcoarse_box_extents = malloc (g_coarse_box_idx[n_rank] * 6 * sizeof(double));

  for(int i = 0; i < n_rank; ++i) {
    g_coarse_box_n  [i] *= 6;
    g_coarse_box_idx[i] *= 6;
  }
  PDM_MPI_Allgatherv (coarse_box_extents , 6 * n_coarse_box, PDM_MPI_DOUBLE,
                      gcoarse_box_extents, g_coarse_box_n  , g_coarse_box_idx, PDM_MPI_DOUBLE,
                      doct->comm);

  free(g_coarse_box_n);
  free(g_coarse_box_idx);
  free(coarse_box_extents);
  free(gcoarse_box_extents);

  PDM_part_to_block_free(ptb);

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
PDM_doctree_free
(
  PDM_doctree_t   *doct
)
{
  free(doct->n_point_cloud    );
  free(doct->pts_g_num        );
  free(doct->pts_coords       );
  free(doct->pts_init_location);

  if(doct->local_octree != NULL) {
    PDM_octree_seq_free(doct->local_octree);
  }

  free(doct);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
