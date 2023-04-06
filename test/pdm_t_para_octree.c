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
#include "pdm_para_octree.h"
#include "pdm_box_gen.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

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
_loop_forward
(
  int  n_nodes,
  int *n_points,
  int *range,
  int *leaf_id,
  int *children_id,
  int *ancestor_id,
  int  n_child,
  int  stack_size
)
{
  PDM_UNUSED(n_nodes);
  PDM_UNUSED(ancestor_id);

  int *stack_id = malloc (sizeof(int) * stack_size);
  int pos_stack = 0;
  stack_id[pos_stack++] = 0;

  int dbg_enabled = 1;
  while (pos_stack > 0) {
    int node_id = stack_id[--pos_stack];

    if (dbg_enabled) {
      log_trace("  node %d : range=%d, n_points=%d, leaf_id=%d\n",
             node_id, range[node_id], n_points[node_id], leaf_id[node_id]);
    }

    if (leaf_id[node_id] >= 0) {
      for (int i = 0; i < n_points[node_id]; i++) {
        int ipt = range[node_id] + i;
        log_trace("\t is leaf ipt = %i \n", ipt);
      }
    }  else {
      for (int i = 0; i < n_child; i++) {
        int child_id = children_id[n_child*node_id + i];
        if (child_id < 0) {
          continue;
        }
        if (dbg_enabled) {
          log_trace("    child %d: id=%d, range=%d, n_points=%d, leaf_id=%d\n",
                 i, child_id, range[child_id], n_points[child_id], leaf_id[child_id]);
        }
        stack_id[pos_stack++] = child_id;
      }
    }
  }
  free(stack_id);
}


static
void
_loop_backward
(
  int  n_nodes,
  int *n_points,
  int *range,
  int *leaf_id,
  int *children_id,
  int *ancestor_id,
  int  n_child,
  int  stack_size
)
{
  PDM_UNUSED(n_nodes);
  PDM_UNUSED(ancestor_id);
  PDM_UNUSED(range);
  PDM_UNUSED(n_child);
  PDM_UNUSED(stack_size);
  PDM_UNUSED(children_id);
  PDM_UNUSED(n_points);

  PDM_log_trace_array_int(leaf_id    , n_nodes, "leaf_id     :");
  PDM_log_trace_array_int(ancestor_id, n_nodes, "ancestor_id :");


  int *node_tag = malloc(n_nodes * sizeof(int));
  int n_ancestor = 0;
  for(int i = 0; i < n_nodes; ++i){
    node_tag[i] = -1;
  }

  /* Recupération des leafs */
  int *extract_leaf_id  = malloc(n_nodes * sizeof(int));
  int *current_ancestor = malloc(n_nodes * sizeof(int));
  int n_leaf = 0;
  for(int i = 0; i < n_nodes; ++i) {
    if(leaf_id[i] != -1){
      extract_leaf_id[n_leaf++] = leaf_id[i];
      int i_ancestor = ancestor_id[i];
      if(i_ancestor >= 0 && node_tag[i_ancestor] == -1){
        node_tag[i_ancestor] = 0;  // 0  :  In stack
        // current_ancestor[n_ancestor++] = i_ancestor;
        current_ancestor[n_ancestor++] = i_ancestor;
      }
    }
  }
  extract_leaf_id = realloc(extract_leaf_id, n_leaf * sizeof(int));

  current_ancestor   = realloc(current_ancestor, n_ancestor * sizeof(int));
  int *next_ancestor = malloc(n_ancestor * sizeof(int));

  // PDM_log_trace_array_int(extract_leaf_id, n_leaf, "extract_leaf_id :");
  // PDM_log_trace_array_int(current_ancestor, n_ancestor, "current_ancestor :");
  // PDM_log_trace_array_int(n_points, n_nodes, "n_points :");

  while(n_ancestor > 0) {

    PDM_log_trace_array_int(current_ancestor, n_ancestor, "current_ancestor :");

    int n_ancestor_next = 0;
    for(int i = 0; i < n_ancestor; ++i) {
      int node_id = current_ancestor[i];
      node_tag[node_id] = 2; // is treated

      int i_ancestor = ancestor_id[node_id];
      if(i_ancestor >= 0 && node_tag[i_ancestor] == -1){
        node_tag[i_ancestor] = 0; // In stack
        next_ancestor[n_ancestor_next++] = i_ancestor;
      }
    }

    /* Swap and update */
    int *tmp_ancestor = current_ancestor;
    current_ancestor  = next_ancestor;
    next_ancestor     = tmp_ancestor;
    n_ancestor        = n_ancestor_next;
  }

  PDM_log_trace_array_int(node_tag, n_nodes, "node_tag :");

  for(int i = 0; i < n_nodes; ++i) {
    if(leaf_id[i] == -1){
      assert(node_tag[i] == 2);
    }
  }

  free(extract_leaf_id);
  free(next_ancestor);
  free(current_ancestor);
  free(node_tag);

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
  double radius = 10.;
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
  PDM_g_num_t *src_g_num = NULL;
  PDM_point_cloud_gen_random (comm,
                              0, // seed
                              0, // geometric_g_num
                              nPts,
                              0., 0., 0.,
                              1., 1., 1.,
                              &n_src,
                              &src_coord,
                              &src_g_num);

  const int depth_max = 31;
  const int points_in_leaf_max = 1;
  PDM_para_octree_t *octree = PDM_para_octree_create(1,
                                                   depth_max,
                                                   points_in_leaf_max,
                                                   1,  //Build neigbor
                                                   comm);

  PDM_para_octree_point_cloud_set(octree,
                                  0,
                                  n_src,
                                  src_coord,
                                  src_g_num);

  PDM_para_octree_build(octree, NULL);
  PDM_para_octree_dump(octree);

  /*
   *  Creation du part_to_part pour transférer d'un nuage de point user <-> octree
   *  Rajouter dans l'API un init_location
   */

  // PDM_para_octree_leaf_get(octree);
  int  n_nodes     = 0;
  int *n_points    = NULL;
  int *range       = NULL;
  int *leaf_id     = NULL;
  int *children_id = NULL;
  int *ancestor_id = NULL;
  int  n_child     = 0;
  int  stack_size  = 0;
  PDM_para_octree_explicit_node_get(octree,
                                    &n_nodes,
                                    &n_points,
                                    &range,
                                    &leaf_id,
                                    &children_id,
                                    &ancestor_id,
                                    &n_child,
                                    &stack_size);

  int          n_pts_octree     = 0;
  double      *pts_coord_octree = NULL;
  PDM_g_num_t *pts_gnum_octree  = NULL;
  PDM_para_octree_points_get(octree,
                             &n_pts_octree,
                             &pts_coord_octree,
                             &pts_gnum_octree);

  if(0 == 1) {
    PDM_log_trace_array_long(pts_gnum_octree   ,     n_pts_octree, "pts_gnum_octree  :");
    PDM_log_trace_array_double(pts_coord_octree, 3 * n_pts_octree, "pts_coord_octree :");
  }

  /* Loop  descending */
  _loop_forward(n_nodes,
                n_points,
                range,
                leaf_id,
                children_id,
                ancestor_id,
                n_child,
                stack_size);

  /* Loop backward */
  _loop_backward(n_nodes,
                n_points,
                range,
                leaf_id,
                children_id,
                ancestor_id,
                n_child,
                stack_size);


  if(1 == 1) {
    PDM_para_octree_export_vtk(octree, "export_octree");
  }

  PDM_para_octree_free(octree);

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
