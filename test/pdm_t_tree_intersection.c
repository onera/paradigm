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
#include "pdm_box_gen.h"
#include "pdm_point_tree_seq.h"
#include "pdm_box_tree.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_unique.h"

#include "pdm_dmesh_nodal.h"
#include "pdm_reader_stl.h"



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
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 PDM_g_num_t   *gn_pts,
 PDM_g_num_t   *gn_box,
 double        *radius,
 int           *tree_type,
 int           *visu
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
        long n = atol(argv[i]);
        *gn_pts = (PDM_g_num_t) n;
      }
    }

    else if (strcmp(argv[i], "-b") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long n = atol(argv[i]);
        *gn_box = (PDM_g_num_t) n;
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

    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *tree_type = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }

    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
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
 int   argc,
 char *argv[]
 )
{
  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  PDM_g_num_t               gn_pts    = 10;
  PDM_g_num_t               gn_box    = 10;
  double                    radius    = 10.;
  PDM_doctree_local_tree_t  tree_type = PDM_DOCTREE_LOCAL_TREE_OCTREE;
  int                       visu      = 0;

  _read_args(argc,
             argv,
             &gn_pts,
             &gn_box,
             &radius,
     (int *) &tree_type,
             &visu);


  double t1, t2;



  /* Random point cloud */
  int          n_pts     = 0;
  double      *pts_coord = NULL;
  PDM_g_num_t *pts_g_num = NULL;
  PDM_point_cloud_gen_random(comm,
                             0, // seed
                             0, // geometric_g_num
                             gn_pts,
                             -radius, -radius, -radius,
                             radius, radius, radius,
                             &n_pts,
                             &pts_coord,
                             &pts_g_num);

  /* Build point tree */
  int depth_max          = 31;
  int points_in_leaf_max = 2;
  const double tolerance = 1e-4;
  PDM_point_tree_seq_t *ptree = PDM_point_tree_seq_create(tree_type,
                                                          depth_max,
                                                          points_in_leaf_max,
                                                          tolerance);

  PDM_point_tree_seq_point_cloud_set(ptree,
                                     n_pts,
                                     pts_coord);

  t1 = PDM_MPI_Wtime();
  PDM_point_tree_seq_build(ptree);
  t2 = PDM_MPI_Wtime();
  double t_point_tree = t2 - t1;
  printf("PDM_point_tree_seq_build        : %12.5es\n", t2 - t1);


  if (visu) {
    char filename2[999];
    sprintf(filename2, "point_tree_%i.vtk", i_rank);
    PDM_point_tree_seq_write_nodes(ptree, filename2);

    sprintf(filename2, "points_%i.vtk", i_rank);
    PDM_vtk_write_point_cloud(filename2,
                              n_pts,
                              pts_coord,
                              pts_g_num,
                              NULL);
  }



  /* (Random) boxes */
  int n_box = 0;
  double      *box_extents = NULL;
  PDM_g_num_t *box_g_num   = NULL;

  double _n = PDM_MAX(2, 1 + pow(gn_box, 1./3.));

  // int n_vtx_x = _n;
  // int n_vtx_y = _n;
  // int n_vtx_z = _n;
  // PDM_box_gen_cartesian(comm,
  //                       n_vtx_x,
  //                       n_vtx_y,
  //                       n_vtx_z,
  //                       -radius, -radius, -radius,
  //                       radius, radius, radius,
  //                       &n_box,
  //                       &box_extents,
  //                       &box_g_num);

  double avg_size = 2*radius/(double) (_n - 1);
  double min_size = 0.5*avg_size;
  double max_size = 1.5*avg_size;
  PDM_box_gen_random(comm,
                     0,
                     0,
                     gn_box,
                     min_size,
                     max_size,
                     -radius, -radius, -radius,
                     radius, radius, radius,
                     &n_box,
                     &box_extents,
                     &box_g_num);

  int *init_location_box = malloc(3 * n_box * sizeof(int));
  for(int i = 0; i < n_box; ++i) {
    init_location_box[3*i  ] = i_rank;
    init_location_box[3*i+1] = 0; // i_part
    init_location_box[3*i+2] = i;
  }

  PDM_box_set_t *box_set = PDM_box_set_create(3,
                                              1,  // No normalization to preserve initial extents
                                              0,  // No projection to preserve initial extents
                                              n_box,
                                              box_g_num,
                                              box_extents,
                                              1,
                                              &n_box,
                                              init_location_box,
                                              comm);

  int   max_boxes_leaf_shared = 10; // Max number of boxes in a leaf for coarse shared BBTree
  int   max_tree_depth_shared = 6; // Max tree depth for coarse shared BBTree
  float max_box_ratio_shared  = 5; // Max ratio for local BBTree (nConnectedBoxe < ratio * nBoxes)

  PDM_box_tree_t *btree = PDM_box_tree_create (max_tree_depth_shared,
                                               max_boxes_leaf_shared,
                                               max_box_ratio_shared);

  t1 = PDM_MPI_Wtime();
  PDM_box_tree_set_boxes (btree,
                          box_set,
                          PDM_BOX_TREE_ASYNC_LEVEL);
  t2 = PDM_MPI_Wtime();
  double t_box_tree = t2 - t1;
  printf("PDM_box_tree_set_boxes          : %12.5es\n", t2 - t1);
  free(init_location_box);


  if (visu) {
    char filename2[999];
    sprintf(filename2, "box_tree_%i.vtk", i_rank);
    PDM_box_tree_write_vtk(filename2,
                           btree,
                           -1,
                           0);

    sprintf(filename2, "boxes_%i.vtk", i_rank);
    PDM_vtk_write_boxes(filename2,
                        n_box,
                        box_extents,
                        box_g_num);
  }


  t1 = PDM_MPI_Wtime();
  int *box_pts_idx = NULL;
  int *box_pts     = NULL;
  PDM_tree_intersection_point_box(btree,
                                  ptree,
                                  &box_pts_idx,
                                  &box_pts);
  t2 = PDM_MPI_Wtime();
  double t_intersection = t2 - t1;
  printf("PDM_tree_intersection_point_box : %12.5es\n", t2 - t1);
  // PDM_log_trace_connectivity_int(box_pts_idx,
  //                                box_pts,
  //                                n_box,
  //                                "box_pts0 : ");

  if (visu) {
    PDM_g_num_t *box_pts_g_num = malloc(sizeof(PDM_g_num_t) * box_pts_idx[n_box]);
    for (int i = 0; i < box_pts_idx[n_box]; i++) {
      box_pts_g_num[i] = pts_g_num[box_pts[i]];
    }

    PDM_log_trace_connectivity_long(box_pts_idx,
                                    box_pts_g_num,
                                    n_box,
                                    "box_pts  : ");
    free(box_pts_g_num);
  }
  free(box_pts_idx);
  free(box_pts);


  t1 = PDM_MPI_Wtime();
  int         *box_pts_idx2   = NULL;
  PDM_g_num_t *box_pts_g_num2 = NULL;
  double      *box_pts_coord2 = NULL;
  PDM_box_tree_points_inside_boxes(btree,
                                   n_pts,
                                   pts_g_num,
                                   pts_coord,
                                   &box_pts_idx2,
                                   &box_pts_g_num2,
                                   &box_pts_coord2);
  t2 = PDM_MPI_Wtime();
  double t_old = t2 - t1;
  printf("PDM_box_tree_points_inside_boxes: %12.5es\n", t2 - t1);

  if (visu) {
    for (int i = 0; i < n_box; i++) {
      PDM_inplace_unique_long(box_pts_g_num2,
                              NULL,
                              box_pts_idx2[i],
                              box_pts_idx2[i+1] - 1);
    }

    PDM_log_trace_connectivity_long(box_pts_idx2,
                                    box_pts_g_num2,
                                    n_box,
                                    "box_pts2 : ");
    free(box_pts_idx2);
    free(box_pts_g_num2);
    free(box_pts_coord2);
  }


  printf("Total intersection : %12.5es\n", t_point_tree + t_box_tree + t_intersection);
  printf("Total old          : %12.5es\n", t_box_tree + t_old);
  /* Free */
  PDM_point_tree_seq_free(ptree);
  PDM_box_tree_destroy(&btree);
  PDM_box_set_destroy (&box_set);

  free(pts_coord);
  free(pts_g_num);

  free(box_extents);
  free(box_g_num);
  PDM_MPI_Finalize ();

  return 0;
}
