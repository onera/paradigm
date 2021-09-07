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
#include "pdm_para_octree.h"
#include "pdm_distrib.h"
#include "pdm_logging.h"

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
 PDM_g_num_t   *gn_box,
 double        *length
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
        long _n = atol(argv[i]);
        *gn_box = (PDM_g_num_t) _n;
      }
    }

    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *length = atof(argv[i]);
      }
    }

    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}



static void
_random_boxes
(
 PDM_MPI_Comm        comm,
 const PDM_g_num_t   gn_box,
 double              length,
 int                *n_box,
 PDM_g_num_t       **box_g_num,
 double            **box_extents
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  *n_box = (int) (gn_box / n_rank);
  if (i_rank < gn_box % n_rank) {
    (*n_box)++;
  }

  PDM_g_num_t* distrib_box = PDM_compute_entity_distribution(comm, (*n_box));
  for(int i = 0; i < 6 * distrib_box[i_rank]; ++i) {
    rand();
  }

  *box_g_num = malloc (sizeof(PDM_g_num_t) * (*n_box));
  for (int i = 0; i < *n_box; i++) {
    (*box_g_num)[i] = distrib_box[i_rank] + i + 1;
  }

  *box_extents = malloc (sizeof(double) * (*n_box) * 6);
  for (int i = 0; i < *n_box; i++) {
    for (int j = 0; j < 3; j++) {
      double x1 = length * (double) rand() / ((double) RAND_MAX);
      double x2 = length * (double) rand() / ((double) RAND_MAX);

      (*box_extents)[6*i + j    ] = PDM_MIN (x1, x2);
      (*box_extents)[6*i + j + 3] = PDM_MAX (x1, x2);
    }
  }

  free (distrib_box);
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
main (int argc, char *argv[])
{
  srand(0);

  PDM_MPI_Init (&argc, &argv);

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);



  PDM_g_num_t gn_box = 10;
  double      length = 3.14;

  _read_args (argc,
              argv,
              &gn_box,
              &length);



  /*
   *  Random boxes
   */
  int n_box;
  PDM_g_num_t *box_g_num   = NULL;
  double      *box_extents = NULL;
  _random_boxes (comm,
                 gn_box,
                 length,
                 &n_box,
                 &box_g_num,
                 &box_extents);

  if (1) {
    for (int i = 0; i < n_box; i++) {
      log_trace("box "PDM_FMT_G_NUM" extents = %f %f %f  %f %f %f\n",
                box_g_num[i],
                box_extents[6*i + 0],
                box_extents[6*i + 1],
                box_extents[6*i + 2],
                box_extents[6*i + 3],
                box_extents[6*i + 4],
                box_extents[6*i + 5]);
    }
  }


  /*
   *  Random point
   */
  int n_pts = 0;
  double      *pts_coord = NULL;
  PDM_g_num_t *pts_g_num = NULL;
  if (i_rank == 0) {
    n_pts = 1;
    pts_coord = malloc (sizeof(double) * n_pts * 3);
    pts_g_num = malloc (sizeof(PDM_g_num_t) * n_pts);

    for (int i = 0; i < 3; i++) {
      pts_coord[i] = 0.5 * length;
    }

    pts_g_num[0] = 1;
  }



  const int octree_depth_max = 15;//?
  const int octree_points_in_leaf_max = 1;
  const int octree_build_leaf_neighbours = 0;

  /* Create octree structure */
  int octree_id = PDM_para_octree_create (1,
                                          octree_depth_max,
                                          octree_points_in_leaf_max,
                                          octree_build_leaf_neighbours,
                                          comm);

  /* Set octree point cloud */
  PDM_para_octree_point_cloud_set (octree_id,
                                   0,
                                   n_pts,
                                   pts_coord,
                                   pts_g_num);

  /* Build parallel octree */
  PDM_para_octree_build (octree_id, NULL);



  /* Find points inside boxes */
  int         *box_pts_idx   = NULL;
  PDM_g_num_t *box_pts_g_num = NULL;
  double      *box_pts_coord = NULL;
  PDM_para_octree_points_inside_boxes_with_copies (octree_id,
                                                   n_box,
                                                   box_extents,
                                                   box_g_num,
                                                   &box_pts_idx,
                                                   &box_pts_g_num,
                                                   &box_pts_coord);

  /* Free octree */
  PDM_para_octree_free (octree_id);



  PDM_log_trace_connectivity_long (box_pts_idx,
                                   box_pts_g_num,
                                   n_box,
                                   "box_pts_g_num :");

  free (box_extents);
  free (box_g_num);
  free (box_pts_idx);
  free (box_pts_g_num);
  free (box_pts_coord);

  if (i_rank == 0) {
    free (pts_coord);
    free (pts_g_num);
  }


  PDM_MPI_Finalize();

  return 0;
}
