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
#include "pdm_closest_points.h"
#include "pdm_version.h"


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
     "  -c       <level> Number of closest points (default : 10).\n\n"
     "  -s       <level> Number of Source points (default : 10).\n\n"
     "  -t       <level> Number of Target points (default : 10).\n\n"
     "  -radius <level>  Radius of domain (default : 10).\n\n"
     "  -local           Number of points is local (default : global).\n\n"
     "  -rand            Random definition of point coordinates (default : false).\n\n"
     "  -clumps          Source points distributed in clumps around target points (default : false).\n\n"
     "  -h               This message.\n\n");

  exit (exit_code);
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

  int i_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  assert (n_rank > 1);

  int n_src = 4;
  PDM_g_num_t src_g_num[4];
  double      src_coord[4*3];

  int n_tgt = 18;
  PDM_g_num_t tgt_g_num[18];
  double      tgt_coord[18*3];

  if (i_rank == 0) {
    int idx = 0;
    for (int j = 0; j < 2; j++) {
      for (int i = 0; i < 2; i++) {
        src_g_num[idx] = idx + 1;
        src_coord[3*idx]     = 2*i + 1;
        src_coord[3*idx + 1] = 2*j + 1;
        src_coord[3*idx + 2] = 1;
        idx++;
      }
    }

    idx = 0;
    for (int k = 0; k < 2; k++) {
      for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 3; i++) {
          tgt_g_num[idx] = idx + 1;
          tgt_coord[3*idx]     = 2*i;
          tgt_coord[3*idx + 1] = 2*j;
          tgt_coord[3*idx + 2] = 2*k;
          idx++;
        }
      }
    }
  }

  else if (i_rank == 1) {
    int idx = 0;
    for (int j = 0; j < 2; j++) {
      for (int i = 0; i < 2; i++) {
        src_g_num[idx] = idx + 5;
        src_coord[3*idx]     = 2*i + 1;
        src_coord[3*idx + 1] = 2*j + 1;
        src_coord[3*idx + 2] = 1;
        idx++;
      }
    }

    idx = 0;
    for (int k = 0; k < 2; k++) {
      for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 3; i++) {
          tgt_g_num[idx] = idx + 10;
          tgt_coord[3*idx]     = 2*i;
          tgt_coord[3*idx + 1] = 2*j;
          tgt_coord[3*idx + 2] = 2*k;
          idx++;
        }
      }
    }
  }

  else {
    n_src = 0;
    n_tgt = 0;
  }



  PDM_closest_point_t *clsp = PDM_closest_points_create (PDM_MPI_COMM_WORLD,
                                                         1,
                                                         PDM_OWNERSHIP_USER);



  PDM_closest_points_n_part_cloud_set (clsp, 1, 1);

  PDM_closest_points_src_cloud_set (clsp,
                                    0,
                                    n_src,
                                    src_coord,
                                    src_g_num);

  PDM_closest_points_tgt_cloud_set (clsp,
                                    0,
                                    n_tgt,
                                    tgt_coord,
                                    tgt_g_num);


  PDM_closest_points_compute (clsp);



  PDM_g_num_t *closest_src_gnum = NULL;
  double      *closest_src_dist = NULL;

  PDM_closest_points_get (clsp,
                          0,
                          &closest_src_gnum,
                          &closest_src_dist);

  for (int i = 0; i < n_tgt; i++) {
    printf("[%d] "PDM_FMT_G_NUM" : "PDM_FMT_G_NUM", %f\n",
           i_rank, tgt_g_num[i],
           closest_src_gnum[i],
           closest_src_dist[i]);
  }



  PDM_MPI_Finalize ();

  return 0;
}
