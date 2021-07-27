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
 int argc,
 char *argv[]
 )
{

  PDM_MPI_Init (&argc, &argv);

  int i_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);

  int numProcs;
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);


  /* Define the numbers of Source/Target points */

  int n_src = 4;
  PDM_g_num_t src_gnum[4];
  double src_coords[12];
  int n_tgt = 18;
  PDM_g_num_t tgt_gnum[18];
  double tgt_coords[54];
  if (i_rank == 0){
    src_gnum[0] = 1;
    src_gnum[1] = 2;
    src_gnum[2] = 3;
    src_gnum[3] = 4;
    src_coords[3*0+0] = 0.25; src_coords[3*0+1] = 0.25; src_coords[3*0+2] = 0.25;
    src_coords[3*1+0] = 0.75; src_coords[3*1+1] = 0.25; src_coords[3*1+2] = 0.25;
    src_coords[3*2+0] = 0.75; src_coords[3*2+1] = 0.75; src_coords[3*2+2] = 0.25;
    src_coords[3*3+0] = 0.25; src_coords[3*3+1] = 0.75; src_coords[3*3+2] = 0.25;

    for (int i=0; i < 18; i ++)
      tgt_gnum[i] = i+1;

    tgt_coords[3* 0+0] = 0.55; tgt_coords[3* 0+1] = 0.05; tgt_coords[3* 0+2] = -0.05;
    tgt_coords[3* 1+0] = 1.05; tgt_coords[3* 1+1] = 0.05; tgt_coords[3* 1+2] = -0.05;
    tgt_coords[3* 2+0] = 1.55; tgt_coords[3* 2+1] = 0.05; tgt_coords[3* 2+2] = -0.05;
    tgt_coords[3* 3+0] = 0.55; tgt_coords[3* 3+1] = 0.55; tgt_coords[3* 3+2] = -0.05;
    tgt_coords[3* 4+0] = 1.05; tgt_coords[3* 4+1] = 0.55; tgt_coords[3* 4+2] = -0.05;
    tgt_coords[3* 5+0] = 1.55; tgt_coords[3* 5+1] = 0.55; tgt_coords[3* 5+2] = -0.05;
    tgt_coords[3* 6+0] = 0.55; tgt_coords[3* 6+1] = 1.05; tgt_coords[3* 6+2] = -0.05;
    tgt_coords[3* 7+0] = 1.05; tgt_coords[3* 7+1] = 1.05; tgt_coords[3* 7+2] = -0.05;
    tgt_coords[3* 8+0] = 1.55; tgt_coords[3* 8+1] = 1.05; tgt_coords[3* 8+2] = -0.05;
    tgt_coords[3* 9+0] = 0.55; tgt_coords[3* 9+1] = 0.05; tgt_coords[3* 9+2] =  0.45;
    tgt_coords[3*10+0] = 1.05; tgt_coords[3*10+1] = 0.05; tgt_coords[3*10+2] =  0.45;
    tgt_coords[3*11+0] = 1.55; tgt_coords[3*11+1] = 0.05; tgt_coords[3*11+2] =  0.45;
    tgt_coords[3*12+0] = 0.55; tgt_coords[3*12+1] = 0.55; tgt_coords[3*12+2] =  0.45;
    tgt_coords[3*13+0] = 1.05; tgt_coords[3*13+1] = 0.55; tgt_coords[3*13+2] =  0.45;
    tgt_coords[3*14+0] = 1.55; tgt_coords[3*14+1] = 0.55; tgt_coords[3*14+2] =  0.45;
    tgt_coords[3*15+0] = 0.55; tgt_coords[3*15+1] = 1.05; tgt_coords[3*15+2] =  0.45;
    tgt_coords[3*16+0] = 1.05; tgt_coords[3*16+1] = 1.05; tgt_coords[3*16+2] =  0.45;
    tgt_coords[3*17+0] = 1.55; tgt_coords[3*17+1] = 1.05; tgt_coords[3*17+2] =  0.45;

  }
  if (i_rank == 1){
    src_gnum[0] = 5;
    src_gnum[1] = 6;
    src_gnum[2] = 7;
    src_gnum[3] = 8;
    src_coords[3*0+0] = 0.25; src_coords[3*0+1] = 0.25; src_coords[3*0+2] = 0.75;
    src_coords[3*1+0] = 0.75; src_coords[3*1+1] = 0.25; src_coords[3*1+2] = 0.75;
    src_coords[3*2+0] = 0.25; src_coords[3*2+1] = 0.75; src_coords[3*2+2] = 0.75;
    src_coords[3*3+0] = 0.75; src_coords[3*3+1] = 0.75; src_coords[3*3+2] = 0.75;

    for (int i=0; i < 9; i ++)
      tgt_gnum[i] = i+10+9;
    for (int i=9; i < 18; i ++)
      tgt_gnum[i] = i+1;

    tgt_coords[3* 0+0] = 0.55; tgt_coords[3* 0+1] = 0.05; tgt_coords[3* 0+2] = 0.95;
    tgt_coords[3* 1+0] = 1.05; tgt_coords[3* 1+1] = 0.05; tgt_coords[3* 1+2] = 0.95;
    tgt_coords[3* 2+0] = 1.55; tgt_coords[3* 2+1] = 0.05; tgt_coords[3* 2+2] = 0.95;
    tgt_coords[3* 3+0] = 0.55; tgt_coords[3* 3+1] = 0.55; tgt_coords[3* 3+2] = 0.95;
    tgt_coords[3* 4+0] = 1.05; tgt_coords[3* 4+1] = 0.55; tgt_coords[3* 4+2] = 0.95;
    tgt_coords[3* 5+0] = 1.55; tgt_coords[3* 5+1] = 0.55; tgt_coords[3* 5+2] = 0.95;
    tgt_coords[3* 6+0] = 0.55; tgt_coords[3* 6+1] = 1.05; tgt_coords[3* 6+2] = 0.95;
    tgt_coords[3* 7+0] = 1.05; tgt_coords[3* 7+1] = 1.05; tgt_coords[3* 7+2] = 0.95;
    tgt_coords[3* 8+0] = 1.55; tgt_coords[3* 8+1] = 1.05; tgt_coords[3* 8+2] = 0.95;
    tgt_coords[3* 9+0] = 0.55; tgt_coords[3* 9+1] = 0.05; tgt_coords[3* 9+2] = 0.45;
    tgt_coords[3*10+0] = 1.05; tgt_coords[3*10+1] = 0.05; tgt_coords[3*10+2] = 0.45;
    tgt_coords[3*11+0] = 1.55; tgt_coords[3*11+1] = 0.05; tgt_coords[3*11+2] = 0.45;
    tgt_coords[3*12+0] = 0.55; tgt_coords[3*12+1] = 0.55; tgt_coords[3*12+2] = 0.45;
    tgt_coords[3*13+0] = 1.05; tgt_coords[3*13+1] = 0.55; tgt_coords[3*13+2] = 0.45;
    tgt_coords[3*14+0] = 1.55; tgt_coords[3*14+1] = 0.55; tgt_coords[3*14+2] = 0.45;
    tgt_coords[3*15+0] = 0.55; tgt_coords[3*15+1] = 1.05; tgt_coords[3*15+2] = 0.45;
    tgt_coords[3*16+0] = 1.05; tgt_coords[3*16+1] = 1.05; tgt_coords[3*16+2] = 0.45;
    tgt_coords[3*17+0] = 1.55; tgt_coords[3*17+1] = 1.05; tgt_coords[3*17+2] = 0.45;
  }

  int all_src_gnum[8] = {1,2,3,4,5,6,7,8};
  double all_src_coords[24];
  all_src_coords[3*0+0] = 0.25; all_src_coords[3*0+1] = 0.25; all_src_coords[3*0+2] = 0.25;
  all_src_coords[3*1+0] = 0.75; all_src_coords[3*1+1] = 0.25; all_src_coords[3*1+2] = 0.25;
  all_src_coords[3*2+0] = 0.75; all_src_coords[3*2+1] = 0.75; all_src_coords[3*2+2] = 0.25;
  all_src_coords[3*3+0] = 0.25; all_src_coords[3*3+1] = 0.75; all_src_coords[3*3+2] = 0.25;
  all_src_coords[3*4+0] = 0.25; all_src_coords[3*4+1] = 0.25; all_src_coords[3*4+2] = 0.75;
  all_src_coords[3*5+0] = 0.75; all_src_coords[3*5+1] = 0.25; all_src_coords[3*5+2] = 0.75;
  all_src_coords[3*6+0] = 0.25; all_src_coords[3*6+1] = 0.75; all_src_coords[3*6+2] = 0.75;
  all_src_coords[3*7+0] = 0.75; all_src_coords[3*7+1] = 0.75; all_src_coords[3*7+2] = 0.75;

  //Brute force expected results
  PDM_g_num_t *expected_closest_gnum = (PDM_g_num_t *) malloc(n_tgt*sizeof(PDM_g_num_t));
  double      *expected_closest_dist = (double *) malloc(n_tgt*sizeof(double));
  for (int i = 0; i < n_tgt; i++) {
    PDM_g_num_t arg_min = -1;
    double min_dist = 1E12;
    for (int j = 0; j < 8; j++) {
      double dist = (tgt_coords[3*i+0] - all_src_coords[3*j+0])*(tgt_coords[3*i+0] - all_src_coords[3*j+0])
                  + (tgt_coords[3*i+1] - all_src_coords[3*j+1])*(tgt_coords[3*i+1] - all_src_coords[3*j+1])
                  + (tgt_coords[3*i+2] - all_src_coords[3*j+2])*(tgt_coords[3*i+2] - all_src_coords[3*j+2]);
      if (dist < min_dist) {
        arg_min = all_src_gnum[j];
        min_dist = dist;
      }
    }
    expected_closest_gnum[i] = arg_min;
    expected_closest_dist[i] = min_dist;
  }


  PDM_closest_point_t* clsp = PDM_closest_points_create (PDM_MPI_COMM_WORLD,
                                                         1,
                                                         PDM_OWNERSHIP_USER);

  PDM_closest_points_n_part_cloud_set (clsp, 1, 1);

  PDM_closest_points_src_cloud_set (clsp,
                                    0,
                                    n_src,
                                    src_coords,
                                    src_gnum);

  PDM_closest_points_tgt_cloud_set (clsp,
                                    0,
                                    n_tgt,
                                    tgt_coords,
                                    tgt_gnum);


  PDM_closest_points_compute (clsp);



  PDM_g_num_t *closest_src_gnum = NULL;
  double      *closest_src_dist = NULL;

  PDM_closest_points_get (clsp,
                          0,
                          &closest_src_gnum,
                          &closest_src_dist);


  PDM_closest_points_free (clsp);


  printf("Results\n");
  if (i_rank == 1)
    sleep(1);
  for (int i = 0; i < n_tgt; i++) {
    printf("[%d] Tgt "PDM_FMT_G_NUM" : closest is "PDM_FMT_G_NUM" with dist %f -- expected : "PDM_FMT_G_NUM" with dist %f\n", i_rank, tgt_gnum[i], closest_src_gnum[i], closest_src_dist[i], expected_closest_gnum[i], expected_closest_dist[i]);
  }


  for (int i = 0; i < n_tgt; i++) {
    /*assert (PDM_ABS(closest_src_dist[i] - expected_closest_dist[i]) < 1E-6);*/
    assert (closest_src_gnum[i] == expected_closest_gnum[i]);
  }


  free(expected_closest_gnum);
  free(expected_closest_dist);
  free(closest_src_gnum);
  free(closest_src_dist);

  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);
  PDM_MPI_Finalize ();

  return 0;
}
