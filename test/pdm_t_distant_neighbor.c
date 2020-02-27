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
#include "pdm_distant_neighbor.h"
#include "pdm_points_merge.h"
#include "pdm_printf.h"
#include "pdm_error.h"

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
  int iRank;
  int nRank;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &iRank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &nRank);

  /*
   * The test case represent two conform interface between 2 partitions, we need to find out connection of each faces
   *    I/  Localisation in the frame of the current/opposite boundary condition
   *    II/ Use distant_neighbor to find out the face number in the partition
   */
  double xyz_j1[45] = {0.5,  2.5,  4. ,  0.5,  1.5,  4. ,  1.5,  2.5,  4. ,  1.5,  1.5,
                       4. ,  2.5,  1.5,  4. ,  3.5,  1.5,  4. ,  1.5,  0.5,  4. ,  2.5,
                       0.5,  4. ,  2.5,  2.5,  4. ,  3.5,  0.5,  4. ,  3.5,  2.5,  4. ,
                       4.5,  0.5,  4. ,  4.5,  1.5,  4. ,  0.5,  0.5,  4. ,  4.5,  2.5,
                       4. };

  double cln_j1[15] = {1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1., 1., 1.};

  double xyz_j2[45] = {0.5,  1.5,  4. ,  1.5,  1.5,  4. ,  1.5,  2.5,  4. ,  2.5,  1.5,
                       4. ,  2.5,  2.5,  4. ,  3.5,  1.5,  4. ,  3.5,  2.5,  4. ,  4.5,
                       1.5,  4. ,  0.5,  0.5,  4. ,  1.5,  0.5,  4. ,  2.5,  0.5,  4. ,
                       3.5,  0.5,  4. ,  4.5,  0.5,  4. ,  4.5,  2.5,  4. ,  0.5,  2.5,
                       4.};

  double cln_j2[15] = {1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1., 1., 1.};

  /*
   * Define i_cloud
   */
  int      n_cloud  = -1;
  int      n_points = 15;
  double** coords;
  double** char_lenght;
  if(nRank == 1){
    n_cloud = 2;
    coords      = (double **) malloc( n_cloud * sizeof(double *));
    char_lenght = (double **) malloc( n_cloud * sizeof(double *));

    coords     [0] = xyz_j1;
    coords     [1] = xyz_j2;
    char_lenght[0] = cln_j1;
    char_lenght[1] = cln_j2;

  } else if ( nRank == 2){
    n_cloud = 1;
    coords      = (double **) malloc( n_cloud * sizeof(double **));
    char_lenght = (double **) malloc( n_cloud * sizeof(double **));
    if(iRank == 0){
      coords     [0] = xyz_j1;
      char_lenght[0] = cln_j1;
    } else if (iRank == 1){
      coords     [0] = xyz_j2;
      char_lenght[0] = cln_j2;
    }

  } else {
    PDM_error(__FILE__, __LINE__, 0, "pdm_t_distant_neighbor error : Bad number of process for test cases \n");
  }


  /*
   * Call pdm_points_merge to find match between the two block
   */
  double tolerance = 0.1;
  int pm_id = PDM_points_merge_create(n_cloud, tolerance, PDM_MPI_COMM_WORLD);

  if(nRank == 1){
    PDM_points_merge_cloud_set(pm_id, 0, n_points, coords[0], char_lenght[0]);
    PDM_points_merge_cloud_set(pm_id, 1, n_points, coords[1], char_lenght[1]);
  } else if(nRank == 1){
    if(iRank == 0){
      PDM_points_merge_cloud_set(pm_id, 0, n_points, coords[0], char_lenght[0]);
    } else if(iRank == 1){
      PDM_points_merge_cloud_set(pm_id, 1, n_points, coords[1], char_lenght[1]);
    }
  }

  PDM_points_merge_process(pm_id);

  /*
   * Get resulting points_merge
   */
  int **candidates_idx  = (int **) malloc( n_cloud * sizeof(int**));
  int **candidates_desc = (int **) malloc( n_cloud * sizeof(int**));
  for(int i_cloud = 0; i_cloud < n_cloud; i_cloud++){
    int n_cloud_points = -1;
    int n_desc   = -1;
    PDM_points_merge_candidates_get(pm_id, i_cloud, &candidates_idx[i_cloud], &candidates_desc[i_cloud]);
    PDM_points_merge_candidates_size_get(pm_id, i_cloud, &n_cloud_points, &n_desc);

    assert(n_desc == candidates_idx[i_cloud][n_cloud_points]);
    if(1 == 1){
      printf("-- n_desc:: %d \n ", n_desc);
      printf("-- candidates_idx[i_cloud][n_cloud_points+1]:: %d \n ", candidates_idx[i_cloud][n_cloud_points]);
      for(int i = 0; i < n_cloud_points; i++){
        printf("-- %d %d ", i_cloud, i);
        for (int j = candidates_idx[i_cloud][i];
                 j <  candidates_idx[i_cloud][i+1]; j++) {
          printf(" : %d", candidates_desc[i_cloud][3*j]);
          printf(" %d", candidates_desc[i_cloud][3*j+1]);
          printf(" %d", candidates_desc[i_cloud][3*j+2]);
        }
        printf("\n");
      }
    }
  }

  /*
   * Free
   */
  PDM_points_merge_free(pm_id);
  free(coords);
  free(char_lenght);
  free(candidates_idx);
  free(candidates_desc);
  PDM_MPI_Finalize ();

  PDM_printf ("\nfin Test\n");

  return 0;

}


  // for(int i_cloud = 0; i_cloud < n_cloud; i_cloud++){
  //   free(coords[i_cloud]);
  //   free(char_lenght[i_cloud]);
  // }
  // for(int i_cloud = 0; i_cloud < n_cloud; i_cloud++){
  //   coords     [i_cloud] = (double *) malloc( 3 * n_points * sizeof(double *));
  //   char_lenght[i_cloud] = (double *) malloc(     n_points * sizeof(double *));
  // }


  // printf("distrib_index : ");
  // for (int i = 0; i < nRank + 1; i++) {
  //   printf(PDM_FMT_G_NUM" ", distrib_index[i]);
  // }
  // printf("\n");

  // PDM_MPI_Allgather (&weight_sum, 1, PDM_MPI_DOUBLE,
  //                    weights_sum_procs, 1, PDM_MPI_DOUBLE,
  //                    PDM_MPI_COMM_WORLD);

  // printf("weights procs :");
  // for (int i = 0; i < nRank; i++) {
  //   printf(" %12.5e", weights_sum_procs[i]);
  // }
  // printf("\n");
