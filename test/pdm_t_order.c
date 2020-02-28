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
#include "pdm_order.h"
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
   * Triplet : (iproc, ipart, index)
   *
   * ----------------------------|           |----------------------------
   *     |                       |           |               |
   *     |       (0,0,2)         |           |               |
   *     |                       |           |  (1, 0, 1)    |
   * ----------------------------|           |               |
   *     |                       |           |               |
   *     |      (0,0,1)          |           |----------------------------
   *     |                       |           |               |
   * ----------------------------|           |               |
   *     |                       |           |  (1, 0, 0)    |
   *     |      (0,0,0)          |           |               |
   *     |                       |           |               |
   * ----------------------------|           |----------------------------
   *                         join_1        join_2
   */

  const int stride  = 3;

  /* Connection de join1 avec join2 */
  // const int n_faces_j1 = 3;
  int connect_triplet_j1[12] = {// Fisrt
                                1, 0, 0,
                                // Second
                                1, 0, 0,
                                1, 0, 1,
                                // Third
                                1, 0, 1};

  /* Connection de join2 avec join1 */
  // const int n_faces_j2 = 3;
  int connect_triplet_j2[12] = {// First
                                0, 0, 0,
                                0, 0, 1,
                                // Second
                                0, 0, 1,
                                0, 0, 2};


  // Ordering
  int nb_ent = 4;
  int order_j1[nb_ent];
  PDM_order_lnum_s(connect_triplet_j1, stride, order_j1, nb_ent);

  printf("order_j1:: ");
  for(int i = 0; i < nb_ent; i++){
    printf("%d ", order_j1[i]);
  }
  printf("\n");

  int order_j2[nb_ent];
  PDM_order_lnum_s(connect_triplet_j2, stride, order_j2, nb_ent);

  printf("order_j2:: ");
  for(int i = 0; i < nb_ent; i++){
    printf("%d ", order_j2[i]);
  }
  printf("\n");

  PDM_printf ("\nfin Test\n");

  return 0;

}
