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

  int n_rank;
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  int n_elt_proc = 10;

  PDM_g_num_t *numabs  = NULL;
  double      *weights = NULL;
  int         *stride  = NULL;
  PDM_malloc(numabs , n_elt_proc, PDM_g_num_t );
  PDM_malloc(weights, n_elt_proc, double      );
  PDM_malloc(stride , n_elt_proc, int         );

  for (int i = 0; i < n_elt_proc; i++) {
    numabs[i] = i_rank * n_elt_proc + i + 1;
    weights[i] = i_rank+1;
    stride[i] = 1;
  }

  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_MERGE,
                                                       1.,
                                                       &numabs,
                                                       &weights,
                                                       &n_elt_proc,
                                                       1,
                                                       PDM_MPI_COMM_WORLD);


  PDM_g_num_t *distrib_index = PDM_part_to_block_distrib_index_get (ptb);

  double *block_weights;
  int *block_stride;

  PDM_part_to_block_exch (ptb,
                          sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          &stride,
                          (void **) &weights,
                           &block_stride,
                          (void **) &block_weights);

  int n_elt_block = PDM_part_to_block_n_elt_block_get (ptb);
  double weight_sum = 0.;

  for (int i = 0; i < n_elt_block; i++) {
    weight_sum += block_weights[i];
  }

  double *weights_sum_procs;
  PDM_malloc(weights_sum_procs, n_rank, double);

  printf("distrib_index : ");
  for (int i = 0; i < n_rank + 1; i++) {
    printf(PDM_FMT_G_NUM" ", distrib_index[i]);
  }
  printf("\n");

  PDM_MPI_Allgather (&weight_sum, 1, PDM_MPI_DOUBLE,
                     weights_sum_procs, 1, PDM_MPI_DOUBLE,
                     PDM_MPI_COMM_WORLD);

  printf("weights procs :");
  for (int i = 0; i < n_rank; i++) {
    printf(" %12.5e", weights_sum_procs[i]);
  }
  printf("\n");

  PDM_part_to_block_free (ptb);

  PDM_free(numabs);
  PDM_free(weights);
  PDM_free(weights_sum_procs);
  PDM_free(stride);
  PDM_free(block_stride);
  PDM_free(block_weights);

  if (i_rank == 0) {
    PDM_printf("-- End\n");
  }
  PDM_MPI_Finalize ();

  return 0;

}
