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
#include "pdm_gnum_from_hash_values.h"
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

  int numProcs;
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);

  int n_part = 1;
  PDM_bool_t equilibrate = PDM_FALSE;

  int n_elt_proc = 10;

  int gnum_fhv_id = PDM_gnum_from_hash_values_create(n_part, equilibrate, PDM_MPI_COMM_WORLD);

  printf("gnum_fhv_id:: %d \n", gnum_fhv_id);


  /*
   * Free
   */
  PDM_gnum_from_hv_free(gnum_fhv_id, 0);

  PDM_MPI_Finalize ();

  PDM_printf ("\nfin Test\n");

  return 0;

}
