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
#include "pdm_gnum.h"

#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_error.h"


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
     "  -n      <level>  Global array size.\n\n"
     "  -l      <level>  Number of nuplet .\n\n"
     "  -h               This message.\n\n");


  exit (exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc       Number of arguments
 * \param [in]      argv       Arguments
 * \param [inout]   n_g_elmts  Number of vertices on the cube side
 * \param [inout]   nuplet     Cube length
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 PDM_g_num_t   *n_g_elmts,
 int           *nuplet
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp (argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp (argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_g_elmts = atol (argv[i]);
        *n_g_elmts = (PDM_g_num_t) _n_g_elmts;
      }
    }
    else if (strcmp (argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else
        *nuplet = atoi (argv[i]);
    }
    else
      _usage (EXIT_FAILURE);
    i++;
  }
}

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
  int           i_rank;
  int           n_rank;

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  /*
   *  Set default values
   */

  PDM_g_num_t   n_g_elmts = 40;
  int           nuplet    = 2;


  /*
   *  Read args
   */
  _read_args (argc,
              argv,
              &n_g_elmts,
              &nuplet);


  /* Generate distribution */
  PDM_g_num_t* distrib = PDM_compute_uniform_entity_distribution(comm, n_g_elmts);
  int dn_elmts = distrib[i_rank+1] - distrib[i_rank];

  PDM_g_num_t* elmts_ln_to_gn = (PDM_g_num_t *) malloc(nuplet * dn_elmts * sizeof(PDM_g_num_t));

  int seed = 0;
  for(int i = 0; i < dn_elmts; ++i) {

    unsigned int _seed = (unsigned int) (distrib[i_rank] + i) + 1;
    srand(_seed + seed);

    for(int k = 0; k < nuplet; ++k) {
      elmts_ln_to_gn[nuplet*i+k] = (PDM_g_num_t ) rand();
    }
  }

  free(distrib);
  free(elmts_ln_to_gn);
  PDM_MPI_Finalize ();

  return 0;
}
