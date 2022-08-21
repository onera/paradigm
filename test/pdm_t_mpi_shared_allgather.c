#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <pdm_mpi.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_part.h"

#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

static
PDM_MPI_Comm setup_numa_graph(PDM_MPI_Comm comm)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // Shared
  PDM_MPI_Comm comm_node;
  PDM_MPI_Comm_split_type(comm, PDM_MPI_SPLIT_NUMA, &comm_node);

  int n_rank_in_node, i_rank_in_node;
  PDM_MPI_Comm_rank (comm_node, &i_rank_in_node);
  PDM_MPI_Comm_size (comm_node, &n_rank_in_node);

  PDM_MPI_Comm comm_same_numa;
  PDM_MPI_Comm_split(comm, i_rank_in_node /* color*/, i_rank /* key */, &comm_same_numa);

  int n_rank_same_numa, i_rank_same_numa;
  PDM_MPI_Comm_rank (comm_same_numa, &i_rank_same_numa);
  PDM_MPI_Comm_size (comm_same_numa, &n_rank_same_numa);

  if(1 == 1) {

    log_trace("i_rank = %i [%i] | i_rank_in_node = %i [%i] | i_rank_same_numa = %i [%i] \n",
              i_rank, n_rank, i_rank_in_node, n_rank_in_node, i_rank_same_numa, n_rank_same_numa);

  }



  PDM_MPI_Comm comm_dist_graph;

  return comm_dist_graph;
}



/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  /*
   *  Init
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // Shared
  PDM_MPI_Comm comm_node;
  PDM_MPI_Comm_split_type(comm, PDM_MPI_SPLIT_NUMA, &comm_node);

  int n_rank_in_node, i_rank_in_node;
  PDM_MPI_Comm_rank (comm_node, &i_rank_in_node);
  PDM_MPI_Comm_size (comm_node, &n_rank_in_node);

  /*
   * Setup graph
   */
  PDM_MPI_Comm comm_dist_graph = setup_numa_graph(comm);


  PDM_MPI_Barrier (comm);
  if (i_rank == 0) {
  PDM_printf ("-- End\n");
  }

  PDM_MPI_Finalize ();

  return 0;
}
