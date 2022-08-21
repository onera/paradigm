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
  PDM_MPI_Comm comm_shared;
  PDM_MPI_Comm_split_type(comm, PDM_MPI_SPLIT_NUMA, &comm_shared);

  int n_rank_in_node, i_rank_in_node;
  PDM_MPI_Comm_rank (comm_shared, &i_rank_in_node);
  PDM_MPI_Comm_size (comm_shared, &n_rank_in_node);

  PDM_MPI_Comm comm_same_numa_id;
  PDM_MPI_Comm_split(comm, i_rank_in_node /* color*/, i_rank /* key */, &comm_same_numa_id);

  int n_rank_same_numa_id, i_rank_same_numa_id;
  PDM_MPI_Comm_rank (comm_same_numa_id, &i_rank_same_numa_id);
  PDM_MPI_Comm_size (comm_same_numa_id, &n_rank_same_numa_id);

  if(1 == 1) {
    log_trace("i_rank = %i [%i] | i_rank_in_node = %i [%i] | i_rank_same_numa_id = %i [%i] \n",
              i_rank, n_rank, i_rank_in_node, n_rank_in_node, i_rank_same_numa_id, n_rank_same_numa_id);

  }

  /*
   * Prepare to save rank_id in shared memory
   */
  PDM_mpi_win_shared_t* wsame_numa_id_core_idx = PDM_mpi_win_shared_create(n_rank_in_node+1, sizeof(int), comm_shared);
  int *same_numa_id_core_idx = PDM_mpi_win_shared_get(wsame_numa_id_core_idx);
  PDM_mpi_win_shared_lock_all (0, wsame_numa_id_core_idx);

  same_numa_id_core_idx[i_rank_in_node+1] = n_rank_same_numa_id;
  PDM_mpi_win_shared_sync(wsame_numa_id_core_idx);
  PDM_MPI_Barrier(comm_shared);

  /*
   * Create index
   */
  if(i_rank_in_node == 0) {
    same_numa_id_core_idx[0] = 0;
    for(int i = 0; i < n_rank_in_node; ++i) {
      same_numa_id_core_idx[i+1] += same_numa_id_core_idx[i];
    }
  }
  PDM_mpi_win_shared_sync(wsame_numa_id_core_idx);
  PDM_MPI_Barrier(comm_shared);

  PDM_log_trace_array_int(same_numa_id_core_idx, n_rank_in_node+1, "wsame_numa_id_core_idx ::");

  /*
   * Get rank number
   */
  PDM_mpi_win_shared_t* wsame_numa_id_rank_gid = PDM_mpi_win_shared_create(same_numa_id_core_idx[n_rank_in_node], sizeof(int), comm_shared);
  int *rank_gid = PDM_mpi_win_shared_get(wsame_numa_id_rank_gid);
  PDM_mpi_win_shared_lock_all (0, wsame_numa_id_rank_gid);

  int *lrank_id = &rank_gid[same_numa_id_core_idx[i_rank_in_node]];

  PDM_MPI_Allgather(&i_rank, 1, PDM_MPI_INT, lrank_id, 1, PDM_MPI_INT, comm_same_numa_id);
  PDM_mpi_win_shared_sync(wsame_numa_id_rank_gid);
  PDM_MPI_Barrier(comm_shared);

  if(1 == 1) {
    PDM_log_trace_array_int(rank_gid, same_numa_id_core_idx[n_rank_in_node], "rank_gid ::");
  }

  /*
   * Il faudrait arriver a ne pas trier lesnieghboor dans le comm shared
   */
  srand(i_rank+11);
  int n_val = (rand() % 8)+1;
  log_trace("n_val = %i \n", n_val);

  int *val = malloc(n_val * sizeof(int));
  for(int i = 0; i < n_val; ++i) {
    val[i] = i_rank;
  }


  /*
   * Build shared array of number of numa
   */
  PDM_mpi_win_shared_t* wval_tot_n = PDM_mpi_win_shared_create(same_numa_id_core_idx[n_rank_in_node], sizeof(int), comm_shared);
  int *val_tot_n  = PDM_mpi_win_shared_get(wval_tot_n);
  int *lval_tot_n = &val_tot_n[same_numa_id_core_idx[i_rank_in_node]];
  PDM_mpi_win_shared_lock_all (0, wval_tot_n);


  PDM_MPI_Allgather(&n_val    , 1, PDM_MPI_INT,
                    lval_tot_n, 1, PDM_MPI_INT, comm_same_numa_id);
  PDM_MPI_Barrier(comm_shared);

  if(1 == 1) {
    PDM_log_trace_array_int(val_tot_n, same_numa_id_core_idx[n_rank_in_node], "val_tot_n ::");
  }

  /*
   * Ce qu'on souhaite c'est  obtenir le même resultat que sur le comm classique (même ordre que avec l'échange via comm_same_numa_id)
   */
  int *recv_count = malloc( n_rank_same_numa_id    * sizeof(int));
  int *recv_shift = malloc((n_rank_same_numa_id+1) * sizeof(int));
  recv_shift[0] = 0;
  for(int i = 0; i < n_rank_same_numa_id; ++i) {
    recv_count[i  ] = lval_tot_n[i];
    recv_shift[i+1] = recv_shift[i] + lval_tot_n[i];
  }

  int *recv_val = malloc(recv_shift[n_rank_same_numa_id] * sizeof(int));


  PDM_MPI_Allgatherv(val     , n_val,                  PDM_MPI_INT,
                     recv_val, recv_count, recv_shift, PDM_MPI_INT, comm_same_numa_id);


  PDM_log_trace_array_int(recv_val, recv_shift[n_rank_same_numa_id], "recv_val ::");

  free(recv_val);
  free(recv_count);
  free(recv_shift);

  PDM_mpi_win_shared_unlock_all(wsame_numa_id_core_idx);
  PDM_mpi_win_shared_unlock_all(wsame_numa_id_rank_gid);
  PDM_mpi_win_shared_unlock_all(wval_tot_n);
  PDM_mpi_win_shared_free(wsame_numa_id_core_idx);
  PDM_mpi_win_shared_free(wsame_numa_id_rank_gid);
  PDM_mpi_win_shared_free(wval_tot_n);

  free(val);

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
  PDM_MPI_Comm comm_shared;
  PDM_MPI_Comm_split_type(comm, PDM_MPI_SPLIT_NUMA, &comm_shared);

  int n_rank_in_node, i_rank_in_node;
  PDM_MPI_Comm_rank (comm_shared, &i_rank_in_node);
  PDM_MPI_Comm_size (comm_shared, &n_rank_in_node);

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
