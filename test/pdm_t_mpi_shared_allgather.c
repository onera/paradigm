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
void
hybrid_exchange_numa(PDM_MPI_Comm comm)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_MPI_Comm comm_shared;
  PDM_MPI_Comm comm_dist_graph;

  int  n_degree_in = 0;
  int *neighbor_in = NULL;

  PDM_MPI_setup_hybrid_dist_comm_graph(comm,
                                       &comm_shared,
                                       &comm_dist_graph,
                                       &n_degree_in,
                                       &neighbor_in);

  int n_rank_in_shm, i_rank_in_shm;
  PDM_MPI_Comm_rank (comm_shared, &i_rank_in_shm);
  PDM_MPI_Comm_size (comm_shared, &n_rank_in_shm);

  /*
   * Allgather to have size + displacement
   */
  srand(i_rank+11);
  int n_val = (rand() % 8)+1;

  int *n_vals_out;
  PDM_malloc(n_vals_out, n_degree_in, int);

  PDM_MPI_Neighbor_allgather(&n_val    , 1, PDM_MPI_INT,
                             n_vals_out, 1, PDM_MPI_INT, comm_dist_graph);

  PDM_mpi_win_shared_t* wshared_vals_out_n   = PDM_mpi_win_shared_create(n_rank  , sizeof(int), comm_shared);
  PDM_mpi_win_shared_t* wshared_vals_out_idx = PDM_mpi_win_shared_create(n_rank+1, sizeof(int), comm_shared);
  int *shared_vals_out_n   = PDM_mpi_win_shared_get(wshared_vals_out_n);
  int *shared_vals_out_idx = PDM_mpi_win_shared_get(wshared_vals_out_idx);
  PDM_mpi_win_shared_lock_all (0, wshared_vals_out_n  );
  PDM_mpi_win_shared_lock_all (0, wshared_vals_out_idx);

  for(int i = 0; i < n_degree_in; ++i) {
    shared_vals_out_n[neighbor_in[i]] = n_vals_out[i];
  }
  PDM_MPI_Barrier(comm_shared);

  /*
   * Tentative allgatherv
   */
  if(i_rank_in_shm == 0) {
    shared_vals_out_idx[0] = 0;
    for(int i = 0; i < n_rank; ++i) {
      shared_vals_out_idx[i+1] = shared_vals_out_idx[i] + shared_vals_out_n[i];
    }
  }
  PDM_MPI_Barrier(comm_shared);

  if(0 == 1) {
    PDM_log_trace_array_int(shared_vals_out_n  , n_rank , "shared_vals_out_n   ::");
    PDM_log_trace_array_int(shared_vals_out_idx, n_rank , "shared_vals_out_idx ::");
  }

  // Hook local recv_shift
  int *recv_shift;
  PDM_malloc(recv_shift, n_degree_in, int);
  for(int i = 0; i < n_degree_in; ++i) {
    recv_shift[i] = shared_vals_out_idx[neighbor_in[i]];
  }

  int *val;
  PDM_malloc(val, n_val, int);
  for(int i = 0; i < n_val; ++i) {
    val[i] = i_rank;
  }

  PDM_mpi_win_shared_t* wshared_vals_out   = PDM_mpi_win_shared_create(shared_vals_out_idx[n_rank]  , sizeof(int), comm_shared);
  int *shared_vals_out   = PDM_mpi_win_shared_get(wshared_vals_out);
  PDM_mpi_win_shared_lock_all (0, wshared_vals_out  );

  PDM_MPI_Neighbor_allgatherv(val            , n_val     , PDM_MPI_INT,
                              shared_vals_out, n_vals_out, recv_shift, PDM_MPI_INT, comm_dist_graph);

  PDM_MPI_Barrier(comm_shared);

  if(0 == 1) {
    PDM_log_trace_array_int(shared_vals_out, shared_vals_out_idx[n_rank], "shared_vals_out ::");
  }

  PDM_mpi_win_shared_unlock_all(wshared_vals_out_n);
  PDM_mpi_win_shared_unlock_all(wshared_vals_out_idx);
  PDM_mpi_win_shared_unlock_all(wshared_vals_out);
  PDM_mpi_win_shared_free(wshared_vals_out_n);
  PDM_mpi_win_shared_free(wshared_vals_out_idx);
  PDM_mpi_win_shared_free(wshared_vals_out);

  PDM_free(n_vals_out);
  PDM_free(recv_shift);
  PDM_free(val);
  PDM_free(neighbor_in);
  PDM_MPI_Comm_free(&comm_shared);
  PDM_MPI_Comm_free(&comm_dist_graph);

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

  int n_rank_in_shm, i_rank_in_shm;
  PDM_MPI_Comm_rank (comm_shared, &i_rank_in_shm);
  PDM_MPI_Comm_size (comm_shared, &n_rank_in_shm);

  /*
   * Setup graph
   */
  hybrid_exchange_numa(comm);


  PDM_MPI_Comm_free(&comm_shared);
  PDM_MPI_Barrier (comm);
  if (i_rank == 0) {
  PDM_printf ("-- End\n");
  }

  PDM_MPI_Finalize ();

  return 0;
}
