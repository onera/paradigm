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

  int n_rank_in_shm, i_rank_in_shm;
  PDM_MPI_Comm_rank (comm_shared, &i_rank_in_shm);
  PDM_MPI_Comm_size (comm_shared, &n_rank_in_shm);

  PDM_MPI_Comm comm_same_numa_id;
  PDM_MPI_Comm_split(comm, i_rank_in_shm /* color*/, i_rank /* key */, &comm_same_numa_id);

  int n_rank_same_numa_id, i_rank_same_numa_id;
  PDM_MPI_Comm_rank (comm_same_numa_id, &i_rank_same_numa_id);
  PDM_MPI_Comm_size (comm_same_numa_id, &n_rank_same_numa_id);

  if(1 == 1) {
    log_trace("i_rank = %i [%i] | i_rank_in_shm = %i [%i] | i_rank_same_numa_id = %i [%i] \n",
              i_rank, n_rank, i_rank_in_shm, n_rank_in_shm, i_rank_same_numa_id, n_rank_same_numa_id);

  }

  /*
   * Prepare to save rank_id in shared memory
   */
  PDM_mpi_win_shared_t* wsame_numa_id_core_idx = PDM_mpi_win_shared_create(n_rank_in_shm+1, sizeof(int), comm_shared);
  int *same_numa_id_core_idx = PDM_mpi_win_shared_get(wsame_numa_id_core_idx);
  PDM_mpi_win_shared_lock_all (0, wsame_numa_id_core_idx);

  same_numa_id_core_idx[i_rank_in_shm+1] = n_rank_same_numa_id;
  PDM_mpi_win_shared_sync(wsame_numa_id_core_idx);
  PDM_MPI_Barrier(comm_shared);

  /*
   * Create index
   */
  if(i_rank_in_shm == 0) {
    same_numa_id_core_idx[0] = 0;
    for(int i = 0; i < n_rank_in_shm; ++i) {
      same_numa_id_core_idx[i+1] += same_numa_id_core_idx[i];
    }
  }
  PDM_mpi_win_shared_sync(wsame_numa_id_core_idx);
  PDM_MPI_Barrier(comm_shared);

  PDM_log_trace_array_int(same_numa_id_core_idx, n_rank_in_shm+1, "wsame_numa_id_core_idx ::");

  /*
   * Get rank number
   */
  PDM_mpi_win_shared_t* wsame_numa_id_rank_gid = PDM_mpi_win_shared_create(same_numa_id_core_idx[n_rank_in_shm], sizeof(int), comm_shared);
  int *rank_gid = PDM_mpi_win_shared_get(wsame_numa_id_rank_gid);
  PDM_mpi_win_shared_lock_all (0, wsame_numa_id_rank_gid);

  int *lrank_id = &rank_gid[same_numa_id_core_idx[i_rank_in_shm]];

  PDM_MPI_Allgather(&i_rank, 1, PDM_MPI_INT, lrank_id, 1, PDM_MPI_INT, comm_same_numa_id);
  PDM_mpi_win_shared_sync(wsame_numa_id_rank_gid);
  PDM_MPI_Barrier(comm_shared);

  if(1 == 1) {
    PDM_log_trace_array_int(rank_gid, same_numa_id_core_idx[n_rank_in_shm], "rank_gid ::");
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
  PDM_mpi_win_shared_t* wval_tot_n = PDM_mpi_win_shared_create(same_numa_id_core_idx[n_rank_in_shm], sizeof(int), comm_shared);
  int *val_tot_n  = PDM_mpi_win_shared_get(wval_tot_n);
  int *lval_tot_n = &val_tot_n[same_numa_id_core_idx[i_rank_in_shm]];
  PDM_mpi_win_shared_lock_all (0, wval_tot_n);


  PDM_MPI_Allgather(&n_val    , 1, PDM_MPI_INT,
                    lval_tot_n, 1, PDM_MPI_INT, comm_same_numa_id);
  PDM_MPI_Barrier(comm_shared);

  if(1 == 1) {
    PDM_log_trace_array_int(val_tot_n, same_numa_id_core_idx[n_rank_in_shm], "val_tot_n ::");
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

static
PDM_MPI_Comm setup_numa_graph2(PDM_MPI_Comm comm)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // Shared
  PDM_MPI_Comm comm_shared;
  PDM_MPI_Comm_split_type(comm, PDM_MPI_SPLIT_NUMA, &comm_shared);

  int n_rank_in_shm, i_rank_in_shm;
  PDM_MPI_Comm_rank (comm_shared, &i_rank_in_shm);
  PDM_MPI_Comm_size (comm_shared, &n_rank_in_shm);

  PDM_MPI_Comm comm_master_of_shm = PDM_MPI_get_group_of_master(comm, comm_shared);

  int i_rank_master_of_shm = -1;
  int n_rank_master_of_shm;
  if(comm_master_of_shm != PDM_MPI_COMM_NULL) {
    PDM_MPI_Comm_rank(comm_master_of_shm, &i_rank_master_of_shm);
    PDM_MPI_Comm_size(comm_master_of_shm, &n_rank_master_of_shm);
  }
  PDM_MPI_Bcast(&n_rank_master_of_shm, 1, PDM_MPI_INT, 0, comm_shared);
  PDM_MPI_Bcast(&i_rank_master_of_shm, 1, PDM_MPI_INT, 0, comm_shared);

  PDM_mpi_win_shared_t* wnuma_by_numa_n = PDM_mpi_win_shared_create(n_rank_master_of_shm, sizeof(int), comm_shared);
  int *numa_by_numa_n  = PDM_mpi_win_shared_get(wnuma_by_numa_n);
  PDM_mpi_win_shared_lock_all (0, wnuma_by_numa_n);


  if(comm_master_of_shm != PDM_MPI_COMM_NULL) {
    PDM_MPI_Allgather(&n_rank_in_shm, 1, PDM_MPI_INT,
                      numa_by_numa_n, 1, PDM_MPI_INT, comm_master_of_shm);
  }
  PDM_mpi_win_shared_sync(wnuma_by_numa_n);
  PDM_MPI_Barrier(comm_shared);

  int n_tot_numa = 0;
  for(int i = 0; i < n_rank_master_of_shm; ++i) {
    n_tot_numa += numa_by_numa_n[i];
  }

  /*
   * Create idx  and  gid of each numa
   */
  PDM_mpi_win_shared_t* wnuma_core_gid    = PDM_mpi_win_shared_create(n_tot_numa               , sizeof(int), comm_shared);
  PDM_mpi_win_shared_t* wnuma_by_numa_idx = PDM_mpi_win_shared_create(n_rank_master_of_shm+1, sizeof(int), comm_shared);
  int *numa_core_gid    = PDM_mpi_win_shared_get(wnuma_core_gid);
  int *numa_by_numa_idx = PDM_mpi_win_shared_get(wnuma_by_numa_idx);
  PDM_mpi_win_shared_lock_all (0, wnuma_core_gid);
  PDM_mpi_win_shared_lock_all (0, wnuma_by_numa_idx);


  if(comm_master_of_shm != PDM_MPI_COMM_NULL) {
    numa_by_numa_idx[0] = 0;
    for(int i = 0; i < n_rank_master_of_shm; ++i) {
      numa_by_numa_idx[i+1] = numa_by_numa_idx[i] + numa_by_numa_n[i];
    }
  }
  PDM_MPI_Barrier(comm_shared);
  PDM_mpi_win_shared_sync(wnuma_by_numa_idx);

  numa_core_gid[numa_by_numa_idx[i_rank_master_of_shm]+i_rank_in_shm] = i_rank;

  PDM_MPI_Barrier(comm_shared);
  PDM_mpi_win_shared_sync(wnuma_core_gid);

  /*
   *  Exchange of the global numbering of rank for each NUMA
   */
  if(comm_master_of_shm != PDM_MPI_COMM_NULL) {
    int *lnuma_core_gid = malloc(n_rank_in_shm * sizeof(int));
    for(int i = 0; i < n_rank_in_shm; ++i) {
      lnuma_core_gid[i] = numa_core_gid[numa_by_numa_idx[i_rank_master_of_shm]+i];
    }
    PDM_MPI_Allgatherv(lnuma_core_gid, n_rank_in_shm, PDM_MPI_INT,
                       numa_core_gid , numa_by_numa_n, numa_by_numa_idx, PDM_MPI_INT, comm_master_of_shm);
    free(lnuma_core_gid);
  }
  PDM_MPI_Barrier(comm_shared);
  PDM_mpi_win_shared_sync(wnuma_core_gid);

  PDM_log_trace_connectivity_int(numa_by_numa_idx, numa_core_gid, n_rank_master_of_shm, "numa_core_gid :: ");

  /*
   * Computation of degree_in
   */
  int n_degrees_in = 0;
  for(int i = 0; i < n_rank_master_of_shm; ++i) {
    for(int j = numa_by_numa_idx[i]; j < numa_by_numa_idx[i+1]; ++j) {
      int lid_rank = (j - numa_by_numa_idx[i]) % n_rank_in_shm; // Donc numero de numa dans le group
      if(lid_rank == i_rank_in_shm){
        n_degrees_in++;
      }
    }
  }

  int* neighbor_in = malloc( (n_degrees_in ) * sizeof(int));
  n_degrees_in = 0;
  for(int i = 0; i < n_rank_master_of_shm; ++i) {
    for(int j = numa_by_numa_idx[i]; j < numa_by_numa_idx[i+1]; ++j) {
      int gid_rank = numa_core_gid[j];
      int lid_rank = (j - numa_by_numa_idx[i]) % n_rank_in_shm;  // Donc numero de numa dans le group
      if(lid_rank == i_rank_in_shm){
        neighbor_in[n_degrees_in++] = gid_rank;
      }
    }
  }

  /*
   * Computation of degree_out : Il faut que tout le monde recupère une liste plein because all gather
   *   Il faut verifier que tous le monde a envoyé sa donné dans chaque numa
   *   Pour chaque range de numa il faut enlever ceux deja envoyer et rajouter les non envoyées
   */
  int *rank_tag = (int *) malloc(n_rank * sizeof(int));

  int n_degrees_out = 0;
  for(int i = 0; i < n_rank_master_of_shm; ++i) {
    for(int j = numa_by_numa_idx[i]; j < numa_by_numa_idx[i+1]; ++j) {
      int lid_rank = (j - numa_by_numa_idx[i]) % numa_by_numa_n[i]; // Donc numero de numa dans le group
      if(lid_rank == i_rank_in_shm){
        n_degrees_out++;
      }
    }
  }

  int* neighbor_out = malloc( (n_degrees_out ) * sizeof(int));
  n_degrees_out = 0;
  for(int i = 0; i < n_rank_master_of_shm; ++i) {
    for(int j = numa_by_numa_idx[i]; j < numa_by_numa_idx[i+1]; ++j) {
      int gid_rank = numa_core_gid[j];
      int lid_rank = (j - numa_by_numa_idx[i]) % numa_by_numa_n[i];  // Donc numero de numa dans le group
      if(lid_rank == i_rank_in_shm){
        neighbor_out[n_degrees_out++] = gid_rank;
      }
    }
  }


  free(rank_tag);


  // n_degrees_out = n_degrees_in;
  // neighbor_out = neighbor_in;

  if(1 == 1) {
    PDM_log_trace_array_int(neighbor_in , n_degrees_in , "neighbor_in  ::");
    PDM_log_trace_array_int(neighbor_out, n_degrees_out, "neighbor_out ::");
  }


  PDM_MPI_Comm comm_dist_graph;
  PDM_MPI_Dist_graph_create_adjacent(comm,
                                     n_degrees_in,
                                     neighbor_in,
                                     n_degrees_out,
                                     neighbor_out,
                                     0,
                                     &comm_dist_graph);


  PDM_mpi_win_shared_unlock_all(wnuma_by_numa_n);
  PDM_mpi_win_shared_unlock_all(wnuma_core_gid);
  PDM_mpi_win_shared_unlock_all(wnuma_by_numa_idx);
  PDM_mpi_win_shared_free(wnuma_by_numa_n);
  PDM_mpi_win_shared_free(wnuma_core_gid);
  PDM_mpi_win_shared_free(wnuma_by_numa_idx);

  free(neighbor_in);
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

  int n_rank_in_shm, i_rank_in_shm;
  PDM_MPI_Comm_rank (comm_shared, &i_rank_in_shm);
  PDM_MPI_Comm_size (comm_shared, &n_rank_in_shm);

  /*
   * Setup graph
   */
  if(0 == 1) {
    PDM_MPI_Comm comm_dist_graph = setup_numa_graph(comm);
  } else {
    PDM_MPI_Comm comm_dist_graph = setup_numa_graph2(comm);
  }


  PDM_MPI_Barrier (comm);
  if (i_rank == 0) {
  PDM_printf ("-- End\n");
  }

  PDM_MPI_Finalize ();

  return 0;
}
