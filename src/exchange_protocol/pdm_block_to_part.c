/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_block_to_part.h"
#include "pdm_array.h"
#include "pdm_binary_search.h"
#include "pdm_block_to_part_priv.h"
#include "pdm_distrib.h"
#include "pdm_io.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include "pdm_priv.h"
#include "pdm_timer.h"


/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \enum _btp_timer_step_t
 *
 */

typedef enum {

  BINARY_SEARCH    = 0, // Binary search step in Block-to-Part creation
  CREATE_EXCHANGE  = 1, // Collective communication in Block-to-Part creation
  DATA_EXCHANGE    = 2  // Collective communication in data exchange

} _btp_timer_step_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

static
void
_comm_graph_statistics
(
 PDM_block_to_part_t* btp
)
{
  /*
   *  Statistic of send --> requested_data_idx
   */
  int min_n_rank_connected  = btp->n_rank+1;
  int max_n_rank_connected  = -1;

  int n_connect_rank = 0;
  for(int i = 0; i < btp->n_rank; ++i) {
    if(btp->requested_data_n[i] > 0) {
      n_connect_rank++;
    }
  }

  double d_n_rank_connected = n_connect_rank;
  double mean_n_rank_connected = 0;
  PDM_MPI_Allreduce(&n_connect_rank    , &max_n_rank_connected , 1, PDM_MPI_INT   , PDM_MPI_MAX, btp->comm);
  PDM_MPI_Allreduce(&n_connect_rank    , &min_n_rank_connected , 1, PDM_MPI_INT   , PDM_MPI_MIN, btp->comm);
  PDM_MPI_Allreduce(&d_n_rank_connected, &mean_n_rank_connected, 1, PDM_MPI_DOUBLE, PDM_MPI_SUM, btp->comm);

  mean_n_rank_connected = mean_n_rank_connected/btp->n_rank;

  if(btp->i_rank == 0) {
    printf("PDM_block_to_part requested statistics : [min/max/mean] = %i / %i / %12.5e \n", min_n_rank_connected, max_n_rank_connected, mean_n_rank_connected);
  }

  n_connect_rank = 0;
  for(int i = 0; i < btp->n_rank; ++i) {
    if(btp->distributed_data_n[i] > 0) {
      n_connect_rank++;
    }
  }

  d_n_rank_connected = n_connect_rank;
  mean_n_rank_connected = 0;
  PDM_MPI_Allreduce(&n_connect_rank    , &max_n_rank_connected , 1, PDM_MPI_INT   , PDM_MPI_MAX, btp->comm);
  PDM_MPI_Allreduce(&n_connect_rank    , &min_n_rank_connected , 1, PDM_MPI_INT   , PDM_MPI_MIN, btp->comm);
  PDM_MPI_Allreduce(&d_n_rank_connected, &mean_n_rank_connected, 1, PDM_MPI_DOUBLE, PDM_MPI_SUM, btp->comm);

  mean_n_rank_connected = mean_n_rank_connected/btp->n_rank;

  if(btp->i_rank == 0) {
    printf("PDM_block_to_part to send statistics : [min/max/mean] = %i / %i / %12.5e \n", min_n_rank_connected, max_n_rank_connected, mean_n_rank_connected);
  }

}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

void
PDM_block_to_part_comm_graph_dump
(
 PDM_block_to_part_t *btp,
 const char          *filename
)
{
  // Write in parallel
  PDM_io_file_t *writer = NULL;
  PDM_l_num_t    ierr;

  PDM_io_open(filename,
              PDM_IO_FMT_BIN,
              PDM_IO_SUFF_MAN,
              "",
              PDM_IO_BACKUP_OFF,
              PDM_IO_KIND_MPI_SIMPLE,
              PDM_IO_MOD_WRITE,
              PDM_IO_NATIVE,
              btp->comm,
              -1.,
              &writer,
              &ierr);

  // Create a node identifier
  PDM_MPI_Comm shared_comm = PDM_MPI_COMM_WORLD;

  PDM_MPI_Comm_split_type(btp->comm, PDM_MPI_SPLIT_SHARED, &shared_comm);

  int i_shared_rank = 0;
  PDM_MPI_Comm_rank(shared_comm, &i_shared_rank);

  int bcast_buffer = 0;
  if (i_shared_rank == 0) {
    bcast_buffer = btp->i_rank;
  }
  PDM_MPI_Bcast(&bcast_buffer, 1, PDM_MPI_INT32_T, 0, shared_comm);

  // Block write i_rank, node and number of send data
  int s_buffer = btp->n_rank * 11 + 40 + 2 + 1; // (10 + 1 space) * n_rank + chaine + space + \n + 1
  char *buffer;
  PDM_malloc(buffer, s_buffer, char);

  for (int i = 0; i < (int) s_buffer; i++) {
    buffer[i] = '\0';
  }

  sprintf(buffer, "i_rank %10d\nnode %10d\nn_send", btp->i_rank, bcast_buffer);

  for (int j_rank = 0; j_rank < btp->n_rank; j_rank++) {
    sprintf(buffer + strlen(buffer), " %10d", btp->distributed_data_n[j_rank]);
  } // end loop on n_rank
  sprintf(buffer + strlen(buffer), " \n");

  PDM_l_num_t one = 1;
  PDM_g_num_t i_rank_gnum = (PDM_g_num_t) (btp->i_rank+1);
  PDM_io_par_interlaced_write(writer,
                              PDM_STRIDE_VAR_INTERLACED,
                              (PDM_l_num_t *) &s_buffer,
                              (PDM_l_num_t) sizeof(char),
                              one,
                              &i_rank_gnum,
                              (const void *) buffer);

  PDM_free(buffer);

  // Finalize parallel write
  PDM_io_close(writer);
  PDM_io_free(writer);
}

PDM_block_to_part_t *
PDM_block_to_part_create_from_sparse_block_and_distrib
(
 const PDM_g_num_t     *block_distrib_idx,
 const PDM_g_num_t     *delt_gnum,  // Should be betwenn [1, N]
 const int              dn_elt,
 const PDM_g_num_t    **gnum_elt,
 const int             *n_elt,
 const int              n_part,
 const PDM_MPI_Comm     comm
)
{
  PDM_block_to_part_t* btp = PDM_block_to_part_create(block_distrib_idx,
                                                      gnum_elt,
                                                      n_elt,
                                                      n_part,
                                                      comm);
  /*
   *  Post traitement du distrib_data
   */
  assert(btp->idx_partial         == NULL);
  assert(btp->n_elt_partial_block == 0);
  PDM_malloc(btp->idx_partial, btp->distributed_data_idx[btp->n_rank], int);


  // PDM_log_trace_array_int(btp->distributed_data_idx, btp->n_rank+1, "distributed_data_idx : ");
  // PDM_log_trace_array_int(btp->distributed_data, btp->distributed_data_idx[btp->n_rank], "distributed_data : ");

  for (int i = 0; i < btp->distributed_data_idx[btp->n_rank]; i++) {
    int lid = btp->distributed_data[i];
    PDM_g_num_t g_num_send = lid + btp->block_distrib_idx[btp->i_rank] + 1;
    if(dn_elt > 0) {
      int idx_in_partial_block = PDM_binary_search_long(g_num_send, delt_gnum, dn_elt);
      btp->idx_partial[i] = idx_in_partial_block;
    } else {
      btp->idx_partial[i] = -1;
    }
  }
  btp->n_elt_partial_block = dn_elt;

  if(0 == 1) {
    PDM_log_trace_array_int(btp->idx_partial, btp->distributed_data_idx[btp->n_rank], "idx_partial : ");
  }

  return btp;
}

PDM_block_to_part_t *
PDM_block_to_part_create_from_sparse_block
(
 const PDM_g_num_t     *delt_gnum,  // Should be betwenn [1, N]
 const int              dn_elt,
 const PDM_g_num_t    **gnum_elt,
 const int             *n_elt,
 const int              n_part,
 const PDM_MPI_Comm     comm
)
{
  int n_rank = -1;
  PDM_MPI_Comm_size (comm, &n_rank);

  PDM_g_num_t *_block_distrib_idx = NULL;
  PDM_malloc(_block_distrib_idx, n_rank+1, PDM_g_num_t);

  PDM_g_num_t max_g_num = 0;

  if(dn_elt > 0) {
    max_g_num = delt_gnum[dn_elt-1];
  }

  PDM_g_num_t max_part_g_num = 0;
  for(int i_part = 0; i_part < n_part; ++i_part) {
    for(int i = 0; i < n_elt[i_part]; ++i) {
      PDM_g_num_t g_num = PDM_ABS(gnum_elt[i_part][i]);
      max_part_g_num = PDM_MAX(max_part_g_num, g_num);
    }
  }

  PDM_MPI_Allgather(&max_g_num,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (&_block_distrib_idx[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    comm);

  _block_distrib_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    _block_distrib_idx[i+1] = PDM_MAX(_block_distrib_idx[i+1], _block_distrib_idx[i]);
  }

  PDM_g_num_t gmax_part_g_num = 0;
  PDM_MPI_Allreduce(&max_part_g_num, &gmax_part_g_num, 1,
                    PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, comm);

  if(_block_distrib_idx[n_rank] == 0) {
    PDM_free(_block_distrib_idx);
    _block_distrib_idx = PDM_compute_uniform_entity_distribution(comm, gmax_part_g_num);
  }

  _block_distrib_idx[n_rank] = PDM_MAX(_block_distrib_idx[n_rank], gmax_part_g_num+1);

  PDM_block_to_part_t* btp = PDM_block_to_part_create_from_sparse_block_and_distrib(_block_distrib_idx,
                                                                                    delt_gnum,
                                                                                    dn_elt,
                                                                                    gnum_elt,
                                                                                    n_elt,
                                                                                    n_part,
                                                                                    comm);
  PDM_free(_block_distrib_idx);
  return btp;
}


PDM_block_to_part_t *
PDM_block_to_part_create
(
 const PDM_g_num_t     *block_distrib_idx,
 const PDM_g_num_t    **gnum_elt,
 const int             *n_elt,
 const int              n_part,
 const PDM_MPI_Comm     comm
)
{
  PDM_block_to_part_t *btp = NULL;
  PDM_malloc(btp, 1, PDM_block_to_part_t);

  btp->comm = comm;

  btp->p2p_factor = 0.25;

  char host[1024];
  gethostname(host, 1023);

  if (!strncmp(host, "sator" , 5)) {
    btp->p2p_factor = -0.1;
  }

  char *env_var = NULL;
  env_var = getenv ("PDM_BLOCK_TO_PART_P2P_FACTOR");
  if (env_var != NULL) {
    btp->p2p_factor = atof (env_var);
  }

  btp->pttopt_comm         = 0;
  btp->n_elt_partial_block = 0;
  btp->idx_partial         = NULL;

  PDM_MPI_Comm_size (comm, &btp->n_rank);
  PDM_MPI_Comm_rank (comm, &btp->i_rank);

  /*
   * Define requested data for each process
   */

  PDM_malloc(btp->block_distrib_idx, btp->n_rank + 1, PDM_g_num_t);
  int max_data_block = -1;
  for (int i = 0; i < btp->n_rank + 1; i++) {
    btp->block_distrib_idx[i] = block_distrib_idx[i];
  }
  for (int i = 0; i < btp->n_rank; i++) {
    max_data_block = PDM_MAX(max_data_block, block_distrib_idx[i+1] - block_distrib_idx[i]) ;
  }

  btp->n_part = n_part;

  PDM_malloc(btp->requested_data_idx, btp->n_rank + 1, int);
  PDM_malloc(btp->requested_data_n  , btp->n_rank    , int);
  for (int i = 0; i < btp->n_rank; i++) {
    btp->requested_data_idx[i] = 0;
    btp->requested_data_n[i] = 0;
  }

  PDM_malloc(btp->n_elt, n_part, int  );
  PDM_malloc(btp->ind  , n_part, int *);

  for (int i = 0; i < n_part; i++) {

    btp->n_elt[i] = n_elt[i];
    PDM_malloc(btp->ind[i], n_elt[i], int);

    const PDM_g_num_t *_gnum_elt = gnum_elt[i];

    for (int j = 0; j < n_elt[i]; j++) {

      int ind = PDM_binary_search_gap_long (PDM_ABS(_gnum_elt[j]) - 1,
                                            block_distrib_idx,
                                            btp->n_rank + 1);
      btp->ind[i][j] = ind; // Temporary use of this array to avoid le PDM_binary_search_gap_long
      // printf(" [%i][%i] --> ind = %i (g_num = %i )\n", i, j, ind, (int) _gnum_elt[j]);
      btp->requested_data_n[ind]++;

    }
  }

  for (int i = 0; i < btp->n_rank; i++) {
    btp->requested_data_idx[i+1] = btp->requested_data_idx[i] + btp->requested_data_n  [i];
  }

  int s_requested_data = btp->requested_data_idx[btp->n_rank - 1]
                       + btp->requested_data_n  [btp->n_rank - 1];

  int *requested_data = NULL;
  PDM_malloc(requested_data, s_requested_data, int);

  for (int i = 0; i < btp->n_rank; i++) {
    btp->requested_data_n[i] = 0;
  }

  for (int i = 0; i < n_part; i++) {

    const PDM_g_num_t *_gnum_elt = gnum_elt[i];

    // printf("n_elt[%i] = %i \n", i, (int) n_elt[i]);
    for (int j = 0; j < n_elt[i]; j++) {

      // int ind = PDM_binary_search_gap_long (_gnum_elt[j] - 1,
      //                                       block_distrib_idx,
      //                                       btp->n_rank + 1);
      int ind = btp->ind[i][j];
      int idx = btp->requested_data_idx[ind] + btp->requested_data_n[ind]++;

      btp->ind[i][j] = idx;

      PDM_g_num_t _requested_data = PDM_ABS(_gnum_elt[j]) - 1 - block_distrib_idx[ind];
      // printf("requested_data[%i] = %i / size_max = %i and gn_m = %i \n", idx, (int) _requested_data, s_requested_data, (int)_gnum_elt[j]);
      requested_data[idx] = (int) _requested_data;
    }
  }

  PDM_malloc(btp->distributed_data_n, btp->n_rank, int);

  PDM_MPI_Alltoall (btp->requested_data_n,   1, PDM_MPI_INT,
                    btp->distributed_data_n, 1, PDM_MPI_INT,
                    comm);

  btp->distributed_data_idx = PDM_array_new_idx_from_sizes_int(btp->distributed_data_n, btp->n_rank);

  PDM_malloc(btp->distributed_data, btp->distributed_data_idx[btp->n_rank], int);

  PDM_MPI_part_of_active_rank(btp->requested_data_n,
                              btp->distributed_data_n,
                              comm,
                              &(btp->part_active_rank));

  if (btp->p2p_factor < btp->part_active_rank) {

    PDM_MPI_Alltoallv (requested_data,
                       btp->requested_data_n,
                       btp->requested_data_idx,
                       PDM_MPI_INT,
                       btp->distributed_data,
                       btp->distributed_data_n,
                       btp->distributed_data_idx,
                       PDM_MPI_INT,
                       comm);
  }

  else {

    PDM_MPI_Alltoallv_p2p (requested_data,
                           btp->requested_data_n,
                           btp->requested_data_idx,
                           PDM_MPI_INT,
                           btp->distributed_data,
                           btp->distributed_data_n,
                           btp->distributed_data_idx,
                           PDM_MPI_INT,
                           comm);

  }

  // For large data

  int coeff = 10;
  if (btp->distributed_data_idx[btp->n_rank] >= coeff * max_data_block) {
    btp->pttopt_comm = 1;
  }

  if(0 == 1) {
    _comm_graph_statistics(btp);
  }

  //PDM_log_trace_array_long(btp->distributed_data_idx, btp->n_rank+1, "block_distrib");

  int tmp;
  PDM_MPI_Allreduce (&(btp->pttopt_comm), &tmp, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);
  btp->pttopt_comm = tmp;

  PDM_free(requested_data);

  int n_rank_recv = 0;
  int n_rank_send = 0;

  for (int i = 0; i < btp->n_rank; i++) {
    if (btp->i_rank != i && btp->distributed_data_n[i] > 0) {
      n_rank_recv += 1;
    }
    if (btp->i_rank != i && btp->requested_data_n[i] > 0) {
      n_rank_send += 1;
    }
  }

  return (PDM_block_to_part_t *) btp;
}


void
PDM_block_to_part_exch_in_place
(
 PDM_block_to_part_t *btp,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                 *block_stride,
 void                *block_data,
 int                **part_stride,
 void               **part_data
)
{
  unsigned char *_block_data = (unsigned char *) block_data;
  unsigned char **_part_data = (unsigned char **) part_data;

  int n_elt_block = btp->block_distrib_idx[btp->i_rank+1] - btp->block_distrib_idx[btp->i_rank];

  size_t *i_send_buffer = NULL;
  size_t *i_recv_buffer = NULL;
  int    *n_send_buffer = NULL;
  int    *n_recv_buffer = NULL;
  PDM_malloc(i_send_buffer, btp->n_rank, size_t);
  PDM_malloc(i_recv_buffer, btp->n_rank, size_t);
  PDM_malloc(n_send_buffer, btp->n_rank, int   );
  PDM_malloc(n_recv_buffer, btp->n_rank, int   );
  int max_n_send_buffer = -1;
  int max_n_recv_buffer = -1;
  int *block_stride_idx = NULL;

  for (int i = 0; i < btp->n_rank; i++) {
    n_send_buffer[i] = 0;
    n_recv_buffer[i] = 0;
    i_send_buffer[i] = 0;
    i_recv_buffer[i] = 0;
  }

  unsigned char **send_buffer = NULL;
  unsigned char *recv_buffer = NULL;

  size_t s_send_buffer = 0;
  size_t s_recv_buffer = 0;

  int n_rank1 = btp->n_rank - 1;

  int s_distributed_data = btp->distributed_data_idx[btp->n_rank];

  int s_data_tot = s_data;
  if(t_stride == PDM_STRIDE_CST_INTERLACED) {
    int cst_stride = *block_stride;
    s_data_tot = s_data * cst_stride;
  }

  PDM_MPI_Datatype mpi_type;
  PDM_MPI_Type_create_contiguous(s_data_tot, PDM_MPI_BYTE, &mpi_type);
  PDM_MPI_Type_commit(&mpi_type);

  /* int step; */

  int rank;
  PDM_MPI_Comm_rank(btp->comm, &rank);

  /*
   * Exchange Stride and build buffer properties
   */

  int *recv_stride = NULL;
  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {

    int s_send_stride = btp->distributed_data_idx[btp->n_rank];

    int s_recv_stride = btp->requested_data_idx[btp->n_rank];

    int *send_stride = NULL;
    PDM_malloc(send_stride, s_send_stride, int);
    PDM_malloc(recv_stride, s_recv_stride, int);

    if(btp->idx_partial == NULL) { // block is full
      for (int i = 0; i < s_send_stride; i++) {
        send_stride[i] = block_stride[btp->distributed_data[i]];
      }
    } else {                       // block is partial and describe by delt_gnum
      for (int i = 0; i < s_send_stride; i++) {
        if(btp->idx_partial[i] != -1) {
          send_stride[i] = block_stride[btp->idx_partial[i]];
        } else {
          send_stride[i] = 0;
        }
      }
    }

    if (btp->p2p_factor < btp->part_active_rank) {

      PDM_MPI_Alltoallv (send_stride,
                         btp->distributed_data_n,
                         btp->distributed_data_idx,
                         PDM_MPI_INT,
                         recv_stride,
                         btp->requested_data_n,
                         btp->requested_data_idx,
                         PDM_MPI_INT,
                         btp->comm);
    }

    else {

      PDM_MPI_Alltoallv_p2p(send_stride,
                             btp->distributed_data_n,
                             btp->distributed_data_idx,
                             PDM_MPI_INT,
                             recv_stride,
                             btp->requested_data_n,
                             btp->requested_data_idx,
                             PDM_MPI_INT,
                             btp->comm);
    }

    for (int i = 0; i < btp->n_part; i++) {
      for (int j = 0; j < btp->n_elt[i]; j++) {
        int ielt = btp->ind[i][j];
        part_stride[i][j] = recv_stride[ielt];
      }
    }

    /*
     * Build buffers
     */

    for (int i = 0; i < btp->n_rank; i++) {
      int ibeg = btp->distributed_data_idx[i];
      int iend = btp->distributed_data_idx[i] +
                 btp->distributed_data_n[i];

      n_send_buffer[i] = 0;
      for (int k = ibeg; k < iend; k++)  {
        n_send_buffer[i] += send_stride[k];
      }

      // n_send_buffer[i] *= (int) s_data;
      max_n_send_buffer = PDM_MAX(max_n_send_buffer, n_send_buffer[i] * s_data_tot);

      if (i > 0) {
        i_send_buffer[i] = i_send_buffer[i-1] + n_send_buffer[i-1];
      } else {
        i_send_buffer[i] = 0;
      }

      ibeg = btp->requested_data_idx[i];
      iend = btp->requested_data_idx[i] + btp->requested_data_n[i];

      n_recv_buffer[i] = 0;
      for (int k = ibeg; k < iend; k++) {
        n_recv_buffer[i] += recv_stride[k];
      }

      // n_recv_buffer[i] *= (int) s_data;
      max_n_recv_buffer = PDM_MAX(max_n_recv_buffer, n_recv_buffer[i]);

      if (i > 0) {
        i_recv_buffer[i] = i_recv_buffer[i-1] + n_recv_buffer[i-1];
      } else {
        i_recv_buffer[i] = 0;
      }
    }

    if(btp->idx_partial == NULL) {
      block_stride_idx = PDM_array_new_idx_from_sizes_int(block_stride, n_elt_block);
    } else {
      block_stride_idx = PDM_array_new_idx_from_sizes_int(block_stride, btp->n_elt_partial_block);
    }
    PDM_free(send_stride);
  }

  else {

    // int cst_stride = *block_stride;
    max_n_send_buffer = 0;
    max_n_recv_buffer = 0;

    for (int i = 0; i < btp->n_rank; i++) {

      i_send_buffer[i] = btp->distributed_data_idx[i]; // * cst_stride; //  * (int) s_data;
      i_recv_buffer[i] = btp->requested_data_idx  [i]; // * cst_stride; //  * (int) s_data;

      n_send_buffer[i] = btp->distributed_data_n[i]; //  * cst_stride; // * (int) s_data;
      n_recv_buffer[i] = btp->requested_data_n  [i]; //  * cst_stride; // * (int) s_data;
      max_n_send_buffer = PDM_MAX(max_n_send_buffer, n_send_buffer[i] * s_data_tot);
      max_n_recv_buffer = PDM_MAX(max_n_recv_buffer, n_recv_buffer[i] * s_data_tot);

    }

    // s_send_buffer = i_send_buffer[n_rank1] + n_send_buffer[n_rank1];
    // s_recv_buffer = i_recv_buffer[n_rank1] + n_recv_buffer[n_rank1];

  }

  s_send_buffer = (i_send_buffer[n_rank1] + n_send_buffer[n_rank1]) * s_data_tot;
  s_recv_buffer = (i_recv_buffer[n_rank1] + n_recv_buffer[n_rank1]) * s_data_tot;

  int n_active_buffer;

  if (btp->pttopt_comm) {
    n_active_buffer = 5;
  }
  else {
    n_active_buffer = 1;
  }

  PDM_malloc(send_buffer, n_active_buffer, unsigned char *);

  if (btp->pttopt_comm) {
    for (int i = 0; i < n_active_buffer; i++) {
      PDM_malloc(send_buffer[i], max_n_send_buffer, unsigned char);
    }
  } else {
    PDM_malloc(send_buffer[0], s_send_buffer, unsigned char);
  }

  PDM_malloc(recv_buffer, s_recv_buffer, unsigned char );

  if (btp->pttopt_comm) {

    PDM_MPI_Request *s_request;
    PDM_MPI_Request *r_request;
    PDM_malloc(s_request, n_active_buffer, PDM_MPI_Request);
    PDM_malloc(r_request, btp->n_rank    , PDM_MPI_Request);

    for (int i = 0; i < btp->n_rank; i++) {
      if (n_recv_buffer[i] > 0) {
        PDM_MPI_Irecv(recv_buffer + i_recv_buffer[i] * s_data_tot,
                      n_recv_buffer[i],
                      mpi_type,
                      i,
                      0,
                      btp->comm,
                      r_request + i);
      }
    }

    int *active_rank = NULL;
    PDM_malloc(active_rank, n_active_buffer, int);
    for (int i = 0; i < n_active_buffer; i++) {
      active_rank[i] = i;
    }

    while (1) {
      int _n_active_buffer = 0;
      for (int i = 0; i < n_active_buffer; i++) {
        if (active_rank[i] < btp->n_rank) {
          _n_active_buffer += 1;
        }
      }

      if (_n_active_buffer == 0) {
        break;
      }

      for (int i = 0; i < _n_active_buffer; i++) {
        if (n_send_buffer[active_rank[i]] > 0) {

          int s_distributed_active_rank = btp->distributed_data_idx[active_rank[i]] +
                                          btp->distributed_data_n [active_rank[i]];

          if (t_stride == PDM_STRIDE_VAR_INTERLACED) {
            int idx1 = 0;

            if(btp->idx_partial == NULL) { // block is full
              for (int j = btp->distributed_data_idx[active_rank[i]];
                       j < s_distributed_active_rank; j++) {

                int ind =  block_stride_idx[btp->distributed_data[j]] * (int) s_data;

                int s_block_unit =  block_stride[btp->distributed_data[j]] * (int) s_data;

                unsigned char *_block_data_deb = _block_data + ind;

                for (int k = 0; k < s_block_unit; k++) {
                  send_buffer[i][idx1++] = _block_data_deb[k];
                }
              }
            } else {  // block is partial and describe by delt_gnum

              for (int j = btp->distributed_data_idx[active_rank[i]];
                       j < s_distributed_active_rank; j++) {

                if(btp->idx_partial[j] != -1) {
                  int ind =  block_stride_idx[btp->idx_partial[j]] * (int) s_data;
                  int s_block_unit =  block_stride[btp->idx_partial[j]] * (int) s_data;
                  unsigned char *_block_data_deb = _block_data + ind;

                  for (int k = 0; k < s_block_unit; k++) {
                    send_buffer[i][idx1++] = _block_data_deb[k];
                  }
                }
              }
            }

          }
          else {
            int cst_stride = *block_stride;
            int s_block_unit = cst_stride * (int) s_data;

            int idx1 = 0;

            if(btp->idx_partial == NULL) { // block is full
              for (int j = btp->distributed_data_idx[active_rank[i]];
                       j < s_distributed_active_rank; j++) {
                int ind = btp->distributed_data[j];
                unsigned char *_block_data_deb = _block_data + ind * cst_stride * (int) s_data;
                for (int k = 0; k < s_block_unit; k++) {
                  send_buffer[i][idx1++] = _block_data_deb[k];
                }
              }
            } else {  // block is partial and describe by delt_gnum
              for (int j = btp->distributed_data_idx[active_rank[i]];
                       j < s_distributed_active_rank; j++) {
                int ind = btp->idx_partial[j];
                if(ind != -1) {
                  unsigned char *_block_data_deb = _block_data + ind * cst_stride * (int) s_data;
                  for (int k = 0; k < s_block_unit; k++) {
                    send_buffer[i][idx1++] = _block_data_deb[k];
                  }
                }
              }
            }
          }

          PDM_MPI_Isend(send_buffer[i],
                        n_send_buffer[active_rank[i]],
                        mpi_type,
                        active_rank[i],
                        0,
                        btp->comm,
                        s_request + i);
        }
      }

      for (int i = 0; i < _n_active_buffer; i++) {
        if (n_send_buffer[active_rank[i]] > 0) {
          PDM_MPI_Wait (s_request + i);
        }
      }

      for (int i = 0; i < n_active_buffer; i++) {
        active_rank[i] += n_active_buffer;
      }

    }

    for (int i = 0; i < btp->n_rank; i++) {
      if (n_recv_buffer[i] > 0) {
        PDM_MPI_Wait (r_request + i);
      }
    }

    PDM_free(s_request);
    PDM_free(r_request);
    PDM_free(active_rank);

  } else {

    if (t_stride == PDM_STRIDE_VAR_INTERLACED) {
      int idx1 = 0;

      if(btp->idx_partial == NULL) { // block is full
        for (int i = 0; i < s_distributed_data; i++) {
          int ind =  block_stride_idx[btp->distributed_data[i]] * (int) s_data;
          int s_block_unit =  block_stride[btp->distributed_data[i]] * (int) s_data;
          unsigned char *_block_data_deb = _block_data + ind;
          for (int k = 0; k < s_block_unit; k++) {
            send_buffer[0][idx1++] = _block_data_deb[k];
          }
        }
      } else { // block is partial and describe by delt_gnum
        for (int i = 0; i < s_distributed_data; i++) {
          if(btp->idx_partial[i] != -1) {
            int ind =  block_stride_idx[btp->idx_partial[i]] * (int) s_data;
            int s_block_unit =  block_stride[btp->idx_partial[i]] * (int) s_data;
            unsigned char *_block_data_deb = _block_data + ind;
            for (int k = 0; k < s_block_unit; k++) {
              send_buffer[0][idx1++] = _block_data_deb[k];
            }
          }
        }
      }
    }
    else {
      int idx1 = 0;
      int cst_stride = *block_stride;
      int s_block_unit = cst_stride * (int) s_data;

      if(btp->idx_partial == NULL) { // block is full
        for (int i = 0; i < s_distributed_data; i++) {
          int ind = btp->distributed_data[i];
          unsigned char *_block_data_deb = _block_data + ind * cst_stride * (int) s_data;
          for (int k = 0; k < s_block_unit; k++) {
            send_buffer[0][idx1++] = _block_data_deb[k];
          }
        }
      } else { // block is partial and describe by delt_gnum
        for (int i = 0; i < s_distributed_data; i++) {
          int ind = btp->idx_partial[i];
          if(ind != -1) {
            unsigned char *_block_data_deb = _block_data + ind * cst_stride * (int) s_data;
            for (int k = 0; k < s_block_unit; k++) {
              send_buffer[0][idx1++] = _block_data_deb[k];
            }
          }
        }
      }
    }


    if (btp->p2p_factor < btp->part_active_rank) {

      int *_i_send_buffer = NULL;
      int *_i_recv_buffer = NULL;
      PDM_malloc(_i_send_buffer, btp->n_rank, int);
      PDM_malloc(_i_recv_buffer, btp->n_rank, int);

      for (int i = 0; i < btp->n_rank; i++) {
        _i_send_buffer[i] = (int) i_send_buffer[i];
        _i_recv_buffer[i] = (int) i_recv_buffer[i];
      }

      PDM_MPI_Alltoallv(send_buffer[0],
                        n_send_buffer,
                        _i_send_buffer,
                        mpi_type,
                        recv_buffer,
                        n_recv_buffer,
                        _i_recv_buffer,
                        mpi_type,
                        btp->comm);

      PDM_free(_i_send_buffer);
      PDM_free(_i_recv_buffer);
    } else {
  
      PDM_MPI_Alltoallv_p2p_l(send_buffer[0],
                              n_send_buffer,
                              i_send_buffer,
                              mpi_type,
                              recv_buffer,
                              n_recv_buffer,
                              i_recv_buffer,
                              mpi_type,
                              btp->comm);
  
    }

    // PDM_MPI_Alltoallv_l(send_buffer[0],
    //                     n_send_buffer,
    //                     i_send_buffer,
    //                     mpi_type,
    //                     recv_buffer,
    //                     n_recv_buffer,
    //                     i_recv_buffer,
    //                     mpi_type,
    //                     btp->comm);
  
  }

  for (int i = 0; i < n_active_buffer; i++) {
    PDM_free(send_buffer[i]);
  }
  PDM_free(send_buffer);
  PDM_free(n_send_buffer);
  PDM_free(i_send_buffer);
  PDM_free(n_recv_buffer);
  PDM_free(i_recv_buffer);

  if (block_stride_idx != NULL) {
    PDM_free(block_stride_idx);
  }

  /*
   * Partitions filling
   */

  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {

    int s_recv_elt = btp->requested_data_idx[n_rank1] +
      btp->requested_data_n[n_rank1];

    int **part_idx = NULL;
    PDM_malloc(part_idx, btp->n_part, int *);
    int  *recv_idx = PDM_array_new_idx_from_sizes_int(recv_stride, s_recv_elt);

    for (int i = 0; i < btp->n_part; i++) {
      part_idx[i] = PDM_array_new_idx_from_sizes_int(part_stride[i], btp->n_elt[i]);
    }

    for (int i = 0; i < btp->n_part; i++) {

      for (int j = 0; j < btp->n_elt[i]; j++) {

        int idx1  = part_idx   [i][j] * (int) s_data;
        int n_elt = part_stride[i][j] * (int) s_data;

        int idx2 = recv_idx[btp->ind[i][j]] * (int) s_data;

        for (int k = 0; k < n_elt; k++) {
          _part_data[i][idx1+k] = recv_buffer[idx2+k];
        }
      }
    }

    for (int i = 0; i < btp->n_part; i++) {
      PDM_free(part_idx[i]);
    }

    PDM_free(recv_idx);
    PDM_free(part_idx);
    PDM_free(recv_stride);
  }

  else if (t_stride == PDM_STRIDE_CST_INTERLACED) {

    const int cst_stride = *block_stride;
    const int s_block_unit = cst_stride * (int) s_data;

    for (int i = 0; i < btp->n_part; i++) {

      for (int j = 0; j < btp->n_elt[i]; j++) {

        int idx1 = j * s_block_unit;
        int idx2 = btp->ind[i][j] * s_block_unit;

        for (int k = 0; k < s_block_unit; k++) {
          _part_data[i][idx1+k] = recv_buffer[idx2+k];
        }
      }
    }
  }

  PDM_MPI_Type_free(&mpi_type);

  PDM_free(recv_buffer);
}

void
PDM_block_to_part_exch
(
 PDM_block_to_part_t *btp,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                 *block_stride,
 void                *block_data,
 int               ***part_stride,
 void              ***part_data
)
{
  int n_elt_block = btp->block_distrib_idx[btp->i_rank+1] - btp->block_distrib_idx[btp->i_rank];

  unsigned char *_block_data = (unsigned char *) block_data;
  unsigned char **_part_data;

  size_t *i_send_buffer = NULL;
  size_t *i_recv_buffer = NULL;
  int    *n_send_buffer = NULL;
  int    *n_recv_buffer = NULL;
  PDM_malloc(i_send_buffer, btp->n_rank, size_t);
  PDM_malloc(i_recv_buffer, btp->n_rank, size_t);
  PDM_malloc(n_send_buffer, btp->n_rank, int   );
  PDM_malloc(n_recv_buffer, btp->n_rank, int   );

  for (int i = 0; i < btp->n_rank; i++) {
    n_send_buffer[i] = 0;
    n_recv_buffer[i] = 0;
    i_send_buffer[i] = 0;
    i_recv_buffer[i] = 0;
  }

  unsigned char *send_buffer = NULL;
  unsigned char *recv_buffer = NULL;

  size_t s_send_buffer = 0;
  size_t s_recv_buffer = 0;

  int n_rank1 = btp->n_rank - 1;

  int s_distributed_data = btp->distributed_data_idx[btp->n_rank];

  int s_data_tot = s_data;
  if(t_stride == PDM_STRIDE_CST_INTERLACED) {
    int cst_stride = *block_stride;
    s_data_tot = s_data * cst_stride;
  }

  PDM_MPI_Datatype mpi_type;
  PDM_MPI_Type_create_contiguous(s_data_tot, PDM_MPI_BYTE, &mpi_type);
  PDM_MPI_Type_commit(&mpi_type);

  /*
   * Exchange Stride and build buffer properties
   */

  int *recv_stride = NULL;
  int **_part_stride = NULL;

  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {

    int s_send_stride = btp->distributed_data_idx[btp->n_rank];

    int s_recv_stride = btp->requested_data_idx[btp->n_rank];

    int *send_stride = NULL;
    PDM_malloc(send_stride, s_send_stride, int);
    PDM_malloc(recv_stride, s_recv_stride, int);

    if(btp->idx_partial == NULL) { // block is full
      for (int i = 0; i < s_send_stride; i++) {
        send_stride[i] = block_stride[btp->distributed_data[i]];
      }
    } else {                       // block is partial and describe by delt_gnum
      for (int i = 0; i < s_send_stride; i++) {
        if(btp->idx_partial[i] != -1) {
          send_stride[i] = block_stride[btp->idx_partial[i]];
        } else {
          send_stride[i] = 0;
        }
      }
    }

    if (btp->p2p_factor < btp->part_active_rank) {

      PDM_MPI_Alltoallv (send_stride,
                         btp->distributed_data_n,
                         btp->distributed_data_idx,
                         PDM_MPI_INT,
                         recv_stride,
                         btp->requested_data_n,
                         btp->requested_data_idx,
                         PDM_MPI_INT,
                         btp->comm);
    }

    else {

      PDM_MPI_Alltoallv_p2p (send_stride,
                             btp->distributed_data_n,
                             btp->distributed_data_idx,
                             PDM_MPI_INT,
                             recv_stride,
                             btp->requested_data_n,
                             btp->requested_data_idx,
                             PDM_MPI_INT,
                             btp->comm);

    }

    PDM_malloc(*part_stride, btp->n_part, int *);
    _part_stride = *part_stride;

    for (int i = 0; i < btp->n_part; i++) {

      PDM_malloc(_part_stride[i], btp->n_elt[i], int);

      for (int j = 0; j < btp->n_elt[i]; j++) {

        int ielt = btp->ind[i][j];
        _part_stride[i][j] = recv_stride[ielt];

      }
    }

    /*
     * Build buffers
     */

    for (int i = 0; i < btp->n_rank; i++) {
      int ibeg = btp->distributed_data_idx[i];
      int iend = btp->distributed_data_idx[i] + btp->distributed_data_n[i];

      n_send_buffer[i] = 0;
      for (int k = ibeg; k < iend; k++)  {
        n_send_buffer[i] += send_stride[k];
      }

      // n_send_buffer[i] *= (int) s_data;

      if (i > 0) {
        i_send_buffer[i] = i_send_buffer[i-1] + n_send_buffer[i-1];
      } else {
        i_send_buffer[i] = 0;
      }

      ibeg = btp->requested_data_idx[i];
      iend = btp->requested_data_idx[i] + btp->requested_data_n[i];

      n_recv_buffer[i] = 0;
      for (int k = ibeg; k < iend; k++) {
        n_recv_buffer[i] += recv_stride[k];
      }

      // n_recv_buffer[i] *= (int) s_data;

      if (i > 0) {
        i_recv_buffer[i] = i_recv_buffer[i-1] + n_recv_buffer[i-1];
      } else {
        i_recv_buffer[i] = 0;
      }

    }

    s_send_buffer = (i_send_buffer[n_rank1] + n_send_buffer[n_rank1]) * s_data_tot;
    s_recv_buffer = (i_recv_buffer[n_rank1] + n_recv_buffer[n_rank1]) * s_data_tot;

    PDM_malloc(send_buffer, s_send_buffer, unsigned char);
    PDM_malloc(recv_buffer, s_recv_buffer, unsigned char);

    // int *send_stride_idx;
    // PDM_malloc(send_stride_idx, s_distributed_data+1, int);
    // send_stride_idx[0] = 0;
    // for (int i = 0; i < s_distributed_data; i++) {
    //   send_stride_idx[i+1] = send_stride_idx[i] + send_stride[i];
    // }

    int idx1 = 0;
    int *block_stride_idx = NULL;
    if(btp->idx_partial == NULL) {
      block_stride_idx = PDM_array_new_idx_from_sizes_int(block_stride, n_elt_block);
    } else {
      // printf("btp->n_elt_partial_block = %i \n", btp->n_elt_partial_block);
      block_stride_idx = PDM_array_new_idx_from_sizes_int(block_stride, btp->n_elt_partial_block);
    }

    if(btp->idx_partial == NULL) { // block is full
      for (int i = 0; i < s_distributed_data; i++) {

        int ind =  block_stride_idx[btp->distributed_data[i]] * (int) s_data;

        int s_block_unit =  block_stride[btp->distributed_data[i]] * (int) s_data;

        unsigned char *_block_data_deb = _block_data + ind;
        for (int k = 0; k < s_block_unit; k++) {
          send_buffer[idx1++] = _block_data_deb[k];
        }
      }
    } 

    else { // block is partial and describe by delt_gnum
      for (int i = 0; i < s_distributed_data; i++) {

        if(btp->idx_partial[i] != -1) {
          int ind =  block_stride_idx[btp->idx_partial[i]] * (int) s_data;

          int s_block_unit =  block_stride[btp->idx_partial[i]] * (int) s_data;

          unsigned char *_block_data_deb = _block_data + ind;
          for (int k = 0; k < s_block_unit; k++) {
            send_buffer[idx1++] = _block_data_deb[k];
          }
        }
      }
    }
    PDM_free(send_stride);
    //PDM_free(send_stride_idx);
    PDM_free(block_stride_idx);

  }

  else if (t_stride == PDM_STRIDE_CST_INTERLACED) {

    int cst_stride = *block_stride;
    int s_block_unit = cst_stride * (int) s_data;

    for (int i = 0; i < btp->n_rank; i++) {

      i_send_buffer[i] = btp->distributed_data_idx[i]; // * cst_stride * (int) s_data;
      i_recv_buffer[i] = btp->requested_data_idx  [i]; // * cst_stride * (int) s_data;

      n_send_buffer[i] = btp->distributed_data_n[i]; // * cst_stride * (int) s_data;
      n_recv_buffer[i] = btp->requested_data_n  [i]; // * cst_stride * (int) s_data;

    }

    s_send_buffer = (i_send_buffer[n_rank1] + n_send_buffer[n_rank1]) * s_data_tot;
    s_recv_buffer = (i_recv_buffer[n_rank1] + n_recv_buffer[n_rank1]) * s_data_tot;

    PDM_malloc(send_buffer, s_send_buffer, unsigned char);
    PDM_malloc(recv_buffer, s_recv_buffer, unsigned char);

    int idx1 = 0;

    if(btp->idx_partial == NULL) { // block is full
      for (int i = 0; i < s_distributed_data; i++) {
        int ind = btp->distributed_data[i];
        unsigned char *_block_data_deb = _block_data + ind * cst_stride * (int) s_data;
        for (int k = 0; k < s_block_unit; k++) {
          send_buffer[idx1++] = _block_data_deb[k];
        }
      }
    }
     else { // block is partial and describe by delt_gnum
      for (int i = 0; i < s_distributed_data; i++) {
        int ind = btp->idx_partial[i];
        if(ind != -1) {
          unsigned char *_block_data_deb = _block_data + ind * cst_stride * (int) s_data;
          for (int k = 0; k < s_block_unit; k++) {
            send_buffer[idx1++] = _block_data_deb[k];
          }
        }
      }
    }
  }

  /*
   * Data exchange
   */
  if (btp->p2p_factor  < btp->part_active_rank) {

    int *_i_send_buffer = NULL;
    int *_i_recv_buffer = NULL;
    PDM_malloc(_i_send_buffer, btp->n_rank, int);
    PDM_malloc(_i_recv_buffer, btp->n_rank, int);

    for (int i = 0; i < btp->n_rank; i++) {
      _i_send_buffer[i] = (int) i_send_buffer[i];
      _i_recv_buffer[i] = (int) i_recv_buffer[i];
    }

    PDM_MPI_Alltoallv(send_buffer,
                      n_send_buffer,
                      _i_send_buffer,
                      mpi_type,
                      recv_buffer,
                      n_recv_buffer,
                      _i_recv_buffer,
                      mpi_type,
                      btp->comm);

    PDM_free(_i_send_buffer);
    PDM_free(_i_recv_buffer);
  } else {

    PDM_MPI_Alltoallv_p2p_l(send_buffer,
                            n_send_buffer,
                            i_send_buffer,
                            mpi_type,
                            recv_buffer,
                            n_recv_buffer,
                            i_recv_buffer,
                            mpi_type,
                            btp->comm);

  }

  // PDM_MPI_Alltoallv_l(send_buffer,
  //                     n_send_buffer,
  //                     i_send_buffer,
  //                     mpi_type,
  //                     recv_buffer,
  //                     n_recv_buffer,
  //                     i_recv_buffer,
  //                     mpi_type,
  //                     btp->comm);

  PDM_free(send_buffer);
  PDM_free(n_send_buffer);
  PDM_free(i_send_buffer);
  PDM_free(n_recv_buffer);
  PDM_free(i_recv_buffer);

  /*
   * Partitions filling
   */

  PDM_malloc(*((unsigned char ***) part_data), btp->n_part, unsigned char *);
  _part_data = (*(unsigned char ***) part_data);

  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {

    int s_recv_elt = btp->requested_data_idx[n_rank1] +
                     btp->requested_data_n[n_rank1];

    int **part_idx;
    PDM_malloc(part_idx, btp->n_part, int *);
    int *recv_idx = PDM_array_new_idx_from_sizes_int(recv_stride, s_recv_elt);

    for (int i = 0; i < btp->n_part; i++) {
      part_idx[i] = PDM_array_new_idx_from_sizes_int(_part_stride[i], btp->n_elt[i]);
    }

    for (int i = 0; i < btp->n_part; i++) {

      int s_part =  part_idx[i][btp->n_elt[i]] * (int) s_data;

      PDM_malloc(_part_data[i], s_part, unsigned char);

      for (int j = 0; j < btp->n_elt[i]; j++) {

        int idx1  = part_idx[i][j] * (int) s_data;
        int n_elt = _part_stride[i][j] * (int) s_data;

        int idx2 = recv_idx[btp->ind[i][j]] * (int) s_data;

        for (int k = 0; k < n_elt; k++) {
           _part_data[i][idx1+k] = recv_buffer[idx2+k];
        }
      }
    }

    for (int i = 0; i < btp->n_part; i++) {
      PDM_free(part_idx[i]);
    }

    PDM_free(recv_idx);
    PDM_free(part_idx);
    PDM_free(recv_stride);
  }

  else if (t_stride == PDM_STRIDE_CST_INTERLACED) {

    const int cst_stride = *block_stride;
    const int s_block_unit = cst_stride * (int) s_data;

    for (int i = 0; i < btp->n_part; i++) {

      PDM_malloc(_part_data[i], s_block_unit * btp->n_elt[i], unsigned char);

      for (int j = 0; j < btp->n_elt[i]; j++) {

        int idx1 = j * s_block_unit;
        int idx2 = btp->ind[i][j] * s_block_unit;

        for (int k = 0; k < s_block_unit; k++) {
           _part_data[i][idx1+k] = recv_buffer[idx2+k];
        }
      }
    }
  }

  PDM_free(recv_buffer);
  PDM_MPI_Type_free(&mpi_type);
}

PDM_block_to_part_t *
PDM_block_to_part_free
(
 PDM_block_to_part_t *btp
)
{

  for (int i = 0; i < btp->n_part; i++) {
    PDM_free(btp->ind[i]);
  }

  PDM_free(btp->ind);
  PDM_free(btp->n_elt);
  PDM_free(btp->block_distrib_idx);
  PDM_free(btp->distributed_data);
  PDM_free(btp->distributed_data_idx);
  PDM_free(btp->distributed_data_n);
  PDM_free(btp->requested_data_idx);
  PDM_free(btp->requested_data_n);

  if(btp->idx_partial != NULL) {
    PDM_free(btp->idx_partial);
  }

  PDM_free(btp);

  return NULL;
}


PDM_l_num_t
PDM_block_to_part_gnum_idx_get
(
 PDM_block_to_part_t *btp,
 PDM_g_num_t gNum
)
{
  return (PDM_l_num_t) (gNum - 1 - btp->block_distrib_idx[btp->i_rank]);
}

int
PDM_block_to_part_n_part_get
(
 PDM_block_to_part_t *btp
 )
{
  assert (btp != NULL);

  return btp->n_part;
}

int
PDM_block_to_part_n_elt_get
(
 PDM_block_to_part_t *btp,
 const int            i_part
 )
{
  assert (btp != NULL);
  assert (i_part < btp->n_part);

  return btp->n_elt[i_part];
}

