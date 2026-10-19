/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_block_to_block.h"
#include "pdm.h"
#include "pdm_block_to_block_priv.h"
#include "pdm_mem_tool.h"
#include "pdm_priv.h"


/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

static inline int _overlap_size(PDM_g_num_t start1, PDM_g_num_t end1, PDM_g_num_t start2, PDM_g_num_t end2)
{
  int diff = PDM_MIN(end1, end2) - PDM_MAX(start1, start2);
  return PDM_MAX(diff, 0);
}

static inline void _compute_displ(int* counts, int* displs, int size) {
  displs[0] = 0;
  for (int i = 1; i < size; ++i) {
    displs[i] = displs[i-1] + counts[i-1];
  }
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_block_to_block_t *
PDM_block_to_block_create
(
 PDM_g_num_t   *block_distrib_ini_idx,
 PDM_g_num_t   *block_distrib_end_idx,
 PDM_MPI_Comm   comm
)
{

  int i_rank = -1;
  int n_rank = 0;
  PDM_MPI_Comm_size (comm, &n_rank);
  PDM_MPI_Comm_rank (comm, &i_rank);

  _pdm_block_to_block_t *btb;
  PDM_malloc(btb, 1, _pdm_block_to_block_t);

  btb->comm = comm;
  btb->i_rank = i_rank;
  btb->n_rank = n_rank;

  /*
   * Define requested data for each process
   */
  PDM_malloc(btb->block_distrib_ini_idx, n_rank + 1, PDM_g_num_t);
  for (int i = 0; i < n_rank + 1; i++) {
    btb->block_distrib_ini_idx[i] = block_distrib_ini_idx[i];
  }

  PDM_malloc(btb->block_distrib_end_idx, n_rank + 1, PDM_g_num_t);
  for (int i = 0; i < n_rank + 1; i++) {
    btb->block_distrib_end_idx[i] = block_distrib_end_idx[i];
  }

  PDM_malloc(btb->n_send_buffer, n_rank, int);
  PDM_malloc(btb->n_recv_buffer, n_rank, int);
  for (int i = 0; i < n_rank; ++i) {
    btb->n_send_buffer[i] = _overlap_size(block_distrib_ini_idx[i_rank], block_distrib_ini_idx[i_rank+1],
                                          block_distrib_end_idx[i], block_distrib_end_idx[i+1]);
    btb->n_recv_buffer[i] = _overlap_size(block_distrib_end_idx[i_rank], block_distrib_end_idx[i_rank+1],
                                          block_distrib_ini_idx[i], block_distrib_ini_idx[i+1]);
  }


  return (PDM_block_to_block_t *) btb;

}

int
PDM_block_to_block_exch
(
 PDM_block_to_block_t *btb,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                  cst_stride,
 int                 *block_stride_ini,
 void                *block_data_ini,
 int                 *block_stride_end,
 void                **block_data_end
)
{
  _pdm_block_to_block_t *_btb = (_pdm_block_to_block_t *) btb;

  int *n_send_buffer = NULL;
  int *n_recv_buffer = NULL;
  int *i_send_buffer = NULL;
  int *i_recv_buffer = NULL;
  PDM_malloc(i_send_buffer, _btb->n_rank, int);
  PDM_malloc(i_recv_buffer, _btb->n_rank, int);

  int s_data_tot    = s_data;

  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {
    /*
    * First exchange for stride
    */
    _compute_displ(_btb->n_send_buffer, i_send_buffer, _btb->n_rank);
    _compute_displ(_btb->n_recv_buffer, i_recv_buffer, _btb->n_rank);

    PDM_MPI_Alltoallv(block_stride_ini,
                      _btb->n_send_buffer,
                      i_send_buffer,
                      PDM_MPI_INT,
                      block_stride_end,
                      _btb->n_recv_buffer,
                      i_recv_buffer,
                      PDM_MPI_INT,
                      _btb->comm);

    // Re count n_send / n_recv knowing stride data
    PDM_malloc(n_send_buffer, _btb->n_rank, int);
    PDM_malloc(n_recv_buffer, _btb->n_rank, int);
    int idx_send = 0;
    int idx_recv = 0;
    for (int i = 0; i < _btb->n_rank; ++i) {
      n_send_buffer[i] = 0;
      n_recv_buffer[i] = 0;
      for (int j = 0; j < _btb->n_send_buffer[i]; ++j) {
        n_send_buffer[i] += block_stride_ini[idx_send++];
      }
      for (int j = 0; j < _btb->n_recv_buffer[i]; ++j) {
        n_recv_buffer[i] += block_stride_end[idx_recv++];
      }
    }
  }

  else if (t_stride == PDM_STRIDE_CST_INTERLACED) {
    n_send_buffer = _btb->n_send_buffer;
    n_recv_buffer = _btb->n_recv_buffer;
    s_data_tot = s_data * cst_stride;
  }

  _compute_displ(n_send_buffer, i_send_buffer, _btb->n_rank);
  _compute_displ(n_recv_buffer, i_recv_buffer, _btb->n_rank);

  int s_recv_buffer = i_recv_buffer[_btb->n_rank-1] + n_recv_buffer[_btb->n_rank-1];

  size_t tot_size = ((size_t) s_data_tot ) * ((size_t) s_recv_buffer);
  unsigned char *recv_buffer = NULL;
  PDM_malloc(recv_buffer, tot_size, unsigned char);


  PDM_MPI_Datatype mpi_type;
  PDM_MPI_Type_create_contiguous(s_data_tot, PDM_MPI_BYTE, &mpi_type);
  PDM_MPI_Type_commit(&mpi_type);

  /*
   * Data exchange
   */
  unsigned char *send_buffer = (unsigned char *) block_data_ini;
  PDM_MPI_Alltoallv(send_buffer,
                    n_send_buffer,
                    i_send_buffer,
                    mpi_type,
                    recv_buffer,
                    n_recv_buffer,
                    i_recv_buffer,
                    mpi_type,
                    _btb->comm);

  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {
    PDM_free(n_send_buffer);
    PDM_free(n_recv_buffer);
  }
  PDM_free(i_send_buffer);
  PDM_free(i_recv_buffer);

  PDM_MPI_Type_free(&mpi_type);

  *block_data_end = recv_buffer;

  return s_recv_buffer/( (int)s_data );

}

int
PDM_block_to_block_exch_with_mpi_type
(
 PDM_block_to_block_t  *btb,
 PDM_stride_t           t_stride,
 PDM_MPI_Datatype       mpi_type,
 int                   *block_stride_ini,
 void                  *block_data_ini,
 int                   *block_stride_end,
 void                 **block_data_end
)
{
  _pdm_block_to_block_t *_btb = (_pdm_block_to_block_t *) btb;

  int *n_send_buffer = NULL;
  int *n_recv_buffer = NULL;
  int *i_send_buffer = NULL;
  int *i_recv_buffer = NULL;
  PDM_malloc(i_send_buffer, _btb->n_rank, int);
  PDM_malloc(i_recv_buffer, _btb->n_rank, int);


  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {
    /*
    * First exchange for stride
    */
    _compute_displ(_btb->n_send_buffer, i_send_buffer, _btb->n_rank);
    _compute_displ(_btb->n_recv_buffer, i_recv_buffer, _btb->n_rank);

    PDM_MPI_Alltoallv(block_stride_ini,
                      _btb->n_send_buffer,
                      i_send_buffer,
                      PDM_MPI_INT,
                      block_stride_end,
                      _btb->n_recv_buffer,
                      i_recv_buffer,
                      PDM_MPI_INT,
                      _btb->comm);

    // Re count n_send / n_recv knowing stride data
    PDM_malloc(n_send_buffer, _btb->n_rank, int);
    PDM_malloc(n_recv_buffer, _btb->n_rank, int);
    int idx_send = 0;
    int idx_recv = 0;
    for (int i = 0; i < _btb->n_rank; ++i) {
      n_send_buffer[i] = 0;
      n_recv_buffer[i] = 0;
      for (int j = 0; j < _btb->n_send_buffer[i]; ++j) {
        n_send_buffer[i] += block_stride_ini[idx_send++];
      }
      for (int j = 0; j < _btb->n_recv_buffer[i]; ++j){
        n_recv_buffer[i] += block_stride_end[idx_recv++];
      }
    }
  } else if (t_stride == PDM_STRIDE_CST_INTERLACED) {
    n_send_buffer = _btb->n_send_buffer;
    n_recv_buffer = _btb->n_recv_buffer;
  }

  _compute_displ(n_send_buffer, i_send_buffer, _btb->n_rank);
  _compute_displ(n_recv_buffer, i_recv_buffer, _btb->n_rank);

  int s_recv_buffer = i_recv_buffer[_btb->n_rank-1] + n_recv_buffer[_btb->n_rank-1];

  int s_type = 0;
  PDM_MPI_Type_size(mpi_type, &s_type);

  size_t tot_size = ((size_t) s_type ) * ((size_t) s_recv_buffer);
  unsigned char *recv_buffer;
  PDM_malloc(recv_buffer, tot_size, unsigned char);

  /*
   * Data exchange
   */
  unsigned char *send_buffer = (unsigned char *) block_data_ini;
  PDM_MPI_Alltoallv(send_buffer,
                    n_send_buffer,
                    i_send_buffer,
                    mpi_type,
                    recv_buffer,
                    n_recv_buffer,
                    i_recv_buffer,
                    mpi_type,
                    _btb->comm);

  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {
    PDM_free(n_send_buffer);
    PDM_free(n_recv_buffer);
  }
  PDM_free(i_send_buffer);
  PDM_free(i_recv_buffer);

  *block_data_end = recv_buffer;

  return s_recv_buffer/( (int)s_type );

}


PDM_block_to_block_t *
PDM_block_to_block_free
(
 PDM_block_to_block_t *btb
)
{
  _pdm_block_to_block_t *_btb = (_pdm_block_to_block_t *) btb;

  PDM_free(_btb->block_distrib_ini_idx);
  PDM_free(_btb->block_distrib_end_idx);
  PDM_free(_btb->n_send_buffer);
  PDM_free(_btb->n_recv_buffer);

  PDM_free(_btb);

  return NULL;
}

