/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_array.h"
#include "pdm_binary_search.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include "pdm_order.h"
#include "pdm_part_comm_graph.h"
#include "pdm_part_comm_graph_priv.h"
#include "pdm_exchange_helper_priv.h"
#include "pdm_priv.h"
#include "pdm_sort.h"
#include "pdm_unique.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

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

static
void
_allocate_send_buffer_strid_cst
(
  PDM_part_comm_graph_t   *pcg,
  size_t                   s_data,
  int                      cst_stride,
  unsigned char          **send_buffer
)
{
  int n_rank;
  PDM_MPI_Comm_size(pcg->comm, &n_rank);

  int send_buff_size = s_data * cst_stride * pcg->send_idx[n_rank];
  PDM_malloc(*send_buffer, send_buff_size, unsigned char);
}

static
void
_fill_send_strid_cst
(
 PDM_part_comm_graph_t   *pcg,
 size_t                   s_data,
 int                      cst_stride,
 void                   **send_entity_data,
 unsigned char           *send_buffer
)
{
  int s_data_tot = s_data * cst_stride;
  unsigned char **_send_entity_data = (unsigned char **) send_entity_data;

  for(int i_part = 0; i_part < pcg->n_part; ++i_part) {
    for(int i = 0; i < pcg->n_entity_graph[i_part]; ++i) {
      for(int k = 0; k < s_data_tot; ++k) {
        int idx_read  = i * s_data_tot + k;
        int idx_write = pcg->part_to_send_buffer[i_part][i] * s_data * cst_stride + k;
        send_buffer[idx_write] = _send_entity_data[i_part][idx_read];
      }
    }
  }
}

static
void
_prepare_send_strid_cst
(
  PDM_part_comm_graph_t   *pcg,
  size_t                   s_data,
  int                      cst_stride,
  void                   **send_entity_data,
  unsigned char          **send_buffer
)
{
  unsigned char *_send_buffer = NULL;
  _allocate_send_buffer_strid_cst(pcg, s_data, cst_stride, &_send_buffer);
  _fill_send_strid_cst(pcg, s_data, cst_stride, send_entity_data, _send_buffer);

  *send_buffer = _send_buffer;
}

static
void
_allocate_recv_buffer_strid_cst
(
  PDM_part_comm_graph_t   *pcg,
  size_t                   s_data,
  int                      cst_stride,
  unsigned char          **recv_buffer
)
{
  int n_rank;
  PDM_MPI_Comm_size(pcg->comm, &n_rank);

  int recv_buff_size = s_data * cst_stride * pcg->recv_idx[n_rank];
  PDM_malloc(*recv_buffer, recv_buff_size, unsigned char);
}


static
void
_allocate_recv_strid_cst
(
  PDM_part_comm_graph_t   *pcg,
  size_t                   s_data,
  int                      cst_stride,
  unsigned char         ***recv_entity_data

)
{
  int s_data_tot = s_data * cst_stride;

  unsigned char **_recv_entity_data = NULL;
  PDM_malloc(_recv_entity_data, pcg->n_part, unsigned char *);
  for(int i_part = 0; i_part < pcg->n_part; ++i_part) {
    _recv_entity_data[i_part] = malloc(pcg->n_entity_graph[i_part] * s_data_tot * sizeof(unsigned char));
  }

  *recv_entity_data = _recv_entity_data;
}

static
void
_fill_recv_strid_cst
(
  PDM_part_comm_graph_t   *pcg,
  size_t                   s_data,
  int                      cst_stride,
  unsigned char          **recv_entity_data,
  unsigned char           *recv_buffer
)
{
  int s_data_tot = s_data * cst_stride;
  for(int i_part = 0; i_part < pcg->n_part; ++i_part) {
    for(int i = 0; i < pcg->n_entity_graph[i_part]; ++i) {
      for(int k = 0; k < s_data_tot; ++k) {
        int idx_write = i * s_data_tot + k;
        int idx_read  = pcg->part_to_recv_buffer[i_part][i] * s_data * cst_stride + k;
        recv_entity_data[i_part][idx_write] = recv_buffer[idx_read];
      }
    }
  }
}

static
void
_post_recv_strid_cst
(
  PDM_part_comm_graph_t   *pcg,
  size_t                   s_data,
  int                      cst_stride,
  unsigned char           *recv_buffer,
  void                  ***recv_entity_data
)
{
  unsigned char **_recv_entity_data = NULL;
  _allocate_recv_strid_cst(pcg, s_data, cst_stride, &_recv_entity_data);
  _fill_recv_strid_cst    (pcg, s_data, cst_stride, _recv_entity_data, recv_buffer);

  *recv_entity_data = (void **) _recv_entity_data;
}

static
void
_exch_strid_cst
(
 PDM_part_comm_graph_t   *pcg,
 size_t                   s_data,
 int                      cst_stride,
 void                   **send_entity_data,
 void                  ***recv_entity_data
)
{
  int n_rank;
  PDM_MPI_Comm_size(pcg->comm, &n_rank);

  int s_data_tot = s_data * cst_stride;

  PDM_MPI_Datatype mpi_type;
  PDM_MPI_Type_create_contiguous(s_data_tot, PDM_MPI_BYTE, &mpi_type);
  PDM_MPI_Type_commit(&mpi_type);

  unsigned char *send_buffer = NULL;
  _prepare_send_strid_cst(pcg,
                          s_data,
                          cst_stride,
                          send_entity_data,
                          &send_buffer);

  unsigned char *recv_buffer = NULL;
  _allocate_recv_buffer_strid_cst(pcg, s_data, cst_stride, &recv_buffer);

  PDM_MPI_Alltoallv(send_buffer,
                    pcg->send_n,
                    pcg->send_idx,
                    mpi_type,
                    recv_buffer,
                    pcg->recv_n,
                    pcg->recv_idx,
                    mpi_type,
                    pcg->comm);
  PDM_free(send_buffer);

  /* Post-traitement */
  _post_recv_strid_cst(pcg,
                       s_data,
                       cst_stride,
                       recv_buffer,
                       recv_entity_data);

  PDM_free(recv_buffer);

  PDM_MPI_Type_free(&mpi_type);
}

static
void
_allocate_send_strid
(
  PDM_part_comm_graph_t   *pcg,
  int                    **send_stride
)
{
  int n_rank;
  PDM_MPI_Comm_size(pcg->comm, &n_rank);
  PDM_malloc(*send_stride, pcg->send_idx[n_rank], int);
}

static
void
_fill_send_strid
(
  PDM_part_comm_graph_t   *pcg,
  int                    **send_entity_stride,
  int                     *send_stride
)
{
  for(int i_part = 0; i_part < pcg->n_part; ++i_part) {
    for(int i = 0; i < pcg->n_entity_graph[i_part]; ++i) {
      int idx_write = pcg->part_to_send_buffer[i_part][i];
      send_stride[idx_write] = send_entity_stride[i_part][i];
    }
  }
}

static
void
_post_send_strid_and_fill_send_buffer
(
  PDM_part_comm_graph_t   *pcg,
  int                      s_data_tot,
  int                     *send_stride,
  int                     *recv_stride,
  int                    **send_entity_stride,
  void                   **send_entity_data,
  int                    **out_recv_stride_idx,
  int                    **out_send_data_idx,
  int                    **out_send_data_n,
  int                    **out_recv_data_idx,
  int                    **out_recv_data_n,
  unsigned char          **out_send_buffer
)
{
  int n_rank;
  PDM_MPI_Comm_size(pcg->comm, &n_rank);

  int *send_stride_idx = NULL;
  PDM_malloc(send_stride_idx, pcg->send_idx[n_rank]+1, int);
  send_stride_idx[0] = 0;
  for(int i = 0; i < pcg->send_idx[n_rank]; ++i) {
    send_stride_idx[i+1] = send_stride_idx[i] + send_stride[i];
  }

  int *recv_stride_idx = NULL;
  PDM_malloc(recv_stride_idx, pcg->recv_idx[n_rank]+1, int);
  recv_stride_idx[0] = 0;
  for(int i = 0; i < pcg->recv_idx[n_rank]; ++i) {
    recv_stride_idx[i+1] = recv_stride_idx[i] + recv_stride[i];
  }

  int *send_data_idx = NULL;
  PDM_malloc(send_data_idx, n_rank+1, int);
  int *send_data_n   = PDM_array_zeros_int(n_rank);
  send_data_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    send_data_idx[i+1] = send_data_idx[i];
    for(int j = pcg->send_idx[i]; j < pcg->send_idx[i+1]; ++j) {
      send_data_idx[i+1] += send_stride[j];
      send_data_n  [i  ] += send_stride[j];
      // send_stride[j] = 0;
    }
  }

  int *recv_data_idx = NULL;
  PDM_malloc(recv_data_idx, n_rank+1, int);
  int *recv_data_n   = PDM_array_zeros_int(n_rank);
  recv_data_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    recv_data_idx[i+1] = recv_data_idx[i];
    for(int j = pcg->recv_idx[i]; j < pcg->recv_idx[i+1]; ++j) {
      recv_data_idx[i+1] += recv_stride[j];
      recv_data_n  [i  ] += recv_stride[j];
    }
  }

  if(0 == 1) {
    PDM_log_trace_array_int(send_data_idx, n_rank+1, "send_data_idx ::");
    PDM_log_trace_array_int(recv_data_idx, n_rank+1, "recv_data_idx ::");
    PDM_log_trace_array_int(send_stride_idx, pcg->send_idx[n_rank]+1, "send_stride_idx ::");
    PDM_log_trace_array_int(recv_stride_idx, pcg->recv_idx[n_rank]+1, "recv_stride_idx ::");
  }

  int send_buff_size = send_data_idx[n_rank] * s_data_tot;
  unsigned char  *send_buffer       = malloc(send_buff_size * sizeof(unsigned char));
  unsigned char **_send_entity_data = (unsigned char **) send_entity_data;

  for(int i_part = 0; i_part < pcg->n_part; ++i_part) {
    int idx_read = 0;
    for(int i = 0; i < pcg->n_entity_graph[i_part]; ++i) {
      for(int j = 0; j < send_entity_stride[i_part][i]; ++j) {
        int idx_buffer = pcg->part_to_send_buffer[i_part][i];
        for(int k = 0; k < s_data_tot; ++k) {
          int idx_write  = (send_stride_idx[idx_buffer] + j) * s_data_tot + k;
          send_buffer[idx_write] = _send_entity_data[i_part][(idx_read+j)*s_data_tot + k];
        }
      }
      idx_read += send_entity_stride[i_part][i];
    }
  }
  PDM_free(send_stride_idx);

  *out_recv_stride_idx = recv_stride_idx;
  *out_send_data_idx   = send_data_idx;
  *out_send_data_n     = send_data_n;
  *out_recv_data_idx   = recv_data_idx;
  *out_recv_data_n     = recv_data_n;
  *out_send_buffer     = send_buffer;

}


static
void
_post_recv_strid_var
(
  PDM_part_comm_graph_t    *pcg,
  int                       s_data_tot,
  int                      *recv_stride,
  int                    ***recv_entity_stride,
  void                   ***recv_entity_data
)
{
  int           **_recv_entity_stride = NULL;
  unsigned char **_recv_entity_data   = NULL;
  PDM_malloc(_recv_entity_stride, pcg->n_part, int           *);
  PDM_malloc(_recv_entity_data  , pcg->n_part, unsigned char *);

  *recv_entity_stride =           _recv_entity_stride;
  *recv_entity_data   = (void **) _recv_entity_data;

  for(int i_part = 0; i_part < pcg->n_part; ++i_part) {
    _recv_entity_stride[i_part] = malloc(pcg->n_entity_graph[i_part] * sizeof(int));
    int recv_buff_size_part = 0;
    for(int i = 0; i < pcg->n_entity_graph[i_part]; ++i) {
      int idx_read  = pcg->part_to_recv_buffer[i_part][i];
      _recv_entity_stride[i_part][i] = recv_stride[idx_read];
      recv_buff_size_part += recv_stride[idx_read];
    }

    PDM_malloc(_recv_entity_data[i_part], recv_buff_size_part * s_data_tot, unsigned char);

    if(1 == 0) {
      PDM_log_trace_array_int(_recv_entity_stride[i_part], pcg->n_entity_graph[i_part], "_recv_entity_stride :");
    }
  }
}

static
void
_post_recv_buffer_strid_var
(
  PDM_part_comm_graph_t    *pcg,
  int                       s_data_tot,
  int                      *recv_stride_idx,
  unsigned char            *recv_buffer,
  int                     **recv_entity_stride,
  unsigned char           **recv_entity_data
)
{
  /*
   * Post-treatment buffer
   */
  for(int i_part = 0; i_part < pcg->n_part; ++i_part) {
    int idx_write = 0;
    for(int i = 0; i < pcg->n_entity_graph[i_part]; ++i) {
      int idx_buffer = pcg->part_to_recv_buffer[i_part][i];
      for(int j = 0; j < recv_entity_stride[i_part][i]; ++j) {
        for(int k = 0; k < s_data_tot; ++k) {
          int idx_read  = (recv_stride_idx[idx_buffer] + j) * s_data_tot + k;
          recv_entity_data[i_part][(idx_write+j)*s_data_tot + k] = recv_buffer[idx_read];
        }
      }
      idx_write += recv_entity_stride[i_part][i];
    }
  }
}


static
void
_exch_strid_var
(
 PDM_part_comm_graph_t   *pcg,
 size_t                   s_data,
 int                    **send_entity_stride,
 void                   **send_entity_data,
 int                   ***recv_entity_stride,
 void                  ***recv_entity_data
)
{
  int n_rank;
  PDM_MPI_Comm_size(pcg->comm, &n_rank);

  int s_data_tot = s_data;

  PDM_MPI_Datatype mpi_type;
  PDM_MPI_Type_create_contiguous(s_data_tot, PDM_MPI_BYTE, &mpi_type);
  PDM_MPI_Type_commit(&mpi_type);

  /* Exchange stride */
  int  *send_stride = NULL;
  _allocate_send_strid(pcg, &send_stride);
  _fill_send_strid    (pcg, send_entity_stride, send_stride);

  int *recv_stride = NULL;
  PDM_malloc(recv_stride, pcg->recv_idx[n_rank], int);
  PDM_MPI_Alltoallv(send_stride,
                    pcg->send_n,
                    pcg->send_idx,
                    PDM_MPI_INT,
                    recv_stride,
                    pcg->recv_n,
                    pcg->recv_idx,
                    PDM_MPI_INT,
                    pcg->comm);

  if(0 == 1) {
    PDM_log_trace_array_int(send_stride, pcg->send_idx[n_rank], "send_stride ::");
    PDM_log_trace_array_int(recv_stride, pcg->recv_idx[n_rank], "recv_stride ::");
  }

  /* Exchange data */
  int           *recv_stride_idx = NULL;
  int           *send_data_idx   = NULL;
  int           *send_data_n     = NULL;
  int           *recv_data_idx   = NULL;
  int           *recv_data_n     = NULL;
  unsigned char *send_buffer     = NULL;
  _post_send_strid_and_fill_send_buffer(pcg,
                                        s_data_tot,
                                        send_stride,
                                        recv_stride,
                                        send_entity_stride,
                                        send_entity_data,
                                        &recv_stride_idx,
                                        &send_data_idx,
                                        &send_data_n,
                                        &recv_data_idx,
                                        &recv_data_n,
                                        &send_buffer);

  int recv_buff_size = recv_data_idx[n_rank] * s_data_tot;
  unsigned char  *recv_buffer = NULL;
  PDM_malloc(recv_buffer, recv_buff_size, unsigned char);

  PDM_MPI_Alltoallv(send_buffer,
                    send_data_n,
                    send_data_idx,
                    mpi_type,
                    recv_buffer,
                    recv_data_n,
                    recv_data_idx,
                    mpi_type,
                    pcg->comm);
  PDM_free(send_buffer);
  PDM_free(send_data_idx);
  PDM_free(recv_data_idx);
  PDM_free(send_data_n);
  PDM_free(recv_data_n);
  PDM_free(send_stride);

  /* Panic verbose */
  // PDM_g_num_t* recv_buffer_dbg = (PDM_g_num_t *) recv_buffer;
  // PDM_log_trace_array_long(recv_buffer_dbg, recv_buff_size/s_data, "recv_buffer_dbg ::");

  /* Post-treatment stride + allocation */
  _post_recv_strid_var(pcg,
                       s_data_tot,
                       recv_stride,
                       recv_entity_stride,
                       recv_entity_data);
  PDM_free(recv_stride);

  /* Post-traitement stride */
  unsigned char **_recv_entity_data = (unsigned char **) *recv_entity_data;
  _post_recv_buffer_strid_var(pcg,
                              s_data_tot,
                              recv_stride_idx,
                              recv_buffer,
                              *recv_entity_stride,
                              _recv_entity_data);

  PDM_free(recv_buffer);
  PDM_free(recv_stride_idx);

  PDM_MPI_Type_free(&mpi_type);
}


/**
 *  \brief Compare lexicographically two nuplets of same size (in absolute value)
 *
 *  \return -1 if a > b, 0 if a == b, 1 if a < b
 */
static inline int
_compare_nuplets_abs
(
  const int size,
  const int a[],
  const int b[]
)
{
  for (int i = 0; i < size; i++) {
    if (PDM_ABS(a[i]) < PDM_ABS(b[i])) {
      return 1;
    }
    if (PDM_ABS(a[i]) > PDM_ABS(b[i])) {
      return -1;
    }
  }
  return 0;
}


static PDM_part_comm_graph_t *
_create
(
  int               n_part,
  int              *pn_entity_graph,
  int             **pentity_graph,
  PDM_ownership_t   owner_graph,
  int               nuplet_size,
  int             **pentity_nuplet,
  PDM_ownership_t   owner_nuplet,
  PDM_bool_t        is_signed,
  PDM_MPI_Comm      comm
)
{
  PDM_part_comm_graph_t *pcg = NULL;
  PDM_malloc(pcg, 1, PDM_part_comm_graph_t);

  pcg->comm   = comm;
  pcg->n_part = n_part;

  pcg->n_active_rank_send = 0;
  pcg->n_active_rank_recv = 0;

  pcg->is_signed = is_signed;

  int n_rank;
  int i_rank;
  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  pcg->owner_graph  = owner_graph;
  pcg->owner_nuplet = owner_nuplet;

  PDM_malloc(pcg->pentity_graph, n_part, int *);
  for (int i_part = 0; i_part < n_part; i_part++) {
    if(owner_graph == PDM_OWNERSHIP_BAD_VALUE) { // We need to copy
      PDM_malloc(pcg->pentity_graph[i_part], 4 * pn_entity_graph[i_part], int);
      for(int i = 0; i < 4 * pn_entity_graph[i_part]; ++i) {
        pcg->pentity_graph[i_part][i] = pentity_graph[i_part][i];
      }
      pcg->owner_graph = PDM_OWNERSHIP_KEEP;
    } else {
      pcg->pentity_graph[i_part] = pentity_graph[i_part];
    }
  }

  pcg->nuplet_size = nuplet_size;
  PDM_malloc(pcg->pentity_nuplet, n_part, int *);
  if(pentity_nuplet != NULL) {
    for (int i_part = 0; i_part < n_part; i_part++) {
      if(owner_nuplet == PDM_OWNERSHIP_BAD_VALUE) { // We need to copy
        PDM_malloc(pcg->pentity_nuplet[i_part], nuplet_size * pn_entity_graph[i_part], int);
        for(int i = 0; i < nuplet_size * pn_entity_graph[i_part]; ++i) {
          pcg->pentity_nuplet[i_part][i] = pentity_nuplet[i_part][i];
        }
        pcg->owner_nuplet = PDM_OWNERSHIP_KEEP;
      } else {
        pcg->pentity_nuplet[i_part] = pentity_nuplet[i_part];
      }
    }
  }

  int n_g_part = 0;
  PDM_MPI_Allreduce(&n_part, &n_g_part, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
  pcg->n_g_part = n_g_part;

  PDM_malloc(pcg->part_to_send_buffer, n_part, int *);
  PDM_malloc(pcg->part_to_recv_buffer, n_part, int *);
  PDM_malloc(pcg->n_entity_graph     , n_part, int  );

  int stride = 3 + nuplet_size;

  PDM_MPI_Datatype mpi_stride_type;
  PDM_MPI_Type_create_contiguous(stride, PDM_MPI_INT, &mpi_stride_type);
  PDM_MPI_Type_commit(&mpi_stride_type);

  int *send_n = PDM_array_zeros_int(n_rank);
  int *recv_n = NULL;
  PDM_malloc(recv_n, n_rank, int);
  for (int i_part = 0; i_part < n_part; i_part++) {
    int n_entity_graph = pn_entity_graph[i_part];
    pcg->n_entity_graph[i_part] = n_entity_graph;

    for(int idx_entity = 0; idx_entity < n_entity_graph; ++idx_entity) {
      int t_rank = pentity_graph[i_part][4*idx_entity+1];
      send_n[t_rank]++;
    }
  }

  int *send_idx = NULL;
  PDM_malloc(send_idx, n_rank+1, int);
  send_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    send_idx[i+1] = send_idx[i] + send_n[i];
    if(send_n[i] > 0) {
      pcg->n_active_rank_send++;
    }
    send_n[i] = 0;
  }
  int *send_buffer = NULL;
  PDM_malloc(send_buffer, stride * send_idx[n_rank], int);

  int sign = (is_signed == PDM_TRUE) ? -1 : 1;

  for (int i_part = 0; i_part < n_part; i_part++) {
    int n_entity_graph = pn_entity_graph[i_part];

    PDM_malloc(pcg->part_to_send_buffer[i_part], n_entity_graph, int);

    for(int idx_entity = 0; idx_entity < n_entity_graph; ++idx_entity) {
      int t_rank = pentity_graph[i_part][4*idx_entity+1];
      int idx_write = send_idx[t_rank] + send_n[t_rank]++;
      send_buffer[stride*idx_write  ] = pentity_graph[i_part][4*idx_entity+2]-1;
      send_buffer[stride*idx_write+1] = i_part;
      send_buffer[stride*idx_write+2] = PDM_ABS(pentity_graph[i_part][4*idx_entity+3])-1;
      for (int i = 0; i < nuplet_size; i++) {
        send_buffer[stride*idx_write+3+i] = sign * pentity_nuplet[i_part][nuplet_size*idx_entity+i];
      }

      pcg->part_to_send_buffer[i_part][idx_entity] = idx_write;
    }
  }

  PDM_MPI_Alltoall(send_n, 1, PDM_MPI_INT,
                   recv_n, 1, PDM_MPI_INT, comm);

  int *recv_idx = NULL;
  PDM_malloc(recv_idx, n_rank+1, int);
  recv_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    recv_idx[i+1] = recv_idx[i] + recv_n[i];
    if(recv_n[i] > 0) {
      pcg->n_active_rank_recv++;
    }
  }

  int *recv_buffer = NULL;
  PDM_malloc(recv_buffer, stride * recv_idx[n_rank], int);

  if(0 == 1) {
    PDM_log_trace_array_int(send_idx, n_rank+1, "send_idx ::");
    PDM_log_trace_array_int(recv_idx, n_rank+1, "recv_idx ::");
  }

  PDM_MPI_Alltoallv(send_buffer,
                    send_n,
                    send_idx,
                    mpi_stride_type,
                    recv_buffer,
                    recv_n,
                    recv_idx,
                    mpi_stride_type,
                    comm);
  PDM_free(send_buffer);

  /* Maintenant on cherche à retrouver la correspondance */
  int **pentity_indices       = NULL;
  int **pentity_indices_order = NULL;
  PDM_malloc(pentity_indices      , n_part, int *);
  PDM_malloc(pentity_indices_order, n_part, int *);

  for (int i_part = 0; i_part < n_part; i_part++) {
    int n_entity_graph = pn_entity_graph[i_part];
    PDM_malloc(pcg->part_to_recv_buffer[i_part], n_entity_graph, int);

    PDM_malloc(pentity_indices      [i_part], stride * n_entity_graph, int);
    PDM_malloc(pentity_indices_order[i_part],          n_entity_graph, int);

    for(int idx_entity = 0; idx_entity < n_entity_graph; ++idx_entity) {
      pentity_indices[i_part][stride*idx_entity  ] = pentity_graph[i_part][4*idx_entity+1];
      pentity_indices[i_part][stride*idx_entity+1] = pentity_graph[i_part][4*idx_entity+2]-1;
      pentity_indices[i_part][stride*idx_entity+2] = pentity_graph[i_part][4*idx_entity]-1; // Indices locaux
      for (int i = 0; i < nuplet_size; i++) {
        pentity_indices[i_part][stride*idx_entity+3+i] = pentity_nuplet[i_part][nuplet_size*idx_entity+i];
      }
    }

    PDM_order_lnum_s(pentity_indices[i_part],
                     stride,
                     pentity_indices_order[i_part],
                     n_entity_graph);

    PDM_order_array(n_entity_graph,
                    stride * sizeof(int),
                    pentity_indices_order[i_part],
                    pentity_indices      [i_part]);
  }

  /* Post buffer */
  int *to_find = NULL;
  PDM_malloc(to_find, stride, int);

  for(int t_rank = 0; t_rank < n_rank; ++t_rank) {
    for(int j = recv_idx[t_rank]; j < recv_idx[t_rank+1]; ++j) {

      int lpart   = recv_buffer[stride*j  ];
      int tpart   = recv_buffer[stride*j+1];
      int lentity = recv_buffer[stride*j+2];
      int n_entity_graph = pn_entity_graph[lpart];

      to_find[0] = t_rank;
      to_find[1] = tpart;
      to_find[2] = lentity;
      for (int i = 0; i < nuplet_size; i++) {
        to_find[3+i] = recv_buffer[stride*j+3+i];
      }

      int pos = PDM_order_binary_search_int(to_find, pentity_indices[lpart], stride, n_entity_graph);
      if (pos == -1) {
        // log_trace("Try to find fail = (%i/%i/%i) --> %i \n", t_rank, lpart, lentity, pos);
        PDM_log_trace_array_int(to_find, stride, "to_find : ");
        log_trace("pentity_indices :\n");
        for (int i = 0; i < n_entity_graph; i++) {
          log_trace("%d : ", i);
          PDM_log_trace_array_int(&pentity_indices[lpart][stride*i], stride, "");
        }
        PDM_error(__FILE__, __LINE__, 0, "Part-comm graph mismatch(see paradigm_*.log)\n");
      }
      pcg->part_to_recv_buffer[lpart][pentity_indices_order[lpart][pos]] = j;

    }
  }
  PDM_free(to_find);

  if(0 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int n_entity_graph = pn_entity_graph[i_part];
      PDM_log_trace_array_int(pcg->part_to_send_buffer[i_part], n_entity_graph, "pcg->part_to_send_buffer ::");
      PDM_log_trace_array_int(pcg->part_to_recv_buffer[i_part], n_entity_graph, "pcg->part_to_recv_buffer ::");
    }
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pentity_indices      [i_part]);
    PDM_free(pentity_indices_order[i_part]);
  }
  PDM_free(pentity_indices       );
  PDM_free(pentity_indices_order );


  PDM_free(recv_buffer);
  PDM_MPI_Type_free(&mpi_stride_type);

  pcg->send_idx = send_idx;
  pcg->recv_idx = recv_idx;
  pcg->send_n   = send_n;
  pcg->recv_n   = recv_n;

  /* Compute p2p array */
  PDM_malloc(pcg->active_rank_send, n_rank, int);
  PDM_malloc(pcg->active_rank_recv, pcg->n_active_rank_recv, int);

  PDM_malloc(pcg->active_send_idx, n_rank, int);
  PDM_malloc(pcg->active_send_n  , n_rank, int);
  PDM_malloc(pcg->active_recv_idx, n_rank, int);
  PDM_malloc(pcg->active_recv_n  , n_rank, int);

  pcg->n_active_rank_send = 0;
  pcg->n_active_rank_recv = 0;

  for(int i = 0; i < n_rank; ++i) {
    if(send_n[i] > 0) {
      pcg->active_send_idx [pcg->n_active_rank_send] = send_idx[i];
      pcg->active_send_n   [pcg->n_active_rank_send] = send_n  [i];
      pcg->active_rank_send[pcg->n_active_rank_send] = i;
      pcg->n_active_rank_send++;
    }
    if(recv_n[i] > 0) {
      pcg->active_recv_idx [pcg->n_active_rank_recv] = recv_idx[i];
      pcg->active_recv_n   [pcg->n_active_rank_recv] = recv_n  [i];
      pcg->active_rank_recv[pcg->n_active_rank_recv] = i;
      pcg->n_active_rank_recv++;
    }
  }

  PDM_realloc(pcg->active_rank_send, pcg->active_rank_send, pcg->n_active_rank_send, int);
  PDM_realloc(pcg->active_rank_recv, pcg->active_rank_recv, pcg->n_active_rank_recv, int);
  PDM_realloc(pcg->active_send_idx , pcg->active_send_idx , pcg->n_active_rank_send, int);
  PDM_realloc(pcg->active_send_n   , pcg->active_send_n   , pcg->n_active_rank_send, int);
  PDM_realloc(pcg->active_recv_idx , pcg->active_recv_idx , pcg->n_active_rank_recv, int);
  PDM_realloc(pcg->active_recv_n   , pcg->active_recv_n   , pcg->n_active_rank_recv, int);

  /*
   * Compute owner
   */
  pcg->bound_owner = NULL;
  PDM_malloc(pcg->bound_owner, n_part, int *);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    int n_entity_graph = pn_entity_graph[i_part];

    PDM_malloc(pcg->bound_owner[i_part], n_entity_graph, int);
    int *lbound_entity = NULL;
    PDM_malloc(lbound_entity, n_entity_graph, int);

    for(int idx_entity = 0; idx_entity < n_entity_graph; ++idx_entity) {
      lbound_entity[idx_entity] = pentity_graph[i_part][4*idx_entity];
      pcg->bound_owner[i_part][idx_entity] = -1;
    }

    int n_unique = PDM_inplace_unique(lbound_entity, 0, n_entity_graph-1);

    int *lowner = PDM_array_const_int(n_unique, -1);
    for(int idx_entity = 0; idx_entity < n_entity_graph; ++idx_entity) {
      int l_entity = pentity_graph[i_part][4*idx_entity];
      int pos = PDM_binary_search_int(l_entity, lbound_entity, n_unique);

      pcg->bound_owner[i_part][idx_entity] = pos; // Stockage temporaire

      int my_location[3] = {i_rank, i_part+1, l_entity};

      if(lowner[pos] != 0) {
        if (_compare_nuplets_abs(3, my_location, &pentity_graph[i_part][4*idx_entity+1]) >= 0) {
          lowner[pos] = 1;
        }
        else {
          lowner[pos] = 0;
        }
      }
    }

    /* Last loop to fill */
    for(int idx_entity = 0; idx_entity < n_entity_graph; ++idx_entity) {
      int pos = pcg->bound_owner[i_part][idx_entity];
      pcg->bound_owner[i_part][idx_entity] = lowner[pos];
    }

    PDM_free(lowner);
    PDM_free(lbound_entity);
  }

  /* Deleguate to exch_helper for Asynchronous and persitent exchange */
  pcg->exch_h = PDM_exchange_helper_create(comm);


  return pcg;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_part_comm_graph_t*
PDM_part_comm_graph_create
(
  int               n_part,
  int              *pn_entity_graph,
  int             **pentity_graph,
  PDM_ownership_t   ownership,
  PDM_MPI_Comm      comm
)
{
  return _create(n_part,
                 pn_entity_graph,
                 pentity_graph,
                 ownership,
                 0,
                 NULL,
                 PDM_OWNERSHIP_BAD_VALUE,
                 PDM_FALSE,
                 comm);
}


PDM_part_comm_graph_t*
PDM_part_comm_graph_with_nuplet_create
(
  int               n_part,
  int              *pn_entity_graph,
  int             **pentity_graph,
  PDM_ownership_t   owner_graph,
  int               nuplet_size,
  int             **pentity_nuplet,
  PDM_ownership_t   owner_nuplet,
  PDM_bool_t        is_signed,
  PDM_MPI_Comm      comm
)
{
  // Check coherence between ranks
  int min_nuplet_size;
  PDM_MPI_Allreduce(&nuplet_size, &min_nuplet_size, 1, PDM_MPI_INT, PDM_MPI_MIN, comm);

  int max_nuplet_size;
  PDM_MPI_Allreduce(&nuplet_size, &max_nuplet_size, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);

  if (min_nuplet_size != max_nuplet_size) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_comm_graph_create_with_nuplet : all ranks must have the same nuplet_size\n");
  }

  return _create(n_part,
                 pn_entity_graph,
                 pentity_graph,
                 owner_graph,
                 nuplet_size,
                 pentity_nuplet,
                 owner_nuplet,
                 is_signed,
                 comm);
}

void
PDM_part_comm_graph_exch
(
  PDM_part_comm_graph_t   *pcg,
  size_t                   s_data,
  PDM_stride_t             t_stride,
  int                      cst_stride,
  int                    **send_entity_stride,
  void                   **send_entity_data,
  int                   ***recv_entity_stride,
  void                  ***recv_entity_data
)
{
  if(t_stride == PDM_STRIDE_CST_INTERLACED) {
    _exch_strid_cst(pcg,
                    s_data,
                    cst_stride,
                    send_entity_data,
                    recv_entity_data);
  } else if (t_stride == PDM_STRIDE_VAR_INTERLACED) {
    _exch_strid_var(pcg,
                    s_data,
                    send_entity_stride,
                    send_entity_data,
                    recv_entity_stride,
                    recv_entity_data);
  } else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_comm_graph_exch, wrong t_stride \n");
  }

}


int
PDM_part_comm_graph_iexch
(
  PDM_part_comm_graph_t   *pcg,
  PDM_mpi_comm_kind_t      kcomm,
  size_t                   s_data,
  PDM_stride_t             t_stride,
  int                      cst_stride,
  int                    **send_entity_stride,
  void                   **send_entity_data,
  int                   ***recv_entity_stride,
  void                  ***recv_entity_data
)
{
  int n_rank;
  PDM_MPI_Comm_size(pcg->comm, &n_rank);

  int s_data_tot = s_data;
  if(t_stride == PDM_STRIDE_CST_INTERLACED) {
    s_data_tot = s_data * cst_stride;
  }

  PDM_MPI_Datatype mpi_type;
  PDM_MPI_Type_create_contiguous(s_data_tot, PDM_MPI_BYTE, &mpi_type);
  PDM_MPI_Type_commit(&mpi_type);

  int request_id = -1;

  if(t_stride == PDM_STRIDE_CST_INTERLACED) {
    unsigned char *send_buffer = NULL;
    _allocate_send_buffer_strid_cst(pcg, s_data, cst_stride, &send_buffer);

    unsigned char *recv_buffer = NULL;
    _allocate_recv_buffer_strid_cst(pcg, s_data, cst_stride, &recv_buffer);

    // Hook internal send_buffer et send_entity_data
    _fill_send_strid_cst(pcg,
                         s_data,
                         cst_stride,
                         send_entity_data,
                         send_buffer);

    request_id = PDM_exchange_helper_iexch(pcg->exch_h,
                                           kcomm,
                                           s_data,
                                           cst_stride,
                                           pcg->send_idx,
                                           pcg->send_n,
                                           send_buffer,
                                           pcg->recv_idx,
                                           pcg->recv_n,
                                           recv_buffer);

    pcg->exch_h->send_buffer  [request_id] = send_buffer;
    pcg->exch_h->recv_buffer  [request_id] = recv_buffer;

    unsigned char **_recv_entity_data = NULL;
    _allocate_recv_strid_cst(pcg, s_data, cst_stride, &_recv_entity_data);
    *recv_entity_data = (void **) _recv_entity_data;

    pcg->exch_h->t_stride     [request_id] = t_stride;
    pcg->exch_h->s_data       [request_id] = s_data;
    pcg->exch_h->cst_stride   [request_id] = cst_stride;
    pcg->exch_h->p_send_stride[request_id] = send_entity_stride;
    pcg->exch_h->p_send_data  [request_id] = send_entity_data;
    if(recv_entity_stride != NULL) {
      pcg->exch_h->p_recv_stride[request_id] = (*recv_entity_stride);
    }
    pcg->exch_h->p_recv_data  [request_id] = (*recv_entity_data);
  } else if (t_stride == PDM_STRIDE_VAR_INTERLACED) {

    int  *send_stride = NULL;
    _allocate_send_strid(pcg, &send_stride);
    _fill_send_strid    (pcg, send_entity_stride, send_stride);

    int *recv_stride = NULL;
    PDM_malloc(recv_stride, pcg->recv_idx[n_rank], int);

    /* Strid exchange : blocking */
    PDM_exchange_helper_mpi_type_exch(pcg->exch_h,
                                      kcomm,
                                      PDM_MPI_INT,
                                      pcg->send_idx,
                                      pcg->send_n,
                                      send_stride,
                                      pcg->recv_idx,
                                      pcg->recv_n,
                                      recv_stride);

    if(0 == 1) {
      PDM_log_trace_array_int(send_stride, pcg->send_idx[n_rank], "send_stride ::");
      PDM_log_trace_array_int(recv_stride, pcg->recv_idx[n_rank], "recv_stride ::");
    }

    /* Exchange data */
    int           *recv_stride_idx = NULL;
    int           *send_data_idx   = NULL;
    int           *send_data_n     = NULL;
    int           *recv_data_idx   = NULL;
    int           *recv_data_n     = NULL;
    unsigned char *send_buffer     = NULL;
    _post_send_strid_and_fill_send_buffer(pcg,
                                          s_data_tot,
                                          send_stride,
                                          recv_stride,
                                          send_entity_stride,
                                          send_entity_data,
                                          &recv_stride_idx,
                                          &send_data_idx,
                                          &send_data_n,
                                          &recv_data_idx,
                                          &recv_data_n,
                                          &send_buffer);
    PDM_free(send_stride);

    int recv_buff_size = recv_data_idx[n_rank] * s_data_tot;
    unsigned char  *recv_buffer = NULL;
    PDM_malloc(recv_buffer, recv_buff_size, unsigned char);

    request_id = PDM_exchange_helper_iexch(pcg->exch_h,
                                           kcomm,
                                           s_data,
                                           1, //
                                           send_data_idx,
                                           send_data_n,
                                           send_buffer,
                                           recv_data_idx,
                                           recv_data_n,
                                           recv_buffer);

    _post_recv_strid_var(pcg,
                         s_data_tot,
                         recv_stride,
                         recv_entity_stride,
                         recv_entity_data);
    PDM_free(recv_stride);

    pcg->exch_h->send_buffer    [request_id] = send_buffer;
    pcg->exch_h->recv_buffer    [request_id] = recv_buffer;

    pcg->exch_h->recv_stride_idx[request_id] = recv_stride_idx;
    pcg->exch_h->send_data_idx  [request_id] = send_data_idx;
    pcg->exch_h->send_data_n    [request_id] = send_data_n;
    pcg->exch_h->recv_data_idx  [request_id] = recv_data_idx;
    pcg->exch_h->recv_data_n    [request_id] = recv_data_n;

    pcg->exch_h->t_stride       [request_id] = t_stride;
    pcg->exch_h->s_data         [request_id] = s_data;
    pcg->exch_h->cst_stride     [request_id] = cst_stride;
    pcg->exch_h->p_send_stride  [request_id] = send_entity_stride;
    pcg->exch_h->p_send_data    [request_id] = send_entity_data;
    pcg->exch_h->p_recv_stride  [request_id] = (*recv_entity_stride);
    pcg->exch_h->p_recv_data    [request_id] = (*recv_entity_data);

  } else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_comm_graph_iexch, wrong t_stride \n");
  }

  PDM_MPI_Type_free(&mpi_type);

  return request_id;
}


int
PDM_part_comm_graph_exch_init
(
  PDM_part_comm_graph_t   *pcg,
  PDM_mpi_comm_kind_t      kcomm,
  size_t                   s_data,
  PDM_stride_t             t_stride,
  int                      cst_stride,
  int                    **send_entity_stride,
  void                   **send_entity_data,
  int                   ***recv_entity_stride,
  void                  ***recv_entity_data
)
{
  PDM_UNUSED(recv_entity_stride);

  int request_id = -1;
  if(t_stride == PDM_STRIDE_CST_INTERLACED) {
    unsigned char *send_buffer = NULL;
    _allocate_send_buffer_strid_cst(pcg, s_data, cst_stride, &send_buffer);

    unsigned char *recv_buffer = NULL;
    _allocate_recv_buffer_strid_cst(pcg, s_data, cst_stride, &recv_buffer);

    request_id = PDM_exchange_helper_exch_init(pcg->exch_h,
                                               kcomm,
                                               s_data,
                                               cst_stride,
                                               pcg->send_idx,
                                               pcg->send_n,
                                               send_buffer,
                                               pcg->recv_idx,
                                               pcg->recv_n,
                                               recv_buffer);

    pcg->exch_h->send_buffer  [request_id] = send_buffer;
    pcg->exch_h->recv_buffer  [request_id] = recv_buffer;

    unsigned char **_recv_entity_data = NULL;
    _allocate_recv_strid_cst(pcg, s_data, cst_stride, &_recv_entity_data);
    *recv_entity_data = (void **) _recv_entity_data;

    pcg->exch_h->t_stride     [request_id] = t_stride;
    pcg->exch_h->s_data       [request_id] = s_data;
    pcg->exch_h->cst_stride   [request_id] = cst_stride;
    pcg->exch_h->p_send_stride[request_id] = send_entity_stride;
    pcg->exch_h->p_send_data  [request_id] = send_entity_data;
    pcg->exch_h->p_recv_data  [request_id] = (*recv_entity_data);
  } else if (t_stride == PDM_STRIDE_VAR_INTERLACED) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_comm_graph_exch_init, not yet implemented for variable stride \n");
  } else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_comm_graph_iexch, wrong t_stride \n");
  }

  return request_id;
}

void
PDM_part_comm_graph_exch_start
(
  PDM_part_comm_graph_t   *pcg,
  int                      request_id
)
{
  if(pcg->exch_h->t_stride[request_id] == PDM_STRIDE_CST_INTERLACED) {
    // Hook internal send_buffer et send_entity_data
    _fill_send_strid_cst(pcg,
                         pcg->exch_h->s_data     [request_id],
                         pcg->exch_h->cst_stride [request_id],
                         pcg->exch_h->p_send_data[request_id],
                         pcg->exch_h->send_buffer[request_id]);
  }
  else if (pcg->exch_h->t_stride[request_id] == PDM_STRIDE_VAR_INTERLACED) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_comm_graph_exch_start, PDM_STRIDE_VAR_INTERLACED not implemented \n");
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_comm_graph_exch_start, wrong t_stride \n");
  }

  PDM_exchange_helper_exch_start(pcg->exch_h, request_id);
}


void
PDM_part_comm_graph_exch_wait
(
  PDM_part_comm_graph_t   *pcg,
  int                      request_id
)
{
  PDM_exchange_helper_exch_wait(pcg->exch_h, request_id);

  unsigned char **_precv_data = (unsigned char **) pcg->exch_h->p_recv_data[request_id];
  if(pcg->exch_h->t_stride[request_id] == PDM_STRIDE_CST_INTERLACED) {
    _fill_recv_strid_cst(pcg,
                         pcg->exch_h->s_data     [request_id],
                         pcg->exch_h->cst_stride [request_id],
                         _precv_data,
                         pcg->exch_h->recv_buffer[request_id]);
  }
  else if (pcg->exch_h->t_stride[request_id] == PDM_STRIDE_VAR_INTERLACED) {

    int **_precv_stri = (int **) pcg->exch_h->p_recv_stride[request_id];
    _post_recv_buffer_strid_var(pcg,
                                pcg->exch_h->s_data         [request_id],
                                pcg->exch_h->recv_stride_idx[request_id],
                                pcg->exch_h->recv_buffer    [request_id],
                                _precv_stri,
                                _precv_data);

    PDM_free(pcg->exch_h->recv_stride_idx[request_id]);
    PDM_free(pcg->exch_h->send_data_idx  [request_id]);
    PDM_free(pcg->exch_h->send_data_n    [request_id]);
    PDM_free(pcg->exch_h->recv_data_idx  [request_id]);
    PDM_free(pcg->exch_h->recv_data_n    [request_id]);

  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_comm_graph_exch_wait, wrong t_stride \n");
  }

  if(pcg->exch_h->is_persistent[request_id] == 0) {
    PDM_part_comm_graph_exch_free(pcg, request_id);
  }

}


void
PDM_part_comm_graph_exch_free
(
  PDM_part_comm_graph_t   *pcg,
  int                      request_id
)
{
  PDM_exchange_helper_exch_free(pcg->exch_h, request_id);

  PDM_free(pcg->exch_h->send_buffer[request_id]);
  PDM_free(pcg->exch_h->recv_buffer[request_id]);

}



int
PDM_part_comm_graph_exch_one_way_raw_init
(
  PDM_part_comm_graph_t      *pcg,
  PDM_exchange_direction_t    direction,
  size_t                      s_data,
  int                         cst_stride,
  int                        *raw_buffer,
  int                         tag
)
{
  int *send_or_recv_idx = NULL;
  int *send_or_recv_n   = NULL;
  int  n_active_rank    = 0;
  int *active_rank      = 0;

  if(direction == PDM_EXCHANGE_DIRECTION_SEND) {
    n_active_rank    = pcg->n_active_rank_send;
    active_rank      = pcg->active_rank_send;
    send_or_recv_idx = pcg->active_send_idx;
    send_or_recv_n   = pcg->active_send_n;
  } else if (direction == PDM_EXCHANGE_DIRECTION_RECV) {
    n_active_rank    = pcg->n_active_rank_recv;
    active_rank      = pcg->active_rank_recv;
    send_or_recv_idx = pcg->active_recv_idx;
    send_or_recv_n   = pcg->active_recv_n;
  } else {
    PDM_error(__FILE__, __LINE__, 0,
              "Error PDM_part_comm_graph_exch_one_way_raw_init not yet implemented with direction = %i\n", direction);
  }

  int request_id = PDM_exchange_helper_exch_one_way_init(pcg->exch_h,
                                                         direction,
                                                         s_data,
                                                         cst_stride,
                                                         n_active_rank,
                                                         active_rank,
                                                         send_or_recv_idx,
                                                         send_or_recv_n,
                                                         tag,
                                                         raw_buffer);

  pcg->exch_h->send_buffer  [request_id] = raw_buffer;
  pcg->exch_h->s_data       [request_id] = s_data;
  pcg->exch_h->cst_stride   [request_id] = cst_stride;

  return request_id;
}


int
PDM_part_comm_graph_iexch_one_way_raw
(
  PDM_part_comm_graph_t      *pcg,
  PDM_exchange_direction_t    direction,
  size_t                      s_data,
  int                         cst_stride,
  int                        *raw_buffer,
  int                         tag
)
{
  int *send_or_recv_idx = NULL;
  int *send_or_recv_n   = NULL;
  int  n_active_rank    = 0;
  int *active_rank      = 0;

  if(direction == PDM_EXCHANGE_DIRECTION_SEND) {
    n_active_rank    = pcg->n_active_rank_send;
    active_rank      = pcg->active_rank_send;
    send_or_recv_idx = pcg->active_send_idx;
    send_or_recv_n   = pcg->active_send_n;
  } else if (direction == PDM_EXCHANGE_DIRECTION_RECV) {
    n_active_rank    = pcg->n_active_rank_recv;
    active_rank      = pcg->active_rank_recv;
    send_or_recv_idx = pcg->active_recv_idx;
    send_or_recv_n   = pcg->active_recv_n;
  } else {
    PDM_error(__FILE__, __LINE__, 0,
              "Error PDM_part_comm_graph_iexch_one_way_raw not yet implemented with direction = %i\n", direction);
  }

  int request_id = PDM_exchange_helper_iexch_one_way(pcg->exch_h,
                                                     direction,
                                                     s_data,
                                                     cst_stride,
                                                     n_active_rank,
                                                     active_rank,
                                                     send_or_recv_idx,
                                                     send_or_recv_n,
                                                     tag,
                                                     raw_buffer);

  pcg->exch_h->send_buffer  [request_id] = raw_buffer;
  pcg->exch_h->s_data       [request_id] = s_data;
  pcg->exch_h->cst_stride   [request_id] = cst_stride;

  return request_id;
}


void
PDM_part_comm_graph_exch_one_way_raw_start
(
  PDM_part_comm_graph_t      *pcg,
  int                         request_id
)
{
  PDM_exchange_helper_exch_start(pcg->exch_h, request_id);
}

void
PDM_part_comm_graph_exch_one_way_raw_wait
(
  PDM_part_comm_graph_t   *pcg,
  int                      request_id
)
{
  PDM_exchange_helper_exch_wait(pcg->exch_h, request_id);
}

void
PDM_part_comm_graph_exch_one_way_raw_free
(
  PDM_part_comm_graph_t      *pcg,
  int                         request_id
)
{
  pcg->exch_h->send_buffer  [request_id] = NULL;

  PDM_exchange_helper_exch_free(pcg->exch_h, request_id);
}


const int*
PDM_part_comm_graph_owner_get
(
  PDM_part_comm_graph_t *pcg,
  int                    i_part
)
{
  return pcg->bound_owner[i_part];
}

void
PDM_part_comm_graph_reorder
(
  PDM_part_comm_graph_t  *pcg,
  int                   **old_to_new
)
{
  int **pentity_graph = pcg->pentity_graph;

  /* Prepare exchange */
  int **send_new_id = NULL;
  PDM_malloc(send_new_id, pcg->n_part, int *);
  for(int i_part = 0; i_part < pcg->n_part; ++i_part) {
    PDM_malloc(send_new_id[i_part], pcg->n_entity_graph[i_part], int);
    for(int i = 0; i < pcg->n_entity_graph[i_part]; ++i) {
      int i_entity = pentity_graph[i_part][4*i  ]-1;
      send_new_id[i_part][i] = old_to_new[i_part][i_entity];
    }
  }

  int **recv_new_id = NULL;
  PDM_part_comm_graph_exch(pcg,
                           sizeof(int),
                           PDM_STRIDE_CST_INTERLACED,
                           1,
                           NULL,
              (void  **)   send_new_id,
                           NULL,
              (void ***)   &recv_new_id);

  for(int i_part = 0; i_part < pcg->n_part; ++i_part) {
    PDM_free(send_new_id[i_part]);
  }
  PDM_free(send_new_id);

  /* Actualisation current and opposite */

  for(int i_part = 0; i_part < pcg->n_part; ++i_part) {
    for(int i = 0; i < pcg->n_entity_graph[i_part]; ++i) {
      int i_entity = pentity_graph[i_part][4*i  ]-1;
      pentity_graph[i_part][4*i  ] = old_to_new[i_part][i_entity]+1; // On suppose que old_to_new commence a 0
      pentity_graph[i_part][4*i+3] = recv_new_id[i_part][i]+1; // On suppose que old_to_new commence a 0
    }
  }

  for(int i_part = 0; i_part < pcg->n_part; ++i_part) {
    PDM_free(recv_new_id[i_part]);
  }
  PDM_free(recv_new_id);
}


int
PDM_part_comm_graph_entity_graph_get
(
  PDM_part_comm_graph_t  *pcg,
  int                     i_part,
  int                   **entity_graph,
  PDM_ownership_t         ownership
)
{
  if (pcg == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_comm_graph_entity_graph_get : Invalid PDM_part_comm_graph_t instance\n");
  }

  if (i_part < 0 || i_part >= pcg->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_comm_graph_entity_graph_get : Invalid i_part (%d / %d)\n", i_part, pcg->n_part);
  }

  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    pcg->owner_graph = ownership;
  }

  *entity_graph = pcg->pentity_graph[i_part];

  return pcg->n_entity_graph[i_part];
}


int
PDM_part_comm_graph_entity_nuplet_get
(
  PDM_part_comm_graph_t  *pcg,
  int                     i_part,
  int                   **entity_nuplet,
  PDM_ownership_t         ownership
)
{
  if (pcg == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_comm_graph_entity_nuplet_get : Invalid PDM_part_comm_graph_t instance\n");
  }

  if (i_part < 0 || i_part >= pcg->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_comm_graph_entity_nuplet_get : Invalid i_part (%d / %d)\n", i_part, pcg->n_part);
  }

  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    pcg->owner_nuplet = ownership;
  }

  if (pcg->pentity_nuplet != NULL) {
    *entity_nuplet = pcg->pentity_nuplet[i_part];
  }
  else {
    *entity_nuplet = NULL;
  }

  return pcg->nuplet_size;
}


void
PDM_part_comm_graph_gather_strided_data
(
  PDM_part_comm_graph_t   *pcg,
  const size_t             size_data,
  PDM_stride_t             t_stride,
  int                     *n_entity,
  int                    **data_stride,
  void                   **data,
  int                   ***out_data_stride,
  void                  ***out_data
)
{
  if(t_stride == PDM_STRIDE_VAR_INTERLACED) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_comm_graph_gather_strided_data: PDM_STRIDE_VAR_INTERLACED not implemented \n");
  }
  else if(t_stride != PDM_STRIDE_CST_INTERLACED) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_comm_graph_gather_strided_data: wrong t_stride\n");
  }

  int n_part = pcg->n_part;

  int           **data_idx        = NULL;
  int            *pn_entity_bound = NULL;
  int           **pentity_bound   = NULL;
  int           **send_data_n     = NULL;
  unsigned char **send_data       = NULL;
  PDM_malloc(data_idx       , n_part, int           *);
  PDM_malloc(pn_entity_bound, n_part, int            );
  PDM_malloc(pentity_bound  , n_part, int           *);
  PDM_malloc(send_data_n    , n_part, int           *);
  PDM_malloc(send_data      , n_part, unsigned char *);

  for (int i_part=0; i_part<n_part; ++i_part) {

    data_idx[i_part] = PDM_array_new_idx_from_sizes_int(data_stride[i_part], n_entity[i_part]);

    pn_entity_bound[i_part] = PDM_part_comm_graph_entity_graph_get(pcg,
                                                                   i_part,
                                                                  &pentity_bound[i_part],
                                                                   PDM_OWNERSHIP_BAD_VALUE);

    PDM_malloc(send_data_n[i_part], pn_entity_bound[i_part], int);
    int send_data_size = 0;
    for(int i = 0; i < pn_entity_bound[i_part]; ++i) {
      int i_entity = pentity_bound[i_part][4*i]-1;
      int n_data = data_stride[i_part][i_entity];
      send_data_n[i_part][i] = n_data;
      send_data_size += send_data_n[i_part][i];
    }

    PDM_malloc(send_data[i_part], send_data_size*size_data, unsigned char);
    unsigned char *_data = (unsigned char* ) data[i_part];

    int i_write = 0;
    for(int i = 0; i < pn_entity_bound[i_part]; ++i) {
      int i_entity = pentity_bound[i_part][4*i]-1;
      for(int k = data_idx[i_part][i_entity]; k < data_idx[i_part][i_entity+1]; ++k) {
        for (int octet = 0; octet <  (int) size_data; ++octet) {
          send_data[i_part][i_write++] = _data[size_data*k + octet];
        }
      }
    }
  }


  int  **recv_data_n = NULL;
  void **recv_data   = NULL;
  PDM_part_comm_graph_exch(pcg,
                           size_data,
                           PDM_STRIDE_VAR_INTERLACED,
                           1,
                           send_data_n,
               (void  **)  send_data,
                          &recv_data_n,
               (void ***) &recv_data);

  for (int i_part=0; i_part<n_part; ++i_part) {
    PDM_free(send_data_n[i_part]);
    PDM_free(send_data  [i_part]);
  }
  PDM_free(send_data_n);
  PDM_free(send_data);


  int           **_out_data_n   = NULL;
  int           **_out_data_idx = NULL;
  unsigned char **_out_data     = NULL;
  PDM_malloc(_out_data_n  , n_part, int           *);
  PDM_malloc(_out_data_idx, n_part, int           *);
  PDM_malloc(_out_data    , n_part, unsigned char *);
  for (int i_part=0; i_part<n_part; ++i_part) {

    /**
     * Count local + rcvd data
     */

    PDM_calloc(_out_data_n  [i_part], n_entity[i_part]  , int);
    PDM_malloc(_out_data_idx[i_part], n_entity[i_part]+1, int); _out_data_idx[i_part][0] = 0;

    int data_size = data_idx[i_part][n_entity[i_part]];
    for(int i = 0; i < pn_entity_bound[i_part]; ++i) {
      int i_entity = pentity_bound[i_part][4*i]-1;
      _out_data_n[i_part][i_entity] += recv_data_n[i_part][i];
      data_size                     += recv_data_n[i_part][i];
    }

    for(int i_entity = 0; i_entity < n_entity[i_part]; ++i_entity) {
      int n_local = data_idx[i_part][i_entity+1]-data_idx[i_part][i_entity];
      _out_data_idx[i_part][i_entity+1] = _out_data_idx[i_part][i_entity] + n_local + _out_data_n[i_part][i_entity];
      _out_data_n  [i_part][i_entity  ] = 0;
    }


    /**
     * Fill out_data with local and received datas
     */
    PDM_malloc(_out_data[i_part], data_size*size_data, unsigned char);
    unsigned char *_data = (unsigned char *) data[i_part];

    int i_write = 0;
    for(int i_entity = 0; i_entity < n_entity[i_part]; ++i_entity) {
      for(int i_read=data_idx[i_part][i_entity];
              i_read<data_idx[i_part][i_entity+1]; ++i_read) {
        for (int octet = 0; octet <  (int) size_data; ++octet) {
          i_write = _out_data_idx[i_part][i_entity] + _out_data_n[i_part][i_entity];
          i_write = size_data*i_write + octet;
          _out_data[i_part][i_write] = _data[size_data*i_read + octet];
        }
        _out_data_n[i_part][i_entity]++;
      }
    }

    int i_readr = 0;
    unsigned char *_recv_data = (unsigned char *) recv_data[i_part];

    for(int i = 0; i < pn_entity_bound[i_part]; ++i) {
      int i_entity = pentity_bound[i_part][4*i]-1;

      for(int k=0; k<recv_data_n[i_part][i]; ++k) {
        for (int octet = 0; octet <  (int) size_data; ++octet) {
          i_write = _out_data_idx[i_part][i_entity] + _out_data_n[i_part][i_entity];
          i_write = size_data*i_write + octet;
          _out_data  [i_part][i_write] = _recv_data[i_readr++];
        }
        _out_data_n[i_part][i_entity]++;
      }
    }

    PDM_free(data_idx     [i_part]);
    PDM_free(recv_data_n  [i_part]);
    PDM_free(recv_data    [i_part]);
    PDM_free(_out_data_idx[i_part]);
  }
  PDM_free(data_idx);
  PDM_free(recv_data_n);
  PDM_free(recv_data);
  PDM_free(_out_data_idx);
  PDM_free(pentity_bound);
  PDM_free(pn_entity_bound);

  *out_data_stride =           _out_data_n;
  *out_data        = (void **) _out_data;
}


void
PDM_part_comm_graph_all_reduce
(
  PDM_part_comm_graph_t   *pcg,
  PDM_MPI_Datatype         datatype,
  int                      stride,
  PDM_MPI_Op               op,
  unsigned char          **pdata
)
{
  int n_part = pcg->n_part;

  if(op != PDM_MPI_SUM &&
     op != PDM_MPI_MIN &&
     op != PDM_MPI_MAX) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_comm_graph_all_reduce only available with op = PDM_MPI_SUM/PDM_MPI_MIN/PDM_MPI_MAX\n");
  }

  unsigned char **send_data = NULL;
  unsigned char **recv_data = NULL;
  PDM_malloc(send_data, n_part, unsigned char *);

  int s_data = 0;
  PDM_MPI_Type_size(datatype, &s_data);

  int _stride = (stride >= 0) ? stride : 1;
  s_data *= _stride;

  for (int i_part = 0; i_part < n_part; ++i_part) {

    PDM_malloc(send_data[i_part], pcg->n_entity_graph[i_part] * s_data, unsigned char);

    for(int i = 0; i < pcg->n_entity_graph[i_part]; ++i) {
      int i_entity = pcg->pentity_graph[i_part][4*i  ]-1;
      for(int k = 0; k < s_data; ++k) {
        send_data[i_part][s_data * i + k] = pdata[i_part][s_data * i_entity + k];
      }
    }
  }

  // Exchange data
  PDM_part_comm_graph_exch(pcg,
                           s_data,
                           PDM_STRIDE_CST_INTERLACED,
                           1,
                           NULL,
              (void **)    send_data,
                           NULL,
              (void ***)   &recv_data);

  // Reduce
  if(datatype == PDM_MPI_DOUBLE) {
    double **_pdata     = (double **) pdata;
    double **_recv_data = (double **) recv_data;
    for (int i_part = 0; i_part < n_part; ++i_part) {
      if(op == PDM_MPI_SUM) {
        // We suppose that current value already init
        for(int i = 0; i < pcg->n_entity_graph[i_part]; ++i) {
          int i_entity = pcg->pentity_graph[i_part][4*i]-1;
          for (int j = 0; j < _stride; j++) {
            _pdata[i_part][_stride*i_entity+j] += _recv_data[i_part][_stride*i+j];
          }
        }
      } else if (op == PDM_MPI_MAX) {
        for(int i = 0; i < pcg->n_entity_graph[i_part]; ++i) {
          int i_entity = pcg->pentity_graph[i_part][4*i]-1;
          for (int j = 0; j < _stride; j++) {
            _pdata[i_part][_stride*i_entity+j] = PDM_MAX(_pdata    [i_part][_stride*i_entity+j],
                                                         _recv_data[i_part][_stride*i       +j]);
          }
        }
      } else if (op == PDM_MPI_MIN) {
        for(int i = 0; i < pcg->n_entity_graph[i_part]; ++i) {
          int i_entity = pcg->pentity_graph[i_part][4*i]-1;
          for (int j = 0; j < _stride; j++) {
            _pdata[i_part][_stride*i_entity+j] = PDM_MIN(_pdata    [i_part][_stride*i_entity+j],
                                                         _recv_data[i_part][_stride*i       +j]);
          }
        }
      }
    }
  }
  else if (datatype == PDM_MPI_INT) {
    int **_pdata     = (int **) pdata;
    int **_recv_data = (int **) recv_data;
    for (int i_part = 0; i_part < n_part; ++i_part) {
      if(op == PDM_MPI_SUM) {
        // We suppose that current value already init
        for(int i = 0; i < pcg->n_entity_graph[i_part]; ++i) {
          int i_entity = pcg->pentity_graph[i_part][4*i]-1;
          for (int j = 0; j < _stride; j++) {
            _pdata[i_part][_stride*i_entity+j] += _recv_data[i_part][_stride*i+j];
          }
        }
      } else if (op == PDM_MPI_MAX) {
        for(int i = 0; i < pcg->n_entity_graph[i_part]; ++i) {
          int i_entity = pcg->pentity_graph[i_part][4*i]-1;
          for (int j = 0; j < _stride; j++) {
            _pdata[i_part][_stride*i_entity+j] = PDM_MAX(_pdata    [i_part][_stride*i_entity+j],
                                                         _recv_data[i_part][_stride*i       +j]);
          }
        }
      } else if (op == PDM_MPI_MIN) {
        for(int i = 0; i < pcg->n_entity_graph[i_part]; ++i) {
          int i_entity = pcg->pentity_graph[i_part][4*i]-1;
          for (int j = 0; j < _stride; j++) {
            _pdata[i_part][_stride*i_entity+j] = PDM_MIN(_pdata    [i_part][_stride*i_entity+j],
                                                         _recv_data[i_part][_stride*i       +j]);
          }
        }
      }
    }
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_comm_graph_all_reduce only available for PDM_MPI_DOUBLE / PDM_MPI_INT \n");
  }

  for (int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(send_data[i_part]);
    PDM_free(recv_data[i_part]);
  }
  PDM_free(send_data);
  PDM_free(recv_data);

}


void
PDM_part_comm_graph_part_to_send_buffer_get
(
  PDM_part_comm_graph_t   *pcg,
  int                   ***out_part_to_send_buffer
)
{
  *out_part_to_send_buffer = pcg->part_to_send_buffer;
}


void
PDM_part_comm_graph_part_to_recv_buffer_get
(
  PDM_part_comm_graph_t   *pcg,
  int                   ***part_to_recv_buffer
)
{
  *part_to_recv_buffer = pcg->part_to_recv_buffer;
}

void
PDM_part_comm_graph_free
(
 PDM_part_comm_graph_t* pcg
)
{
  if (pcg == NULL) {
    return;
  }

  for(int i_part = 0; i_part < pcg->n_part; ++i_part) {
    PDM_free(pcg->part_to_send_buffer[i_part]);
    PDM_free(pcg->part_to_recv_buffer[i_part]);
    PDM_free(pcg->bound_owner        [i_part]);
  }
  PDM_free(pcg->part_to_send_buffer);
  PDM_free(pcg->part_to_recv_buffer);
  PDM_free(pcg->n_entity_graph);
  PDM_free(pcg->bound_owner);

  PDM_free(pcg->send_idx);
  PDM_free(pcg->recv_idx);
  PDM_free(pcg->send_n);
  PDM_free(pcg->recv_n);
  PDM_free(pcg->active_rank_send);
  PDM_free(pcg->active_rank_recv);
  PDM_free(pcg->active_send_idx);
  PDM_free(pcg->active_send_n  );
  PDM_free(pcg->active_recv_idx);
  PDM_free(pcg->active_recv_n  );

  if (pcg->owner_graph == PDM_OWNERSHIP_KEEP) {
    for (int i_part = 0; i_part < pcg->n_part; i_part++) {
      PDM_free(pcg->pentity_graph[i_part]);
    }
  }
  PDM_free(pcg->pentity_graph);

  if (pcg->owner_nuplet == PDM_OWNERSHIP_KEEP) {
    for (int i_part = 0; i_part < pcg->n_part; i_part++) {
      PDM_free(pcg->pentity_nuplet[i_part]);
    }
  }
  PDM_free(pcg->pentity_nuplet);
  PDM_exchange_helper_free(pcg->exch_h);

  PDM_free(pcg);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
