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
_prepare_send_strid_cst
(
 PDM_part_comm_graph_t   *pcg,
 size_t                   s_data,
 int                      cst_stride,
 void                   **send_entity_data,
 unsigned char          **send_buffer
)
{

  int n_rank;
  PDM_MPI_Comm_size(pcg->comm, &n_rank);

  int s_data_tot = s_data * cst_stride;

  int send_buff_size = s_data * cst_stride * pcg->send_idx[n_rank];
  unsigned char *_send_buffer = NULL;
  PDM_malloc(_send_buffer, send_buff_size, unsigned char);

  unsigned char **_send_entity_data = (unsigned char **) send_entity_data;

  for(int i_part = 0; i_part < pcg->n_part; ++i_part) {
    for(int i = 0; i < pcg->n_entity_graph[i_part]; ++i) {
      for(int k = 0; k < s_data_tot; ++k) {
        int idx_read  = i * s_data_tot + k;
        int idx_write = pcg->part_to_send_buffer[i_part][i] * s_data * cst_stride + k;
        _send_buffer[idx_write] = _send_entity_data[i_part][idx_read];
      }
    }
  }
  *send_buffer = _send_buffer;
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
  int s_data_tot = s_data * cst_stride;

  unsigned char **_recv_entity_data = NULL;
  PDM_malloc(_recv_entity_data, pcg->n_part, unsigned char *);
  *recv_entity_data = (void **) _recv_entity_data;
  for(int i_part = 0; i_part < pcg->n_part; ++i_part) {
    _recv_entity_data[i_part] = malloc(pcg->n_entity_graph[i_part] * s_data_tot * sizeof(unsigned char));
    for(int i = 0; i < pcg->n_entity_graph[i_part]; ++i) {
      for(int k = 0; k < s_data_tot; ++k) {
        int idx_write = i * s_data_tot + k;
        int idx_read  = pcg->part_to_recv_buffer[i_part][i] * s_data * cst_stride + k;
        _recv_entity_data[i_part][idx_write] = recv_buffer[idx_read];
      }
    }
  }
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

  int recv_buff_size = s_data * cst_stride * pcg->recv_idx[n_rank];
  unsigned char *recv_buffer = NULL;
  PDM_malloc(recv_buffer, recv_buff_size, unsigned char);
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
  PDM_malloc(send_stride, pcg->send_idx[n_rank], int);
  for(int i_part = 0; i_part < pcg->n_part; ++i_part) {
    for(int i = 0; i < pcg->n_entity_graph[i_part]; ++i) {
      int idx_write = pcg->part_to_send_buffer[i_part][i];
      send_stride[idx_write] = send_entity_stride[i_part][i];
    }
  }

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

  /* Panic verbose */
  // PDM_g_num_t* recv_buffer_dbg = (PDM_g_num_t *) recv_buffer;
  // PDM_log_trace_array_long(recv_buffer_dbg, recv_buff_size/s_data, "recv_buffer_dbg ::");

  /* Post-traitement stride */
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

  /*
   * Post-treatment buffer
   */
  for(int i_part = 0; i_part < pcg->n_part; ++i_part) {
    int idx_write = 0;
    for(int i = 0; i < pcg->n_entity_graph[i_part]; ++i) {
      int idx_buffer = pcg->part_to_recv_buffer[i_part][i];
      for(int j = 0; j < _recv_entity_stride[i_part][i]; ++j) {
        for(int k = 0; k < s_data_tot; ++k) {
          int idx_read  = (recv_stride_idx[idx_buffer] + j) * s_data_tot + k;
          _recv_entity_data[i_part][(idx_write+j)*s_data_tot + k] = recv_buffer[idx_read];
        }
      }
      idx_write += _recv_entity_stride[i_part][i];
    }
  }

  PDM_free(recv_buffer);
  PDM_free(recv_stride_idx);

  PDM_free(send_data_idx);
  PDM_free(recv_data_idx);
  PDM_free(send_data_n);
  PDM_free(recv_data_n);

  PDM_free(send_stride);
  PDM_free(recv_stride);

  PDM_MPI_Type_free(&mpi_type);
}

/**
 *  \brief Compare lexicographically two nuplets of same size
 *
 *  \return -1 if a > b, 0 if a == b, 1 if a < b
 */
static inline int
_compare_nuplets
(
  const int size,
  const int a[],
  const int b[]
)
{
  for (int i = 0; i < size; i++) {
    if (a[i] < b[i]) {
      return 1;
    }
    if (a[i] > b[i]) {
      return -1;
    }
  }
  return 0;
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

  pcg->nuplet_size    = nuplet_size;
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
  PDM_malloc(pcg->active_rank_send, pcg->n_active_rank_send, int);
  PDM_malloc(pcg->active_rank_recv, pcg->n_active_rank_recv, int);

  pcg->n_active_rank_send = 0;
  pcg->n_active_rank_recv = 0;

  for(int i = 0; i < n_rank; ++i) {
    if(send_n[i] > 0) {
      pcg->active_rank_send[pcg->n_active_rank_send++] = i;
    }
    if(recv_n[i] > 0) {
      pcg->active_rank_recv[pcg->n_active_rank_recv++] = i;
    }
  }

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


void
PDM_part_comm_graph_entity1_to_part_comm_graph_entity2
(
  PDM_part_comm_graph_t   *ptpgc_entity1,
  int                     *pn_entity1,
  int                     *pn_entity2,
  int                    **entity2_entity1_idx,
  int                    **entity2_entity1,
  PDM_part_comm_graph_t  **out_ptpgc_entity2
)
{

  int  *pn_entity2_graph = NULL;
  int **pentity2_graph   = NULL;
  int **pentity2_nuplet  = NULL;

  PDM_part_comm_graph_entity1_to_entity2(ptpgc_entity1->comm,
                                         ptpgc_entity1->n_part,
                                         ptpgc_entity1->n_entity_graph,
                                         ptpgc_entity1->pentity_graph,
                                         ptpgc_entity1->nuplet_size,
                                         ptpgc_entity1->pentity_nuplet,
                                         pn_entity1,
                                         pn_entity2,
                                         entity2_entity1_idx,
                                         entity2_entity1,
                                         &pn_entity2_graph,
                                         &pentity2_graph,
                                         &pentity2_nuplet);

  PDM_part_comm_graph_t* ptpgc_entity2 = NULL;
  if(ptpgc_entity1->nuplet_size == 0) {
    ptpgc_entity2 = PDM_part_comm_graph_create(ptpgc_entity1->n_part,
                                               pn_entity2_graph,
                                               pentity2_graph,
                                               PDM_OWNERSHIP_KEEP,
                                               ptpgc_entity1->comm);
  } else {
    ptpgc_entity2 = PDM_part_comm_graph_with_nuplet_create(ptpgc_entity1->n_part,
                                                           pn_entity2_graph,
                                                           pentity2_graph,
                                                           PDM_OWNERSHIP_KEEP,
                                                           ptpgc_entity1->nuplet_size,
                                                           pentity2_nuplet,
                                                           PDM_OWNERSHIP_KEEP,
                                                           PDM_TRUE, // is_signed
                                                           ptpgc_entity1->comm);
  }

  PDM_free(pentity2_nuplet);
  PDM_free(pentity2_graph);
  PDM_free(pn_entity2_graph);

  *out_ptpgc_entity2 = ptpgc_entity2;
}


void
PDM_part_comm_graph_entity1_to_entity2
(
  PDM_MPI_Comm             comm,
  int                      n_part,
  int                     *pn_entity1_graph,
  int                    **pentity1_graph,
  int                      nuplet_size,
  int                    **pentity1_nuplet,
  int                     *pn_entity1,
  int                     *pn_entity2,
  int                    **entity2_entity1_idx,
  int                    **entity2_entity1,
  int                    **out_pn_entity2_graph,
  int                   ***out_pentity2_graph,
  int                   ***out_pentity2_nuplet
)
{

  if(0 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_log_trace_array_int(pentity1_graph     [i_part], 4 * pn_entity1_graph[i_part]                   , "pentity1_graph      ::");
      PDM_log_trace_array_int(entity2_entity1_idx[i_part], pn_entity2[i_part]+1                           , "entity2_entity1_idx ::");
      PDM_log_trace_array_int(entity2_entity1    [i_part], entity2_entity1_idx[i_part][pn_entity2[i_part]], "entity2_entity1     ::");
    }
  }

  int n_rank;
  int i_rank;
  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  int i_have_nuplet = (nuplet_size != 0);
  int have_nuplet;
  PDM_MPI_Allreduce(&i_have_nuplet, &have_nuplet, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);

  if (i_have_nuplet != have_nuplet) {
    PDM_error(__FILE__, __LINE__, 0, "Error : inconsistent 'i_have_nuplet'\n");
  }

  /* Create transpose graph information - More pratical */
  int **entity1_to_graph_comm_idx = NULL;
  int **entity1_to_graph_comm     = NULL;

  PDM_malloc(entity1_to_graph_comm_idx, n_part, int *);
  PDM_malloc(entity1_to_graph_comm    , n_part, int *);

  for(int i_part = 0; i_part < n_part; ++i_part) {

    PDM_malloc(entity1_to_graph_comm_idx[i_part], pn_entity1[i_part]+1, int);
    int *_entity1_to_graph_comm_idx = entity1_to_graph_comm_idx[i_part];

    int *entity1_to_graph_comm_n = PDM_array_zeros_int(pn_entity1[i_part]);
    for(int i = 0; i < pn_entity1_graph[i_part]; ++i) {
      int i_entity = pentity1_graph[i_part][4*i]-1;
      entity1_to_graph_comm_n[i_entity]++;
    }

    _entity1_to_graph_comm_idx[0] = 0;
    for(int i = 0; i < pn_entity1[i_part]; ++i) {
      _entity1_to_graph_comm_idx[i+1] = _entity1_to_graph_comm_idx[i] + entity1_to_graph_comm_n[i];
      entity1_to_graph_comm_n[i] = 0;
    }

    PDM_malloc(entity1_to_graph_comm[i_part], _entity1_to_graph_comm_idx[pn_entity1[i_part]], int);
    int *_entity1_to_graph_comm = entity1_to_graph_comm[i_part];

    for(int i = 0; i < pn_entity1_graph[i_part]; ++i) {
      int i_entity = pentity1_graph[i_part][4*i  ]-1;
      int idx_write = _entity1_to_graph_comm_idx[i_entity] + entity1_to_graph_comm_n[i_entity]++;
      _entity1_to_graph_comm[idx_write] = i; // Forcement croissant !!!
    }

    PDM_free(entity1_to_graph_comm_n);

  }

  /* Extract all entity2 in graph_comm */
  int *send_n     = PDM_array_zeros_int(n_rank);
  int *send_key_n = PDM_array_zeros_int(n_rank);

  int  *pn_entity2_graph = NULL;
  int **pentity2_graph   = NULL;
  int **pentity2_nuplet  = NULL;
  PDM_malloc(pn_entity2_graph, n_part, int  );
  PDM_malloc(pentity2_graph  , n_part, int *);
  if (have_nuplet) {
    PDM_malloc(pentity2_nuplet, n_part, int *);
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {

    int *_entity1_to_graph_comm_idx = entity1_to_graph_comm_idx[i_part];
    int *_entity1_to_graph_comm     = entity1_to_graph_comm    [i_part];

    /* Count */
    pn_entity2_graph[i_part] = 0;
    int n_data = 0;
    for(int i_entity2 = 0; i_entity2 < pn_entity2[i_part]; ++i_entity2) {

      int is_on_comm_graph = 1;
      for(int idx_entity2 = entity2_entity1_idx[i_part][i_entity2]; idx_entity2 < entity2_entity1_idx[i_part][i_entity2+1]; ++idx_entity2) {
        int i_entity1 = PDM_ABS(entity2_entity1[i_part][idx_entity2])-1;
        if(_entity1_to_graph_comm_idx[i_entity1+1] - _entity1_to_graph_comm_idx[i_entity1] == 0) {
          is_on_comm_graph = 0;
        }
      }

      if(is_on_comm_graph == 1) {
        for(int idx_entity2 = entity2_entity1_idx[i_part][i_entity2]; idx_entity2 < entity2_entity1_idx[i_part][i_entity2+1]; ++idx_entity2) {
          int i_entity1 = PDM_ABS(entity2_entity1[i_part][idx_entity2])-1;
          n_data += _entity1_to_graph_comm_idx[i_entity1+1] - _entity1_to_graph_comm_idx[i_entity1];
        }
      }
    }

    PDM_malloc(pentity2_graph[i_part], 3 * n_data, int);
    int *_pentity2_graph = pentity2_graph[i_part];

    int *_pentity2_nuplet = NULL;
    if (have_nuplet) {
      PDM_malloc(pentity2_nuplet[i_part], nuplet_size * n_data, int);
      _pentity2_nuplet = pentity2_nuplet[i_part];
    }

    for(int i_entity2 = 0; i_entity2 < pn_entity2[i_part]; ++i_entity2) {

      int is_on_comm_graph = 1;
      for(int idx_entity2 = entity2_entity1_idx[i_part][i_entity2]; idx_entity2 < entity2_entity1_idx[i_part][i_entity2+1]; ++idx_entity2) {
        int i_entity1 = PDM_ABS(entity2_entity1[i_part][idx_entity2])-1;
        if(_entity1_to_graph_comm_idx[i_entity1+1] - _entity1_to_graph_comm_idx[i_entity1] == 0) {
          is_on_comm_graph = 0;
        }
      }

      if(is_on_comm_graph == 1) {

        for(int idx_entity2 = entity2_entity1_idx[i_part][i_entity2]; idx_entity2 < entity2_entity1_idx[i_part][i_entity2+1]; ++idx_entity2) {
          int i_entity1 = PDM_ABS(entity2_entity1[i_part][idx_entity2])-1;

          for(int idx_graph = _entity1_to_graph_comm_idx[i_entity1]; idx_graph < _entity1_to_graph_comm_idx[i_entity1+1]; ++idx_graph) {
            int idx_bound = _entity1_to_graph_comm[idx_graph];
            _pentity2_graph[3*pn_entity2_graph[i_part]  ] = pentity1_graph[i_part][4*idx_bound+1];
            _pentity2_graph[3*pn_entity2_graph[i_part]+1] = pentity1_graph[i_part][4*idx_bound+2]-1;
            _pentity2_graph[3*pn_entity2_graph[i_part]+2] = i_entity2+1;
            if (have_nuplet) {
              memcpy(&_pentity2_nuplet       [nuplet_size*pn_entity2_graph[i_part]],
                     &pentity1_nuplet[i_part][nuplet_size*idx_bound],
                     sizeof(int) * nuplet_size);
            }
            pn_entity2_graph[i_part]++;
          }
        }
      }
    }


    /**
     * At this stage we have the raw graph of entity2 :
     *     - Some info inside can be wrong due to connectivity
     *     - We filter then all
     *
     * On filtre pour ne garder que les entités qui ont toutes leurs entity1 qui existent sur une même partition distante
     */

    // Concatenate graph and nuplet
    int stride = 3;

    int *entity2_graph_nuplet = _pentity2_graph;
    if (have_nuplet) {
      stride += nuplet_size;
      PDM_malloc(entity2_graph_nuplet, pn_entity2_graph[i_part] * stride, int);

      for (int i_entity2 = 0; i_entity2 < pn_entity2_graph[i_part]; i_entity2++) {
        for (int i = 0; i < 3; i++) {
          entity2_graph_nuplet[stride*i_entity2+i] = _pentity2_graph[3*i_entity2+i];
        }
        for (int i = 0; i < nuplet_size; i++) {
          entity2_graph_nuplet[stride*i_entity2+3+i] = _pentity2_nuplet[nuplet_size*i_entity2+i];
        }
      }
    }

    int *tmp_order = NULL;
    PDM_malloc(tmp_order, pn_entity2_graph[i_part], int);
    pn_entity2_graph[i_part] = PDM_order_inplace_unique_int(pn_entity2_graph[i_part], stride, entity2_graph_nuplet, tmp_order);

    PDM_free(tmp_order);

    if (have_nuplet) {
      for (int i_entity2 = 0; i_entity2 < pn_entity2_graph[i_part]; i_entity2++) {
        for (int i = 0; i < 3; i++) {
          _pentity2_graph[3*i_entity2+i] = entity2_graph_nuplet[stride*i_entity2+i];
        }
        for (int i = 0; i < nuplet_size; i++) {
          _pentity2_nuplet[nuplet_size*i_entity2+i] = entity2_graph_nuplet[stride*i_entity2+3+i];
        }
      }

      PDM_free(entity2_graph_nuplet);
      PDM_realloc(_pentity2_nuplet, _pentity2_nuplet, pn_entity2_graph[i_part] * nuplet_size, int);
    }

    PDM_realloc(_pentity2_graph, _pentity2_graph, pn_entity2_graph[i_part] * 3, int);

    // PDM_log_trace_array_int(_pentity2_graph, 3 * pn_entity2_graph[i_part], "_pentity2_graph ::");

    int n_valid = 0;
    for(int idx = 0; idx < pn_entity2_graph[i_part]; ++idx) {

      int i_proc_opp = _pentity2_graph[3*idx  ];
      int i_part_opp = _pentity2_graph[3*idx+1];
      int i_entity2  = _pentity2_graph[3*idx+2]-1;

      /* On check si tous les entités sous jacentes sont valides */
      int is_valid = 1;
      for(int idx_entity2 = entity2_entity1_idx[i_part][i_entity2]; idx_entity2 < entity2_entity1_idx[i_part][i_entity2+1]; ++idx_entity2) {

        int i_entity1 = PDM_ABS(entity2_entity1[i_part][idx_entity2])-1;
        int found = 0;
        for(int idx_graph = _entity1_to_graph_comm_idx[i_entity1]; idx_graph < _entity1_to_graph_comm_idx[i_entity1+1]; ++idx_graph) {
          int idx_bound = _entity1_to_graph_comm[idx_graph];

          int t_proc = pentity1_graph[i_part][4*idx_bound+1];
          int t_part = pentity1_graph[i_part][4*idx_bound+2]-1;
          int same_nuplet = 1;
          if (have_nuplet) {
            same_nuplet = (_compare_nuplets(nuplet_size,
                                            &pentity1_nuplet[i_part][nuplet_size*idx_bound],
                                            &_pentity2_nuplet       [nuplet_size*idx]) == 0);
          }

          if(t_proc == i_proc_opp && t_part == i_part_opp && same_nuplet) {
            found = 1;
          }
        }

        if(found == 0) {
          is_valid = 0;
        }

      }

      if(is_valid == 1) { // Compress inplace
        _pentity2_graph[3*n_valid  ] = i_proc_opp;
        _pentity2_graph[3*n_valid+1] = i_part_opp;
        _pentity2_graph[3*n_valid+2] = i_entity2+1;
        if (have_nuplet) {
          for (int i = 0; i < nuplet_size; i++) {
            _pentity2_nuplet[nuplet_size*n_valid+i] = _pentity2_nuplet[nuplet_size*idx+i];
          }
        }
        assert(n_valid <= idx);
        n_valid++;
      }

      // log_trace("idx = %i - (%i/%i) - i_entity2 = %i / is_valid = %i \n\n", idx, i_proc_opp, i_part_opp, i_entity2, is_valid);
    }

    PDM_realloc(_pentity2_graph, _pentity2_graph, 3 * n_valid, int);
    pentity2_graph  [i_part] = _pentity2_graph;
    pn_entity2_graph[i_part] = n_valid;
    if (have_nuplet) {
      PDM_realloc(_pentity2_nuplet, _pentity2_nuplet, 3 * n_valid, int);
      pentity2_nuplet[i_part] = _pentity2_nuplet;
    }

    // Compute send
    for(int idx = 0; idx < pn_entity2_graph[i_part]; ++idx) {

      int i_proc_opp = _pentity2_graph[3*idx  ];
      int i_part_opp = _pentity2_graph[3*idx+1];
      int i_entity2  = _pentity2_graph[3*idx+2]-1;

      send_n[i_proc_opp] += 6;

      for(int idx_entity2 = entity2_entity1_idx[i_part][i_entity2]; idx_entity2 < entity2_entity1_idx[i_part][i_entity2+1]; ++idx_entity2) {
        int i_entity1 = PDM_ABS(entity2_entity1[i_part][idx_entity2])-1;
        for(int idx_graph = _entity1_to_graph_comm_idx[i_entity1]; idx_graph < _entity1_to_graph_comm_idx[i_entity1+1]; ++idx_graph) {
          int idx_bound = _entity1_to_graph_comm[idx_graph];

          int t_proc = pentity1_graph[i_part][4*idx_bound+1];
          int t_part = pentity1_graph[i_part][4*idx_bound+2]-1;

          int same_nuplet = 1;
          if (have_nuplet) {
            same_nuplet = (_compare_nuplets(nuplet_size,
                                            &pentity1_nuplet[i_part][nuplet_size*idx_bound],
                                            &pentity2_nuplet[i_part][nuplet_size*idx]) == 0);
          }

          if(t_proc == i_proc_opp && t_part == i_part_opp && same_nuplet) {
            send_n[i_proc_opp] += 2;
          }
        }
      }

      send_key_n[i_proc_opp] += 1;
    }
  }

  int *send_idx = NULL;
  PDM_malloc(send_idx, n_rank+1, int);

  send_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    send_idx[i+1] = send_idx[i] + send_n[i];
    send_n[i] = 0;
  }

  int *send_data = NULL;
  PDM_malloc(send_data, send_idx[n_rank], int);

  /* Prepare send */
  for(int i_part = 0; i_part < n_part; ++i_part) {
    int *_entity1_to_graph_comm_idx = entity1_to_graph_comm_idx[i_part];
    int *_entity1_to_graph_comm     = entity1_to_graph_comm    [i_part];
    int *_pentity2_graph            = pentity2_graph           [i_part];

    for(int idx = 0; idx < pn_entity2_graph[i_part]; ++idx) {

      int i_proc_opp = _pentity2_graph[3*idx  ];
      int i_part_opp = _pentity2_graph[3*idx+1];
      int i_entity2  = _pentity2_graph[3*idx+2]-1;

      int idx_write = send_idx[i_proc_opp] + send_n[i_proc_opp];
      send_data[idx_write++] = entity2_entity1_idx[i_part][i_entity2+1] - entity2_entity1_idx[i_part][i_entity2];
      send_data[idx_write++] = i_entity2+1;
      send_data[idx_write++] = i_rank;
      send_data[idx_write++] = i_proc_opp;
      send_data[idx_write++] = i_part;
      send_data[idx_write++] = i_part_opp;
      send_n[i_proc_opp] += 6;

      for(int idx_entity2 = entity2_entity1_idx[i_part][i_entity2]; idx_entity2 < entity2_entity1_idx[i_part][i_entity2+1]; ++idx_entity2) {

        int i_entity1 = PDM_ABS(entity2_entity1[i_part][idx_entity2])-1;
        int found = 0;
        for(int idx_graph = _entity1_to_graph_comm_idx[i_entity1]; idx_graph < _entity1_to_graph_comm_idx[i_entity1+1]; ++idx_graph) {
          int idx_bound = _entity1_to_graph_comm[idx_graph];

          int t_proc = pentity1_graph[i_part][4*idx_bound+1];
          int t_part = pentity1_graph[i_part][4*idx_bound+2]-1;

          int same_nuplet = 1;
          if (have_nuplet) {
            same_nuplet = (_compare_nuplets(nuplet_size,
                                            &pentity1_nuplet[i_part][nuplet_size*idx_bound],
                                            &pentity2_nuplet[i_part][nuplet_size*idx]) == 0);
          }

          if(t_proc == i_proc_opp && t_part == i_part_opp && found == 0 && same_nuplet) {

            send_data[idx_write++] = pentity1_graph[i_part][4*idx_bound  ];
            send_data[idx_write++] = pentity1_graph[i_part][4*idx_bound+3];
            send_n[i_proc_opp] += 2;
            found = 1;

          }
        }
      }

    }
  }

  int *recv_n     = NULL;
  int *recv_key_n = NULL;
  int *recv_idx   = NULL;
  int *recv_data  = NULL;
  PDM_malloc(recv_n    , n_rank  , int);
  PDM_malloc(recv_key_n, n_rank  , int);
  PDM_malloc(recv_idx  , n_rank+1, int);

  PDM_MPI_Alltoall(send_n, 1, PDM_MPI_INT,
                   recv_n, 1, PDM_MPI_INT, comm);

  PDM_MPI_Alltoall(send_key_n, 1, PDM_MPI_INT,
                   recv_key_n, 1, PDM_MPI_INT, comm);

  recv_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    recv_idx[i+1] = recv_idx[i] + recv_n[i];
  }

  PDM_malloc(recv_data, recv_idx[n_rank], int);

  PDM_MPI_Alltoallv(send_data,
                    send_n,
                    send_idx,
                    PDM_MPI_INT,
                    recv_data,
                    recv_n,
                    recv_idx,
                    PDM_MPI_INT,
                    comm);

  if(0 == 1) {
    PDM_log_trace_array_int(send_data, send_idx[n_rank], "send_data ::");
    PDM_log_trace_array_int(recv_data, recv_idx[n_rank], "recv_data ::");
  }

  int n_key_send_tot = 0;
  int n_key_recv_tot = 0;
  for(int i = 0; i < n_rank; ++i) {
    n_key_send_tot += send_key_n[i];
    n_key_recv_tot += recv_key_n[i];
  }


  int ln_entity2_max = 0;
  for(int i_part = 0; i_part < n_part; ++i_part) {
    ln_entity2_max = PDM_MAX(ln_entity2_max, pn_entity2[i_part]);
  }

  int pn_entity2_max = 0;
  PDM_MPI_Allreduce(&ln_entity2_max, &pn_entity2_max, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);

  /**
   * Post-treatment :
   *   - Hash table
   *   - Compute keys and manage conflict directly on buffer
   */
  int idx_read = 0;
  int *send_data_idx = NULL;
  int *send_data_key = NULL;
  PDM_malloc(send_data_idx, n_key_send_tot+1, int);
  PDM_malloc(send_data_key, n_key_send_tot  , int);
  send_data_idx[0] = 0;
  int n_max_connect = 0;
  for(int i = 0; i < n_key_send_tot; ++i) {
    int n_connec = send_data[idx_read++];
    int n_data   = 6+2*n_connec;
    send_data_idx[i+1] = send_data_idx[i] + n_data;

    n_max_connect = PDM_MAX(n_max_connect, n_data);

    idx_read++; // Ignore entity2
    send_data_key[i] = 0;
    for(int i_data = 0; i_data < 4 + 2 * n_connec; ++i_data) {
      send_data_key[i] += send_data[idx_read++];
    }
    send_data_key[i] = send_data_key[i] % pn_entity2_max;
  }

  if(0 == 1) {
    PDM_log_trace_array_int(send_data_key, n_key_send_tot, "send_data_key = ");
  }

  int n_g_max_connect = 0;
  PDM_MPI_Allreduce(&n_max_connect, &n_g_max_connect, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);
  n_max_connect = n_g_max_connect;

  int *order     = NULL;
  int *is_solved = NULL;
  PDM_malloc(order    , n_key_send_tot, int);
  PDM_malloc(is_solved, n_key_send_tot, int);

  for (int i = 0; i < n_key_send_tot; ++i) {
    order    [i] = i;
    is_solved[i] = 0;
  }
  PDM_sort_int(send_data_key, order, n_key_send_tot);

  /* Identify keys in conflict */
  int n_conflit_to_solve = 0;
  int last_key = -1;

  int *key_conflict_idx = NULL;
  int *key_to_conflict  = NULL;
  PDM_malloc(key_conflict_idx, n_key_send_tot+1, int);
  PDM_malloc(key_to_conflict , n_key_send_tot  , int);

  key_conflict_idx[0] = 0;
  for (int i = 0; i < n_key_send_tot; ++i) {
    if (send_data_key[i] != last_key){
      key_conflict_idx[n_conflit_to_solve+1] = key_conflict_idx[n_conflit_to_solve]+1;
      n_conflit_to_solve++;
      last_key = send_data_key[i];
    } else {
      key_conflict_idx[n_conflit_to_solve]++;
    }
    key_to_conflict[i] = n_conflit_to_solve-1;
  }

  int *tmp_connec = NULL;
  PDM_malloc(tmp_connec, n_max_connect, int);

  int **tmp_pentity2_nuplet = NULL;
  if (have_nuplet) {
    tmp_pentity2_nuplet = pentity2_nuplet;
    PDM_malloc(pentity2_nuplet, n_part, int *);
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_realloc(pentity2_graph[i_part], pentity2_graph[i_part], 4 * pn_entity2_graph[i_part], int);
    if (have_nuplet) {
      PDM_malloc(pentity2_nuplet[i_part], pn_entity2_graph[i_part] * nuplet_size, int);
    }
    pn_entity2_graph[i_part] = 0;
  }

  /* Link with receive */
  idx_read = 0;
  for(int i_key = 0; i_key < n_key_recv_tot; ++i_key) {

    int n_connec_opp  = recv_data[idx_read++];
    int i_entity2_opp = recv_data[idx_read++];

    int key_opp = 0;
    for(int i_data = 0; i_data < 4 + 2 * n_connec_opp; ++i_data) {
      key_opp += recv_data[idx_read+i_data];
    }
    key_opp = key_opp % pn_entity2_max;

    int i_proc_cur = recv_data[idx_read++];
    int i_proc_opp = recv_data[idx_read++];
    int i_part_cur = recv_data[idx_read++];
    int i_part_opp = recv_data[idx_read++];

    for(int i_data = 0; i_data < 2 * n_connec_opp; ++i_data) {
      tmp_connec[i_data] = recv_data[idx_read++];
    }

    int pos_key = -1;
    if(n_key_send_tot > 0) {
      pos_key = PDM_binary_search_int(key_opp, send_data_key, n_key_send_tot);
    }
    if(pos_key == -1) { // Completly wrong opposite
      continue;
    }
    int i_conflict = key_to_conflict[pos_key];

    // log_trace("key_opp = %i / i_conflict = %i \n", key_opp, i_conflict);
    // PDM_log_trace_array_int(tmp_connec, 2 * n_connec_opp, "tmp_connec ::");

    /* Brute force */
    int n_conflict_entitys = key_conflict_idx[i_conflict+1] - key_conflict_idx[i_conflict];

    for (int idx_entity = 0; idx_entity < n_conflict_entitys; ++idx_entity) {

      int idx_decompose_entity1 = order[key_conflict_idx[i_conflict]+idx_entity];
      // log_trace("\t idx_decompose_entity1 = %i \n", idx_decompose_entity1);

      int beg_data = send_data_idx[idx_decompose_entity1];

      int n_connec_cur  = send_data[beg_data  ];
      int i_entity2_cur = send_data[beg_data+1];

      int i_proc_cur2 = send_data[beg_data+2];
      int i_proc_opp2 = send_data[beg_data+3];
      int i_part_cur2 = send_data[beg_data+4];
      int i_part_opp2 = send_data[beg_data+5];

      if(i_proc_opp == i_proc_cur2 && i_part_opp == i_part_cur2 &&
         i_proc_cur == i_proc_opp2 && i_part_cur == i_part_opp2 && n_connec_cur == n_connec_opp && is_solved[idx_decompose_entity1] == 0) {

        int* connect_cur = &send_data[beg_data+6];

        // PDM_log_trace_array_int(connect_cur, 2 * n_connec_cur, "connect_cur ::");

        // On doit chercher par couple
        int i_entity1_cur1 = tmp_connec[0];
        int i_entity1_cur2 = tmp_connec[1];

        int first_idx = 0;
        for(int idx = 0; idx < n_connec_cur; ++idx) {
          int t_entity1_cur1 = connect_cur[2*idx  ];
          int t_entity1_cur2 = connect_cur[2*idx+1];

          if(t_entity1_cur1 == i_entity1_cur2 && t_entity1_cur2 == i_entity1_cur1 ) {
            first_idx = idx;
            break;
          }
        }
        int first_idx_save = first_idx;

        // log_trace("Match proc/part -> first_idx = %i \n", first_idx);

        // Cas particulier pour les edges
        int is_same = 1;
        int sens    = 1;
        if(n_connec_cur == 2) {

          int c1_entity1_cur1 = tmp_connec[0];
          int c1_entity1_cur2 = tmp_connec[1];
          int t1_entity1_cur1 = connect_cur[0];
          int t1_entity1_cur2 = connect_cur[1];

          int c2_entity1_cur1 = tmp_connec[2];
          int c2_entity1_cur2 = tmp_connec[3];
          int t2_entity1_cur1 = connect_cur[2];
          int t2_entity1_cur2 = connect_cur[3];

          if(t1_entity1_cur1 == c1_entity1_cur2 &&
             t1_entity1_cur2 == c1_entity1_cur1 &&
             t2_entity1_cur1 == c2_entity1_cur2 &&
             t2_entity1_cur2 == c2_entity1_cur1) {
            is_same = 1;
            sens = 1;
          } else if (t1_entity1_cur1 == c2_entity1_cur2 &&
                     t1_entity1_cur2 == c2_entity1_cur1 &&
                     t2_entity1_cur1 == c1_entity1_cur2 &&
                     t2_entity1_cur2 == c1_entity1_cur1) {
            is_same = 1;
            sens    = -1;
          } else {
            is_same = 0;
          }

        } else {

          for(int idx = 0; idx < n_connec_cur; ++idx) {

            i_entity1_cur1 = tmp_connec[2*idx  ];
            i_entity1_cur2 = tmp_connec[2*idx+1];
            int t_entity1_cur1 = connect_cur[2*first_idx  ];
            int t_entity1_cur2 = connect_cur[2*first_idx+1];

            if(t_entity1_cur1 != i_entity1_cur2 && t_entity1_cur2 != i_entity1_cur1 ) {
              is_same = 0;
            }

            first_idx++;
            if(first_idx == n_connec_cur) {
              first_idx = first_idx % n_connec_cur;
            }
          }

          // Si ca echoue on tente le revert
          if(is_same == 0) {
            is_same = 1;
            first_idx = first_idx_save;
            for(int idx = 0; idx < n_connec_cur; ++idx) {

              i_entity1_cur1 = tmp_connec[2*idx  ];
              i_entity1_cur2 = tmp_connec[2*idx+1];
              int t_entity1_cur1 = connect_cur[2*first_idx  ];
              int t_entity1_cur2 = connect_cur[2*first_idx+1];

              if(t_entity1_cur1 != i_entity1_cur2 && t_entity1_cur2 != i_entity1_cur1 ) {
                is_same = 0;
              }

              if(first_idx == 0) {
                first_idx = n_connec_cur;
              }
              first_idx--;
            }
            if(is_same == 1) {
              sens = -1;
            }
          }
        }

        // log_trace("Match proc/part -> is_same = %i \n", is_same);

        if(is_same == 1) {

          // Rebuild graph
          pentity2_graph[i_part_cur2][4*pn_entity2_graph[i_part_cur2]  ] = i_entity2_cur;
          pentity2_graph[i_part_cur2][4*pn_entity2_graph[i_part_cur2]+1] = i_proc_cur;
          pentity2_graph[i_part_cur2][4*pn_entity2_graph[i_part_cur2]+2] = i_part_cur+1; // Because i_part start at 1 in paradigm
          pentity2_graph[i_part_cur2][4*pn_entity2_graph[i_part_cur2]+3] = i_entity2_opp * sens; // Pour l'instant on le met la le sens

          // Rebuild nuplet
          if (have_nuplet) {
            for (int i = 0; i < nuplet_size; i++) {
              pentity2_nuplet[i_part_cur2][nuplet_size*pn_entity2_graph[i_part_cur2]+i] = tmp_pentity2_nuplet[i_part_cur2][nuplet_size*idx_decompose_entity1+i];
            }
          }

          is_solved[idx_decompose_entity1] = 1;
          pn_entity2_graph[i_part_cur2]++;
        }
      }
    }
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_realloc(pentity2_graph[i_part], pentity2_graph[i_part], 4 * pn_entity2_graph[i_part], int);
    if (have_nuplet) {
      PDM_realloc(pentity2_nuplet[i_part], pentity2_nuplet[i_part], nuplet_size * pn_entity2_graph[i_part], int);
      PDM_free(tmp_pentity2_nuplet[i_part]);
    }
  }
  PDM_free(tmp_pentity2_nuplet);


  if(0 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_log_trace_array_int(pentity2_graph[i_part], 4 * pn_entity2_graph[i_part], "pentity2_graph ::");
    }
  }

  PDM_free(tmp_connec);
  PDM_free(key_conflict_idx);
  PDM_free(key_to_conflict);
  PDM_free(is_solved);
  PDM_free(order);
  PDM_free(send_data_key);
  PDM_free(send_data_idx);

  PDM_free(send_n);
  PDM_free(send_idx);
  PDM_free(send_data);
  PDM_free(send_key_n);

  PDM_free(recv_n);
  PDM_free(recv_idx);
  PDM_free(recv_data);
  PDM_free(recv_key_n);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(entity1_to_graph_comm_idx[i_part]);
    PDM_free(entity1_to_graph_comm    [i_part]);
  }
  PDM_free(entity1_to_graph_comm_idx);
  PDM_free(entity1_to_graph_comm    );
  PDM_free(send_n                   );

  *out_pn_entity2_graph = pn_entity2_graph;
  *out_pentity2_graph   = pentity2_graph;
  if (have_nuplet) {
    *out_pentity2_nuplet = pentity2_nuplet;
  }

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

  PDM_free(pcg);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
