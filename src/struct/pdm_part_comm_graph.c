/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <unistd.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_part_comm_graph.h"
#include "pdm_part_comm_graph_priv.h"
#include "pdm.h"
#include "pdm_timer.h"
#include "pdm_priv.h"
#include "pdm_binary_search.h"
#include "pdm_sort.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_array.h"
#include "pdm_order.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
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
 PDM_part_comm_graph_t   *ptpgc,
 size_t                   s_data,
 int                      cst_stride,
 void                   **send_entity_data,
 unsigned char          **send_buffer
)
{

  int n_rank;
  PDM_MPI_Comm_size(ptpgc->comm, &n_rank);

  int s_data_tot = s_data * cst_stride;

  int send_buff_size = s_data * cst_stride * ptpgc->send_idx[n_rank];
  unsigned char *_send_buffer = NULL;
  PDM_malloc(_send_buffer, send_buff_size, unsigned char);

  unsigned char **_send_entity_data = (unsigned char **) send_entity_data;

  for(int i_part = 0; i_part < ptpgc->n_part; ++i_part) {
    for(int i = 0; i < ptpgc->n_entity_graph[i_part]; ++i) {
      for(int k = 0; k < s_data_tot; ++k) {
        int idx_read  = i * s_data_tot + k;
        int idx_write = ptpgc->part_to_send_buffer[i_part][i] * s_data * cst_stride + k;
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
 PDM_part_comm_graph_t   *ptpgc,
 size_t                   s_data,
 int                      cst_stride,
 unsigned char           *recv_buffer,
 void                  ***recv_entity_data
)
{
  int s_data_tot = s_data * cst_stride;

  unsigned char **_recv_entity_data = NULL;
  PDM_malloc(_recv_entity_data, ptpgc->n_part, unsigned char *);
  *recv_entity_data = (void **) _recv_entity_data;
  for(int i_part = 0; i_part < ptpgc->n_part; ++i_part) {
    _recv_entity_data[i_part] = malloc(ptpgc->n_entity_graph[i_part] * s_data_tot * sizeof(unsigned char));
    for(int i = 0; i < ptpgc->n_entity_graph[i_part]; ++i) {
      for(int k = 0; k < s_data_tot; ++k) {
        int idx_write = i * s_data_tot + k;
        int idx_read  = ptpgc->part_to_recv_buffer[i_part][i] * s_data * cst_stride + k;
        _recv_entity_data[i_part][idx_write] = recv_buffer[idx_read];
      }
    }
  }

}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_part_comm_graph_t*
PDM_part_comm_graph_create
(
  int            n_part,
  int           *pn_entity_graph,
  int          **pentity_graph,
  PDM_MPI_Comm   comm
)
{
  PDM_part_comm_graph_t *ptpgc = NULL;
  PDM_malloc(ptpgc, 1 ,PDM_part_comm_graph_t);

  ptpgc->comm   = comm;
  ptpgc->n_part = n_part;

  ptpgc->n_active_rank_send = 0;
  ptpgc->n_active_rank_recv = 0;

  int n_rank;
  int i_rank;
  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  ptpgc->pentity_graph = pentity_graph;

  int n_g_part = 0;
  PDM_MPI_Allreduce(&n_part, &n_g_part, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
  ptpgc->n_g_part = n_g_part;

  PDM_malloc(ptpgc->part_to_send_buffer, n_part, int *);
  PDM_malloc(ptpgc->part_to_recv_buffer, n_part, int *);
  PDM_malloc(ptpgc->n_entity_graph     , n_part, int  );

  PDM_MPI_Datatype mpi_doublet_type;
  PDM_MPI_Type_create_contiguous(2, PDM_MPI_INT, &mpi_doublet_type);
  PDM_MPI_Type_commit(&mpi_doublet_type);

  int *send_n = PDM_array_zeros_int(n_rank);
  int *recv_n = NULL;
  PDM_malloc(recv_n, n_rank, int);
  for (int i_part = 0; i_part < n_part; i_part++) {
    int n_entity_graph = pn_entity_graph[i_part];
    ptpgc->n_entity_graph[i_part] = n_entity_graph;

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
      ptpgc->n_active_rank_send++;
    }
    send_n[i] = 0;
  }
  int *send_doublet = NULL;
  PDM_malloc(send_doublet, 2 * send_idx[n_rank], int);

  for (int i_part = 0; i_part < n_part; i_part++) {
    int n_entity_graph = pn_entity_graph[i_part];

    PDM_malloc(ptpgc->part_to_send_buffer[i_part], n_entity_graph, int);

    for(int idx_entity = 0; idx_entity < n_entity_graph; ++idx_entity) {
      int t_rank = pentity_graph[i_part][4*idx_entity+1];
      int idx_write = send_idx[t_rank] + send_n[t_rank]++;
      send_doublet[2*idx_write  ] = pentity_graph[i_part][4*idx_entity+2]-1;
      send_doublet[2*idx_write+1] = pentity_graph[i_part][4*idx_entity+3]-1;

      ptpgc->part_to_send_buffer[i_part][idx_entity] = idx_write;
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
      ptpgc->n_active_rank_recv++;
    }
  }

  int *recv_doublet = NULL;
  PDM_malloc(recv_doublet, 2 * recv_idx[n_rank], int);

  if(0 == 1) {
    PDM_log_trace_array_int(send_idx, n_rank+1, "send_idx ::");
    PDM_log_trace_array_int(recv_idx, n_rank+1, "recv_idx ::");
  }

  PDM_MPI_Alltoallv(send_doublet,
                    send_n,
                    send_idx,
                    mpi_doublet_type,
                    recv_doublet,
                    recv_n,
                    recv_idx,
                    mpi_doublet_type,
                    comm);

  /* Maintenant on cherche à retrouver la correspondance */
  int **pentity_indices       = NULL;
  int **pentity_indices_order = NULL;
  PDM_malloc(pentity_indices      , n_part, int *);
  PDM_malloc(pentity_indices_order, n_part, int *);

  for (int i_part = 0; i_part < n_part; i_part++) {
    int n_entity_graph = pn_entity_graph[i_part];
    PDM_malloc(ptpgc->part_to_recv_buffer[i_part], n_entity_graph, int);

    PDM_malloc(pentity_indices      [i_part], 2 * n_entity_graph, int);
    PDM_malloc(pentity_indices_order[i_part],     n_entity_graph, int);

    for(int idx_entity = 0; idx_entity < n_entity_graph; ++idx_entity) {
      pentity_indices[i_part][2*idx_entity  ] = pentity_graph[i_part][4*idx_entity+1];
      pentity_indices[i_part][2*idx_entity+1] = pentity_graph[i_part][4*idx_entity]-1; // Indices locaux
    }

    PDM_order_lnum_s(pentity_indices[i_part],
                     2,
                     pentity_indices_order[i_part],
                     n_entity_graph);

    PDM_order_array(n_entity_graph,
                    2 * sizeof(int),
                    pentity_indices_order[i_part],
                    pentity_indices      [i_part]);

  }

  /* Post buffer */
  for(int t_rank = 0; t_rank < n_rank; ++t_rank) {
    for(int j = recv_idx[t_rank]; j < recv_idx[t_rank+1]; ++j) {

      int lpart   = recv_doublet[2*j  ];
      int lentity = recv_doublet[2*j+1];
      int n_entity_graph = pn_entity_graph[lpart];

      int to_find[2] = {t_rank, lentity};

      int pos = PDM_order_binary_search_int(to_find, pentity_indices[lpart], 2, n_entity_graph);
      // log_trace("Try to find = (%i/%i) --> %i \n", t_rank, lentity, pos);

      ptpgc->part_to_recv_buffer[lpart][pentity_indices_order[lpart][pos]] = j;

    }
  }

  if(0 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int n_entity_graph = pn_entity_graph[i_part];
      PDM_log_trace_array_int(ptpgc->part_to_send_buffer[i_part], n_entity_graph, "ptpgc->part_to_send_buffer ::");
      PDM_log_trace_array_int(ptpgc->part_to_recv_buffer[i_part], n_entity_graph, "ptpgc->part_to_recv_buffer ::");
    }
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pentity_indices      [i_part]);
    PDM_free(pentity_indices_order[i_part]);
  }
  PDM_free(pentity_indices       );
  PDM_free(pentity_indices_order );


  PDM_free(send_doublet);
  PDM_free(recv_doublet);
  PDM_MPI_Type_free(&mpi_doublet_type);

  ptpgc->send_idx = send_idx;
  ptpgc->recv_idx = recv_idx;
  ptpgc->send_n   = send_n;
  ptpgc->recv_n   = recv_n;

  /* Compute p2p array */
  PDM_malloc(ptpgc->active_rank_send, ptpgc->n_active_rank_send, int);
  PDM_malloc(ptpgc->active_rank_recv, ptpgc->n_active_rank_recv, int);

  ptpgc->n_active_rank_send = 0;
  ptpgc->n_active_rank_recv = 0;

  for(int i = 0; i < n_rank; ++i) {
    if(send_n[i] > 0) {
      ptpgc->active_rank_send[ptpgc->n_active_rank_send++] = i;
    }
    if(recv_n[i] > 0) {
      ptpgc->active_rank_recv[ptpgc->n_active_rank_recv++] = i;
    }
  }

  /*
   * Compute owner
   */
  ptpgc->bound_owner = NULL;
  PDM_malloc(ptpgc->bound_owner, n_part, int *);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    int n_entity_graph = pn_entity_graph[i_part];

    PDM_malloc(ptpgc->bound_owner[i_part], n_entity_graph, int);
    int *lbound_entity       = NULL;
    PDM_malloc(lbound_entity, n_entity_graph, int);

    for(int idx_entity = 0; idx_entity < n_entity_graph; ++idx_entity) {
      lbound_entity      [idx_entity] = pentity_graph[i_part][4*idx_entity];
      ptpgc->bound_owner[i_part][idx_entity] = -1;
    }

    int n_unique = PDM_inplace_unique(lbound_entity, 0, n_entity_graph-1);

    int *lowner = PDM_array_const_int(n_unique, -1);
    for(int idx_entity = 0; idx_entity < n_entity_graph; ++idx_entity) {
      int l_entity = pentity_graph[i_part][4*idx_entity];
      int t_rank   = pentity_graph[i_part][4*idx_entity+1];
      int t_part   = pentity_graph[i_part][4*idx_entity+2]-1;
      int pos = PDM_binary_search_int(l_entity, lbound_entity, n_unique);

      ptpgc->bound_owner[i_part][idx_entity] = pos; // Stockage temporaire

      if(lowner[pos] != 0) {
        if(i_rank < t_rank || (i_rank == t_rank && i_part < t_part)) {
          lowner[pos] = 1;
        } else {
          lowner[pos] = 0;
        }
      }
    }

    /* Last loop to fill */
    for(int idx_entity = 0; idx_entity < n_entity_graph; ++idx_entity) {
      int pos = ptpgc->bound_owner[i_part][idx_entity];
      ptpgc->bound_owner[i_part][idx_entity] = lowner[pos];
    }

    PDM_free(lowner);
    PDM_free(lbound_entity);
  }

  return ptpgc;
}

void
PDM_part_comm_graph_exch
(
 PDM_part_comm_graph_t   *ptpgc,
 size_t                   s_data,
 PDM_stride_t             t_stride,
 int                      cst_stride,
 int                    **send_entity_stride,
 void                   **send_entity_data,
 int                   ***recv_entity_stride,
 void                  ***recv_entity_data
)
{
  int s_data_tot = s_data;
  if(t_stride == PDM_STRIDE_CST_INTERLACED) {
    s_data_tot = s_data * cst_stride;
  }

  PDM_MPI_Datatype mpi_type;
  PDM_MPI_Type_create_contiguous(s_data_tot, PDM_MPI_BYTE, &mpi_type);
  PDM_MPI_Type_commit(&mpi_type);

  int n_rank;
  PDM_MPI_Comm_size(ptpgc->comm, &n_rank);

  if(t_stride == PDM_STRIDE_CST_INTERLACED) {

    unsigned char *send_buffer = NULL;
    _prepare_send_strid_cst(ptpgc,
                            s_data,
                            cst_stride,
                            send_entity_data,
                            &send_buffer);

    int recv_buff_size = s_data * cst_stride * ptpgc->recv_idx[n_rank];
    unsigned char *recv_buffer = NULL;
    PDM_malloc(recv_buffer, recv_buff_size, unsigned char);
    PDM_MPI_Alltoallv(send_buffer,
                      ptpgc->send_n,
                      ptpgc->send_idx,
                      mpi_type,
                      recv_buffer,
                      ptpgc->recv_n,
                      ptpgc->recv_idx,
                      mpi_type,
                      ptpgc->comm);
    PDM_free(send_buffer);

    /* Post-traitement */
    _post_recv_strid_cst(ptpgc,
                         s_data,
                         cst_stride,
                         recv_buffer,
                         recv_entity_data);

    PDM_free(recv_buffer);

  } else if (t_stride == PDM_STRIDE_VAR_INTERLACED) {

    /* Exchange stride */
    int  *send_stride = NULL;
    PDM_malloc(send_stride, ptpgc->send_idx[n_rank], int);
    for(int i_part = 0; i_part < ptpgc->n_part; ++i_part) {
      for(int i = 0; i < ptpgc->n_entity_graph[i_part]; ++i) {
        int idx_write = ptpgc->part_to_send_buffer[i_part][i];
        send_stride[idx_write] = send_entity_stride[i_part][i];
      }
    }

    int *recv_stride = NULL;
    PDM_malloc(recv_stride, ptpgc->recv_idx[n_rank], int);
    PDM_MPI_Alltoallv(send_stride,
                      ptpgc->send_n,
                      ptpgc->send_idx,
                      PDM_MPI_INT,
                      recv_stride,
                      ptpgc->recv_n,
                      ptpgc->recv_idx,
                      PDM_MPI_INT,
                      ptpgc->comm);

    // PDM_log_trace_array_int(send_stride, ptpgc->send_idx[n_rank], "send_stride ::");
    // PDM_log_trace_array_int(recv_stride, ptpgc->recv_idx[n_rank], "recv_stride ::");


    /* Exchange data */
    int *send_stride_idx = NULL;
    PDM_malloc(send_stride_idx, ptpgc->send_idx[n_rank]+1, int);
    send_stride_idx[0] = 0;
    for(int i = 0; i < ptpgc->send_idx[n_rank]; ++i) {
      send_stride_idx[i+1] = send_stride_idx[i] + send_stride[i];
    }

    int *recv_stride_idx = NULL;
    PDM_malloc(recv_stride_idx, ptpgc->recv_idx[n_rank]+1, int);
    recv_stride_idx[0] = 0;
    for(int i = 0; i < ptpgc->recv_idx[n_rank]; ++i) {
      recv_stride_idx[i+1] = recv_stride_idx[i] + recv_stride[i];
    }

    int *send_data_idx = NULL;
    PDM_malloc(send_data_idx, n_rank+1, int);
    int *send_data_n   = PDM_array_zeros_int(n_rank);
    send_data_idx[0] = 0;
    for(int i = 0; i < n_rank; ++i) {
      send_data_idx[i+1] = send_data_idx[i];
      for(int j = ptpgc->send_idx[i]; j < ptpgc->send_idx[i+1]; ++j) {
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
      for(int j = ptpgc->recv_idx[i]; j < ptpgc->recv_idx[i+1]; ++j) {
        recv_data_idx[i+1] += recv_stride[j];
        recv_data_n  [i  ] += recv_stride[j];
      }
    }

    // PDM_log_trace_array_int(send_data_idx, n_rank+1, "send_data_idx ::");
    // PDM_log_trace_array_int(recv_data_idx, n_rank+1, "recv_data_idx ::");
    // PDM_log_trace_array_int(send_stride_idx, ptpgc->send_idx[n_rank]+1, "send_stride_idx ::");
    // PDM_log_trace_array_int(recv_stride_idx, ptpgc->recv_idx[n_rank]+1, "recv_stride_idx ::");

    int send_buff_size = send_data_idx[n_rank] * s_data_tot;
    unsigned char  *send_buffer       = malloc(send_buff_size * sizeof(unsigned char));
    unsigned char **_send_entity_data = (unsigned char **) send_entity_data;

    for(int i_part = 0; i_part < ptpgc->n_part; ++i_part) {
      int idx_read = 0;
      for(int i = 0; i < ptpgc->n_entity_graph[i_part]; ++i) {
        for(int j = 0; j < send_entity_stride[i_part][i]; ++j) {
          int idx_buffer = ptpgc->part_to_send_buffer[i_part][i];
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

    // log_trace("send_buff_size : %i \n", send_buff_size);
    // log_trace("recv_buff_size : %i \n", recv_buff_size);
    // PDM_log_trace_array_int(recv_data_idx, ptpgc->recv_idx[n_rank], "recv_data_idx ::");
    PDM_MPI_Alltoallv(send_buffer,
                      send_data_n,
                      send_data_idx,
                      mpi_type,
                      recv_buffer,
                      recv_data_n,
                      recv_data_idx,
                      mpi_type,
                      ptpgc->comm);
    PDM_free(send_buffer);

    /* Panic verbose */
    // PDM_g_num_t* recv_buffer_dbg = (PDM_g_num_t *) recv_buffer;
    // PDM_log_trace_array_long(recv_buffer_dbg, recv_buff_size/s_data, "recv_buffer_dbg ::");

    /* Post-traitement stride */
    int           **_recv_entity_stride = NULL;
    unsigned char **_recv_entity_data   = NULL;
    PDM_malloc(_recv_entity_stride, ptpgc->n_part, int           *);
    PDM_malloc(_recv_entity_data  , ptpgc->n_part, unsigned char *);
    *recv_entity_stride =           _recv_entity_stride;
    *recv_entity_data   = (void **) _recv_entity_data;

    for(int i_part = 0; i_part < ptpgc->n_part; ++i_part) {
      _recv_entity_stride[i_part] = malloc(ptpgc->n_entity_graph[i_part] * sizeof(int));
      int recv_buff_size_part = 0;
      for(int i = 0; i < ptpgc->n_entity_graph[i_part]; ++i) {
        int idx_read  = ptpgc->part_to_recv_buffer[i_part][i];
        _recv_entity_stride[i_part][i] = recv_stride[idx_read];
        recv_buff_size_part += recv_stride[idx_read];
      }

      PDM_malloc(_recv_entity_data[i_part], recv_buff_size_part * s_data_tot, unsigned char);

      // PDM_log_trace_array_int(_recv_entity_stride[i_part], ptpgc->n_entity_graph[i_part], "_recv_entity_stride :");
    }

    /*
     * Post-treatment buffer
     */
    for(int i_part = 0; i_part < ptpgc->n_part; ++i_part) {
      int idx_write = 0;
      for(int i = 0; i < ptpgc->n_entity_graph[i_part]; ++i) {
        int idx_buffer = ptpgc->part_to_recv_buffer[i_part][i];
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

  } else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_comm_graph_exch, wrong t_stride \n");
  }

  PDM_MPI_Type_free(&mpi_type);

}

const int*
PDM_part_comm_graph_owner_get
(
 PDM_part_comm_graph_t *ptpgc,
 int                    i_part
)
{
  return ptpgc->bound_owner[i_part];
}

void
PDM_part_comm_graph_reorder
(
  PDM_part_comm_graph_t  *ptpgc,
  int                   **pentity_graph,
  int                   **old_to_new
)
{

  /* Prepare exchange */
  int **send_new_id = NULL;
  PDM_malloc(send_new_id, ptpgc->n_part, int         *);
  for(int i_part = 0; i_part < ptpgc->n_part; ++i_part) {
    PDM_malloc(send_new_id[i_part], ptpgc->n_entity_graph[i_part], int        );
    for(int i = 0; i < ptpgc->n_entity_graph[i_part]; ++i) {
      int i_entity = pentity_graph[i_part][4*i  ]-1;
      send_new_id[i_part][i] = old_to_new[i_part][i_entity];
    }
  }

  int **recv_new_id = NULL;
  PDM_part_comm_graph_exch(ptpgc,
                                   sizeof(int),
                                   PDM_STRIDE_CST_INTERLACED,
                                   1,
                                   NULL,
                      (void  **)   send_new_id,
                                   NULL,
                      (void ***)   &recv_new_id);

  for(int i_part = 0; i_part < ptpgc->n_part; ++i_part) {
    PDM_free(send_new_id[i_part]);
  }
  PDM_free(send_new_id);

  /* Actualisation current and opposite */

  for(int i_part = 0; i_part < ptpgc->n_part; ++i_part) {
    for(int i = 0; i < ptpgc->n_entity_graph[i_part]; ++i) {
      int i_entity = pentity_graph[i_part][4*i  ]-1;
      pentity_graph[i_part][4*i  ] = old_to_new[i_part][i_entity]+1; // On suppose que old_to_new commence a 0
      pentity_graph[i_part][4*i+3] = recv_new_id[i_part][i]+1; // On suppose que old_to_new commence a 0
    }
  }

  for(int i_part = 0; i_part < ptpgc->n_part; ++i_part) {
    PDM_free(recv_new_id[i_part]);
  }
  PDM_free(recv_new_id);
}

void
PDM_part_comm_graph_free
(
 PDM_part_comm_graph_t* ptpgc
)
{

  for(int i_part = 0; i_part < ptpgc->n_part; ++i_part) {
    PDM_free(ptpgc->part_to_send_buffer[i_part]);
    PDM_free(ptpgc->part_to_recv_buffer[i_part]);
    PDM_free(ptpgc->bound_owner        [i_part]);
  }
  PDM_free(ptpgc->part_to_send_buffer);
  PDM_free(ptpgc->part_to_recv_buffer);
  PDM_free(ptpgc->n_entity_graph);
  PDM_free(ptpgc->bound_owner);

  PDM_free(ptpgc->send_idx);
  PDM_free(ptpgc->recv_idx);
  PDM_free(ptpgc->send_n);
  PDM_free(ptpgc->recv_n);
  PDM_free(ptpgc->active_rank_send);
  PDM_free(ptpgc->active_rank_recv);

  PDM_free(ptpgc);
}




#ifdef __cplusplus
}
#endif /* __cplusplus */
