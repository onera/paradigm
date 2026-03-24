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
#include "pdm_part_comm_graph_algorithm.h"
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

/*=============================================================================
 * Public function definitions
 *============================================================================*/

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
            send_n[i_proc_opp] += 2;
            found = 1;
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


void
PDM_part_comm_graph_selected_entity1_to_selected_entity2
(
  int                     *n_selected_entity1,
  int                    **selected_entity1,
  int                    **entity1_entity2_idx,
  int                    **entity1_entity2,
  PDM_part_comm_graph_t   *pcg_entity2,
  int                    **out_n_selected_entity2,
  int                   ***out_selected_entity2
)
{
  if (pcg_entity2 == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Part Comm Graph is NULL\n");
  }

  int n_part = pcg_entity2->n_part;

  PDM_malloc(*out_n_selected_entity2, n_part, int  );
  PDM_malloc(*out_selected_entity2,   n_part, int *);

  int  *n_selected_entity2 = *out_n_selected_entity2;
  int **selected_entity2   = *out_selected_entity2;


  int **entity2_flag = NULL;
  PDM_malloc(entity2_flag, n_part, int *);
  for (int i_part = 0; i_part < n_part; i_part++) {

    /* Compute upper bound for n_entity2 */
    int n_entity2_ub = 0;
    if (selected_entity1 == NULL) {
      // Account for all entities1
      for (int idx_entity2 = 0; idx_entity2 < entity1_entity2_idx[i_part][n_selected_entity1[i_part]]; idx_entity2++) {
        n_entity2_ub = PDM_MAX(n_entity2_ub, PDM_ABS(entity1_entity2[i_part][idx_entity2]));
      }
    }
    else {
      // Only account for subset of entities1
      for (int idx_entity1 = 0; idx_entity1 < n_selected_entity1[i_part]; idx_entity1++) {
        int i_entity1 = selected_entity1[i_part][idx_entity1] - 1;
        for (int idx_entity2 = entity1_entity2_idx[i_part][i_entity1]; idx_entity2 < entity1_entity2_idx[i_part][i_entity1+1]; idx_entity2++) {
          n_entity2_ub = PDM_MAX(n_entity2_ub, PDM_ABS(entity1_entity2[i_part][idx_entity2]));
        }
      }
    }

    int *graph_entity2 = NULL;
    int graph_entity2_n = PDM_part_comm_graph_entity_graph_get(pcg_entity2,
                                                               i_part,
                                                              &graph_entity2,
                                                               PDM_OWNERSHIP_BAD_VALUE);
    for (int idx_entity2 = 0; idx_entity2 < graph_entity2_n; idx_entity2++) {
      n_entity2_ub = PDM_MAX(n_entity2_ub, graph_entity2[4*idx_entity2]);
    }



    /* Flag entities2 incident to *local* entities1 */
    entity2_flag[i_part] = PDM_array_zeros_int(n_entity2_ub);

    n_selected_entity2[i_part] = 0;
    PDM_malloc(selected_entity2[i_part], n_entity2_ub, int);

    for (int idx_entity1 = 0; idx_entity1 < n_selected_entity1[i_part]; idx_entity1++) {

      int i_entity1 = (selected_entity1 == NULL) ? idx_entity1 : selected_entity1[i_part][idx_entity1] - 1;

      for (int idx_entity2 = entity1_entity2_idx[i_part][i_entity1]; idx_entity2 < entity1_entity2_idx[i_part][i_entity1+1]; idx_entity2++) {
        int i_entity2 = entity1_entity2[i_part][idx_entity2] - 1;

        if (entity2_flag[i_part][i_entity2] == 0) {
          entity2_flag    [i_part][i_entity2] = 1;
          selected_entity2[i_part][n_selected_entity2[i_part]++] = i_entity2 + 1;
        }
      }
    } // End loop on entities1

  } // End loop on parts


  /* Synchronize part boundaries */
  /* Possible optimizations :
   *  - create "sub" Part Comm Graph restricted to vertices incident to current pmne
   *  - exchange in PDM_STRIDE_VAR_INTERLACED to filter out vertices not incident to current pmne
   */
  int **send_flag = NULL;
  PDM_malloc(send_flag, n_part, int *);
  for (int i_part = 0; i_part < n_part; i_part++) {

    int *graph_entity2 = NULL;
    int graph_entity2_n = PDM_part_comm_graph_entity_graph_get(pcg_entity2,
                                                               i_part,
                                                              &graph_entity2,
                                                               PDM_OWNERSHIP_BAD_VALUE);

    PDM_malloc(send_flag[i_part], graph_entity2_n, int);
    for (int idx_entity2 = 0; idx_entity2 < graph_entity2_n; idx_entity2++) {
      int i_entity2 = graph_entity2[4*idx_entity2] - 1;
      send_flag[i_part][idx_entity2] = entity2_flag[i_part][i_entity2];
    }
  } // End loop on parts

  int **recv_flag = NULL;
  PDM_part_comm_graph_exch(pcg_entity2,
                           sizeof(int),
                           PDM_STRIDE_CST_INTERLACED,
                           1,
                           NULL,
                (void  **) send_flag,
                           NULL,
                (void ***) &recv_flag);


  for (int i_part = 0; i_part < n_part; i_part++) {
    int *graph_entity2 = NULL;
    int graph_entity2_n = PDM_part_comm_graph_entity_graph_get(pcg_entity2,
                                                               i_part,
                                                              &graph_entity2,
                                                               PDM_OWNERSHIP_BAD_VALUE);

    for (int idx_entity2 = 0; idx_entity2 < graph_entity2_n; idx_entity2++) {
      int i_entity2 = graph_entity2[4*idx_entity2] - 1;
      if (recv_flag[i_part][idx_entity2] && !entity2_flag[i_part][i_entity2]) {
        entity2_flag[i_part][i_entity2] = 1;
        selected_entity2[i_part][n_selected_entity2[i_part]++] = i_entity2 + 1;
      }
    }
    PDM_free(send_flag   [i_part]);
    PDM_free(recv_flag   [i_part]);
    PDM_free(entity2_flag[i_part]);

    PDM_realloc(selected_entity2[i_part], (*out_selected_entity2)[i_part], n_selected_entity2[i_part], int);

  } // End loop on parts
  PDM_free(send_flag);
  PDM_free(recv_flag);
  PDM_free(entity2_flag);
}


PDM_part_comm_graph_t *
PDM_part_comm_graph_concatenate
(
  PDM_MPI_Comm            comm,
  int                     n_pcg,
  PDM_part_comm_graph_t **pcgs
)
{
  /**
   * Check that all pcg have same partition number while getting nuplet size
   */
  int n_part             = -1;
  int concat_is_signed   =  0;
  int concat_nuplet_size =  0;
  for (int i_pcg=0; i_pcg<n_pcg; ++i_pcg) {

    if (i_pcg == 0) {
      n_part = pcgs[i_pcg]->n_part;
    }
    else {
      if (n_part != pcgs[i_pcg]->n_part) {
        PDM_error(__FILE__, __LINE__, 0, "pcg %d has not same n_part (=%d) as others (=%d)", i_pcg, pcgs[i_pcg]->n_part, n_part);
      }
    }
    concat_is_signed   = PDM_MAX(concat_is_signed  , pcgs[i_pcg]->is_signed);
    concat_nuplet_size = PDM_MAX(concat_nuplet_size, pcgs[i_pcg]->nuplet_size);
  }


  /**
   * Count number of concatenated graph entities
   */
  int  *concat_n_entity_graph = NULL;
  int **concat_pentity_graph  = NULL;
  int **concat_pentity_nuplet = NULL;
  PDM_malloc(concat_n_entity_graph, n_part, int  );
  PDM_malloc(concat_pentity_graph , n_part, int *);
  if (concat_nuplet_size>0) {
    PDM_malloc(concat_pentity_nuplet, n_part, int *);
  }

  for (int i_part=0; i_part<n_part; ++i_part) {

    concat_n_entity_graph[i_part] = 0;

    for (int i_pcg=0; i_pcg<n_pcg; ++i_pcg) {
      int *entity_graph  = NULL;
      int n_entity_graph = PDM_part_comm_graph_entity_graph_get(pcgs[i_pcg],
                                                                i_part,
                                                               &entity_graph,
                                                                PDM_OWNERSHIP_BAD_VALUE);

      concat_n_entity_graph[i_part] += n_entity_graph;
    }
  }


  /**
   * Create concatenate pcg arrays
   */
  for (int i_part=0; i_part<n_part; ++i_part) {

    PDM_malloc(concat_pentity_graph[i_part], 4*concat_n_entity_graph[i_part], int);
    if (concat_nuplet_size>0) {
      PDM_malloc(concat_pentity_nuplet[i_part], concat_nuplet_size*concat_n_entity_graph[i_part], int);
    }
    concat_n_entity_graph[i_part] = 0;

    for (int i_pcg=0; i_pcg<n_pcg; ++i_pcg) {

      int *entity_graph  = NULL;
      int *entity_nuplet = NULL;
      int n_entity_graph = PDM_part_comm_graph_entity_graph_get(pcgs[i_pcg],
                                                                i_part,
                                                               &entity_graph,
                                                                PDM_OWNERSHIP_BAD_VALUE);
      PDM_part_comm_graph_entity_nuplet_get(pcgs[i_pcg],
                                            i_part,
                                           &entity_nuplet,
                                            PDM_OWNERSHIP_BAD_VALUE);

      memcpy(&concat_pentity_graph[i_part][4*concat_n_entity_graph[i_part]], entity_graph, 4*n_entity_graph*sizeof(int));

      if (concat_nuplet_size>0) {

        int i_write = concat_n_entity_graph[i_part];

        for (int i_entity=0; i_entity<n_entity_graph; i_entity++) {

          for (int i_nuplet=0; i_nuplet<pcgs[i_pcg]->nuplet_size; i_nuplet++) {
            concat_pentity_nuplet[i_part][concat_nuplet_size*i_write+i_nuplet] = entity_nuplet[pcgs[i_pcg]->nuplet_size*i_entity + i_nuplet];
          }
          for (int i_nuplet=pcgs[i_pcg]->nuplet_size; i_nuplet<concat_nuplet_size; i_nuplet++) {
            concat_pentity_nuplet[i_part][concat_nuplet_size*i_write+i_nuplet] = 0;
          }

          i_write++;
        }
      }
      concat_n_entity_graph[i_part] += n_entity_graph;
    }
  }


  /**
   * Create pcg object
   */
  PDM_part_comm_graph_t *pcg = NULL;
  if (concat_nuplet_size>0) {
    pcg = PDM_part_comm_graph_with_nuplet_create(n_part,
                                                 concat_n_entity_graph,
                                                 concat_pentity_graph,
                                                 PDM_OWNERSHIP_KEEP,
                                                 concat_nuplet_size,
                                                 concat_pentity_nuplet,
                                                 PDM_OWNERSHIP_KEEP,
                                                 concat_is_signed,
                                                 comm);
  }
  else {
    pcg = PDM_part_comm_graph_create(n_part,
                                     concat_n_entity_graph,
                                     concat_pentity_graph,
                                     PDM_OWNERSHIP_KEEP,
                                     comm);
  }


  PDM_free(concat_n_entity_graph);
  PDM_free(concat_pentity_graph);
  if (concat_nuplet_size>0) {
    PDM_free(concat_pentity_nuplet);
  }

  return pcg;
}


void
PDM_part_comm_graph_split
(
  PDM_part_comm_graph_t  *pcg,
  const int               n_color,
  const int             **entity_color,
  PDM_part_comm_graph_t **split_pcgs
)
{
  PDM_MPI_Comm comm        = PDM_part_comm_graph_comm_get(pcg);
  int          n_part      = PDM_part_comm_graph_n_part_get(pcg);
  int          nuplet_size = PDM_part_comm_graph_nuplet_size(pcg);
  int          is_signed   = PDM_part_comm_graph_is_signed(pcg);

  int  **split_n_entity_graph = NULL;
  int ***split_entity_graph = NULL;
  int ***split_entity_nuplt = NULL;
  PDM_calloc(split_n_entity_graph, n_color, int  *);
  PDM_malloc(split_entity_graph  , n_color, int **);
  if (nuplet_size>0) {
    PDM_malloc(split_entity_nuplt, n_color, int **);
  }
  for (int i_color = 0; i_color < n_color; ++i_color) {
    PDM_calloc(split_n_entity_graph[i_color], n_part, int  );
    PDM_malloc(split_entity_graph  [i_color], n_part, int *);
    if (nuplet_size>0) {
      PDM_malloc(split_entity_nuplt[i_color], n_part, int *);
    }
  }

  for (int i_part=0; i_part<n_part; ++i_part) {
    int *entity_graph = NULL;
    int *entity_nuplt = NULL;
    int n_entity = PDM_part_comm_graph_entity_graph_get(pcg,
                                                        i_part,
                                                        &entity_graph,
                                                        PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_comm_graph_entity_nuplet_get(pcg,
                                          i_part,
                                          &entity_nuplt,
                                          PDM_OWNERSHIP_BAD_VALUE);

    for (int i_color = 0; i_color < n_color; ++i_color) {
      PDM_malloc(split_entity_graph[i_color][i_part], 4*n_entity, int);
      if (nuplet_size>0) {
        PDM_malloc(split_entity_nuplt[i_color][i_part], n_entity, int);
      }
    }

    for (int i_entity=0; i_entity<n_entity; ++i_entity) {
      int color = entity_color[i_part][i_entity];
      log_trace("i_entity = %d, color = %d\n", i_entity, color);
      if (color<0) {
        PDM_error(__FILE__, __LINE__, 0, "PDM_part_comm_graph_split: color[i_part=%d][i_entity=%d] = %d, but should be >0", i_part, i_entity, color);
      }
      if (color>=n_color) {
        PDM_error(__FILE__, __LINE__, 0, "PDM_part_comm_graph_split: color[i_part=%d][i_entity=%d] = %d, but should be < n_color (= %d)", i_part, i_entity, color, n_color);
      }
      int i_write = split_n_entity_graph[color][i_part];

      split_entity_graph[color][i_part][4*i_write  ] = entity_graph[4*i_entity  ];
      split_entity_graph[color][i_part][4*i_write+1] = entity_graph[4*i_entity+1];
      split_entity_graph[color][i_part][4*i_write+2] = entity_graph[4*i_entity+2];
      split_entity_graph[color][i_part][4*i_write+3] = entity_graph[4*i_entity+3];

      if (nuplet_size>0) {
        for (int i_nuplet=0; i_nuplet<nuplet_size; ++i_nuplet) {
          split_entity_nuplt[color][i_part][nuplet_size*i_write+i_nuplet] = entity_nuplt[nuplet_size*i_entity+i_nuplet];
        }
      }

      split_n_entity_graph[color][i_part]++;

    }
  }

  for (int color = 0; color < n_color; ++color) {
    if (nuplet_size>0) {
      split_pcgs[color] = PDM_part_comm_graph_with_nuplet_create(n_part,
                                                                 split_n_entity_graph[color],
                                                                 split_entity_graph[color],
                                                                 PDM_OWNERSHIP_KEEP,
                                                                 nuplet_size,
                                                                 split_entity_nuplt[color],
                                                                 PDM_OWNERSHIP_KEEP,
                                                                 is_signed,
                                                                 comm);
    }
    else {
      split_pcgs[color] = PDM_part_comm_graph_create(n_part,
                                                     split_n_entity_graph[color],
                                                     split_entity_graph[color],
                                                     PDM_OWNERSHIP_KEEP,
                                                     comm);
    }
  }

  for (int i_color = 0; i_color < n_color; ++i_color) {
    PDM_free(split_n_entity_graph[i_color]);
    PDM_free(split_entity_graph  [i_color]);
    if (nuplet_size>0) {
      PDM_free(split_entity_nuplt[i_color]);
    }
  }
  PDM_free(split_n_entity_graph);
  PDM_free(split_entity_graph);
  PDM_free(split_entity_nuplt);

}


#ifdef __cplusplus
}
#endif /* __cplusplus */
