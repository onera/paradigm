/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_array.h"
#include "pdm_distrib.h"
#include "pdm_gnum.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include "pdm_part_graph_dual.h"
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


/*=============================================================================
 * Public function definitions
 *============================================================================*/


void
PDM_part_assembly_dual_graph
(
  PDM_MPI_Comm            comm,
  int                     n_part,
  int                    *n_node,
  int                    *n_arc,
  int                   **select_node,
  int                   **node_arc_idx,
  int                   **node_arc,
  int                   **arc_node_idx,
  int                   **arc_node,
  int                   **node_weight,
  int                   **arc_weight,
  PDM_part_comm_graph_t  *pcg_node,
  PDM_part_comm_graph_t  *pcg_arc,
  int                    *out_n_tot_node,
  PDM_g_num_t           **out_gnode_node_idx,
  PDM_g_num_t           **out_gnode_node,
  int                   **out_gnode_weight,
  int                   **out_garc_weight,
  PDM_g_num_t           **out_distrib_node,
  int                  ***out_part_to_graph
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int internal_pcg_node = 0;
  if(pcg_node == NULL) {
    internal_pcg_node = 1;
    int  *pn_node_graph = NULL;
    int **pnode_graph   = NULL;
    PDM_malloc(pn_node_graph, n_part, int  );
    PDM_malloc(pnode_graph  , n_part, int *);
    for(int i_part = 0; i_part < n_part; ++i_part) {
      pn_node_graph[i_part] = 0;
      PDM_malloc(pnode_graph[i_part], 4 * pn_node_graph[i_part], int);
    }

    pcg_node = PDM_part_comm_graph_create(n_part,
                                          pn_node_graph,
                                          pnode_graph,
                                          PDM_OWNERSHIP_KEEP,
                                          comm);

    PDM_free(pn_node_graph);
    PDM_free(pnode_graph);
  }

  /*
   * Manage empty cases / special cases
   */
  int  *pn_node_graph = NULL;
  int  *pn_arc_graph  = NULL;
  int **pnode_graph   = NULL;
  int **parc_graph    = NULL;

  PDM_malloc(pn_node_graph, n_part, int  );
  PDM_malloc(pn_arc_graph , n_part, int  );
  PDM_malloc(pnode_graph  , n_part, int *);
  PDM_malloc(parc_graph   , n_part, int *);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    if(pcg_node != NULL) {
      pn_node_graph[i_part] = PDM_part_comm_graph_entity_graph_get(pcg_node,
                                                                   i_part,
                                                                   &pnode_graph[i_part],
                                                                   PDM_OWNERSHIP_BAD_VALUE);
    } else {
      pn_node_graph[i_part] = 0;
      pnode_graph  [i_part] = NULL;
    }

    if(pcg_arc != NULL) {
      pn_arc_graph[i_part] = PDM_part_comm_graph_entity_graph_get(pcg_arc,
                                                                  i_part,
                                                                  &parc_graph[i_part],
                                                                  PDM_OWNERSHIP_BAD_VALUE);
    } else {
      pn_arc_graph[i_part] = 0;
      parc_graph  [i_part] = NULL;
    }
  }

  /*
   * En noeuds centrés, pas besoin de synchronisés les arcs (car il sont deja commun a chaque partition connectés)
   * On peu eviter les échanges sur les arc en mettant pcg_arc == NULL
   * De la même manière en cellules centrés, à priori pas besoin de syncho les celluls, on peut faire pcg_node == NULL
   */
  PDM_part_comm_graph_t *pcg_subnode = NULL;
  PDM_g_num_t** pnode_ln_to_gn = NULL;
  PDM_malloc(pnode_ln_to_gn, n_part, PDM_g_num_t *);
  if(select_node != NULL) {

    int  *pn_select_node   = NULL;
    int **pnode_old_to_new = NULL;
    PDM_malloc(pn_select_node  , n_part, int  );
    PDM_malloc(pnode_old_to_new, n_part, int *);
    for(int i_part = 0; i_part < n_part; ++i_part) {
      pn_select_node[i_part] = 0;
      PDM_malloc(pnode_old_to_new[i_part], n_node[i_part], int);
      for(int i = 0; i < n_node[i_part]; ++i) {
        pnode_old_to_new[i_part][i] = -1; // Unused if all thing are going well
        if(select_node[i_part][i] == 1) {
          pnode_old_to_new[i_part][i] = pn_select_node[i_part];
          pn_select_node[i_part]++;
        }
      }
    }

    // We need to recreate part_comm_graph but for subset :
    //   1/ On echange les selected avant
    //   2/ On filtre localement et on refait le part_comm_graph
    PDM_part_comm_graph_all_reduce(pcg_node,
                                   PDM_MPI_INT,
                                   PDM_MPI_MAX,
               (unsigned char **)  select_node);

    int  *pn_sub_node_graph = NULL;
    int **psub_node_graph   = NULL;

    PDM_malloc(pn_sub_node_graph, n_part, int  );
    PDM_malloc(psub_node_graph  , n_part, int *);

    for(int i_part = 0; i_part < n_part; ++i_part) {

      int *node_graph = NULL;
      int n_node_graph = PDM_part_comm_graph_entity_graph_get(pcg_node,
                                                              i_part,
                                                              &node_graph,
                                                              PDM_OWNERSHIP_BAD_VALUE);

      pn_sub_node_graph[i_part] = 0;
      PDM_malloc(psub_node_graph[i_part], 4 * n_node_graph, int);
      for(int i = 0; i < n_node_graph; ++i) {
        int i_node = node_graph[4*i]-1;
        if(select_node[i_part][i_node] == 1) {
          psub_node_graph[i_part][4*pn_sub_node_graph[i_part]  ] = node_graph[4*i  ];
          psub_node_graph[i_part][4*pn_sub_node_graph[i_part]+1] = node_graph[4*i+1];
          psub_node_graph[i_part][4*pn_sub_node_graph[i_part]+2] = node_graph[4*i+2];
          psub_node_graph[i_part][4*pn_sub_node_graph[i_part]+3] = node_graph[4*i+3];
          pn_sub_node_graph[i_part]++;
        }
      }
    }

    pcg_subnode = PDM_part_comm_graph_create(n_part,
                                             pn_sub_node_graph,
                                             psub_node_graph,
                                             PDM_OWNERSHIP_KEEP,
                                             comm);

    /* Update sub numbering */
    PDM_part_comm_graph_reorder(pcg_subnode,
                                pnode_old_to_new);

    /* Generate gnum only on subset */
    PDM_gen_gnum_t *gen_gnum_node = PDM_gnum_create(3, 1, PDM_TRUE, 1e-6, comm, PDM_OWNERSHIP_KEEP);

    PDM_gnum_set_from_part_comm_graph(gen_gnum_node,
                                      pn_select_node,
                                      pcg_subnode);
    PDM_gnum_compute(gen_gnum_node);
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_malloc(pnode_ln_to_gn[i_part], n_node[i_part], PDM_g_num_t);

      PDM_g_num_t* psub_node_ln_to_gn = PDM_gnum_get(gen_gnum_node, i_part);

      /* Hack here - We put -1 on non selected entities */
      int idx_read = 0;
      for(int i = 0; i < n_node[i_part]; ++i) {
        if(select_node[i_part][i] == 0) {
          pnode_ln_to_gn[i_part][i] = -1;
        } else {
          pnode_ln_to_gn[i_part][i] = psub_node_ln_to_gn[idx_read++];
        }
      }
      PDM_free(pnode_old_to_new[i_part]);
      if(0 == 1) {
        PDM_log_trace_array_long(pnode_ln_to_gn[i_part], n_node[i_part], "pnode_ln_to_gn ::");
      }
    }
    PDM_free(pnode_old_to_new);
    PDM_gnum_free(gen_gnum_node);

    PDM_part_comm_graph_free(pcg_subnode);

    PDM_free(pn_sub_node_graph);
    PDM_free(psub_node_graph  );
    PDM_free(pn_select_node);
  } else {

    /*
     * En cellule centré, le pcg_node == NULL => A gerer
     */
    PDM_gen_gnum_t *gen_gnum_node = PDM_gnum_create(3, 1, PDM_TRUE, 1e-6, comm, PDM_OWNERSHIP_USER);

    PDM_gnum_set_from_part_comm_graph(gen_gnum_node,
                                      n_node,
                                      pcg_node);
    PDM_gnum_compute(gen_gnum_node);
    for(int i_part = 0; i_part < n_part; ++i_part) {
      pnode_ln_to_gn[i_part] = PDM_gnum_get(gen_gnum_node, i_part);
      PDM_log_trace_array_long(pnode_ln_to_gn[i_part], n_node[i_part], "pnode_ln_to_gn ::");
    }
    PDM_gnum_free(gen_gnum_node);
  }

  /*
   *
   */
  int **is_owner_node = NULL;
  PDM_malloc(is_owner_node, n_part, int *);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    const int *is_owner = PDM_part_comm_graph_owner_get(pcg_node, i_part);
    is_owner_node[i_part] = PDM_array_const_int(n_node[i_part], 1);
    int  n_node_graph  = pn_node_graph[i_part];
    int *_node_graph   = pnode_graph  [i_part];
    for(int i_graph_node = 0; i_graph_node < n_node_graph; ++i_graph_node) {
      int i_node = _node_graph[4*i_graph_node]-1;
      is_owner_node[i_part][i_node] = is_owner[i_graph_node];
    }
  }

  int n_tot_node = 0;
  if(select_node == NULL) {
    // Hack here : we use is_owner as old_to_new (Caution, this contains shift for partition also)
    for(int i_part = 0; i_part < n_part; ++i_part) {
      for(int i = 0; i < n_node[i_part]; ++i) {
        if(is_owner_node[i_part][i] == 1) {
          is_owner_node[i_part][i] = n_tot_node++;
        } else {
          is_owner_node[i_part][i] = -1;
        }
      }
    }
  } else {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      for(int i = 0; i < n_node[i_part]; ++i) {
        if(is_owner_node[i_part][i] == 1 && select_node[i_part][i] == 1) {
          is_owner_node[i_part][i] = n_tot_node++;
        } else {
          is_owner_node[i_part][i] = -1;
        }
      }
    }
  }

  PDM_g_num_t* distrib_node = PDM_compute_entity_distribution(comm, n_tot_node);

  int         **send_node_n     = NULL;
  PDM_g_num_t **send_node       = NULL;
  int         **send_weight     = NULL;
  int         **send_arc_node_n = NULL;
  PDM_g_num_t **send_arc_node   = NULL;
  int         **send_arc_weight = NULL;
  PDM_malloc(send_node_n    , n_part, int         *);
  PDM_malloc(send_node      , n_part, PDM_g_num_t *);
  PDM_malloc(send_weight    , n_part, int         *);
  PDM_malloc(send_arc_node_n, n_part, int         *);
  PDM_malloc(send_arc_node  , n_part, PDM_g_num_t *);
  PDM_malloc(send_arc_weight, n_part, int         *);
  for(int i_part = 0; i_part < n_part; ++i_part) {

    int  n_arc_graph   = pn_arc_graph [i_part];
    int  n_node_graph  = pn_node_graph[i_part];
    int *_arc_graph    = parc_graph   [i_part];
    int *_node_graph   = pnode_graph  [i_part];
    int *_node_arc_idx = node_arc_idx [i_part];
    int *_node_arc     = node_arc     [i_part];
    int *_arc_node_idx = arc_node_idx [i_part];
    int *_arc_node     = arc_node     [i_part];

    const int *is_owner = PDM_part_comm_graph_owner_get(pcg_node, i_part);

    /*
     * Count
     */
    PDM_malloc(send_node_n[i_part], n_node_graph, int);
    for(int i_graph_node = 0; i_graph_node < n_node_graph; ++i_graph_node) {
      int i_node = _node_graph[4*i_graph_node]-1;
      send_node_n[i_part][i_graph_node] = 0;
      if(is_owner[i_graph_node] == 1) {
        continue;
      }
      if(pnode_ln_to_gn[i_part][i_node] == -1) { // Donc pas selectioné
        continue;
      }

      for(int idx_arc = _node_arc_idx[i_node]; idx_arc < _node_arc_idx[i_node+1]; ++idx_arc) {
        int i_arc = PDM_ABS(_node_arc[idx_arc]) - 1;
        send_node_n[i_part][i_graph_node] += _arc_node_idx[i_arc+1] - _arc_node_idx[i_arc];
      }
    }

    PDM_malloc(send_arc_node_n[i_part], n_arc_graph, int);
    for(int i_graph_arc = 0; i_graph_arc < n_arc_graph; ++i_graph_arc) {
      int i_arc = _arc_graph[4*i_graph_arc]-1;
      send_arc_node_n[i_part][i_graph_arc] = _arc_node_idx[i_arc+1] - _arc_node_idx[i_arc];
    }

    /*
     * Allocate
     */
    int *send_node_idx = NULL;
    PDM_malloc(send_node_idx, n_node_graph+1, int);
    send_node_idx[0] = 0;
    for(int i_graph_node = 0; i_graph_node < n_node_graph; ++i_graph_node) {
      send_node_idx[i_graph_node+1] = send_node_idx[i_graph_node] + send_node_n[i_part][i_graph_node];
      send_node_n  [i_part][i_graph_node] = 0;
    }
    PDM_malloc(send_node  [i_part], send_node_idx[n_node_graph], PDM_g_num_t);
    PDM_malloc(send_weight[i_part], send_node_idx[n_node_graph], int        );

    int *send_arc_node_idx = NULL;
    PDM_malloc(send_arc_node_idx, n_arc_graph+1, int);
    send_arc_node_idx[0] = 0;
    for(int i_graph_arc = 0; i_graph_arc < n_arc_graph; ++i_graph_arc) {
      send_arc_node_idx[i_graph_arc+1] = send_arc_node_idx[i_graph_arc] + send_arc_node_n[i_part][i_graph_arc];
      send_arc_node_n  [i_part][i_graph_arc] = 0;
    }
    PDM_malloc(send_arc_node  [i_part], send_arc_node_idx[n_arc_graph], PDM_g_num_t);
    PDM_malloc(send_arc_weight[i_part], send_arc_node_idx[n_arc_graph], int        );

    /*
     * Fill
     */
    for(int i_graph_node = 0; i_graph_node < n_node_graph; ++i_graph_node) {
      int i_node = _node_graph[4*i_graph_node]-1;
      if(is_owner[i_graph_node] == 1) {
        continue;
      }
      if(pnode_ln_to_gn[i_part][i_node] == -1) { // Donc pas selectioné
        continue;
      }
      for(int idx_arc = _node_arc_idx[i_node]; idx_arc < _node_arc_idx[i_node+1]; ++idx_arc) {
        int i_arc = PDM_ABS(_node_arc[idx_arc]) - 1;
        for(int idx_node = _arc_node_idx[i_arc]; idx_node < _arc_node_idx[i_arc+1]; ++idx_node) {
          int i_node2 = PDM_ABS(_arc_node[idx_node])-1;
          int idx_write = send_node_idx[i_graph_node] + send_node_n[i_part][i_graph_node]++;
          send_node  [i_part][idx_write] = pnode_ln_to_gn[i_part][i_node2];
          send_weight[i_part][idx_write] = arc_weight    [i_part][i_arc];
        }
      }
    }

    for(int i_graph_arc = 0; i_graph_arc < n_arc_graph; ++i_graph_arc) {
      int i_arc = _arc_graph[4*i_graph_arc]-1;
      for(int idx_node = _arc_node_idx[i_arc]; idx_node < _arc_node_idx[i_arc+1]; ++idx_node) {
        int i_node = PDM_ABS(_arc_node[idx_node])-1;
        int idx_write = send_arc_node_idx[i_graph_arc] + send_arc_node_n[i_part][i_graph_arc]++;
        log_trace("idx_write = %i (%i) \n", idx_write, send_arc_node_idx[n_arc_graph]);
        log_trace("pnode_ln_to_gn[%i] = %i \n", i_node, i_node);
        send_arc_node  [i_part][idx_write] = pnode_ln_to_gn[i_part][i_node];
        send_arc_weight[i_part][idx_write] = arc_weight    [i_part][i_arc];
      }
    }

    PDM_free(send_arc_node_idx);
    PDM_free(send_node_idx);
  }

  /*
   * Exchange
   */
  int         **recv_node_n = NULL;
  PDM_g_num_t **recv_node   = NULL;
  PDM_part_comm_graph_exch(pcg_node,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_VAR_INTERLACED,
                           1,
                           send_node_n,
                 (void **) send_node,
                           &recv_node_n,
                (void ***) &recv_node);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(recv_node_n[i_part]);
  }
  PDM_free(recv_node_n);

  int **recv_weight = NULL;
  PDM_part_comm_graph_exch(pcg_node,
                           sizeof(int),
                           PDM_STRIDE_VAR_INTERLACED,
                           1,
                           send_node_n,
                 (void **) send_weight,
                           &recv_node_n,
                (void ***) &recv_weight);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(send_node  [i_part]);
    PDM_free(send_node_n[i_part]);
    PDM_free(send_weight[i_part]);
  }
  PDM_free(send_node  );
  PDM_free(send_node_n);
  PDM_free(send_weight);


  int         **recv_arc_node_n = NULL;
  PDM_g_num_t **recv_arc_node   = NULL;
  PDM_part_comm_graph_exch(pcg_arc,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_VAR_INTERLACED,
                           1,
                           send_arc_node_n,
                 (void **) send_arc_node,
                           &recv_arc_node_n,
                (void ***) &recv_arc_node);


  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(recv_arc_node_n[i_part]);
  }
  PDM_free(recv_arc_node_n);

  int **recv_arc_weight = NULL;
  PDM_part_comm_graph_exch(pcg_arc,
                           sizeof(int),
                           PDM_STRIDE_VAR_INTERLACED,
                           1,
                           send_arc_node_n,
                 (void **) send_arc_weight,
                           &recv_arc_node_n,
                (void ***) &recv_arc_weight);


  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(send_arc_node  [i_part]);
    PDM_free(send_arc_node_n[i_part]);
    PDM_free(send_arc_weight[i_part]);
  }
  PDM_free(send_arc_node  );
  PDM_free(send_arc_node_n);
  PDM_free(send_arc_weight);

  /*
   * Node are by construction ordered by increasing gnum numbering
   *   - We compute dual graph and shift properly
   *   - Add also weight if any
   */
  int max_size   = 0;
  int *node_node_n = PDM_array_zeros_int(n_tot_node);
  for(int i_part = 0; i_part < n_part; ++i_part) {

    int  n_arc_graph   = pn_arc_graph [i_part];
    int  n_node_graph  = pn_node_graph[i_part];
    int *_arc_graph    = parc_graph   [i_part];
    int *_node_graph   = pnode_graph  [i_part];
    int *_node_arc_idx = node_arc_idx [i_part];
    int *_node_arc     = node_arc     [i_part];
    int *_arc_node_idx = arc_node_idx [i_part];
    int *_arc_node     = arc_node     [i_part];

    // Rebuild dual graph and update into gnum
    for(int i_node = 0; i_node < n_node[i_part]; ++i_node) {
      int l_node = is_owner_node[i_part][i_node];
      if(l_node == -1) {
        continue;
      }
      if(pnode_ln_to_gn[i_part][i_node] == -1) { // Donc pas selectioné
        continue;
      }
      for(int idx_arc = _node_arc_idx[i_node]; idx_arc < _node_arc_idx[i_node+1]; ++idx_arc) {
        int i_arc = PDM_ABS(_node_arc[idx_arc])-1;
        for(int idx_node = _arc_node_idx[i_arc]; idx_node < _arc_node_idx[i_arc+1]; ++idx_node) {
          int i_node2 = PDM_ABS(_arc_node[idx_node])-1;
          if(pnode_ln_to_gn[i_part][i_node2] == -1) { // Donc pas selectioné
            continue;
          }
          node_node_n[l_node] += 1;
        }
        max_size += _arc_node_idx[i_arc+1] - _arc_node_idx[i_arc];
      }
    }

    for(int i_graph_node = 0; i_graph_node < n_node_graph; ++i_graph_node) {
      int i_node = _node_graph[4*i_graph_node]-1;
      int l_node = is_owner_node[i_part][i_node];
      if(l_node == -1) {
        continue;
      }
      node_node_n[l_node] += recv_node_n[i_part][i_graph_node];
      max_size += recv_node_n[i_part][i_graph_node];
    }

    for(int i_graph_arc = 0; i_graph_arc < n_arc_graph; ++i_graph_arc) {
      int i_arc = _arc_graph[4*i_graph_arc]-1;
      for(int idx_node = _arc_node_idx[i_arc]; idx_node < _arc_node_idx[i_arc+1]; ++idx_node) {
        int i_node = PDM_ABS(_arc_node[idx_node])-1;
        int l_node = is_owner_node[i_part][i_node];
        if(l_node == -1) {
          continue;
        }
        node_node_n[l_node] += recv_arc_node_n[i_part][i_graph_arc];
        max_size += recv_arc_node_n[i_part][i_graph_arc];
      }
    }
  }

  /*
   * Count
   */
  PDM_g_num_t *node_node_idx = NULL;
  PDM_malloc(node_node_idx, n_tot_node+1, PDM_g_num_t);

  int max_node_node = 0;
  node_node_idx[0] = 0;
  for(int i = 0; i < n_tot_node; ++i) {
    node_node_idx[i+1] = node_node_idx[i] + node_node_n[i];
    max_node_node = PDM_MAX(max_node_node, node_node_n[i]);
    node_node_n[i] = 0;
  }

  /*
   * Manage node_weight
   */
  int *gnode_weight = NULL;
  if(node_weight != NULL) {
    PDM_malloc(gnode_weight, n_tot_node, int);
    for(int i_part = 0; i_part < n_part; ++i_part) {
      for(int i_node = 0; i_node < n_node[i_part]; ++i_node) {
        int l_node = is_owner_node[i_part][i_node];
        if(l_node == -1) {
          continue;
        }
        if(pnode_ln_to_gn[i_part][i_node] == -1) { // Donc pas selectioné
          continue;
        }
        gnode_weight[l_node] = node_weight[i_part][i_node];
      }
    }
  }

  /*
   * Fill
   */
  PDM_g_num_t *gnode_node   = NULL;
  int         *garc_weight  = NULL;
  PDM_malloc(gnode_node  , max_size  , PDM_g_num_t);
  PDM_malloc(garc_weight , max_size  , int        );
  for(int i_part = 0; i_part < n_part; ++i_part) {

    int  n_arc_graph   = pn_arc_graph [i_part];
    int  n_node_graph  = pn_node_graph[i_part];
    int *_arc_graph    = parc_graph   [i_part];
    int *_node_graph   = pnode_graph  [i_part];
    int *_node_arc_idx = node_arc_idx [i_part];
    int *_node_arc     = node_arc     [i_part];
    int *_arc_node_idx = arc_node_idx [i_part];
    int *_arc_node     = arc_node     [i_part];

    // Rebuild dual graph and update into gnum
    for(int i_node = 0; i_node < n_node[i_part]; ++i_node) {
      int l_node = is_owner_node[i_part][i_node];
      if(l_node == -1) {
        continue;
      }
      if(pnode_ln_to_gn[i_part][i_node] == -1) { // Donc pas selectioné
        continue;
      }
      for(int idx_arc = _node_arc_idx[i_node]; idx_arc < _node_arc_idx[i_node+1]; ++idx_arc) {
        int i_arc = PDM_ABS(_node_arc[idx_arc])-1;
        for(int idx_node = _arc_node_idx[i_arc]; idx_node < _arc_node_idx[i_arc+1]; ++idx_node) {
          int i_node2 = PDM_ABS(_arc_node[idx_node])-1;
          if(pnode_ln_to_gn[i_part][i_node2] == -1) { // Donc pas selectioné
            continue;
          }
          int idx_write = node_node_idx[l_node] + node_node_n[l_node]++;
          gnode_node [idx_write] = pnode_ln_to_gn[i_part][i_node2];
          garc_weight[idx_write] = arc_weight    [i_part][i_arc];
        }
      }
    }

    int idx_read = 0;
    for(int i_graph_node = 0; i_graph_node < n_node_graph; ++i_graph_node) {
      int i_node = _node_graph[4*i_graph_node]-1;
      int l_node = is_owner_node[i_part][i_node];
      if(l_node == -1) {
        idx_read += recv_node_n[i_part][i_graph_node];
        continue;
      }
      for(int j = 0; j < recv_node_n[i_part][i_graph_node]; ++j) {
        int idx_write = node_node_idx[l_node] + node_node_n[l_node]++;
        gnode_node [idx_write] = recv_node  [i_part][idx_read];
        garc_weight[idx_write] = recv_weight[i_part][idx_read];
        idx_read++;
      }
    }

    /* From arc */
    idx_read = 0;
    for(int i_graph_arc = 0; i_graph_arc < n_arc_graph; ++i_graph_arc) {
      int i_arc = _arc_graph[4*i_graph_arc]-1;
      for(int idx_node = _arc_node_idx[i_arc]; idx_node < _arc_node_idx[i_arc+1]; ++idx_node) {
        int i_node = PDM_ABS(_arc_node[idx_node])-1;
        int l_node = is_owner_node[i_part][i_node];
        if(l_node == -1) {
          continue;
        }
        for(int j = 0; j < recv_arc_node_n[i_part][i_graph_arc]; ++j) {
          int idx_write = node_node_idx[l_node] + node_node_n[l_node]++;
          gnode_node [idx_write] = recv_arc_node  [i_part][idx_read+j];
          garc_weight[idx_write] = recv_arc_weight[i_part][idx_read+j];
        }
      }
      idx_read += recv_arc_node_n[i_part][i_graph_arc];
    }
  }

  PDM_free(node_node_n);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(recv_node      [i_part]);
    PDM_free(recv_node_n    [i_part]);
    PDM_free(recv_weight    [i_part]);
    PDM_free(recv_arc_node  [i_part]);
    PDM_free(recv_arc_node_n[i_part]);
    PDM_free(recv_arc_weight[i_part]);
  }
  PDM_free(recv_node      );
  PDM_free(recv_node_n    );
  PDM_free(recv_weight    );
  PDM_free(recv_arc_node_n);
  PDM_free(recv_arc_node  );
  PDM_free(recv_arc_weight);

  /*
   * At this stage we have dual graph but :
   *   - Not unique
   *   - Not sorted
   *   - Not compacted (with owner / not owner)
   */
  int *lorder  = NULL;
  int *lweight = NULL;
  PDM_malloc(lorder , max_node_node, int);
  PDM_malloc(lweight, max_node_node, int);

  int idx_read  = 0;
  int idx_write = 0;
  for(int i_node = 0; i_node < n_tot_node; ++i_node) {

    int beg   = idx_read;
    int end   = node_node_idx[i_node+1];
    int n_adj = end - beg;

    for(int j = 0; j < n_adj; ++j) {
      lorder[j] = j;
    }

    int n_unique = PDM_inplace_unique_long_and_order(&gnode_node[beg], lorder, 0, n_adj-1);

    // Copy weight (mandatory because we copy in place)
    for(int i = 0; i < n_adj; ++i) {
      lweight[i] = garc_weight[beg+i];
    }

    // Tassage + move weight
    PDM_g_num_t gnum = i_node + distrib_node[i_rank] + 1;
    for(int i = 0; i < n_unique; ++i) {
      if(gnode_node[beg+i] != gnum && gnode_node[beg+i] != -1) {
        gnode_node [idx_write] = gnode_node[beg+i];
        garc_weight[idx_write] = lweight[lorder[i]];
        idx_write++;
      }
    }

    idx_read = node_node_idx[i_node+1];

    node_node_idx[i_node+1] = idx_write;

  }

  PDM_realloc(gnode_node , gnode_node , node_node_idx[n_tot_node], PDM_g_num_t);
  PDM_realloc(garc_weight, garc_weight, node_node_idx[n_tot_node], int        );

  if(0 == 1) {
    log_trace("gnode_node ----- \n");
    for(int i = 0; i < n_tot_node; ++i) {
      log_trace("ln_to_gn = "PDM_FMT_G_NUM" \n", distrib_node[i_rank]+i+1);
      for(int j = node_node_idx[i]; j < node_node_idx[i+1]; ++j) {
        // log_trace(""PDM_FMT_G_NUM" (%i) ", gnode_node[j], garc_weight[j]);
        log_trace(""PDM_FMT_G_NUM" ", gnode_node[j]);
      }
      log_trace("\n");
    }
    PDM_log_trace_connectivity_long(node_node_idx, gnode_node, n_tot_node, "node_node ::");
  }

  /*
   * Shift to zero
   */
  for(int i = 0; i < n_tot_node; ++i) {
    for(int j = node_node_idx[i]; j < node_node_idx[i+1]; ++j) {
      gnode_node[j] -= 1;
    }
  }

  /*
   * Free
   */
  PDM_free(lorder);
  PDM_free(lweight);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pnode_ln_to_gn[i_part]);
  }
  PDM_free(pnode_ln_to_gn);
  PDM_free(pn_node_graph);
  PDM_free(pnode_graph  );
  PDM_free(pn_arc_graph );
  PDM_free(parc_graph   );

  if(internal_pcg_node == 1) {
    PDM_part_comm_graph_free(pcg_node);
  }

  /*
   * Fix output
   */
  *out_n_tot_node     = n_tot_node;
  *out_gnode_node_idx = node_node_idx;
  *out_gnode_node     = gnode_node;
  *out_gnode_weight   = gnode_weight;
  *out_garc_weight    = garc_weight;
  *out_distrib_node   = distrib_node;
  *out_part_to_graph  = is_owner_node;
}


void
PDM_transfer_entity1_part_id_to_entity2_part_id
(
  int    n_part,
  int   *pn_entity1,
  int  **entity1_part_id,
  int   *pn_entity2,
  int  **pentity2_entity1_idx,
  int  **pentity2_entity1,
  int ***out_entity2_part_id
)
{
  int max_connectivity = 0;
  for(int i_part = 0; i_part < n_part; ++i_part) {
    for(int i_entity2 = 0; i_entity2 < pn_entity1[i_part]; ++i_entity2) {
      max_connectivity = PDM_MAX(max_connectivity, pentity2_entity1_idx[i_part][i_entity2+1] - pentity2_entity1_idx[i_part][i_entity2]);
    }
  }

  int *lpart_id = NULL;
  PDM_malloc(lpart_id, max_connectivity, int);

  int **entity2_part_id = NULL;
  PDM_malloc(entity2_part_id, n_part, int *);
  for(int i_part = 0; i_part < n_part; ++i_part) {

    PDM_malloc(entity2_part_id[i_part], pn_entity2[i_part], int);
    for(int i_entity2 = 0; i_entity2 < pn_entity2[i_part]; ++i_entity2) {

      int beg = pentity2_entity1_idx[i_part][i_entity2];
      int end = pentity2_entity1_idx[i_part][i_entity2+1];
      int n_adj = end - beg;

      int idx_write = 0;
      for(int idx_entity2 = beg; idx_entity2 < end; ++idx_entity2) {
        int i_entity1 = PDM_ABS(pentity2_entity1[i_part][idx_entity2])-1;
        lpart_id[idx_write++] = entity1_part_id[i_part][i_entity1];
      }
      PDM_sort_int(lpart_id, NULL, n_adj);

      int current_id    = lpart_id[0];
      int winner_id     = lpart_id[0];
      int current_count = 0;
      int max_count     = 0;
      for(int i = 0; i < n_adj; ++i) {
        if(lpart_id[i] == current_id) {
          current_count++;
        } else {
          if (current_count > max_count) {
            max_count = current_count;
            winner_id = current_id;
          }
          current_id = lpart_id[i];
          current_count = 1;
        }
      }
      // Last
      if (current_count > max_count) {
        winner_id = current_id;
      }
      entity2_part_id[i_part][i_entity2] = winner_id;
    }
  }

  PDM_free(lpart_id);

  *out_entity2_part_id = entity2_part_id;
}



#ifdef __cplusplus
}
#endif /* __cplusplus */
