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
  int                   **out_gnode_node_idx,
  PDM_g_num_t           **out_gnode_node,
  int                   **out_garc_weight,
  PDM_g_num_t           **out_distrib_node
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);
  /*
   * To avoid to much computation, we directly pre-alloc the node-node :
   *   - We need to know size with shared node (only owner will take graph)
   *   - Only compute node global numbering ?
   */

  /*
   * En noeuds centrés, pas besoin de synchronisés les arcs (car il sont deja commun a chaque partition connectés)
   * On peu eviter les échanges sur les arc en mettant pcg_arc == NULL
   * De la même manière en cellules centrés, à priori pas besoin de syncho les celluls, on peut faire pcg_node == NULL
   */
  PDM_part_comm_graph_t *pcg_subnode = NULL;
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
    PDM_free(pn_select_node);

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


    PDM_free(pn_sub_node_graph);
    PDM_free(psub_node_graph  );
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
   * En cellule centré, le pcg_node == NULL => A gerer
   */
  PDM_gen_gnum_t *gen_gnum_node = PDM_gnum_create(3, 1, PDM_TRUE, 1e-6, comm, PDM_OWNERSHIP_USER);

  // if(pcg_node == NULL) {
  // } else
  PDM_gnum_set_from_part_comm_graph(gen_gnum_node,
                                    n_node,
                                    pcg_node);
  PDM_gnum_compute(gen_gnum_node);
  PDM_g_num_t** pnode_ln_to_gn = NULL;
  PDM_malloc(pnode_ln_to_gn, n_part, PDM_g_num_t *);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pnode_ln_to_gn[i_part] = PDM_gnum_get(gen_gnum_node, i_part);
  }
  PDM_gnum_free(gen_gnum_node);

  PDM_g_num_t* distrib_node = PDM_compute_uniform_entity_distribution_from_partition(comm,
                                                                                     n_part,
                                                                                     n_node,
                                                              (const PDM_g_num_t **) pnode_ln_to_gn);


  PDM_g_num_t **psub_node_ln_to_gn = NULL;

  if(select_node != NULL) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int idx_read = 0;
      for(int i = 0; i < n_node[i_part]; ++i) {
        if(select_node[i_part][i] == 0) {
          pnode_ln_to_gn[i_part][i] = -1;
        } else {
          pnode_ln_to_gn[i_part][i] = psub_node_ln_to_gn[i_part][idx_read++];
        }
      }
    }
  }


  /*
   *
   */
  int         **send_node_n     = NULL;
  PDM_g_num_t **send_node       = NULL;
  int         **send_weight     = NULL;
  int         **send_arc_node_n = NULL;
  PDM_g_num_t **send_arc_node   = NULL;
  int         **send_arc_weight = NULL;
  int         **is_owner_node   = NULL;
  PDM_malloc(send_node_n    , n_part, int         *);
  PDM_malloc(send_node      , n_part, PDM_g_num_t *);
  PDM_malloc(send_weight    , n_part, int         *);
  PDM_malloc(send_arc_node_n, n_part, int         *);
  PDM_malloc(send_arc_node  , n_part, PDM_g_num_t *);
  PDM_malloc(send_arc_weight, n_part, int         *);
  PDM_malloc(is_owner_node  , n_part, int         *);


  int n_tot_node = 0;
  if(select_node == NULL) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      n_tot_node += n_node[i_part];
    }
  } else {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      for(int i = 0; i < n_node[i_part]; ++i) {
        n_tot_node += select_node[i_part][i];
      }
    }
  }


  for(int i_part = 0; i_part < n_part; ++i_part) {

    int  n_arc_graph   = pn_arc_graph [i_part];
    int  n_node_graph  = pn_node_graph[i_part];
    int *_arc_graph    = parc_graph   [i_part];
    int *_node_graph   = pnode_graph  [i_part];
    int *_node_arc_idx = node_arc_idx [i_part];
    int *_node_arc     = node_arc     [i_part];
    int *_arc_node_idx = arc_node_idx [i_part];
    int *_arc_node     = arc_node     [i_part];

    /*
     * Shift computation
     */
    const int *is_owner = PDM_part_comm_graph_owner_get(pcg_node, i_part);
    is_owner_node[i_part] = PDM_array_const_int(n_node[i_part], 1);
    for(int i_graph_node = 0; i_graph_node < n_node_graph; ++i_graph_node) {
      int i_node = _node_graph[4*i_graph_node]-1;
      is_owner_node[i_part][i_node] = is_owner[i_graph_node];
    }

    /*
     * Count
     */
    PDM_malloc(send_node_n[i_part], n_node_graph, int);
    for(int i_graph_node = 0; i_graph_node < n_node_graph; ++i_graph_node) {
      int i_node = _node_graph[4*i_graph_node]-1;
      send_node_n[i_part][i_graph_node] = 0;
      if(is_owner[i_graph_node] == 0) {
        continue;
      }
      // if(pnode_ln_to_gn[i_part][i_node] == -1) { // Donc pas selectioné
      //   continue;
      // }

      for(int idx_arc = _node_arc_idx[i_node]; idx_arc < _node_arc_idx[i_node+1]; ++idx_arc) {
        int i_arc = _node_arc[idx_arc] - 1;
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
      if(is_owner[i_graph_node] == 0) {
        continue;
      }
      for(int idx_arc = _node_arc_idx[i_node]; idx_arc < _node_arc_idx[i_node+1]; ++idx_arc) {
        int i_arc = _node_arc[idx_arc] - 1;
        for(int idx_node = _arc_node_idx[i_arc]; idx_node < _arc_node_idx[i_arc+1]; ++idx_node) {
          int i_node2 = _arc_node[idx_node]-1;
          int idx_write = send_node_idx[i_graph_node] + send_node_n[i_part][i_graph_node]++;
          send_node  [i_part][idx_write] = pnode_ln_to_gn[i_part][i_node2];
          send_weight[i_part][idx_write] = arc_weight    [i_part][i_arc];
        }
      }
    }

    for(int i_graph_arc = 0; i_graph_arc < n_arc_graph; ++i_graph_arc) {
      int i_arc = _arc_graph[4*i_graph_arc]-1;
      for(int idx_node = _arc_node_idx[i_arc]; idx_node < _arc_node_idx[i_arc+1]; ++idx_node) {
        int i_node = _arc_node[idx_node]-1;
        int idx_write = send_arc_node_idx[i_graph_arc] + send_arc_node_n[i_part][i_graph_arc]++;
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
  int shift_part = 0;
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
      if(is_owner_node[i_part][i_node] == 0) {
        continue;
      }
      if(pnode_ln_to_gn[i_part][i_node] == -1) { // Donc pas selectioné
        continue;
      }
      for(int idx_arc = _node_arc_idx[i_node]; idx_arc < _node_arc_idx[i_node+1]; ++idx_arc) {
        int i_arc = _node_arc[idx_arc]-1;
        node_node_n[shift_part+i_node] += _arc_node_idx[i_arc+1] - _arc_node_idx[i_arc];
        max_size += _arc_node_idx[i_arc+1] - _arc_node_idx[i_arc];
      }
    }

    for(int i_graph_node = 0; i_graph_node < n_node_graph; ++i_graph_node) {
      int i_node = _node_graph[4*i_graph_node]-1;
      node_node_n[shift_part+i_node] += recv_node_n[i_part][i_graph_node];
      max_size += recv_node_n[i_part][i_graph_node];
    }

    for(int i_graph_arc = 0; i_graph_arc < n_arc_graph; ++i_graph_arc) {
      int i_arc = _arc_graph[4*i_graph_arc]-1;
      for(int idx_node = _arc_node_idx[i_arc]; idx_node < _arc_node_idx[i_arc+1]; ++idx_node) {
        int i_node = _arc_node[idx_node]-1;
        node_node_n[shift_part+i_node] += recv_arc_node_n[i_part][i_graph_arc];
        max_size += recv_arc_node_n[i_part][i_graph_arc];
      }
    }

    shift_part += n_node[i_part];
  }

  /*
   * Count
   */
  int *node_node_idx = NULL;
  PDM_malloc(node_node_idx, n_tot_node+1, int);

  int max_node_node = 0;
  node_node_idx[0] = 0;
  for(int i = 0; i < n_tot_node; ++i) {
    node_node_idx[i+1] = node_node_idx[i] + node_node_n[i];
    max_node_node = PDM_MAX(max_node_node, node_node_n[i]);
    node_node_n[i] = 0;
  }

  /*
   * Fill
   */
  PDM_g_num_t *gnode_node  = NULL;
  int         *garc_weight = NULL;
  PDM_malloc(gnode_node , max_size, PDM_g_num_t);
  PDM_malloc(garc_weight, max_size, int        );
  shift_part = 0;
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
      if(is_owner_node[i_part][i_node] == 0) {
        continue;
      }
      if(pnode_ln_to_gn[i_part][i_node] == -1) { // Donc pas selectioné
        continue;
      }
      for(int idx_arc = _node_arc_idx[i_node]; idx_arc < _node_arc_idx[i_node+1]; ++idx_arc) {
        int i_arc = _node_arc[idx_arc]-1;
        for(int idx_node = _arc_node_idx[i_arc]; idx_node < _arc_node_idx[i_arc+1]; ++idx_node) {
          int i_node2 = _arc_node[idx_node]-1;
          int idx_write = node_node_idx[shift_part+i_node] + node_node_n[shift_part+i_node]++;
          gnode_node [idx_write] = pnode_ln_to_gn[i_part][i_node2];
          garc_weight[idx_write] = arc_weight    [i_part][i_arc];
        }
      }
    }

    int idx_read = 0;
    for(int i_graph_node = 0; i_graph_node < n_node_graph; ++i_graph_node) {
      int i_node = _node_graph[4*i_graph_node]-1;
      for(int j = 0; j < recv_node_n[i_part][i_graph_node]; ++j) {
        int idx_write = node_node_idx[shift_part+i_node] + node_node_n[shift_part+i_node]++;
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
        int i_node = _arc_node[idx_node]-1;
        for(int j = 0; j < recv_arc_node_n[i_part][i_graph_arc]; ++j) {
          int idx_write = node_node_idx[shift_part+i_node] + node_node_n[shift_part+i_node]++;
          gnode_node [idx_write] = recv_arc_node  [i_part][idx_read+j];
          garc_weight[idx_write] = recv_arc_weight[i_part][idx_read+j];
        }
      }
      idx_read += recv_arc_node_n[i_part][i_graph_arc];
    }

    shift_part += n_node[i_part];
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

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(is_owner_node[i_part]);
  }
  PDM_free(is_owner_node);

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
      if(gnode_node[beg+i] != gnum || gnode_node[beg+i] == -1) {
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

  if(1 == 1) {
    log_trace("gnode_node ----- \n");
    for(int i = 0; i < n_tot_node; ++i) {
      log_trace("ln_to_gn = "PDM_FMT_G_NUM" \n", distrib_node[i_rank]+i+1);
      for(int j = node_node_idx[i]; j < node_node_idx[i+1]; ++j) {
        log_trace(""PDM_FMT_G_NUM" (%i) ", gnode_node[j], garc_weight[j]);
      }
      log_trace("\n");
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

  /*
   * Fix output
   */
  *out_gnode_node_idx = node_node_idx;
  *out_gnode_node     = gnode_node;
  *out_garc_weight    = garc_weight;
  *out_distrib_node   = distrib_node;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
