#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pdm.h"
#include "pdm_array.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_distrib.h"
#include "pdm_gnum.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include "pdm_mpi.h"
#include "pdm_multipart.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_algorithm.h"
#include "pdm_part_mesh_nodal_to_part_mesh.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_printf.h"
#include "pdm_priv.h"
#include "pdm_sort.h"
#include "pdm_unique.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/


/**
 *
 * \brief  Usage
 *
 */

static void
_usage
(
int exit_code
)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number vtx in side of mesh A (default : 10).\n\n"
     "  -n_part <level>  Number of partition                        .\n\n"
     "  -post            Ensight outputs (only if n_part == 1). \n\n"
     "  -t               Element kind .\n\n"
     "  -h               This message.\n\n");
  exit (exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]    argc   Number of arguments
 * \param [in]    argv   Arguments
 * \param [inout] nPts   Number of points
 * \param [inout] ls     Low scalability
 * \param [inout] length Length of domains
 *
 */

static void
_read_args
(
 int                    argc,
 char                 **argv,
 PDM_g_num_t           *n_vtx_a,
 int                   *post
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _n_vtx_a = atol(argv[i]);
        *n_vtx_a = (PDM_g_num_t) _n_vtx_a;
      }
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}

static
PDM_part_mesh_nodal_t*
_generate_mesh
(
  PDM_MPI_Comm    pdm_comm,
  int             n_vtx_seg
)
{
  int              n_part       = 1;
  PDM_split_dual_t split_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;

  /* Warmup */
  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create(pdm_comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        2.,
                                                        -1,
                                                        -1,
                                                        -1,
                                                        PDM_MESH_NODAL_TRIA3,
                                                        1,
                                                        PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);
  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

  int n_domain = 1;
  int n_part_domains = n_part;
  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                &n_part_domains,
                                                PDM_FALSE,
                                                split_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                pdm_comm,
                                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart,
                                       -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");

  PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);

  if(0 == 1) {
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  }
  PDM_multipart_compute(mpart);

  PDM_part_mesh_nodal_t* pmesh_nodal = NULL;
  PDM_multipart_get_part_mesh_nodal(mpart, 0, &pmesh_nodal, PDM_OWNERSHIP_USER);

  PDM_part_mesh_nodal_complete_part_comm_graph(pmesh_nodal);

  PDM_multipart_free(mpart);
  PDM_dcube_nodal_gen_free(dcube);
  return pmesh_nodal;
}

static
void
_part_split
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
  PDM_part_comm_graph_t  *pcg_arc
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

  if(select_node != NULL) {
    int *pn_select_node = NULL;
    PDM_malloc(pn_select_node, n_part, int  );
    for(int i_part = 0; i_part < n_part; ++i_part) {
      pn_select_node[i_part] = 0;
      for(int i = 0; i < n_node[i_part]; ++i) {
        if(select_node[i_part][i] == 1) {
          pn_select_node[i_part]++;
        }
      }
    }
    PDM_free(pn_select_node);
  }

  PDM_part_comm_graph_t *pcg_subnode = NULL;
  if(select_node != NULL) {
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

  PDM_g_num_t l_shift    = 0;
  int         n_tot_node = 0;
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
    l_shift    += n_node[i_part];
    n_tot_node += n_node[i_part];
    for(int i_graph_node = 0; i_graph_node < n_node_graph; ++i_graph_node) {
      if(is_owner[i_graph_node] == 0) {
        l_shift--;
      }
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
      if(is_owner[i_graph_node] == 1) {
        continue;
      }
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
      if(is_owner[i_graph_node] == 1) {
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
   * Compute global shift
   */
  // PDM_g_num_t g_shift = 0;
  // PDM_MPI_Exscan(&l_shift, &g_shift, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);

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
      if(gnode_node[beg+i] != gnum) {
        gnode_node [idx_write] = gnode_node[beg+i];
        garc_weight[idx_write] = lweight[lorder[i]];
        idx_write++;
      }
    }

    idx_read = node_node_idx[i_node+1];

    node_node_idx[i_node+1] = node_node_idx[i_node] + n_unique - 1;

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


  // A renvoyer
  PDM_free(node_node_idx);
  PDM_free(gnode_node);
  PDM_free(garc_weight);
  PDM_free(distrib_node );

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
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief  Main
 *
 */

int
main
(
  int   argc,
  char *argv[]
)
{
  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_g_num_t n_vtx_a = 10;

  int post = 0;

  _read_args(argc,
             argv,
             &n_vtx_a,
             &post);

  PDM_part_mesh_nodal_t* pmn = _generate_mesh(comm, n_vtx_a);

  if(post) {
    PDM_part_mesh_nodal_dump_vtk(pmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  }

  /*
   * Differents modes de splitting
   *   - Split les sommets (idéal pour l'adaptation de maillage )
   *   - Split les cellules (le plus classiques )
   *   - On peut partir d'un pmn ou d'un pm
   */

  /* Warm-up */
  int n_part = PDM_part_mesh_nodal_n_part_get(pmn);
  int  *pn_elt       = NULL;
  int  *pn_vtx       = NULL;
  int **pelt_vtx_idx = NULL;
  int **pelt_vtx     = NULL;
  int **pvtx_vtx_idx = NULL;
  int **pvtx_vtx     = NULL;
  PDM_malloc(pn_elt      , n_part, int  );
  PDM_malloc(pn_vtx      , n_part, int  );
  PDM_malloc(pelt_vtx_idx, n_part, int *);
  PDM_malloc(pelt_vtx    , n_part, int *);
  PDM_malloc(pvtx_vtx_idx, n_part, int *);
  PDM_malloc(pvtx_vtx    , n_part, int *);
  for(int i_part = 0; i_part < n_part; ++i_part) {

    pn_vtx[i_part] = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
    PDM_geometry_kind_t principal_kind = PDM_part_mesh_nodal_principal_geom_kind_get(pmn);
    pn_elt[i_part] = PDM_part_mesh_nodal_cell_vtx_connect_get(pmn,
                                                              principal_kind,
                                                              i_part,
                                                              &pelt_vtx_idx[i_part],
                                                              &pelt_vtx    [i_part]);

    // Vtx-vtx connectivity
    int *pvtx_elt_idx = NULL;
    int *pvtx_elt     = NULL;
    PDM_connectivity_transpose(pn_elt       [i_part],
                               pn_vtx       [i_part],
                               pelt_vtx_idx [i_part],
                               pelt_vtx     [i_part],
                               &pvtx_elt_idx,
                               &pvtx_elt);

    PDM_combine_connectivity(pn_vtx       [i_part],
                             pvtx_elt_idx,
                             pvtx_elt,
                             pelt_vtx_idx [i_part],
                             pelt_vtx     [i_part],
                             &pvtx_vtx_idx[i_part],
                             &pvtx_vtx    [i_part]);

    PDM_free(pvtx_elt_idx);
    PDM_free(pvtx_elt    );
  }

  /*
   * Graph assembly
   */
  PDM_part_mesh_nodal_to_part_mesh_t* pmn_to_pm = PDM_part_mesh_nodal_to_part_mesh_create(pmn,
                                                                                          PDM_FALSE,
                                                                                          PDM_OWNERSHIP_USER);

  // Connectivities
  int dim = 2;
  PDM_part_mesh_t *pm = NULL;
  PDM_part_mesh_nodal_to_part_mesh_g_nums_enable(pmn_to_pm, PDM_MESH_ENTITY_VTX );
  PDM_part_mesh_nodal_to_part_mesh_g_nums_enable(pmn_to_pm, PDM_MESH_ENTITY_EDGE);
  if(dim == 3) {
    PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm,
                                                         PDM_CONNECTIVITY_TYPE_CELL_FACE);
    PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm,
                                                         PDM_CONNECTIVITY_TYPE_FACE_VTX);
  } else {
    PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm,
                                                         PDM_CONNECTIVITY_TYPE_FACE_EDGE);
    PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm,
                                                         PDM_CONNECTIVITY_TYPE_EDGE_VTX);
  }

  PDM_part_mesh_nodal_to_part_mesh_compute(pmn_to_pm);

  PDM_part_mesh_nodal_to_part_mesh_part_mesh_get(pmn_to_pm,
                                                 &pm,
                                                 PDM_OWNERSHIP_USER);

  PDM_part_mesh_nodal_to_part_mesh_free(pmn_to_pm);

  // Hook
  PDM_part_mesh_part_comm_graph_compute_from_gnum(pm, PDM_MESH_ENTITY_VTX );
  PDM_part_mesh_part_comm_graph_compute_from_gnum(pm, PDM_MESH_ENTITY_EDGE);


  int  *pn_node       = NULL;
  int  *pn_arc        = NULL;
  int **pselect_node  = NULL;
  int **pnode_arc_idx = NULL;
  int **pnode_arc     = NULL;
  int **parc_node_idx = NULL;
  int **parc_node     = NULL;
  int **pnode_weight  = NULL;
  int **parc_weight   = NULL;

  PDM_malloc(pn_node      , n_part, int  );
  PDM_malloc(pn_arc       , n_part, int  );
  // PDM_malloc(pselect_node , n_part, int *);
  // PDM_malloc(pnode_arc_idx, n_part, int *);
  // PDM_malloc(pnode_arc    , n_part, int *);
  PDM_malloc(parc_node_idx, n_part, int *);
  PDM_malloc(parc_node    , n_part, int *);
  PDM_malloc(pnode_weight , n_part, int *);
  PDM_malloc(parc_weight  , n_part, int *);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    pn_node[i_part] = PDM_part_mesh_n_entity_get(pm, i_part, PDM_MESH_ENTITY_VTX );
    pn_arc [i_part] = PDM_part_mesh_n_entity_get(pm, i_part, PDM_MESH_ENTITY_EDGE);

    PDM_part_mesh_connectivity_get(pm,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                   &parc_node    [i_part],
                                   &parc_node_idx[i_part],
                                   PDM_OWNERSHIP_KEEP);

    double *vtx_coords = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part, PDM_OWNERSHIP_BAD_VALUE);

    if(parc_node_idx[i_part] == NULL) {
      PDM_malloc(parc_node_idx[i_part], pn_arc[i_part] + 1, int);
      for(int i = 0; i < pn_arc[i_part]+1; ++i) {
        parc_node_idx[i_part][i] = 2*i;
      }
    }

    PDM_malloc(pnode_weight[i_part], pn_node[i_part], int);
    PDM_malloc(parc_weight [i_part], pn_arc [i_part], int);

    for(int i = 0; i < pn_node[i_part]; ++i) {
      pnode_weight[i_part][i] = 1;
    }
    for(int i = 0; i < pn_arc[i_part]; ++i) {
      parc_weight[i_part][i] = 1;
    }

    // This is stupid but why not
    double vdir[3] = {1., 0., 0.};
    for(int i = 0; i < pn_arc[i_part]; ++i) {
      int i_vtx1 = parc_node[i_part][2*i  ]-1;
      int i_vtx2 = parc_node[i_part][2*i+1]-1;

      double v[3] = {vtx_coords[3*i_vtx2  ] - vtx_coords[3*i_vtx1  ],
                     vtx_coords[3*i_vtx2+1] - vtx_coords[3*i_vtx1+1],
                     vtx_coords[3*i_vtx2+2] - vtx_coords[3*i_vtx1+2]};
      double mod = PDM_MODULE(v);
      v[0] = v[0] / mod;
      v[1] = v[1] / mod;
      v[2] = v[2] / mod;

      double vdot = PDM_DOT_PRODUCT(vdir, v);
      int idot = (int) ( PDM_ABS(vdot) * 10) ;

      log_trace("i_arc = %i / i_vtx1 = %i / i_vtx2 = %i --> idot = %i (%12.5e) \n ", i, i_vtx1, i_vtx2, idot, vdot);

      parc_weight[i_part][i] += idot;
    }

    PDM_log_trace_array_int(parc_weight[i_part], pn_arc[i_part], "parc_weight ::");
  }

  // Transpose
  PDM_part_connectivity_transpose(n_part,
                                  pn_arc,
                                  pn_node,
                                  parc_node_idx,
                                  parc_node,
                                  &pnode_arc_idx,
                                  &pnode_arc);


  PDM_part_comm_graph_t *pcg_node = NULL;
  PDM_part_mesh_part_comm_graph_get(pm,
                                    PDM_MESH_ENTITY_VTX,
                                    &pcg_node,
                                    PDM_OWNERSHIP_KEEP);

  PDM_part_comm_graph_t *pcg_arc = NULL;
  PDM_part_mesh_part_comm_graph_get(pm,
                                    PDM_MESH_ENTITY_EDGE,
                                    &pcg_arc,
                                    PDM_OWNERSHIP_KEEP);

  /*
   * Tester avec cell_vtx + vtx_cell aussi -> Shortcut for mesh adaptation + quality
   */
  _part_split(comm,
              n_part,
              pn_node,
              pn_arc,
              pselect_node,
              pnode_arc_idx,
              pnode_arc,
              parc_node_idx,
              parc_node,
              pnode_weight,
              parc_weight,
              pcg_node,
              pcg_arc);


  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pnode_arc_idx[i_part]);
    PDM_free(pnode_arc    [i_part]);
    PDM_free(pnode_weight [i_part]);
    PDM_free(parc_weight  [i_part]);
  }

  /* Free all */
  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pelt_vtx_idx[i_part]);
    PDM_free(pelt_vtx    [i_part]);
    PDM_free(pvtx_vtx_idx[i_part]);
    PDM_free(pvtx_vtx    [i_part]);
  }
  PDM_free(pelt_vtx_idx);
  PDM_free(pelt_vtx    );
  PDM_free(pvtx_vtx_idx);
  PDM_free(pvtx_vtx    );
  PDM_free(pn_elt);
  PDM_free(pn_vtx);

  PDM_free(pn_node      );
  PDM_free(pn_arc       );
  // PDM_free(pselect_node );
  PDM_free(pnode_arc_idx);
  PDM_free(pnode_arc    );
  PDM_free(parc_node_idx);
  PDM_free(parc_node    );
  PDM_free(pnode_weight );
  PDM_free(parc_weight  );


  PDM_part_mesh_free(pm);
  PDM_part_mesh_nodal_free(pmn);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }
  PDM_MPI_Finalize();

  return EXIT_SUCCESS;
}
