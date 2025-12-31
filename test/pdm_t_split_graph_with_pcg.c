#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pdm.h"
#include "pdm_array.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_gnum.h"
#include "pdm_mem_tool.h"
#include "pdm_mpi.h"
#include "pdm_multipart.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_algorithm.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_printf.h"
#include "pdm_priv.h"

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
  PDM_part_comm_graph_t  *pcg_arc,
  PDM_split_dual_t        split_method,
  int                   **ppart_id_idx,
  int                   **ppart_id
)
{
  /*
   * To avoid to much computation, we directly pre-alloc the node-node :
   *   - We need to know size with shared node (only owner will take graph)
   *   - Only compute node global numbering ?
   */

  PDM_gen_gnum_t *gen_gnum_node = PDM_gnum_create(3, 1, PDM_TRUE, 1e-6, comm, PDM_OWNERSHIP_USER);
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

  /*
   *
   */
  int         **send_node_n = NULL;
  PDM_g_num_t **send_node   = NULL;
  PDM_malloc(send_node_n, n_part, int         *);
  PDM_malloc(send_node  , n_part, PDM_g_num_t *);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    int *node_graph = NULL;
    int n_node_graph = PDM_part_comm_graph_entity_graph_get(pcg_node,
                                                            i_part,
                                                            &node_graph,
                                                            PDM_OWNERSHIP_BAD_VALUE);

    const int *is_owner = PDM_part_comm_graph_owner_get(pcg_node, i_part);

    int *_node_arc_idx = node_arc_idx[i_part];
    int *_node_arc     = node_arc    [i_part];
    int *_arc_node_idx = arc_node_idx[i_part];
    int *_arc_node     = arc_node    [i_part];

    /*
     * Count
     */
    PDM_malloc(send_node_n[i_part], n_node_graph, int);
    for(int i_graph_node = 0; i_graph_node < n_node_graph; ++i_graph_node) {
      int i_node = node_graph[4*i_graph_node]-1;
      send_node_n[i_part][i_graph_node] = 0;
      if(is_owner[i_graph_node] == 1) {
        continue;
      }
      for(int idx_arc = _node_arc_idx[i_node]; idx_arc < _node_arc_idx[i_node+1]; ++idx_arc) {
        int i_arc = _node_arc[idx_arc] - 1;
        send_node_n[i_part][i_graph_node] += _arc_node_idx[i_arc+1] - _arc_node_idx[i_arc];
      }
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
    PDM_malloc(send_node[i_part], send_node_idx[n_node_graph], PDM_g_num_t);

    /*
     * Fill
     */
    for(int i_graph_node = 0; i_graph_node < n_node_graph; ++i_graph_node) {
      int i_node = node_graph[4*i_graph_node]-1;
      if(is_owner[i_graph_node] == 1) {
        continue;
      }
      for(int idx_arc = _node_arc_idx[i_node]; idx_arc < _node_arc_idx[i_node+1]; ++idx_arc) {
        int i_arc = _node_arc[idx_arc] - 1;
        for(int idx_node = _arc_node_idx[i_arc]; idx_node < _arc_node_idx[i_arc+1]; ++idx_node) {
          int i_node2 = _arc_node[idx_node]-1;
          int idx_write = send_node_idx[i_graph_node] + send_node_n[i_part][i_graph_node]++;
          send_node[i_part][idx_write] = pnode_ln_to_gn[i_part][i_node2];
        }
      }
    }


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
                           &send_node_n,
                 (void **) &send_node,
                           &recv_node_n,
                (void ***) &recv_node);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(send_node  [i_part]);
    PDM_free(send_node_n[i_part]);
  }
  PDM_free(send_node_n);
  PDM_free(send_node  );

  /*
   * Node are by construction ordered by increasing gnum numbering
   */




  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(recv_node  [i_part]);
    PDM_free(recv_node_n[i_part]);
  }
  PDM_free(recv_node_n);
  PDM_free(recv_node  );

  // for(int i_part = 0; i_part < n_part; ++i_part) {

  //   for(int i_arc = 0; i_arc < n_arc[i_part]; ++i_arc) {
  //     int i_node1 = arc_node[2*i_arc  ]-1;
  //     int i_node2 = arc_node[2*i_arc+1]-1;

  //     // We need to have symetric weight on graph
  //     int pos11 = PDM_binary_search_int(i_node1+1, &node_node[node_node_idx[i_node1]], n_adj_node1);
  //     int pos12 = PDM_binary_search_int(i_node1+1, &node_node[node_node_idx[i_node2]], n_adj_node2);
  //     int pos21 = PDM_binary_search_int(i_node2+1, &node_node[node_node_idx[i_node1]], n_adj_node1);
  //     int pos22 = PDM_binary_search_int(i_node2+1, &node_node[node_node_idx[i_node2]], n_adj_node2);

  //     node_node_weight[beg1+pos1] += arc_weight[i_arc];
  //     node_node_weight[beg2+pos2] += arc_weight[i_arc];

  //   }

  // }


  /*
   * Synchronise node weight
   */

  /*
   * Synchronise arc weight
   */


  // int *part_id = NULL;
  // PDM_malloc(part_id, n_vtx_owner, int);
  // double *part_fraction = NULL;
  // PDM_para_graph_split(split_method,
  //                      vtx_distrib,
  //                      gvtx_vtx_idx,
  //                      gvtx_vtx,
  //                      vtx_weight,
  //                      gvtx_vtx_weight,
  //                      n_rank, // Number of partition
  //                      part_fraction,
  //                      part_id,
  //                      mawr->comm);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pnode_ln_to_gn[i_part]);
  }
  PDM_free(pnode_ln_to_gn);
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


  PDM_part_mesh_nodal_free(pmn);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }
  PDM_MPI_Finalize();

  return EXIT_SUCCESS;
}
