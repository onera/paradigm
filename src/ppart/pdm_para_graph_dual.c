
/*============================================================================
 * Parallel partitioning
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_timer.h"
#include "pdm_mpi.h"
#include "pdm_mpi_ext_dependencies.h"

#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_sort.h"
#include "pdm_binary_search.h"
#include "pdm_handles.h"
#include "pdm_part_to_block.h"

#include "pdm_para_graph_dual.h"


/*----------------------------------------------------------------------------
 *  Optional headers
 *----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

void
PDM_compress_connectivity
(
 PDM_g_num_t *dual_graph,
 int         *dual_graph_idx,
 int         *dual_graph_n,
 int          dn_elt
)
{
  int idx_comp  = 0; /* Compressed index use to fill the buffer */
  int idx_block = 0; /* Index in the block to post-treat        */
  dual_graph_idx[0] = 0;
  int need_shift = 0;
  for(int i = 0; i < dn_elt; ++i){
    int n_cell_connect = dual_graph_n[i];
    /* Reshift next value in compressed block to avoid create a new shift */
    if(need_shift) {
      int idx_new = idx_comp;
      for(int j = idx_block; j < idx_block+n_cell_connect; ++j){
        dual_graph[idx_new++] = dual_graph[j];
      }
    }

    int end_connect         = idx_comp + n_cell_connect - 1;
    // printf(" idx_comp:: %d | end_connect:: %d \n", idx_comp, end_connect);

    int n_cell_connect_comp = PDM_inplace_unique_long(dual_graph, idx_comp, end_connect);
    // printf(" n_cell_connect:: %d | n_cell_connect_comp:: %d \n", n_cell_connect, n_cell_connect_comp);

    if(n_cell_connect_comp < n_cell_connect) {
      need_shift = 1;
    } else {
      need_shift = 0;
    }

    dual_graph_idx[i+1] = dual_graph_idx[i] + n_cell_connect_comp;
    idx_comp  += n_cell_connect_comp;
    idx_block += n_cell_connect;

    // printf("idx_comp : %d | idx_block : %d \n", idx_comp, idx_block);
  }
}


/**
 *
 * \brief Compute in parallel the dual graph of an unstructured graph represented
 *        by its arc (edges of the graph) to node (vertices of the graph) connectivity.
 *        Arc and edge terminology is employed to avoid confusion with geometric entities
 *        such as vertices, edges, etc.
 *        Usually for a CFD mesh, the nodes of the graph are the cells of the mesh
 *        and the arcs of the graph are thus the faces of the mesh.
 *
 *        The dual graph computed by this function is the node to node connectivity
 *        as know as adjacency list, requested by graph partitioners.
 *
 *        Additively, this function computes the node to arc connectivity if
 *        compute_node_to_arc is true.
 *
 * \param [in]   comm               PDM_MPI communicator
 * \param [in]   graph_node_distrib distribution of nodes over the procs (size=n_rank+1)
 * \param [in]   graph_node_distrib distribution of arcs  over the procs (size=n_rank+1)
 * \param [in]   darc_to_node       Arc to node connectivity (size=2*dn_arc)
 * \param [out]  dual_graph_idx     Node to node connectivity indexes (size=dn_node+1)
 * \param [out]  dual_graph         Node to node connectivity (size=dual_graph_idx[dn_node])
 * \param [in]   compute_dnode_to_arc Compute or not node to arc connectivity
 * \param [out]  dnode_to_arc_idx   Node to arc connectivity indexes (size=dn_node+1)
 * \param [out]  dnode_to_arc       Node to arc connectivity (size=dnode_to_arc_idx[dn_node])
 *
 */
void
PDM_para_graph_dual_from_arc2node
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *graph_node_distrib,
 const PDM_g_num_t     *graph_arc_distrib,
 const PDM_g_num_t     *darc_to_node,
       int            **dual_graph_idx,
       PDM_g_num_t    **dual_graph,
 const int              compute_dnode_to_arc,
       int            **dnode_to_arc_idx,
       PDM_g_num_t    **dnode_to_arc
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_arc  = graph_arc_distrib[i_rank+1]  -  graph_arc_distrib[i_rank];
  int dn_node = graph_node_distrib[i_rank+1] -  graph_node_distrib[i_rank];
  // printf("dn_node : %d\n", dn_node);
  // printf("dn_arc  : %d\n", dn_arc);

  /*
   * We need for each node the connectivity with other nodes ( in global numbering )
   *    -> We can use part to block on arc2node to setup the correct connectivity
   *       because for each nodes we receive the contribution of each arc connectivity
   */
  PDM_g_num_t* dnode_ln_to_gn = (PDM_g_num_t *) malloc( 2 * dn_arc * sizeof(PDM_g_num_t));
  PDM_g_num_t* dopposite_node = (PDM_g_num_t *) malloc( 2 * dn_arc * sizeof(PDM_g_num_t));
  PDM_g_num_t* darc_g;

  int* arc_strid  = (int *) malloc(sizeof(int) * 2 * dn_arc);
  int* node_strid = (int *) malloc(sizeof(int) * 2 * dn_arc);


  if(compute_dnode_to_arc){
    darc_g = (PDM_g_num_t *) malloc( 2 * dn_arc * sizeof(PDM_g_num_t));
  }

  int shift_arc_g   = graph_arc_distrib[i_rank]; // Entre 1 et N
  int dn_arc_int    = 0;
  int idx_data_arc  = 0;
  int idx_data_node = 0;
  for(int i_arc = 0; i_arc < dn_arc; ++i_arc){

    PDM_g_num_t g_node1 = darc_to_node[2*i_arc  ];
    PDM_g_num_t g_node2 = darc_to_node[2*i_arc+1];

    dnode_ln_to_gn[dn_arc_int] = g_node1;

    if(compute_dnode_to_arc){
      arc_strid[dn_arc_int]  = 1;
      darc_g[idx_data_arc++] = shift_arc_g + i_arc ;
    }

    if(g_node2 > 0){
      node_strid[dn_arc_int++]        = 1;
      dopposite_node[idx_data_node++] = g_node2;

      if(compute_dnode_to_arc){
        arc_strid[dn_arc_int]  = 1;
        darc_g[idx_data_arc++] = shift_arc_g + i_arc;
      }

      node_strid[dn_arc_int]          = 1;
      dnode_ln_to_gn[dn_arc_int++]    = g_node2;
      dopposite_node[idx_data_node++] = g_node1;

    } else {
      node_strid[dn_arc_int++] = 0;
    }

  }

  if( 0 == 1){
    printf("dnode_ln_to_gn::");
    for(int i = 0; i < dn_arc_int; ++i){
      printf("%d ", dnode_ln_to_gn[i]);
    }
    printf("\n");

    printf("node_strid::");
    for(int i = 0; i < dn_arc_int; ++i){
      printf("%d ", node_strid[i]);
    }
    printf("\n");

    printf("dopposite_node::");
    for(int i = 0; i < idx_data_node ; ++i){
      printf("%d ", dopposite_node[i]);
    }
    printf("\n");

    printf("idx_data_node ::%d\n", idx_data_node );
    printf("dn_arc_int    ::%d\n", dn_arc_int);
    printf("dn_arc        ::%d\n", dn_arc);
  }


  dnode_ln_to_gn = realloc(dnode_ln_to_gn, dn_arc_int    * sizeof(PDM_g_num_t) );
  dopposite_node = realloc(dopposite_node, idx_data_node * sizeof(PDM_g_num_t) );

  PDM_g_num_t* graph_node_distrib_ptb = (PDM_g_num_t * ) malloc( sizeof(PDM_g_num_t) * (n_rank+1) );
  for(int i = 0; i < n_rank+1; ++i){
    graph_node_distrib_ptb[i] = graph_node_distrib[i] - 1;
  }

  /*
   * Initialize part_to_block for the computation of node_node
   *    -> Si on force une distribution utilisateur on devra passer le cell_distribution
   *           --> Semble necessaire pour parMetis mais pas scotch
   */
  PDM_part_to_block_t *ptb_dual =
   PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                             PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM,
                             1.,
                             &dnode_ln_to_gn,
                             graph_node_distrib_ptb,
                             &dn_arc_int,
                             1,
                             comm);

  const int n_node_block = PDM_part_to_block_n_elt_block_get (ptb_dual);
  const PDM_g_num_t* blk_gnum = PDM_part_to_block_block_gnum_get(ptb_dual);

  /*
   * We exchange the dnode_ln_to_gn (almost asynchrone - send strid is synchrone )
   */
  int req_id_dual = PDM_part_to_block_async_exch(ptb_dual,
                                                 sizeof(PDM_g_num_t),
                                                 PDM_STRIDE_VAR,
                                                 1,
                                                 &node_strid,
                                       (void **) &dopposite_node);



  int req_id_node_arc = -1;
  if(compute_dnode_to_arc){

    req_id_node_arc = PDM_part_to_block_async_exch(ptb_dual,
                                                   sizeof(PDM_g_num_t),
                                                   PDM_STRIDE_VAR,
                                                   1,
                                                   &arc_strid,
                                         (void **) &darc_g);
  }

  /*
   * Reception asynchrone
   */
  // int*            node_node_n = NULL;
  int*             recv_strid = NULL;
  PDM_g_num_t* recv_node_node = NULL;
  PDM_part_to_block_async_wait(ptb_dual, req_id_dual);
  int n_data_recv = PDM_part_to_block_asyn_get_raw(ptb_dual,
                                                   req_id_dual,
                                                  &recv_strid,
                                        (void **) &recv_node_node);

  /*
   * The data is recv in raw format - We need to post-treat them as int
   */
  *dual_graph = (PDM_g_num_t *) malloc( n_data_recv * sizeof(PDM_g_num_t));
  PDM_g_num_t* _dual_graph     = (PDM_g_num_t  *) *dual_graph;

  /*
   * Allocate and setup convenient pointeur
   */
  *dual_graph_idx      = (int* ) malloc( sizeof(int) * (n_node_block+1));
  int* _dual_graph_idx = *dual_graph_idx;

  int* node_node_n   = (int *) malloc( (n_node_block+1) * sizeof(int)); /* Suralloc */

  /*
   * Count - In our case we know that recv_strid == 1 or 0
   */
  for(int i_recv = 0; i_recv < n_node_block+1; ++i_recv) {
    _dual_graph_idx[i_recv] = 0;
    node_node_n    [i_recv] = 0;
  }

  /* Stride can be 0 or 1 */
  for(int i_recv = 0; i_recv < n_data_recv; ++i_recv) {
    // if(recv_strid[i_recv] != 0){
      int ielmt = blk_gnum[i_recv] - graph_node_distrib[i_rank];
      // printf("ielmt::[%d] --> [%d]\n",blk_gnum[i_recv], ielmt );
      _dual_graph_idx[ielmt+1] += recv_strid[i_recv];
    // }
  }

  /* Index computation */
  for(int i_recv = 1; i_recv < n_node_block; ++i_recv) {
    _dual_graph_idx[i_recv+1] += _dual_graph_idx[i_recv];
  }

  /* Panic verbose */
  if( 0 == 1 ){
    printf("node_node_idx::");
    for(int i_recv = 0; i_recv < n_node_block+1; ++i_recv) {
      printf(" %d", _dual_graph_idx[i_recv]);
    }
    printf("\n");
  }

  /*
   * Fill buffer
   */
  int idx_recv = 0;
  for(int i_recv = 0; i_recv < n_data_recv; ++i_recv) {
    if(recv_strid[i_recv] != 0){
      int ielmt = blk_gnum[i_recv] - graph_node_distrib[i_rank];
      _dual_graph[_dual_graph_idx[ielmt] + node_node_n[ielmt]++] = recv_node_node[idx_recv++];
    }
  }

  free(recv_strid);
  free(recv_node_node);

  /*
   * Each block can have multiple same cell, we need to compress them
   *   We do it inplace cause unique will always have inferior size of complete array
   *
   */
  PDM_compress_connectivity(_dual_graph, _dual_graph_idx, node_node_n, n_node_block);

  free(node_node_n);

  /*
   * Realloc
   */
  *dual_graph     = (PDM_g_num_t*) realloc(*dual_graph, sizeof(PDM_g_num_t) * _dual_graph_idx[n_node_block] );

  // For now we can change it later
  printf("n_node_block::%d \n", n_node_block);
  printf("dn_node      ::%d \n", dn_node );
  assert(n_node_block == dn_node );

  /*
   * Panic verbose
   */
  if( 0 == 1 ){
    printf("n_node_block:: %d \n", n_node_block);
    for(int i = 0; i < n_node_block; ++i){
      printf(" _dual_graph_idx = %d ---> \n", _dual_graph_idx[i]);
      for(int i_data = _dual_graph_idx[i]; i_data < _dual_graph_idx[i+1]; ++i_data){
        // printf("%d ", _dual_graph[i_data]);
        printf("\t _dual_graph[%d] = %d \n", i_data, _dual_graph[i_data]);
      }
      printf("\n");
    }
  }

  /* Scoth begin at 0 even if we put base value to 1 */
  for(int i = 0; i < n_node_block; ++i){
    for(int i_data = _dual_graph_idx[i]; i_data < _dual_graph_idx[i+1]; ++i_data){
      _dual_graph[i_data] = _dual_graph[i_data] - 1;
    }
  }

  /*
   * Async recv for cell_face
   */
  if(compute_dnode_to_arc){

    int*   recv_node2arc_strid = NULL;
    PDM_g_num_t* recv_node2arc = NULL;
    PDM_part_to_block_async_wait(ptb_dual, req_id_node_arc);
    int n_data_cf_recv = PDM_part_to_block_asyn_get_raw(ptb_dual,
                                                        req_id_node_arc,
                                                       &recv_node2arc_strid,
                                             (void **) &recv_node2arc);
    free(recv_node2arc_strid); // Always 1

    /*
     * Post treatment
     */
    *dnode_to_arc      = (PDM_g_num_t *) malloc(  n_data_cf_recv  * sizeof(PDM_g_num_t));
    *dnode_to_arc_idx  = (int*         ) malloc( (n_node_block+1) * sizeof(int        ));
    int* node_to_arc_n = (int*         ) malloc( (n_node_block+1) * sizeof(int        ));

    /* Short-cut */
    int*         _dnode_to_arc_idx = *dnode_to_arc_idx;
    PDM_g_num_t* _dnode_to_arc     = *dnode_to_arc;


    for(int i = 0; i < n_node_block+1; ++i) {
      _dnode_to_arc_idx[i] = 0;
    }

    for(int i_recv = 0; i_recv < n_data_cf_recv; ++i_recv) {
      int ielmt = blk_gnum[i_recv] - graph_node_distrib[i_rank];
      // printf("ielmt::[%d] --> [%d] - [%d] \n",blk_gnum[i_recv], ielmt, recv_node2arc_strid[i_recv]);
      _dnode_to_arc_idx[ielmt+1] += 1;
    }

    /* Index computation */
    for(int i = 1; i < n_node_block; ++i) {
      _dnode_to_arc_idx[i+1] += _dnode_to_arc_idx[i];
    }

    for(int i = 0; i < n_node_block+1; ++i){
      node_to_arc_n[i] = 0;
    }

    /*
     * Fill buffer -  Cas particulier ou recv_stride == 1
     */
    for(int i_recv = 0; i_recv < n_data_cf_recv; ++i_recv) {
      int ielmt = blk_gnum[i_recv] - graph_node_distrib[i_rank];
      // printf(" ielmt :: %d | i_recv :: %d | _dnode_to_arc_idx :: %d | node_to_arc_n :: %d \n", ielmt, i_recv, _dnode_to_arc_idx[ielmt], node_to_arc_n[ielmt]);
      _dnode_to_arc[_dnode_to_arc_idx[ielmt] + node_to_arc_n[ielmt]++] = recv_node2arc[i_recv];
    }

    free(recv_node2arc);

    if( 0 == 1 ){
      printf("n_node_block:: %d \n", n_node_block);
      int idx_block = 0;
      for(int i = 0; i < n_node_block; ++i){
        printf(" node_to_arc_n = %d ---> ", node_to_arc_n[i]);
        for(int i_data = 0; i_data < node_to_arc_n[i]; ++i_data){
          printf("%d ", _dnode_to_arc[idx_block]);
          idx_block++;
        }
        printf("\n");
      }
    }

    // _dnode_to_arc_idx[0] = 0;
    // for(int i_node = 0; i_node < n_node_block; ++i_node){
    //   _dnode_to_arc_idx[i_node+1] = _dnode_to_arc_idx[i_node] + node_to_arc_n[i_node];
    // }
    free(node_to_arc_n);

  }

  // abort();

  /*
   * Exchange is done we can free direclty memory
   */
  free(dnode_ln_to_gn);
  free(arc_strid);
  free(node_strid);
  free(dopposite_node);
  free(graph_node_distrib_ptb);
  if(compute_dnode_to_arc){
    free(darc_g);
  }

  PDM_part_to_block_free (ptb_dual);
}

/**
 *
 * \brief Compute the dual graph in parallel for a cell face connectivity
 */
void
PDM_para_graph_dual_from_cell_face
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *cell_distribution,
 const PDM_g_num_t     *face_distribution,
 const PDM_g_num_t     *dcell_face,
       int            **dual_graph_idx,
       PDM_g_num_t    **dual_graph

)
{
}


/**
 * \brief   Split graph
 */
void
PDM_split_graph
(
 const PDM_MPI_Comm  comm,
 int                *dual_graph_idx,
 PDM_g_num_t        *dual_graph,
 int                *delmt_weight,
 int                *cell_part,
 int                 dn_elmt,
 int                 n_part
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  for (int i = 0; i < dn_elmt; i++) {
    cell_part[i] = 0;
  }

  /*
   *  Uniquement scotch pour l'instant car l'interface est simple
   */
  int check = 1;
  int *edgeWeight = NULL;

  printf("PDM_SCOTCH_dpart %d | %d \n", n_part, dn_elmt);
  PDM_SCOTCH_dpart (dn_elmt,
                    dual_graph_idx,
                    dual_graph,
                    delmt_weight,
                    edgeWeight,
                    check,
                    comm,
                    n_part,
                    cell_part);
  printf("PDM_SCOTCH_dpart end \n");
}





#ifdef __cplusplus
}
#endif /* __cplusplus */
