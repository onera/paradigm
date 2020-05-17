
/*============================================================================
 * Parallel partitioning
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

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

/**
 * \def _PDM_part_MIN(a,b)
 * Computes the minimum of \a x and \a y.
 *
 */

#define _PDM_part_MIN(a,b) ((a) > (b) ? (b) : (a))

/**
 * \def _PDM_part_MAX(a,b)
 * Computes the maximum of \a x and \a y.
 *
 */

#define _PDM_part_MAX(a,b) ((a) < (b) ? (b) : (a))

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

// void
// PDM_para_graph_dual_from_face_cell
// (
//  const PDM_MPI_Comm     comm,
//        PDM_g_num_t     *cell_distribution,
//        PDM_g_num_t     *face_distribution,
//        PDM_g_num_t     *dface_cell,
//        PDM_g_num_t    **dual_graph,
//        PDM_g_num_t    **dual_graph_idx
// )
// {

//   int dn_face = face_distribution[i_rank+1] -  face_distribution[i_rank];
//   int dn_cell = cell_distribution[i_rank+1] -  cell_distribution[i_rank];

// }

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

    int n_cell_connect_comp = PDM_inpace_unique_long(dual_graph, idx_comp, end_connect);
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
 * \brief Compute the dual graph in parallel for a face cell connectivity
 *
 */
void
PDM_para_graph_dual_from_face_cell
(
 const PDM_MPI_Comm     comm,
 const int              dn_cell,
 const int              dn_face,
       PDM_g_num_t     *dface_cell,
       PDM_g_num_t    **dual_graph,
       PDM_g_num_t    **dual_graph_idx
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  printf("dn_cell : %d\n", dn_cell);
  printf("dn_face : %d\n", dn_face);

  /*
   * We need for each cell the connectivity with other cells ( in global numbering )
   *    -> We can use part to block on face_cell to setup the correct connectivity
   *       because for each cells we receive the contribution of each face connectivity
   */
  PDM_g_num_t* dcell_ln_to_gn = (PDM_g_num_t *) malloc( 2 * dn_face * sizeof(PDM_g_num_t));
  PDM_g_num_t* dcell_opp      = (PDM_g_num_t *) malloc( 2 * dn_face * sizeof(PDM_g_num_t));

  int dn_face_int = 0;
  for(int i_face = 0; i_face < dn_face; ++i_face){
    PDM_g_num_t g_cell1 = dface_cell[2*i_face  ];
    PDM_g_num_t g_cell2 = dface_cell[2*i_face+1];

    if(g_cell2 > 0){
      dcell_ln_to_gn[dn_face_int]   = g_cell1;
      dcell_opp[dn_face_int++]      = g_cell2;

      dcell_ln_to_gn[dn_face_int]   = g_cell2;
      dcell_opp[dn_face_int++]      = g_cell1;

    }

  }

  // printf("dcell_ln_to_gn::");
  // for(int i = 0; i < dn_face_int; ++i){
  //   printf("%d ", dcell_ln_to_gn[i]);
  // }
  // printf("\n");

  dcell_ln_to_gn = realloc(dcell_ln_to_gn, dn_face_int * sizeof(PDM_g_num_t) );
  dcell_opp      = realloc(dcell_opp     , dn_face_int * sizeof(PDM_g_num_t) );

  /*
   * Initialize part_to_block for the computation of cell_cell
   */
  PDM_part_to_block_t *ptb_dual =
   PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                             PDM_PART_TO_BLOCK_POST_MERGE,
                             1.,
                             &dcell_ln_to_gn,
                              NULL,
                             &dn_face_int,
                             1,
                             comm);


  /*
   * We exchange the dcell_ln_to_gn
   *    NB : Eric, on pourrai faire un part_to_block stride cst --> Stride variable ?
   */
  int* face_strid = (int *) malloc(sizeof(int) * dn_face_int);
  for(int i_face = 0; i_face < dn_face_int; ++i_face ){
    face_strid[i_face] = 1;
  }

  int* cell_cell_n = NULL;

  PDM_part_to_block_exch (ptb_dual,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          1,
                          &face_strid,
                (void **) &dcell_opp,
                          &cell_cell_n,
                (void **) &*dual_graph);

  //
  const int n_cell_block = PDM_part_to_block_n_elt_block_get (ptb_dual);

  if( 0 == 1){
    printf("n_cell_block:: %d \n", n_cell_block);
    int idx_block = 0;
    for(int i = 0; i < n_cell_block; ++i){
      printf(" cell_cell_n = %d ---> ", cell_cell_n[i]);
      for(int i_data = 0; i_data < cell_cell_n[i]; ++i_data){
        printf("%d ", dual_graph[0][idx_block]);
        idx_block++;
      }
      printf("\n");
    }
  }

  /*
   * Exchange is done we can free direclty memory
   */
  free(dcell_ln_to_gn);
  free(face_strid);
  free(dcell_opp);

  /*
   * Allocate and setup convenient pointeur
   */
  *dual_graph_idx = (int*        ) malloc( sizeof(int        ) * (n_cell_block+1));

  int*         _dual_graph_idx = (int          *) *dual_graph_idx;
  PDM_g_num_t* _dual_graph     = (PDM_g_num_t  *) *dual_graph;

  /*
   * Each block can have multiple same cell, we need to compress them
   *   We do it inplace cause unique will always have inferior size of complete array
   *
   */
  PDM_compress_connectivity(_dual_graph, _dual_graph_idx, cell_cell_n, n_cell_block);

  /*
   * Realloc
   */
  *dual_graph     = (PDM_g_num_t*) realloc(*dual_graph, sizeof(PDM_g_num_t) * _dual_graph_idx[n_cell_block] );

  /*
   * Panic verbose
   */
  if( 0 == 1 ){
    printf("n_cell_block:: %d \n", n_cell_block);
    for(int i = 0; i < n_cell_block; ++i){
      printf(" _dual_graph_idx = %d ---> \n", _dual_graph_idx[i]);
      for(int i_data = _dual_graph_idx[i]; i_data < _dual_graph_idx[i+1]; ++i_data){
        // printf("%d ", _dual_graph[i_data]);
        printf("\t _dual_graph[%d] = %d \n", i_data, _dual_graph[i_data]);
      }
      printf("\n");
    }
  }

  PDM_part_to_block_free (ptb_dual);
  free(cell_cell_n);

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
       PDM_g_num_t    **dual_graph,
       PDM_g_num_t    **dual_graph_idx

)
{
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
