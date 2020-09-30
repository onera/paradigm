
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
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_block_to_part.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_error.h"
#include "pdm_timer.h"
#include "pdm_logging.h"

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
 * Static global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/


/*=============================================================================
 * Public function definitions
 *============================================================================*/

void
PDM_deduce_combine_connectivity
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *cell_distrib,
 const PDM_g_num_t     *face_distrib,
 const int             *dcell_face_idx,
 const PDM_g_num_t     *dcell_face,
 const int             *dface_vtx_idx,
 const PDM_g_num_t     *dface_vtx,
       int            **dcell_vtx_idx,
       PDM_g_num_t    **dcell_vtx
)
{

  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_cell = cell_distrib[i_rank+1] - cell_distrib[i_rank];
  int dn_face = face_distrib[i_rank+1] - face_distrib[i_rank];

  int* dface_vtx_n = (int * ) malloc( dn_face * sizeof(int));
  for(int i = 0; i < dn_face; ++i) {
    dface_vtx_n[i] = dface_vtx_idx[i+1] - dface_vtx_idx[i];
  }

  PDM_g_num_t* entity_distribution_ptb = (PDM_g_num_t * ) malloc( sizeof(PDM_g_num_t) * (n_rank+1) );
  for(int i = 0; i < n_rank+1; ++i){
    entity_distribution_ptb[i] = face_distrib[i] - 1;
  }

  PDM_g_num_t* dcell_face_unsigned = (PDM_g_num_t* ) malloc( dcell_face_idx[dn_cell] * sizeof(PDM_g_num_t));
  for(int i = 0; i < dcell_face_idx[dn_cell]; ++i) {
    dcell_face_unsigned[i] = PDM_ABS(dcell_face[i]);
  }

  PDM_log_trace_array_int(dcell_face_idx, dn_cell+1, "dcell_face_idx::");
  PDM_log_trace_array_int(dcell_face_unsigned, dcell_face_idx[dn_cell], "dcell_face::");

  /*
   * First compute the dcell_vtx connectivity
   *      -  We use ln_to_gn = dcell_face in order to have for each partition the faces
   *      -  We exhange the face_vtx, face_vtx_idx
   * So, for each face describe by dcell_face, we receive all vtx
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(entity_distribution_ptb,
                               (const PDM_g_num_t **) &dcell_face_unsigned,
                                                      &dcell_face_idx[dn_cell],
                                                      1,
                                                      comm);

  /*
   * Exchange
   */
  int**         pface_vtx_n;
  PDM_g_num_t** pface_vtx;
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          dface_vtx_n,
             (void *  )   dface_vtx,
             (int  ***)  &pface_vtx_n,
             (void ***)  &pface_vtx);
  free(dface_vtx_n);


  int* pface_vtx_idx = (int*) malloc( (dcell_face_idx[dn_cell] + 1) * sizeof(int));
  pface_vtx_idx[0] = 0;
  for(int i = 0; i < dcell_face_idx[dn_cell]; ++i) {
    pface_vtx_idx[i+1] = pface_vtx_idx[i] + pface_vtx_n[0][i];
  }

  PDM_log_trace_array_int(pface_vtx_n[0], dcell_face_idx[dn_cell], "pface_vtx_n::");
  PDM_log_trace_array_int(pface_vtx_idx, dcell_face_idx[dn_cell]+1, "pface_vtx_idx::");
  PDM_log_trace_array_long(pface_vtx[0], pface_vtx_idx[dcell_face_idx[dn_cell]], "pface_vtx::");

  /*
   * Free
   */
  free(entity_distribution_ptb);
  free(dcell_face_unsigned);
  PDM_block_to_part_free(btp);

  /*
   * Post-treat
   */
  *dcell_vtx_idx = (int*) malloc( (dn_cell + 1) * sizeof(int));
  int* _dcell_vtx_idx = *dcell_vtx_idx;
  int* dcell_vtx_n   = (int*) malloc( (dn_cell    ) * sizeof(int));

  int idx = 0;
  _dcell_vtx_idx[0] = 0;
  for(int i_cell = 0; i_cell < dn_cell; ++i_cell) {

    dcell_vtx_n  [i_cell  ] = 0;
    _dcell_vtx_idx[i_cell+1] = _dcell_vtx_idx[i_cell];

    int n_face_per_cell = dcell_face_idx[i_cell+1] - dcell_face_idx[i_cell];
    printf("n_face_per_cell::%i\n", n_face_per_cell);

    for(int i_face = 0; i_face < n_face_per_cell; ++i_face) {
      _dcell_vtx_idx[i_cell+1] += pface_vtx_n[0][idx];
      dcell_vtx_n[i_cell]     += pface_vtx_n[0][idx++];
    }
  }
  printf("idx::%i\n", idx);
  printf("dcell_face_idx[dn_cell]::%i\n", dcell_face_idx[dn_cell]);

  assert(idx == dcell_face_idx[dn_cell]);

  *dcell_vtx = pface_vtx[0]; // Car 1 seule partition

  PDM_log_trace_array_int(_dcell_vtx_idx, dn_cell+1, "_dcell_vtx_idx::");
  PDM_log_trace_array_int(dcell_vtx_n  , dn_cell  , "dcell_vtx_n::");

  PDM_para_graph_compress_connectivity2(dn_cell, dcell_vtx_idx, dcell_vtx_n, dcell_vtx);

  /*
   * Realloc
   */
  *dcell_vtx = (PDM_g_num_t *) realloc(dcell_vtx, sizeof(PDM_g_num_t) * _dcell_vtx_idx[dn_cell] );
  PDM_g_num_t* _dcell_vtx = *dcell_vtx;

  PDM_log_trace_array_int (_dcell_vtx_idx, dn_cell+1              , "after -> _dcell_vtx_idx::");
  PDM_log_trace_array_int (dcell_vtx_n   , dn_cell                , "after -> dcell_vtx_n::");
  PDM_log_trace_array_long(_dcell_vtx    , _dcell_vtx_idx[dn_cell], "after -> dcell_vtx::");

  /*
   * Free
   */
  free(dcell_vtx_n);
  free(pface_vtx_idx);
  free(pface_vtx);

}


#ifdef  __cplusplus
}
#endif
