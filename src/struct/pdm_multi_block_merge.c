
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_multi_block_merge.h"
#include "pdm_multi_block_merge_priv.h"
#include "pdm_block_to_part.h"
#include "pdm_multi_block_to_part.h"
#include "pdm_array.h"
#include "pdm_part_geom.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_distrib.h"
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

/*----------------------------------------------------------------------------
 * Maximum number of sections depending of section type
 *----------------------------------------------------------------------------*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/


PDM_multi_block_merge_t*
PDM_multi_block_merge_create
(
       PDM_g_num_t  **block_distrib_idx,
 const int            n_block,
       int           *n_selected,
       PDM_g_num_t  **selected_g_num,
       int          **dmerge_idx,
       int          **dmerge_block_id,
       PDM_g_num_t  **dmerge_g_num,
       PDM_MPI_Comm   comm
)
{
  PDM_UNUSED(dmerge_idx);
  PDM_UNUSED(dmerge_block_id);
  PDM_UNUSED(dmerge_g_num);

  PDM_multi_block_merge_t *mbm = (PDM_multi_block_merge_t *) malloc (sizeof(PDM_multi_block_merge_t));

  PDM_MPI_Comm_size (comm, &mbm->n_rank);
  PDM_MPI_Comm_rank (comm, &mbm->i_rank);

  mbm->comm        = comm;
  mbm->n_block     = n_block;

  mbm->multi_block_distrib = malloc( (n_block + 1) * sizeof(PDM_g_num_t));
  mbm->multi_block_distrib[0] = 0;
  for(int i_block = 0; i_block < n_block; ++i_block) {
    mbm->multi_block_distrib[i_block+1] = mbm->multi_block_distrib[i_block] + block_distrib_idx[i_block][mbm->n_rank];
  }

  /*
   * Shift des données d'entrés
   */
  PDM_g_num_t **_selected_g_num  = malloc( (n_block * 2) * sizeof(PDM_g_num_t *));
  PDM_g_num_t **_send_orig_g_num = malloc( (n_block * 2) * sizeof(PDM_g_num_t *));
  int         **_select_kind     = malloc( (n_block * 2) * sizeof(int         *));
  int         **_stride_one      = malloc( (n_block * 2) * sizeof(int         *));
  PDM_g_num_t  *_n_selected      = malloc( (n_block * 2) * sizeof(PDM_g_num_t  ));
  for(int i_block = 0; i_block < n_block; ++i_block) {

    _n_selected     [i_block] = n_selected[i_block];
    _selected_g_num [i_block] = malloc(n_selected[i_block] * sizeof(PDM_g_num_t));
    _select_kind    [i_block] = malloc(n_selected[i_block] * sizeof(int        ));
    _stride_one     [i_block] = malloc(n_selected[i_block] * sizeof(int        ));
    _send_orig_g_num[i_block] = _selected_g_num[i_block];

    for(int i = 0; i < n_selected[i_block]; ++i) {
      _selected_g_num[i_block][i] = mbm->multi_block_distrib[i_block] + selected_g_num[i_block][i];
      _select_kind   [i_block][i] = 1;
      _stride_one    [i_block][i] = 1;
    }
  }


  for(int i_block = 0; i_block < n_block; ++i_block) {
    int dn_size = block_distrib_idx[i_block][mbm->i_rank+1] - block_distrib_idx[i_block][mbm->i_rank];
    _n_selected     [n_block+i_block] = dmerge_idx[i_block][dn_size];
    _selected_g_num [n_block+i_block] = malloc(dmerge_idx[i_block][dn_size] * sizeof(PDM_g_num_t));
    _send_orig_g_num[n_block+i_block] = malloc(dmerge_idx[i_block][dn_size] * sizeof(PDM_g_num_t));
    _select_kind    [n_block+i_block] = malloc(dmerge_idx[i_block][dn_size] * sizeof(int        ));
    _stride_one     [n_block+i_block] = malloc(dmerge_idx[i_block][dn_size] * sizeof(int        ));
    for(int i = 0; i < dn_size; ++i) {
      for(int j = dmerge_idx[i_block][i]; j < dmerge_idx[i_block][i+1]; ++j) {
        int i_block_opp = dmerge_block_id[i_block][j];
        _selected_g_num [n_block+i_block][j] = mbm->multi_block_distrib[i_block_opp] + dmerge_g_num[i_block][j];
        _send_orig_g_num[n_block+i_block][j] = mbm->multi_block_distrib[i_block    ] + i + 1;
        if(i_block_opp < i_block) {
          _select_kind   [n_block+i_block][j] = 2;
        } else {
          _select_kind   [n_block+i_block][j] = 3;
        }
        _stride_one    [n_block+i_block][j] = 1;
      }
    }
  }

  /*
   * Use ptb to equilibrate AND to reforms a block and keep link with the concatenate global numbering
   */
  PDM_part_to_block_t* ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      _selected_g_num,
                                                      NULL,
                                                      _n_selected,
                                                      2 * n_block,
                                                      comm);


  int *blk_select_kind_n = NULL;
  int *blk_select_kind   = NULL;
  int dn_debug = PDM_part_to_block_exch(ptb,
                                        sizeof(int),
                                        PDM_STRIDE_VAR,
                                        1,
                                        _stride_one,
                              (void **) _select_kind,
                                        &blk_select_kind_n,
                              (void **) &blk_select_kind);
  free(blk_select_kind_n);
  PDM_g_num_t *blk_send_orig_g_num = NULL;
  PDM_part_to_block_exch(ptb,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR,
                         1,
                         _stride_one,
               (void **) _send_orig_g_num,
                         &blk_select_kind_n,
               (void **) &blk_send_orig_g_num);


  if(1 == 1) {
    PDM_log_trace_array_int (blk_select_kind    , dn_debug, "blk_select_kind :: ");
    PDM_log_trace_array_long(blk_send_orig_g_num, dn_debug, "blk_send_orig_g_num :: ");
  }

  int dn_merge  = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_g_num_t* blk_gnum = PDM_part_to_block_block_gnum_get (ptb);

  int* blk_select_kind_idx = malloc( (dn_merge + 1) * sizeof(int));
  blk_select_kind_idx[0] = 0;
  for(int i = 0; i < dn_merge; ++i) {
    blk_select_kind_idx[i+1] = blk_select_kind_idx[i] + blk_select_kind_n[i];
  }


  int* dnew_to_old_idx = (int *) malloc( (dn_merge+1) * sizeof(int));
  dnew_to_old_idx[0] = 0;
  int idx_read = 0;
  int dn_new_block = 0;
  for(int i = 0; i < dn_merge; ++i ) {

    int dn_loc = blk_select_kind_idx[i+1] - blk_select_kind_idx[i];
    PDM_log_trace_array_long(&blk_send_orig_g_num[blk_select_kind_idx[i]], dn_loc, "blk_send_orig_g_num :: ");
    PDM_log_trace_array_int (&blk_select_kind    [blk_select_kind_idx[i]], dn_loc, "blk_select_kind     :: ");

    /*
     *  Normal / Join /
     */
    int is_normal = 1;
    int is_merge  = 0;
    int is_remove = 0;
    dnew_to_old_idx[i+1] = dnew_to_old_idx[i];
    for(int j = blk_select_kind_idx[i]; j < blk_select_kind_idx[i+1]; ++j) {
      if(blk_select_kind[j] == 2) {
        is_merge  = 1;
        is_normal = 0;
      } else if(blk_select_kind[j] == 3) {
        is_remove = 1;
        is_normal = 0;
      }
      dnew_to_old_idx[i+1]++;
    }

    if(is_remove) {
      assert(is_merge == 0);
    }
    if(is_merge) {
      assert(is_remove == 0);
    }

    log_trace("is_normal = %i | is_merge = %i | is_remove = %i \n", is_normal, is_merge, is_remove);

  }
  // assert(dn_debug == idx_read);



  for(int i_block = 0; i_block < 2 * n_block; ++i_block) {
    free(_selected_g_num[i_block]);
  }

  for(int i_block = n_block; i_block < 2 * n_block; ++i_block) {
    free(_send_orig_g_num[i_block]);
  }
  free(_selected_g_num);
  free(_send_orig_g_num);
  free(_n_selected);
  free(blk_select_kind_n);
  free(blk_select_kind);
  free(blk_select_kind_idx);



  PDM_log_trace_array_long(blk_gnum, dn_merge, "blk_gnum :: ");
  PDM_log_trace_array_int (dnew_to_old_idx, dn_merge+1, "dnew_to_old_idx :: ");

  mbm->distrib_merge = PDM_compute_entity_distribution(comm, dn_merge);


  for(int i_block = 0; i_block < 2 * n_block; ++i_block) {
    free(_select_kind[i_block]);
    free(_stride_one[i_block]);
  }
  free(_select_kind);
  free(_stride_one);

  /*
   *  blk_gnum is partial because we implicitly remove all elements not refered
   *     -> The blk_gnum is exactly the new_to_old
   *     -> The information we need is the reverse one -> old_to_new to update afterwards connecitivty
   */
  // for(int i = 0; i < dn_merge; ++i) {
  //   dnew_to_old_idx[i+1] = dnew_to_old_idx[i] + 1;
  // }

  PDM_log_trace_array_long(mbm->multi_block_distrib, n_block+1, "mbm->multi_block_distrib_idx");

  // Compute the old global distribution
  PDM_g_num_t n_g_old      = mbm->multi_block_distrib[n_block];
  PDM_g_num_t* old_distrib = PDM_compute_uniform_entity_distribution(mbm->comm, n_g_old);

  PDM_dconnectivity_transpose(comm,
                              mbm->distrib_merge,
                              old_distrib,
                              dnew_to_old_idx,
                              blk_send_orig_g_num,
                              0,
                              &mbm->dold_to_new_idx,
                              &mbm->dold_to_new);

  free(dnew_to_old_idx);

  /*
   * If Hilbert we done it here to find a permutation to reorder the block
   */
  // if(reorder_block == PDM_HILBERT) {
  // }

  int dn_orig = old_distrib[mbm->i_rank+1] - old_distrib[mbm->i_rank];
  PDM_log_trace_array_long(blk_send_orig_g_num, mbm->dold_to_new_idx[dn_orig], "new_to_old_g_num");

  free(old_distrib);
  mbm->mbtp = PDM_multi_block_to_part_create(mbm->multi_block_distrib,
                                             mbm->n_block,
                       (const PDM_g_num_t**) block_distrib_idx,
                       (const PDM_g_num_t**)&blk_send_orig_g_num,
                                            &dn_merge,
                                             1,
                                             mbm->comm);
  PDM_part_to_block_free(ptb);
  free(blk_send_orig_g_num);

  return mbm;
}


void
PDM_multi_block_merge_free
(
 PDM_multi_block_merge_t* mbm
)
{

  free(mbm->multi_block_distrib);
  free(mbm->dold_to_new_idx);
  free(mbm->dold_to_new);
  free(mbm->distrib_merge);

  PDM_multi_block_to_part_free(mbm->mbtp);
  free(mbm);
}


// PDM_multi_block_merge_renum_from

#ifdef __cplusplus
}
#endif /* __cplusplus */
