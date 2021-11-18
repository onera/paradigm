
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
 const PDM_g_num_t  **block_distrib_idx,
 const int            n_block,
       int           *n_selected,
       PDM_g_num_t  **selected_g_num,
       int          **dmerge_idx,
       PDM_g_num_t  **dmerge_block_id,
       PDM_g_num_t  **dmerge_g_num,
       PDM_MPI_Comm   comm
)
{

  /*
   *  Le merge courant est doublement implicite :
   *     --> Pour chaque bloc on connait un voisins potentiel, can be NULL prt ?
   */



  PDM_multi_block_merge_t *mbm =
    (PDM_multi_block_merge_t *) malloc (sizeof(PDM_multi_block_merge_t));

  PDM_MPI_Comm_size (comm, &mbm->n_rank);
  PDM_MPI_Comm_rank (comm, &mbm->i_rank);
  mbm->comm        = comm;
  mbm->n_block     = n_block;

  // Ecrire l'algo dans distrib qui fait le multi_block from block

  // A/ Shift des blocs (selected_g_num)
  //    Shift des neighbor

  /*
   * Use ptb to equilibrate AND to reforms a block and keep link with the concatenate global numbering
   */
  PDM_part_to_block_t* ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      selected_g_num,
                                                      NULL,
                                                      n_selected,
                                                      n_block,
                                                      comm);

  int dn_merge  = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_g_num_t* blk_gnum = PDM_part_to_block_block_gnum_get (ptb);

  mbm->distrib_merge = PDM_compute_entity_distribution(comm, dn_merge);

  /*
   *  blk_gnum is partial because we implicitly remove all elements not refered
   *     -> The blk_gnum is exactly the new_to_old
   *     -> The information we need is the reverse one -> old_to_new to update afterwards connecitivty
   */
  int* dnew_to_old_idx = (int *) malloc( (dn_merge+1) * sizeof(int));
  dnew_to_old_idx[0] = 0;
  for(int i = 0; i < dn_merge; ++i) {
    dnew_to_old_idx[i+1] = dnew_to_old_idx[i] + 1;
  }

  /*
   *
   */
  PDM_g_num_t** multi_block_distrib_idx = (PDM_g_num_t **) malloc (mbm->n_block * sizeof(PDM_g_num_t *));
  PDM_g_num_t shift = 0;
  for(int i_block = 0; i_block < mbm->n_block; ++i_block) {
    multi_block_distrib_idx[i_block] = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (mbm->n_rank + 1));
    for(int i = 0; i < mbm->n_rank + 1; ++i) {
      multi_block_distrib_idx[i_block][i] = shift + block_distrib_idx[i_block][i];
    }
    shift += block_distrib_idx[i_block][mbm->n_rank];
    // PDM_log_trace_array_long(multi_block_distrib_idx[i_block], mbm->n_rank + 1, "multi_block_distrib_idx:: ");
  }

  // Compute the old global distribution
  PDM_g_num_t n_g_old = block_distrib_idx[n_block][mbm->n_rank];
  PDM_g_num_t* old_distrib = PDM_compute_uniform_entity_distribution(mbm->comm, n_g_old);

  PDM_dconnectivity_transpose(comm,
                              mbm->distrib_merge,
                              old_distrib,
                              dnew_to_old_idx,
                              blk_gnum,
                              0,
                              &mbm->dold_to_new_idx,
                              &mbm->dold_to_new);

  free(dnew_to_old_idx);
  free(old_distrib);
  for(int i_block = 0; i_block < mbm->n_block; ++i_block) {
    free(multi_block_distrib_idx[i_block]);
  }
  free(multi_block_distrib_idx);


  /*
   * If Hilbert we done it here to find a permutation to reorder the block
   */
  // if(reorder_block == PDM_HILBERT) {
  // }


  PDM_g_num_t* new_to_old_g_num = blk_gnum;

  mbm->mbtp = PDM_multi_block_to_part_create(mbm->multi_distrib_idx,
                                             mbm->n_block,
                       (const PDM_g_num_t**) mbm->block_distrib_idx,
                       (const PDM_g_num_t**)&new_to_old_g_num,
                                            &dn_merge,
                                             1,
                                             mbm->comm);
  PDM_part_to_block_free(ptb);

  // PDM_UNUSED(multi_distrib_idx);
  PDM_UNUSED(n_block);
  PDM_UNUSED(block_distrib_idx);
  // PDM_UNUSED(n_part);
  PDM_UNUSED(selected_g_num);

  return mbm;
}

// PDM_multi_block_merge_renum_from


#ifdef __cplusplus
}
#endif /* __cplusplus */
