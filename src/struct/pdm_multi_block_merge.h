#ifndef __PDM_MULTI_BLOCK_MERGE_H__
#define __PDM_MULTI_BLOCK_MERGE_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_io.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Types definition
 *============================================================================*/

typedef struct _pdm_multi_block_merge_t      PDM_multi_block_merge_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/**
 * \brief Create a redistribute structure
 *
 * \param [in]   comm     MPI communicator
 *
 * \return       New mesh nodal handle
 *
 */
PDM_multi_block_merge_t*
PDM_multi_block_merge_create
(
 PDM_g_num_t  *distrib_block_init,
 PDM_g_num_t  *distri_block_add,
 int           n_update,
 int           n_remove,
 int           n_add,
 PDM_g_num_t  *ln_to_gn_update,
 PDM_g_num_t  *ln_to_gn_remove,
 PDM_g_num_t  *ln_to_gn_add_in,
 PDM_MPI_Comm  comm
);



void
PDM_multi_block_merge_reorder_block_hilbert
(
 PDM_multi_block_merge_t  *mbm,
 double                   *block_coord,
 double                   *padd_coord
);

void
PDM_multi_block_merge_exch
(
 PDM_multi_block_merge_t  *mbm,
 size_t                    s_data,
 PDM_stride_t              t_stride,
 int                       cst_stride,
 int                      *blk_merge_stride,
 void                     *blk_merge_data,
 int                      *padd_stride,
 void                     *padd_data,
 int                     **block_strid_out,
 void                    **block_data_out
);


PDM_g_num_t*
PDM_multi_block_merge_compute_old_to_new
(
 PDM_multi_block_merge_t  *mbm
);

void
PDM_multi_block_merge_get_old_to_new
(
 PDM_multi_block_merge_t  *mbm,
 int                **dold_to_new_idx,
 PDM_g_num_t        **dold_to_new
);

int
PDM_multi_block_merge_get_n_block
(
 PDM_multi_block_merge_t  *mbm
);

PDM_g_num_t*
PDM_multi_block_merge_get_parent_blk_g_num
(
 PDM_multi_block_merge_t  *mbm
);

/**
 * \brief Create a redistribute structure
 *
 * \param [in]   comm     MPI communicator
 *
 * \return       New mesh nodal handle
 *
 */
void
PDM_multi_block_merge_free
(
PDM_multi_block_merge_t   *mbm
);

void
PDM_multi_block_merge_exch_and_update_child_g_num
(
 PDM_multi_block_merge_t  *mbm,
 PDM_g_num_t              *child_old_distrib,
 int                      *dchild_old_to_new_idx,
 PDM_g_num_t              *dchild_old_to_new,
 PDM_stride_t              t_stride,
 int                       cst_stride,
 int                      *blk_merge_stride,
 void                     *blk_merge_data,
 int                      *padd_stride,
 void                     *padd_data,
 int                     **block_strid_out,
 void                    **block_data_out
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MULTI_BLOCK_MERGE_H__ */
