/*
 * \file
 */

#ifndef __PDM_SIZE_IDX_FROM_STRIDE__
#define __PDM_SIZE_IDX_FROM_STRIDE__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_part_geom.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


/**
 *
 * \brief Return the integer size required to store an idx array from stride array
 *
 * \param [in]   stride    Stride array
 * \param [in]   s_sride   Size of the stride array
 * \param [in]   comm      MPI communicator
 *
 * \return the size of int required
 *
 */

int
PDM_size_idx_from_stride
(
 int          *stride,
 int          s_stride,
 PDM_MPI_Comm comm
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_SIZE_IDX_FROM_STRIDE__ */
