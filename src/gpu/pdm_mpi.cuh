#ifndef __PDM_MPI_CUH__
#define __PDM_MPI_CUH__

/*============================================================================
 * Bibliotheque de messagerie
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_rank (wrapping de la fonction MPI_Comm_rank)
 *
 *----------------------------------------------------------------------------*/

__device__ int PDM_MPI_Comm_rank_GPU(PDM_MPI_Comm comm, int *rank);

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_size (wrapping de la fonction MPI_Comm_size)
 *
 *----------------------------------------------------------------------------*/

__device__ int PDM_MPI_Comm_size_GPU(PDM_MPI_Comm comm, int *size);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MPI_CUH__ */