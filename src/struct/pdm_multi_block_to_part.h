#ifndef PDM_MULTI_BLOCK_TO_PART_H
#define PDM_MULTI_BLOCK_TO_PART_H

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#if !defined (__hpux) && !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

#ifdef  __cplusplus
extern "C" {
#endif

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \struct PDM_multi_block_to_part_t
 * \brief  Block to partition redistribution
 *
 */

typedef struct _pdm_multi_block_to_part_t PDM_multi_block_to_part_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Create a block to partitions redistribution
 *
 * \param [in]   blockDistribIdx Block distribution (size : \ref size of \ref comm + 1)
 * \param [in]   gnum_elt        Element global number (size : \ref n_part)
 * \param [in]   n_elt           Local number of elements (size : \ref n_part)
 * \param [in]   n_part          Number of partition
 * \param [in]   comm            MPI communicator
 *
 * \return   Initialized \ref PDM_block_to_part instance
 *
 */

PDM_multi_block_to_part_t *
PDM_multi_block_to_part_create
(
 const PDM_g_num_t   *blockDistribIdx,
 const PDM_g_num_t  **gnum_elt,
 const int           *n_elt,
 const int            n_part,
 const PDM_MPI_Comm   comm
);

#ifdef  __cplusplus
}
#endif

#endif  /* PDM_BLOCK_TO_PART_H */
