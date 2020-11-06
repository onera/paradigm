/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_multi_block_to_part.h"
#include "pdm_multi_block_to_part_priv.h"
#include "pdm_binary_search.h"
#include "pdm_priv.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/


/*=============================================================================
 * Public function definitions
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
)
{
  printf(" ola que tal PDM_multi_block_to_part_t \n");

  _pdm_multi_block_to_part_t *mbtp =
    (_pdm_multi_block_to_part_t *) malloc (sizeof(_pdm_multi_block_to_part_t));

  mbtp->comm = comm;
  PDM_MPI_Comm_size (comm, &mbtp->n_rank);
  PDM_MPI_Comm_rank (comm, &mbtp->i_rank);

  return (PDM_multi_block_to_part_t *) mbtp;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
