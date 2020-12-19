/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_part_extension.h"
#include "pdm_part_extension_priv.h"

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
 * \brief Create a part extension structure
 *
 * \param [in]   comm           Communicator
 *
 */

PDM_part_extension_t*
PDM_part_extension_create
(
 const int              n_block,
 const int             *n_part,
 const PDM_MPI_Comm     comm,
 const PDM_ownership_t  owner
)
{
 PDM_part_extension_t *part_ext = (PDM_part_extension_t *) malloc(sizeof(PDM_part_extension_t));

  part_ext->n_block = n_block;
  part_ext->n_part  = n_part;
  part_ext->comm    = comm;
  part_ext->owner   = owner;

  return part_ext;
}


/**
 *
 * \brief Free a part extension structure
 *
 * \param [in]   part_ext          PDM_part_extension_t
 *
 */
void
PDM_part_extension_free
(
 PDM_part_extension_t *part_ext
)
{
  part_ext->n_part = NULL;
  free(part_ext);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
