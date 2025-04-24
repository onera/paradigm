/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_exchange_helper.h"
#include "pdm.h"
#include "pdm_array.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include "pdm_mpi.h"
#include "pdm_exchange_helper_priv.h"
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

PDM_exchange_helper_t *
PDM_exchange_helper_create
(
 const PDM_MPI_Comm    comm,
       int             n_request_init
)
{
  PDM_exchange_helper_t *exch_helper = NULL;
  PDM_malloc(exch_helper, 1, PDM_exchange_helper_t);

  exch_helper->comm      = comm;
  exch_helper->n_request = n_request_init;

  return exch_helper;
}

void
PDM_exchange_helper_free
(
  PDM_exchange_helper_t *exch_helper
)
{
  PDM_free(exch_helper);
}






#ifdef __cplusplus
}
#endif
