
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
 PDM_g_num_t  *distrib_block_init,
 PDM_g_num_t  *distri_block_add,
 int           n_update,
 int           n_remove,
 int           n_add,
 PDM_g_num_t  *ln_to_gn_update,
 PDM_g_num_t  *ln_to_gn_remove,
 PDM_g_num_t  *ln_to_gn_add_in,
 PDM_MPI_Comm  comm
)
{
  PDM_UNUSED(distrib_block_init);
  PDM_UNUSED(distri_block_add);
  PDM_UNUSED(n_update);
  PDM_UNUSED(n_remove);
  PDM_UNUSED(n_add);
  PDM_UNUSED(ln_to_gn_update);
  PDM_UNUSED(ln_to_gn_remove);
  PDM_UNUSED(ln_to_gn_add_in);
  PDM_UNUSED(comm);

  return NULL;
}




#ifdef __cplusplus
}
#endif /* __cplusplus */
