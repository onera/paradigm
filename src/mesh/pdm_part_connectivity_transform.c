
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_block_to_part.h"
#include "pdm_part_to_block.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_para_graph_dual.h"
#include "pdm_error.h"
#include "pdm_timer.h"
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

/*============================================================================
 * Static global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/



void
PDM_part_combine_connectivity
(
 int   n_entity1,
 int  *entity1_entity2_idx,
 int  *entity1_entity2,
 int  *entity2_entity3_idx,
 int  *entity2_entity3,
 int **entity1_entity3_idx,
 int **entity1_entity3
)
{
  PDM_UNUSED(n_entity1);
  PDM_UNUSED(entity1_entity2_idx);
  PDM_UNUSED(entity1_entity2);
  PDM_UNUSED(entity2_entity3_idx);
  PDM_UNUSED(entity2_entity3);
  PDM_UNUSED(entity1_entity3_idx);
  PDM_UNUSED(entity1_entity3);
}



void
PDM_part_connectivity_transpose
(
 int   n_entity1,
 int  *entity1_entity2_idx,
 int  *entity1_entity2,
 int **entity2_entity1_idx,
 int **entity2_entity1
)
{
  PDM_UNUSED(n_entity1);
  PDM_UNUSED(entity1_entity2_idx);
  PDM_UNUSED(entity1_entity2);
  PDM_UNUSED(entity2_entity1_idx);
  PDM_UNUSED(entity2_entity1);

}
