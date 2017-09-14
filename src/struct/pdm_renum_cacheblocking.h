#ifndef __PDM_RENUM_CACHE_BLOCKING_H__
#define __PDM_RENUM_CACHE_BLOCKING_H__

/*============================================================================
 * Hilbert space-filling curve construction for coordinates.
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

#include "pdm_part.h"
#include "pdm_part_priv.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Macro and type definitions
 *============================================================================*/


/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Perform a cells renumbering with cache blocking (Synchrone)
 *
 * parameters:
 *   part       --> Mesh Partition 
 *---------------------------------------------------------------------------*/

void 
PDM_renum_cacheblocking
(
 _part_t     *part,
int          split_method, 
int          nCellPerCacheWanted, 
int          isAsynchrone,  
int          isVectorisation 
);


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_RENUM_CACHE_BLOCKING_H__ */
