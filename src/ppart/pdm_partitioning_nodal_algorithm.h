#ifndef __PDM_PARTITIONING_NODAL_ALGORITHM_H__
#define __PDM_PARTITIONING_NODAL_ALGORITHM_H__

#include <stdio.h>
#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_dmesh_nodal.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

void
PDM_dmesh_nodal_to_pmesh_nodal
(
 PDM_dmesh_nodal_t* dmn
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_PARTITIONING_NODAL_ALGORITHM_H__ */
