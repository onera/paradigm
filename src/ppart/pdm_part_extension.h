#ifndef __PDM_PART_EXTENSION_H__
#define __PDM_PART_EXTENSION_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"

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

typedef struct _pdm_part_extension_t PDM_part_extension_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


PDM_part_extension_t*
PDM_part_extension_create
(
 const int              n_block,
 const int             *n_part,
 const PDM_MPI_Comm     comm,
 const PDM_ownership_t  owner
);


void
PDM_part_extension_free
(
 PDM_part_extension_t *part_ext
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_EXTENSION_H__ */
