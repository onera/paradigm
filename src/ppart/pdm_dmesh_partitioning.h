#ifndef __PDM_PART_H__
#define __PDM_PART_H__

#include <stdio.h>
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

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \enum PDM_part_split_t
 * \brief Split method
 *
 */

typedef enum {
  PDM_PART_SPLIT_PARMETIS = 1,
  PDM_PART_SPLIT_PTSCOTCH = 2,
  PDM_PART_SPLIT_HILBERT  = 3
} PDM_part_split_t;

typedef struct _PDM_part_t PDM_part_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_part_H__ */
