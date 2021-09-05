#ifndef __PDM_PREDICATE_H__
#define __PDM_PREDICATE_H__


/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

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


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*
 *  from Triangle (https://github.com/libigl/triangle)
 */


#define REAL double
#define vertex double *

#define INEXACT /* Nothing */
/* #define INEXACT volatile */
void PDM_predicate_exactinit(void);

REAL PDM_predicate_incircle(vertex pa, vertex pb, vertex pc, vertex pd);

REAL PDM_predicate_insphere
(
 REAL *pa,
 REAL *pb,
 REAL *pc,
 REAL *pd,
 REAL *pe
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PREDICATE_H__ */
