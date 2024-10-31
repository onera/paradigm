#ifndef __PDM_QUATERNION_H__
#define __PDM_QUATERNION_H__

/*----------------------------------------------------------------------------*/

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
 * Type definitions
 *============================================================================*/

typedef struct _pdm_quaternion_t PDM_quaternion_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


PDM_quaternion_t*
PDM_quaternion_create
(
  double w,
  double v1,
  double v2,
  double v3
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_QUATERNION_H__ */
