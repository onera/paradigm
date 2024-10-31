#ifndef __PDM_QUATERNION_PRIV_H__
#define __PDM_QUATERNION_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation csback to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

struct _pdm_quaternion_t {
  double w;       /**< Scalar part */
  double v[3];    /**< Vector part */
};


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_QUATERNION_PRIV_H__ */
