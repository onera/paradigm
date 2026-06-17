#ifndef __PDM_QUATERNION_PRIV_H__
#define __PDM_QUATERNION_PRIV_H__

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

struct _pdm_quaternion_t {

  double q[4];
  double q_squared[4];

};

/*=============================================================================
 * Static global variables
 *============================================================================*/


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif	/* __PDM_QUATERNION_PRIV_H__ */
