#ifndef __PDM_ROTATION_PRIV_H__
#define __PDM_ROTATION_PRIV_H__

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif
/*============================================================================
 * Type
 *============================================================================*/

typedef struct {

    double q[4];
    double q_squared[4];

} _pdm_quaternion_t;

/**
 * \struct _pdm_twist_t
 *
 * \brief  Structuring encoding a twist (in french vissage) motion
 *
 */

typedef struct {

  _pdm_quaternion_t quaternion;
  double            translation[3];

} _pdm_twist_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/


#ifdef	__cplusplus
}
#endif

#endif	/* __PDM_ROTATION_PRIV_H__ */