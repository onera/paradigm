#ifndef __PDM_POINT_LOCATION_H__
#define __PDM_POINT_LOCATION_H__

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

void
PDM_point_location_nodal
(
 const int            mesh_nodal_id,
 const int            n_pts,
 const int           *pts_idx,
 const double        *pts_coords,
 const PDM_g_num_t   *pts_g_num,//
 const double         tolerance,
 int                  base_element_num,//?
 float              **distance,
 double             **projected_coords,
 int                **bar_coord_idx,
 double             **bar_coord
 );



//-->>
void
PDM_locate_points_in_tetra (const double       vtx_xyz[12],
                            const int          n_pts,
                            const double       pts_xyz[],
                            const PDM_g_num_t  pts_g_num[],
                            float             *distance,
                            double            *bary_coords);
//<<--
#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_POINT_LOCATION_H__ */
