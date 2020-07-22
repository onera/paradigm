#ifndef __PDM_POINT_LOCATION_H__
#define __PDM_POINT_LOCATION_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
//-->>
#include "pdm_mesh_nodal.h"
//<<--

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

PDM_bool_t PDM_point_location_uvw (const PDM_Mesh_nodal_elt_t elt_type,
                                   const double               vtx_coord[],
                                   const double               tolerance,
                                   const double               pt_coord[3],
                                   double                     uvw[3]);

void
PDM_point_location_nodal
(
 const int            mesh_nodal_id,
 const int            n_pts,
 const int           *pts_idx,
 const double        *pts_coords,
 const double         tolerance,
 float              **distance,
 double             **projected_coords,
 int                **bar_coord_idx,
 double             **bar_coord
 );


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_POINT_LOCATION_H__ */
