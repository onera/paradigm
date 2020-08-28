#ifndef __PDM_POINT_LOCATION_H__
#define __PDM_POINT_LOCATION_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mesh_nodal.h"

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
//-->> TO REMOVE
PDM_bool_t PDM_point_location_uvw (const PDM_Mesh_nodal_elt_t elt_type,
                                   const double               vtx_coord[],
                                   const double               tolerance,
                                   const double               pt_coord[3],
                                   double                     uvw[3]);

void PDM_locate_points_in_cell
(
 const PDM_Mesh_nodal_elt_t  elt_type,
 const PDM_l_num_t           cell_vtx[],
 const double                vtx_coord[],
 const int                   n_pts,
 const double                pts_coord[],
 float                      *distance,
 double                     *bar_coord
 );
//<<--



void
PDM_point_location_nodal
(
 const int           type_idx[],
 const PDM_g_num_t   elt_g_num[],
 const int           elt_vtx_idx[],
 const double        elt_vtx_coord[],
 const PDM_l_num_t   poly3d_face_idx[],
 const PDM_l_num_t   face_vtx_idx[],
 const PDM_l_num_t   face_vtx[],
 const int           face_orientation[],
 const double        poly3d_char_length[],
 const int           pts_idx[],
 const double        pts_coord[],
 const double        tolerance,
 float             **distance,
 //double            **projected_coord,//high-order
 int               **bar_coord_idx,
 double            **bar_coord
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_POINT_LOCATION_H__ */
