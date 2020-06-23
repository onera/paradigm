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
PDM_locate_points_on_triangles (const int          n_tri,
                                const PDM_l_num_t  tri_vtx[],
                                const double       vtx_coord[],
                                const int          n_pts,
                                const double       pts_coord[],
                                const PDM_g_num_t  pts_g_num[],//debug only
                                int                location[],
                                float              distance[],
                                double             bar_coord[]);

void
PDM_locate_points_on_quad (const double       vtx_xyz[12],
                           const int          n_pts,
                           const double       pts_xyz[],
                           const PDM_g_num_t  pts_g_num[],
                           float             *distance,
                           double            *bary_coords);

void
PDM_locate_points_in_tetra (const double       vtx_xyz[12],
                            const int          n_pts,
                            const double       pts_xyz[],
                            const PDM_g_num_t  pts_g_num[],
                            float             *distance,
                            double            *bary_coords);

void
PDM_locate_points_in_cell (const PDM_Mesh_nodal_elt_t  elt_type,
                           const PDM_l_num_t           cell_vtx[],
                           const double                vtx_coord[],
                           const int                   n_pts,
                           const double                pts_coord[],
                           const PDM_g_num_t           pts_g_num[],//debug only
                           float                      *distance,
                           double                     *bar_coord);

#if 1
void PDM_point_location_distance
(
 const PDM_Mesh_nodal_elt_t elt_type,
 const int                  n_pts,
 const double               uvw[],
 double                     shapef[],
 double                     distance[]
 );
#endif
//<<--

void
PDM_locate_in_polyhedron
(
 const PDM_l_num_t n_vtx,
 const double      vtx_coord[],
 const PDM_l_num_t n_face,
 const PDM_l_num_t face_vtx_idx[],
 const PDM_l_num_t face_vtx[],
 const int         face_orientation[],
 const int         n_pts,
 const double      pts_coord[],
 const double      char_length,
 float             distance[],
 double            bar_coord[]
 );


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_POINT_LOCATION_H__ */
