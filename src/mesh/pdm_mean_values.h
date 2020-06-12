#ifndef __PDM_MEAN_VALUES_H__
#define __PDM_MEAN_VALUES_H__

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

/**
 * \brief Compute mean values of a list of points in a polygon
 *
 * \param [in]    n_pts            Number of points to locate
 * \param [in]    pts              xyz-Coordinates of points locate
 * \param [in]    n_vtx            Number of polygon vertices
 * \param [in]    poly_vtx         Polygon connectivity
 * \param [in]    mesh_vtx_coords  Coordinates of mesh vertices
 * \param [inout] mesh_vtx_coords  Mean value coordinates of points to locate
 *
 *
 */

void
PDM_mean_values_polygon_compute
(
 const int     n_pts,
 const double *pts,
 const int     n_vtx,
 const int    *poly_vtx,
 const double *mesh_vtx_coords,
 double       *mean_values
 );


void
PDM_mean_values_polygon_compute2
(
 const int     n_pts,
 const double *pts,
 const int     n_vtx,
 const int    *poly_vtx,
 const double *mesh_vtx_coords,
 double       *mean_values
 );

void
PDM_mean_values_polygon_compute3
(
 const int     n_pts,
 const double *pts_xyz,
 const int     n_vtx,
 const double *vtx_xyz,
 double       *mean_value_coords
 );







void
PDM_mean_values_polyhedron
(
 const int     n_pts,
 const double  pts_xyz[],
 const int     n_vtx,
 const double  vtx_xyz[],
 const int     n_faces,
 const int     face_vtx_idx[],
 const int     face_vtx[],
 double       *mean_value_coords
 );







void
PDM_mean_value_coordinates_polygon_2d
(
 const int    n_vtx,
 const double vtx_coord[],
 const int    n_pts,
 const double pts_coord[],
 double       mean_value_coord[]
 );

void
PDM_mean_value_coordinates_polygon_3d
(
 const int    n_vtx,
 const double vtx_coord[],
 const int    n_pts,
 const double pts_coord[],
 double       mean_value_coord[]
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MEAN_VALUES_H__ */
