/*
 * \file
 */

 #ifndef __PDM_CONVEX_H__
 #define __PDM_CONVEX_H__

 /*
   This file is part of the ParaDiGM library.

   Copyright (C) 2025       ONERA

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

 /*----------------------------------------------------------------------------
  *  Local headers
  *----------------------------------------------------------------------------*/

 #include "pdm.h"

 /*----------------------------------------------------------------------------*/

 #ifdef __cplusplus
 extern "C" {
 #if 0
 } /* Fake brace to force Emacs auto-indentation back to column 0 */
 #endif
 #endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

/**
 *
 * \brief Classify target points with respect to the convex hull of given source points.
 *
 * \warning Coordinates must always be defined in dimension 3, even if \p dim is lower
 *
 * \param [in]    dim         Spatial dimension (<= 3)
 * \param [in]    n_src       Number of source points
 * \param [in]    src_coord   Coordinates of source points (size = 3 * \p n_src)
 * \param [in]    n_tgt       Number of target points
 * \param [in]    tgt_coord   Coordinates of target points (size = 3 * \p n_tgt)
 * \param [inout] tgt_status  Status of each target point (0 = outside, 1 = inside) (size = \p n_tgt)
 *
 */
void
PDM_points_inside_convex_hull
(
  const int     dim,
  const int     n_src,
  const double *src_coord,
  const int     n_tgt,
  const double *tgt_coord,
        int    *tgt_status
);


/**
 *
 * \brief Determine if the current box intersects a given convex volume bounded by planes
 *
 * \param [in]   n_plane          Number of planes in the current volume
 * \param [in]   plane_origin     Coordinates of a point on each plane
 * \param [in]   plane_normal     Normal vector of each plane
 * \param [in]   box_extents      Extents of the box (x_min, y_min, z_min, x_max, y_max, z_max)
 *
 */
int
PDM_intersect_convex_volume_box
(
  const int  n_plane,
  double    *plane_origin,
  double    *plane_normal,
  double    *box_extents
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_CONVEX_H__ */
