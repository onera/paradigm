#ifndef __FVM_POINT_LOCATION_H__
#define __FVM_POINT_LOCATION_H__

/*============================================================================
 * Locate local points in a nodal representation associated with a mesh
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011-2018  ONERA

  Copyright (C) 2007-2009  EDF

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

/*----------------------------------------------------------------------------*/

/* #include "fvmc_config.h" */
/* #include "fvmc_config_defs.h" */

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_defs.h"
/* #include "fvmc_nodal.h" */

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

#define FVMC_POLYGON_FAILURE -1
#define FVMC_POLYGON_OUTSIDE 0
#define FVMC_POLYGON_INSIDE 1
#define FVMC_POLYGON_INTERSECTION 2
#define FVMC_POLYGON_ON_LINE 3

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute distance to polygons
 *
 * parameters:
 *   dim               <-- dimension
 *   n_poly            <-- number of polygon
 *   connectivity_idx  <-- polygon connectivity index
 *   connectivity      <-- polygon connectivity
 *   vertex_coords     <-- polygon connectivity
 *   n_points          <-- polygon connectivity
 *   point_coords      <-- polygon connectivity
 *   distance          --> 3d surf : distance from point to element indicated by
 *                         location[]: < 0 if unlocated, or absolute
 *                         distance to a surface element (size: n_points)
 *                         2d or 3d : distance from point to element indicated by
 *                         location[]: < 0 if unlocated, 0 - 1 if inside,
 *                          > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

void
fvmc_point_dist_closest_polygon(const int            dim,
                                const fvmc_lnum_t    n_poly,
                                const fvmc_lnum_t    connectivity_idx[],
                                const fvmc_lnum_t    connectivity[],
                                const fvmc_coord_t   vertex_coords[],
                                const fvmc_lnum_t    n_points,
                                const fvmc_lnum_t    point_ids[],
                                const fvmc_coord_t   point_coords[],
                                fvmc_lnum_t          location[],
                                float                distance[]);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVM_POINT_LOCATION_H__ */
