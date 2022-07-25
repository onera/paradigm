#ifndef __PDM_HO_BEZIER_H__
#define __PDM_HO_BEZIER_H__

/*============================================================================
 * Functions about high order meshes
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2018       ONERA

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

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/


void
PDM_ho_bezier_de_casteljau_tria
(
 const int     dim,
 const int     order,
 const double  u,
 const double  v,
 double       *b,
 double       *val,
 double       *atr,
 double       *ars,
 double       *ast
 );


void
PDM_ho_bezier_triangle_derivatives
(
 const int  dim,
 const int  order,
 double    *b,
 double    *bu,
 double    *bv
 );

double
PDM_ho_bezier_tria_location
(
 const int     order,
 const int     n_node,
       double *node_coord,
       double *point_coord,
       double *projected_coord,
       double *uvw
 );

#endif /* __PDM_HO_BEZIER_H__ */
