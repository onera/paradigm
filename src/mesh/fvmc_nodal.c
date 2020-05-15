/*============================================================================
 * Main structure for a nodal representation associated with a mesh
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

  Copyright (C) 2004-2009  EDF

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

/* #include <bftc_mem.h> */
/* #include <bftc_printf.h> */

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/* #include "fvmc_config_defs.h" */
#include "fvmc_defs.h"
/* #include "fvmc_io_num.h" */
/* #include "fvmc_parall.h" */
/* #include "fvmc_tesselation.h" */

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_nodal.h"
/* #include "fvmc_nodal_priv.h" */

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get number of vertices.
 *
 * returns:
 *   Number of vertices
 *----------------------------------------------------------------------------*/

int
fvmc_nodal_n_vertices_element (fvmc_element_t type, int order)
{
 int n_vtx = 0;
 int _order = order;
 if (order == -1) {
   _order = 1;
 }

 switch(type) {
 case FVMC_EDGE:               /* Edge */
   n_vtx = (_order+1);
   break;
 case FVMC_FACE_TRIA:          /* Triangle */
   n_vtx = (_order+1)*(_order+2)/2;
   break;
 case FVMC_FACE_QUAD:          /* Quadrangle */
   n_vtx = (_order+1)*(_order+1);
   break;
 case FVMC_FACE_POLY:          /* Simple Polygon */
   n_vtx = -1;
   break;
 case FVMC_CELL_TETRA:         /* Tetrahedron */
   n_vtx = (_order+1)*(_order+2)*(_order+3)/6;
   break;
 case FVMC_CELL_PYRAM:         /* Pyramid */
   n_vtx = (_order+1)*(_order+2)*(2*_order+3)/6;
   break;
 case FVMC_CELL_PRISM:         /* Prism (pentahedron) */
   n_vtx = (_order+1)*(_order+1)*(_order+2)/2;
   break;
 case FVMC_CELL_HEXA:         /* Hexahedron (brick) */
   n_vtx = (_order+1)*(_order+1)*(_order+1);
   break;
 case FVMC_CELL_POLY:          /* Simple Polyhedron (convex or quasi-convex) */
   n_vtx = -1;
   break;
 default:
   n_vtx = -1;
 }

 return n_vtx;
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
