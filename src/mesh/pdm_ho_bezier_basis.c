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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_error.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_vtk.h"
#include "pdm_ho_location.h"
#include "pdm_box_tree.h"
#include "pdm_mesh_nodal.h"
#include "pdm_ho_seg_intersect.h"
#include "pdm_line.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Edge Bézier basis
 *
 * \param [in]  order           Order
 * \param [in]  n_pts           Number of points
 * \param [in]  u               Parametric coordinates (size = \ref n_pts)
 * \param [out] weights         Weights (size = n_nodes * \ref n_pts)
 *
 */

static void
_basis_bezier_edge
(
 const int              order,
 const int              n_pts,
 const double *restrict u,
 double       *restrict weights
)
{

  const int n_nodes = order + 1;
}

/**
 *
 * \brief Triangle Bézier basis
 *
 * \param [in]  order           Order
 * \param [in]  n_pts           Number of points
 * \param [in]  u               Parametric coordinates (size = \ref n_pts)
 * \param [out] weights         Weights (size = n_nodes * \ref n_pts)
 *
 */

static void
_basis_bezier_tria
(
 const int              order,
 const int              n_pts,
 const double *restrict u,
 double       *restrict weights
)
{
  const int n_nodes = order + 1;

  if (order == 1) {
  }

  if (order == 2) {
  }

  if (order == 3) {
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Evaluate high-order basis Bézier functions
 *
 *
 * \param [in]  type      Element type structure
 * \param [in]  order     Element order
 * \param [in]  n_nodes   Number of nodes
 * \param [in]  n_pts     Number of points
 * \param [in]  uvw       Parametric coordinates of the points (size = elt_dim * \ref n_pts)
 * \param [out] weights   Weights (size = \ref n_pts * \ref n_nodes)
 *
 */

void
PDM_ho_bezier_basis
(
 const PDM_Mesh_nodal_elt_t  type,
 const int                   order,
 const int                   n_nodes,
 const int                   n_pts,
 const double               *uvw,
 double                     *weights
)
{
 switch (type) {

  case PDM_MESH_NODAL_BAR2:
  case PDM_MESH_NODAL_BARHO:
    _basis_bezier_edge(order, n_pts, uvw, weights);
    break;

  case PDM_MESH_NODAL_TRIA3:
  case PDM_MESH_NODAL_TRIAHO:
    _basis_bezier_tria(order, n_pts, uvw, weights);
    break;
  default:
    break;
  }
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
