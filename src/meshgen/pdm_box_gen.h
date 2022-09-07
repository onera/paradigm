#ifndef __PDM_BOX_GEN_H__
#define __PDM_BOX_GEN_H__

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2019       ONERA

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*----------------------------------------------------------------------------*/

#ifdef  __cplusplus
extern "C" {
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

void
PDM_box_gen_cartesian
(
  int           n_vtx_x,
  int           n_vtx_y,
  int           n_vtx_z,
  double        length,
  int          *n_box_out,
  double      **box_coord_out,
  PDM_g_num_t **box_gnum_out
);


/**
 *
 * \brief Generate a random set of boxes
 *
 * \param [in]   comm                   MPI Communicator id
 * \param [in]   seed                   Random seed
 * \param [in]   geometric_g_num        Compute global ids from coordinates
 * \param [in]   gn_box                 Global number of boxes
 * \param [in]   min_size               Minimal box size
 * \param [in]   max_size               Maximal box size
 * \param [in]   x_min                  Minimal X-coordinate for box centers
 * \param [in]   y_min                  Minimal Y-coordinate for box centers
 * \param [in]   z_min                  Minimal Z-coordinate for box centers
 * \param [in]   x_max                  Maximal X-coordinate for box centers
 * \param [in]   y_max                  Maximal Y-coordinate for box centers
 * \param [in]   z_max                  Maximal Z-coordinate for box centers
 * \param [out]  n_box                  Local number of boxes
 * \param [out]  box_extents            Extents of the local boxes
 * \param [out]  box_ln_to_gn           Global ids of the local boxes
 *
 */

void
PDM_box_gen_random
(
 PDM_MPI_Comm   comm,
 int            seed,
 int            geometric_g_num,
 PDM_g_num_t    gn_box,
 double         min_size,
 double         max_size,
 double         x_min,
 double         y_min,
 double         z_min,
 double         x_max,
 double         y_max,
 double         z_max,
 int           *n_box,
 double       **box_extents,
 PDM_g_num_t  **box_ln_to_gn
 );


#ifdef __cplusplus
}
#endif

#endif // PDM_BOX_GEN_H
