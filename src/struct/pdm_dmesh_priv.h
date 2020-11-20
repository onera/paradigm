#ifndef __PDM_DMESH_PRIV_H__
#define __PDM_DMESH_PRIV_H__

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2017       ONERA

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
#include "pdm_mpi.h"

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

/**
 * \struct _pdm_dmesh_t
 * \brief  Define a distributed mesh. Arrays are shared: this structure does not
 *         holds the data.
 *
 */

struct _pdm_dmesh_t
{
  int               dn_cell;          /*!< Number of distributed cells         */
  int               dn_face;          /*!< Number of distributed faces         */
  int               dn_vtx;           /*!< Number of distributed vertices      */
  int               n_bnd;            /*!< Number of boundaries                */
  int               n_join;           /*!< Number of interfaces with other zone*/
  const PDM_g_num_t *_dface_cell;     /*!< Face-cell connectivity of distributed
                                        faces (size = 2 * dn_face)
                                        if iface is a boundary face,
                                        _dface_cell[2*iface + 1] = 0           */
  const int         *_dface_vtx_idx;   /*!< Face-vertex connectivity index of
                                        distributed faces (size = dn_face + 1) */
  const PDM_g_num_t *_dface_vtx;      /*!< Face-vertex connectivity of
                                        distributed faces (size = dface_vtx_idx[
                                        dn_face])                              */
  const double      *_dvtx_coord;     /*!< Coordinates of ditributed vertices
                                        (size = 3 * dn_vtx)                    */
  const int         *_dface_bound_idx; /*!< Index of distributed faces list of
                                        each boundary (size = n_bnd + 1)
                                        or NULL                               */
  const PDM_g_num_t *_dface_bound;    /*!< Distributed faces list of each
                                       boundary (size = dface_bound_idx[n_bnd])
                                        or NULL                               */
  const int         *_joins_glob_id;  /*!< Global id of each joi (size=n_join)
                                           or NULL. Same data for all procs   */
  const int         *_dface_join_idx;  /*!< Index of distributed faces list of
                                        each join (size = n_join + 1)
                                        or NULL                               */
  const PDM_g_num_t *_dface_join;     /*!< Distributed faces list of each
                                       join (size = dface_join_idx[n_join])
                                        or NULL                               */
};


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DMESH_H__ */
