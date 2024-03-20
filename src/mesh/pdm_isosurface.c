/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2023       ONERA

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
#include "pdm_mesh_nodal.h"
#include "pdm_isosurface.h"
#include "pdm_isosurface_priv.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

PDM_isosurface_t *
PDM_isosurface_create
(
 PDM_MPI_Comm             comm,
 int                      mesh_dimension,
 PDM_Mesh_nodal_elt_t     elt_type
)
{
  PDM_isosurface_t *isos = (PDM_isosurface_t *) malloc(sizeof(PDM_isosurface_t));
  isos->is_dist_or_part  = -1; 
  return isos;
}


int
PDM_isosurface_add
(
 PDM_isosurface_t       *isos,
 PDM_iso_surface_kind_t  kind,
 int                     n_isovalues,
 double                 *isovalues
)
{

}


void
PDM_isosurface_equation_set
(
 PDM_isosurface_t *isos,
 int               id_isosurface,
 double           *coeff,
 int               use_gradient
)
{

}


void
PDM_isosurface_field_function_set
(
 PDM_isosurface_t                *isos,
 int                              id_isosurface,
 PDM_isosurface_field_function_t  func
)
{

}


void
PDM_isosurface_reset
(
 PDM_isosurface_t *isos,
 int               id_isosurface
)
{

}


void
PDM_isosurface_compute
(
 PDM_isosurface_t *isos,
 int               id_isosurface
)
{

}


void
PDM_isosurface_dump_times
(
 PDM_isosurface_t *isos
)
{

}


void
PDM_isosurface_free
(
  PDM_isosurface_t  *isos
)
{

}


#ifdef  __cplusplus
}
#endif