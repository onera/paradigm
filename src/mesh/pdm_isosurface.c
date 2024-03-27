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


void
_check_entry_mesh_coherence
(
 PDM_isosurface_t *isos,
 int               entry_mesh_type
)
{
  if (isos->entry_mesh_type==0) {
    isos->entry_mesh_type=entry_mesh_type;
  }
  else if (isos->entry_mesh_type!=entry_mesh_type) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t:entry_mesh_type already set to %d.\n", isos->entry_mesh_type);
  }
}


PDM_isosurface_t *
PDM_isosurface_create
(
 PDM_MPI_Comm             comm,
 int                      mesh_dimension,
 PDM_Mesh_nodal_elt_t     elt_type
)
{
  PDM_isosurface_t *isos = (PDM_isosurface_t *) malloc(sizeof(PDM_isosurface_t));

  // > Entry mesh information
  isos->is_dist_or_part  = -1; 
  isos->entry_mesh_type  =  0; 
  isos->entry_mesh_dim   =  mesh_dimension;

  // > Isosurface mesh information
  isos->iso_elt_type     =  elt_type; 

  // > Isovalues
  isos->n_isovalues = 0;

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
  isos->kind = kind;
  isos->n_isovalues = n_isovalues;
  isos->isovalues = malloc(sizeof(double) * n_isovalues);
  for (int i=0; i<n_isovalues; ++i) {
    isos->isovalues[i] = isovalues[i];
  }
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
  if (isos->n_isovalues>0) {
    free(isos->isovalues);
  }

  free(isos);
}


#ifdef  __cplusplus
}
#endif