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
#include "pdm_error.h"
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


static void
_check_is_not_part
(
 PDM_isosurface_t *isos
)
{
  if (isos->is_dist_or_part==-1) {
    isos->is_dist_or_part=0;
  }
  else if (isos->is_dist_or_part==1) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t already set as partitioned.\n");
  }
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/


void
PDM_isosurface_dconnectivity_set
(
 PDM_isosurface_t        *isos,
 PDM_connectivity_type_t  connectivity_type,
 int                     *dconnect_idx,
 PDM_g_num_t             *dconnect
)
{
  _check_is_not_part(isos);

  switch (connectivity_type) {
    case PDM_CONNECTIVITY_TYPE_CELL_FACE:
      isos->dcell_face     = dconnect;
      isos->dcell_face_idx = dconnect_idx;
      break;
    case PDM_CONNECTIVITY_TYPE_FACE_EDGE:
      isos->dface_edge     = dconnect;
      isos->dface_edge_idx = dconnect_idx;
      break;
    case PDM_CONNECTIVITY_TYPE_FACE_VTX:
      isos->dface_vtx     = dconnect;
      isos->dface_vtx_idx = dconnect_idx;
      break;
    case PDM_CONNECTIVITY_TYPE_EDGE_VTX:
      isos->dedge_vtx     = dconnect;
      break;
    default:
      PDM_error(__FILE__, __LINE__, 0, "invalid connectivity_type for iso_surface %d\n", connectivity_type);
      break;
  }
}


void
PDM_isosurface_dvtx_coord_set
(
 PDM_isosurface_t *isos,
 double           *dvtx_coord
)
{
  _check_is_not_part(isos);
  
  isos->dvtx_coord = dvtx_coord;
}


void
PDM_isosurface_distrib_set
(
 PDM_isosurface_t    *isos,
 PDM_mesh_entities_t  entity_type,
 PDM_g_num_t         *distrib
)
{
  _check_is_not_part(isos);

  switch (entity_type) {
    case PDM_MESH_ENTITY_CELL:
      isos->distrib_cell = distrib;
      break;
    case PDM_MESH_ENTITY_FACE:
      isos->distrib_face = distrib;
      break;
    case PDM_MESH_ENTITY_EDGE:
      isos->distrib_edge = distrib;
      break;
    case PDM_MESH_ENTITY_VTX:
      isos->distrib_vtx  = distrib;
      break;
    default:
      PDM_error(__FILE__, __LINE__, 0, "invalid entity_type for iso_surface %d\n", entity_type);
      break;
  }
}


void
PDM_isosurface_dgroup_set
(
 PDM_isosurface_t    *isos,
 PDM_mesh_entities_t  entity_type,
 int                  n_group,
 int                 *dgroup_entity_idx,
 PDM_g_num_t         *dgroup_entity
)
{
  _check_is_not_part(isos);
  
  switch (entity_type) {
    case PDM_MESH_ENTITY_FACE:
      isos->n_group_face    = n_group;
      isos->dgroup_face_idx = dgroup_entity_idx;
      isos->dgroup_face     = dgroup_entity;
      break;
    case PDM_MESH_ENTITY_EDGE:
      isos->n_group_edge    = n_group;
      isos->dgroup_edge_idx = dgroup_entity_idx;
      isos->dgroup_edge     = dgroup_entity;
      break;
    case PDM_MESH_ENTITY_VTX:
      isos->n_group_vtx     = n_group;
      isos->dgroup_vtx_idx  = dgroup_entity_idx;
      isos->dgroup_vtx      = dgroup_entity;
      break;
    default:
      PDM_error(__FILE__, __LINE__, 0, "invalid entity_type for iso_surface boundary %d\n", entity_type);
      break;
  }
}


void
PDM_isosurface_dmesh_set
(
 PDM_isosurface_t *isos,
 PDM_dmesh_t      *dmesh
)
{
  _check_is_not_part(isos);
}

void
PDM_isosurface_dmesh_nodal_set
(
 PDM_isosurface_t  *isos,
 PDM_dmesh_nodal_t *dmn
)
{
  _check_is_not_part(isos);
}


void
PDM_isosurface_dfield_set
(
 PDM_isosurface_t *isos,
 int               id_isosurface,
 double           *dfield
)
{
  _check_is_not_part(isos);

  isos->dfield[id_isosurface] = dfield;
}


void
PDM_isosurface_dgradient_set
(
 PDM_isosurface_t *isos,
 int               id_isosurface,
 double           *dgradient
)
{
  _check_is_not_part(isos);
 
  isos->dgradient[id_isosurface] = dgradient;
}

int
PDM_isosurface_dconnectivity_get
(
 PDM_isosurface_t         *isos,
 int                       id_isosurface,
 PDM_connectivity_type_t   connectivity_type,
 int                     **dconnect_idx,
 PDM_g_num_t             **dconnect,
 PDM_ownership_t           ownership
)
{
  _check_is_not_part(isos);

  return 0;
}


int
PDM_isosurface_dvtx_coord_get
(
 PDM_isosurface_t  *isos,
 int                id_isosurface,
 double           **dvtx_coord,
 PDM_ownership_t    ownership
)
{
  _check_is_not_part(isos);

  return 0;
}


int
PDM_isosurface_dgroup_get
(
 PDM_isosurface_t     *isos,
 int                   id_isosurface,
 PDM_mesh_entities_t   entity_type,
 int                  *n_group,
 int                 **dgroup_entity_idx,
 PDM_g_num_t         **dgroup_entity,
 PDM_ownership_t       ownership
)
{
  _check_is_not_part(isos);

  return 0;
}


#ifdef  __cplusplus
}
#endif