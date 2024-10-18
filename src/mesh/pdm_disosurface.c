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
#include <stdlib.h>

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_error.h"
#include "pdm_distrib.h"
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
 * Macro definitions
 *============================================================================*/

#define CHECK_IS_NOT_PART(isos) \
  if ((isos)->entry_is_part==-1) { \
    (isos)->entry_is_part=0; \
  } \
  else if ((isos)->entry_is_part==1) { \
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t already set as partitioned.\n"); \
  }

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
  /*
   * TODO: transform connectivity as dmesh or not ? 
   */

  CHECK_IS_NOT_PART(isos);
  PDM_ISOSURFACE_CHECK_ENTRY_MESH_COHERENCE(isos, 1);

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
      PDM_error(__FILE__, __LINE__, 0, "invalid connectivity_type (%d) for isosurface.\n", connectivity_type);
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
  CHECK_IS_NOT_PART(isos);
  PDM_ISOSURFACE_CHECK_ENTRY_MESH_COHERENCE(isos, 1);
  
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
  CHECK_IS_NOT_PART(isos);
  PDM_ISOSURFACE_CHECK_ENTRY_MESH_COHERENCE(isos, 1);

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
      PDM_error(__FILE__, __LINE__, 0, "invalid entity_type (%d) for isosurface.\n", entity_type);
      break;
  }
}


void
PDM_isosurface_dgroup_set
(
 PDM_isosurface_t    *isos,
 PDM_mesh_entities_t  entity_type,
 int                 *dgroup_entity_idx,
 PDM_g_num_t         *dgroup_entity
)
{
  CHECK_IS_NOT_PART(isos);
  PDM_ISOSURFACE_CHECK_ENTRY_MESH_COHERENCE(isos, 1);
  
  switch (entity_type) {
    case PDM_MESH_ENTITY_FACE:
      isos->dgroup_face_idx = dgroup_entity_idx;
      isos->dgroup_face     = dgroup_entity;
      break;
    // case PDM_MESH_ENTITY_EDGE:
    //   isos->dgroup_edge_idx = dgroup_entity_idx;
    //   isos->dgroup_edge     = dgroup_entity;
    //   break;
    // case PDM_MESH_ENTITY_VTX:
    //   isos->dgroup_vtx_idx  = dgroup_entity_idx;
    //   isos->dgroup_vtx      = dgroup_entity;
    //   break;
    default:
      PDM_error(__FILE__, __LINE__, 0, "invalid entity_type (%d) for isosurface boundary.\n", entity_type);
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
  CHECK_IS_NOT_PART(isos);
  PDM_ISOSURFACE_CHECK_ENTRY_MESH_COHERENCE(isos, 2);
 
  /**
   * TODO: add check on entry dmesh to be sure that all required data are in. 
   */

  isos->dmesh = dmesh;

  /* Unpack dmesh */
  isos->entry_mesh_type = 1; // héhé
  for (int i_entity = PDM_MESH_ENTITY_CELL; i_entity < PDM_MESH_ENTITY_MAX; i_entity++) {
    PDM_g_num_t *distrib = NULL;
    PDM_dmesh_distrib_get     (dmesh, i_entity, &distrib);
    PDM_isosurface_distrib_set(isos,  i_entity,  distrib);
  }

  // Cells
  int         *dcell_face_idx = NULL;
  PDM_g_num_t *dcell_face     = NULL;
  PDM_dmesh_connectivity_get(dmesh,
                             PDM_CONNECTIVITY_TYPE_CELL_FACE,
                             &dcell_face,
                             &dcell_face_idx,
                             PDM_OWNERSHIP_BAD_VALUE);
  PDM_isosurface_dconnectivity_set(isos,
                                   PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                   dcell_face_idx,
                                   dcell_face);

  // Faces
  int         *dface_edge_idx = NULL;
  PDM_g_num_t *dface_edge     = NULL;
  PDM_dmesh_connectivity_get(dmesh,
                             PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                             &dface_edge,
                             &dface_edge_idx,
                             PDM_OWNERSHIP_BAD_VALUE);
  PDM_isosurface_dconnectivity_set(isos,
                                   PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                   dface_edge_idx,
                                   dface_edge);

  int         *dface_vtx_idx = NULL;
  PDM_g_num_t *dface_vtx     = NULL;
  PDM_dmesh_connectivity_get(dmesh,
                             PDM_CONNECTIVITY_TYPE_FACE_VTX,
                             &dface_vtx,
                             &dface_vtx_idx,
                             PDM_OWNERSHIP_BAD_VALUE);
  PDM_isosurface_dconnectivity_set(isos,
                                   PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                   dface_vtx_idx,
                                   dface_vtx);

  // Edges
  int         *dedge_vtx_idx = NULL;
  PDM_g_num_t *dedge_vtx     = NULL;
  PDM_dmesh_connectivity_get(dmesh,
                             PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                             &dedge_vtx,
                             &dedge_vtx_idx,
                             PDM_OWNERSHIP_BAD_VALUE);

  PDM_isosurface_dconnectivity_set(isos,
                                   PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                   NULL,
                                   dedge_vtx);

  // Vertices
  double *dvtx_coord = NULL;
  PDM_dmesh_vtx_coord_get(dmesh, &dvtx_coord, PDM_OWNERSHIP_BAD_VALUE);
  PDM_isosurface_dvtx_coord_set(isos, dvtx_coord);

  // Groups (surfaces)
  int         *dsurface_face_idx = NULL;
  PDM_g_num_t *dsurface_face     = NULL;
  int n_surface = PDM_dmesh_bound_get(dmesh,
                                      PDM_BOUND_TYPE_FACE,
                                      &dsurface_face,
                                      &dsurface_face_idx,
                                      PDM_OWNERSHIP_BAD_VALUE);

  PDM_isosurface_n_group_set(isos,
                             PDM_MESH_ENTITY_FACE,
                             n_surface);

  PDM_isosurface_dgroup_set(isos,
                            PDM_MESH_ENTITY_FACE,
                            dsurface_face_idx,
                            dsurface_face);
}

void
PDM_isosurface_dmesh_nodal_set
(
 PDM_isosurface_t  *isos,
 PDM_dmesh_nodal_t *dmn
)
{
  CHECK_IS_NOT_PART(isos);
  PDM_ISOSURFACE_CHECK_ENTRY_MESH_COHERENCE(isos, 3);
  
  /**
   * TODO: add check on entry dmesh_nodal to be sure that all required data are in. 
   */

  isos->dmesh_nodal = dmn;
}


void
PDM_isosurface_dfield_set
(
 PDM_isosurface_t *isos,
 int               id_isosurface,
 double           *dfield
)
{
  CHECK_IS_NOT_PART(isos);

  PDM_ISOSURFACE_CHECK_ID(isos, id_isosurface);

  _isosurface_t *_iso = &isos->isosurfaces[id_isosurface];
  _iso->dfield = dfield;
}


int
PDM_isosurface_dconnectivity_get
(
  PDM_isosurface_t         *isos,
  int                       id_iso,
  PDM_connectivity_type_t   connectivity_type,
  int                     **dconnect_idx,
  PDM_g_num_t             **dconnect,
  PDM_ownership_t           ownership
)
{
  CHECK_IS_NOT_PART(isos);

  PDM_ISOSURFACE_CHECK_ID      (isos, id_iso);
  PDM_ISOSURFACE_CHECK_COMPUTED(isos, id_iso);

  _isosurface_t *_iso = &isos->isosurfaces[id_iso];

  int n_entity = 0;

  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    _iso->iso_owner_dconnec[connectivity_type] = ownership;
  }

  if (connectivity_type==PDM_CONNECTIVITY_TYPE_EDGE_VTX) {
    n_entity      = _iso->iso_dn_entity[PDM_MESH_ENTITY_EDGE];
    *dconnect_idx = NULL;
    *dconnect     = _iso->iso_dconnec[connectivity_type];
  }
  else if (connectivity_type==PDM_CONNECTIVITY_TYPE_FACE_VTX) {
    n_entity      = _iso->iso_dn_entity[PDM_MESH_ENTITY_FACE];
    *dconnect_idx = _iso->iso_dconnec_idx[connectivity_type];
    *dconnect     = _iso->iso_dconnec    [connectivity_type];
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: has no connectivity of type %d.\n", connectivity_type);
  }

  return n_entity;
}


int
PDM_isosurface_disovalue_entity_get
(
  PDM_isosurface_t     *isos,
  int                   id_isosurface,
  PDM_mesh_entities_t   entity_type,
  int                 **disovalue_entity_idx,
  PDM_g_num_t         **disovalue_entity,
  PDM_ownership_t       ownership
)
{
  CHECK_IS_NOT_PART(isos);

  PDM_ISOSURFACE_CHECK_ID      (isos, id_isosurface);
  PDM_ISOSURFACE_CHECK_COMPUTED(isos, id_isosurface);

  _isosurface_t *_iso = &isos->isosurfaces[id_isosurface];

  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    _iso->iso_owner_disovalue_entity[entity_type] = ownership;
  }

  *disovalue_entity_idx = _iso->disovalue_entity_idx[entity_type];
  *disovalue_entity     = _iso->disovalue_entity    [entity_type];

  return _iso->n_isovalues;
}


int
PDM_isosurface_dvtx_parent_weight_get
(
  PDM_isosurface_t     *isos,
  int                   id_iso,
  int                 **dvtx_parent_idx,
  double              **dvtx_parent_weight,
  PDM_ownership_t       ownership
)
{
  CHECK_IS_NOT_PART(isos);

  PDM_ISOSURFACE_CHECK_ID      (isos, id_iso);
  PDM_ISOSURFACE_CHECK_COMPUTED(isos, id_iso);

  _isosurface_t *_iso = &isos->isosurfaces[id_iso];

  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    _iso->iso_owner_dparent_idx[PDM_MESH_ENTITY_VTX] = ownership;
    _iso->iso_owner_dvtx_parent_weight               = ownership;
  }

  *dvtx_parent_idx    = _iso->iso_dentity_parent_idx[PDM_MESH_ENTITY_VTX];
  *dvtx_parent_weight = _iso->iso_dvtx_parent_weight                     ;

  return _iso->iso_dn_entity[PDM_MESH_ENTITY_VTX];
}


int
PDM_isosurface_dvtx_coord_get
(
  PDM_isosurface_t  *isos,
  int                id_iso,
  double           **dvtx_coord,
  PDM_ownership_t    ownership
)
{
  CHECK_IS_NOT_PART(isos);

  PDM_ISOSURFACE_CHECK_ID      (isos, id_iso);
  PDM_ISOSURFACE_CHECK_COMPUTED(isos, id_iso);

  _isosurface_t *_iso = &isos->isosurfaces[id_iso];

  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    _iso->iso_owner_dvtx_coord = ownership;
  }

  *dvtx_coord = _iso->iso_dvtx_coord;

  return _iso->iso_dn_entity[PDM_MESH_ENTITY_VTX];
}


int
PDM_isosurface_distrib_get
(
  PDM_isosurface_t     *isos,
  int                   id_iso,
  PDM_mesh_entities_t   entity_type,
  PDM_g_num_t         **distribution
)
{
  CHECK_IS_NOT_PART(isos);

  PDM_ISOSURFACE_CHECK_ID      (isos, id_iso);
  PDM_ISOSURFACE_CHECK_COMPUTED(isos, id_iso);

  PDM_ISOSURFACE_CHECK_ENTITY_TYPE(entity_type);

  _isosurface_t *_iso = &isos->isosurfaces[id_iso];
  *distribution = PDM_compute_entity_distribution(isos->comm, _iso->iso_dn_entity[entity_type]);

  int i_rank;
  PDM_MPI_Comm_rank(isos->comm, &i_rank);
  return distribution[i_rank+1]-distribution[i_rank];
}


int
PDM_isosurface_dgroup_get
(
  PDM_isosurface_t     *isos,
  int                   id_iso,
  PDM_mesh_entities_t   entity_type,
  int                 **dgroup_entity_idx,
  PDM_g_num_t         **dgroup_entity,
  PDM_ownership_t       ownership
)
{
  CHECK_IS_NOT_PART(isos);

  PDM_ISOSURFACE_CHECK_ID      (isos, id_iso);
  PDM_ISOSURFACE_CHECK_COMPUTED(isos, id_iso);

  _isosurface_t *_iso = &isos->isosurfaces[id_iso];

  int n_group = 0;

  if (entity_type==PDM_MESH_ENTITY_EDGE) {

    if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
      _iso->iso_owner_dedge_bnd = ownership;
    }

    n_group            = _iso->iso_n_edge_group;
    *dgroup_entity_idx = _iso->iso_dedge_group_idx ;
    *dgroup_entity     = _iso->iso_dedge_group_gnum;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: has no group for entity_type %d.\n", entity_type);
  }

  return n_group;
}


#undef CHECK_IS_NOT_PART


#ifdef  __cplusplus
}
#endif
