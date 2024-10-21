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

#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/
#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_error.h"
#include "pdm_mesh_nodal.h"
#include "pdm_array.h"
#include "pdm_isosurface.h"
#include "pdm_isosurface_priv.h"
#include "pdm_priv.h"

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

#define CHECK_IS_NOT_DIST(isos) \
  if ((isos)->entry_is_part==-1) { \
    (isos)->entry_is_part=1; \
  } \
  else if ((isos)->entry_is_part==0) { \
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t already set as distributed.\n"); \
  }

#define CHECK_I_PART_SET(isos, i_part) \
  if (i_part >= (isos)->n_part) { \
    PDM_error(__FILE__, __LINE__, 0, "Invalid i_part (%d / %d).\n", \
              i_part, (isos)->n_part); \
  }

#define CHECK_I_PART_GET(isos, i_part) \
  if (i_part >= (isos)->iso_n_part) { \
    PDM_error(__FILE__, __LINE__, 0, "Invalid i_part (%d / %d).\n", \
              i_part, (isos)->iso_n_part); \
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
PDM_isosurface_n_part_set
(
 PDM_isosurface_t *isos,
 int               n_part
)
{
  CHECK_IS_NOT_DIST(isos);
  PDM_ISOSURFACE_CHECK_ENTRY_MESH_COHERENCE(isos, -1);

  isos->n_part = n_part;
  isos->n_cell = PDM_array_const_int(n_part, -1);
  isos->n_face = PDM_array_const_int(n_part, -1);
  isos->n_edge = PDM_array_const_int(n_part, -1);
  isos->n_vtx  = PDM_array_const_int(n_part, -1);

  PDM_malloc(isos->cell_face    , n_part, int         *);
  PDM_malloc(isos->cell_face_idx, n_part, int         *);
  PDM_malloc(isos->face_edge    , n_part, int         *);
  PDM_malloc(isos->face_edge_idx, n_part, int         *);
  PDM_malloc(isos->face_vtx     , n_part, int         *);
  PDM_malloc(isos->face_vtx_idx , n_part, int         *);
  PDM_malloc(isos->edge_vtx     , n_part, int         *);
  PDM_malloc(isos->vtx_coord    , n_part, double      *);
  PDM_malloc(isos->cell_gnum    , n_part, PDM_g_num_t *);
  PDM_malloc(isos->face_gnum    , n_part, PDM_g_num_t *);
  PDM_malloc(isos->edge_gnum    , n_part, PDM_g_num_t *);
  PDM_malloc(isos->vtx_gnum     , n_part, PDM_g_num_t *);

  PDM_malloc(isos->group_face_idx , n_part, int         *);
  PDM_malloc(isos->group_face     , n_part, int         *);
  PDM_malloc(isos->group_face_gnum, n_part, PDM_g_num_t *);

  for (int i_part = 0; i_part < n_part; i_part++) {
    isos->cell_face    [i_part] = NULL;
    isos->cell_face_idx[i_part] = NULL;
    isos->face_edge    [i_part] = NULL;
    isos->face_edge_idx[i_part] = NULL;
    isos->face_vtx     [i_part] = NULL;
    isos->face_vtx_idx [i_part] = NULL;
    isos->edge_vtx     [i_part] = NULL;
    isos->vtx_coord    [i_part] = NULL;
    isos->cell_gnum    [i_part] = NULL;
    isos->face_gnum    [i_part] = NULL;
    isos->edge_gnum    [i_part] = NULL;
    isos->vtx_gnum     [i_part] = NULL;

    isos->group_face_idx [i_part] = NULL;
    isos->group_face     [i_part] = NULL;
    isos->group_face_gnum[i_part] = NULL;
  }
}

void
PDM_isosurface_connectivity_set
(
 PDM_isosurface_t        *isos,
 int                      i_part,
 PDM_connectivity_type_t  connectivity_type,
 int                      n_entity,
 int                     *connect_idx,
 int                     *connect
)
{
  CHECK_IS_NOT_DIST(isos);
  PDM_ISOSURFACE_CHECK_ENTRY_MESH_COHERENCE(isos, -1);
  CHECK_I_PART_SET(isos, i_part);

  switch (connectivity_type) {
    case PDM_CONNECTIVITY_TYPE_CELL_FACE:
      isos->n_cell       [i_part] = n_entity;
      isos->cell_face    [i_part] = connect;
      isos->cell_face_idx[i_part] = connect_idx;
      break;
    case PDM_CONNECTIVITY_TYPE_FACE_EDGE:
      isos->n_face       [i_part] = n_entity;
      isos->face_edge    [i_part] = connect;
      isos->face_edge_idx[i_part] = connect_idx;
      break;
    case PDM_CONNECTIVITY_TYPE_FACE_VTX:
      isos->n_face       [i_part] = n_entity;
      isos->face_vtx     [i_part] = connect;
      isos->face_vtx_idx [i_part] = connect_idx;
      break;
    case PDM_CONNECTIVITY_TYPE_EDGE_VTX:
      isos->n_edge       [i_part] = n_entity;
      isos->edge_vtx     [i_part] = connect;
      break;
    default:
      PDM_error(__FILE__, __LINE__, 0, "invalid connectivity_type (%d) for isosurface.\n", connectivity_type);
      break;
  }
}


void
PDM_isosurface_vtx_coord_set
(
 PDM_isosurface_t *isos,
 int               i_part,
 int               n_vtx,
 double           *vtx_coord
)
{
  CHECK_IS_NOT_DIST(isos);
  PDM_ISOSURFACE_CHECK_ENTRY_MESH_COHERENCE(isos, -1);
  CHECK_I_PART_SET(isos, i_part);
  
  isos->n_vtx    [i_part] = n_vtx;
  isos->vtx_coord[i_part] = vtx_coord;
}


void
PDM_isosurface_ln_to_gn_set
(
 PDM_isosurface_t    *isos,
 int                  i_part,
 PDM_mesh_entities_t  entity_type,
 PDM_g_num_t         *ln_to_gn
)
{
  CHECK_IS_NOT_DIST(isos);
  PDM_ISOSURFACE_CHECK_ENTRY_MESH_COHERENCE(isos, -1);
  CHECK_I_PART_SET(isos, i_part);
  
  switch (entity_type) {
    case PDM_MESH_ENTITY_CELL:
      isos->cell_gnum[i_part] = ln_to_gn;
      break;
    case PDM_MESH_ENTITY_FACE:
      isos->face_gnum[i_part] = ln_to_gn;
      break;
    case PDM_MESH_ENTITY_EDGE:
      isos->edge_gnum[i_part] = ln_to_gn;
      break;
    case PDM_MESH_ENTITY_VTX:
      isos->vtx_gnum[i_part]  = ln_to_gn;
      break;
    default:
      PDM_error(__FILE__, __LINE__, 0, "invalid entity_type (%d) for isosurface.\n", entity_type);
      break;
  }
}


void
PDM_isosurface_group_set
(
 PDM_isosurface_t    *isos,
 int                  i_part,
 PDM_mesh_entities_t  entity_type,
 int                 *group_entity_idx,
 int                 *group_entity,
 PDM_g_num_t         *group_entity_ln_to_gn
)
{
  CHECK_IS_NOT_DIST(isos);
  PDM_ISOSURFACE_CHECK_ENTRY_MESH_COHERENCE(isos, -1);
  CHECK_I_PART_SET(isos, i_part);

  switch (entity_type) {
    case PDM_MESH_ENTITY_FACE:
      isos->group_face_idx [i_part] = group_entity_idx;
      isos->group_face     [i_part] = group_entity;
      isos->group_face_gnum[i_part] = group_entity_ln_to_gn;
      break;
    // case PDM_MESH_ENTITY_EDGE:
    //   isos->group_edge_idx[i_part] = group_entity_idx;
    //   isos->group_edge    [i_part] = group_entity;
    //   break;
    // case PDM_MESH_ENTITY_VTX:
    //   isos->group_vtx_idx[i_part]  = group_entity_idx;
    //   isos->group_vtx    [i_part]  = group_entity;
    //   break;
    default:
      PDM_error(__FILE__, __LINE__, 0, "invalid entity_type (%d) for isosurface boundary.\n", entity_type);
      break;
  }
}


void
PDM_isosurface_part_mesh_set
(
 PDM_isosurface_t *isos,
 PDM_part_mesh_t  *pmesh
)
{
  CHECK_IS_NOT_DIST(isos);
  PDM_ISOSURFACE_CHECK_ENTRY_MESH_COHERENCE(isos, -2);
  
  /**
   * TODO: add check on entry pmesh to be sure that all required data are in. 
   */

  isos->pmesh = pmesh;

  /* Unpack part_mesh */
  isos->entry_mesh_type = -1; // héhé
  PDM_isosurface_n_part_set(isos, PDM_part_mesh_n_part_get(pmesh));

  isos->n_group_face = PDM_part_mesh_n_bound_get(pmesh, PDM_BOUND_TYPE_FACE);

  for (int i_part = 0; i_part < isos->n_part; i_part++) {
    // Number of entities
    isos->n_cell[i_part] = PDM_part_mesh_n_entity_get(pmesh,
                                                      i_part,
                                                      PDM_MESH_ENTITY_CELL);

    isos->n_face[i_part] = PDM_part_mesh_n_entity_get(pmesh,
                                                      i_part,
                                                      PDM_MESH_ENTITY_FACE);

    isos->n_edge[i_part] = PDM_part_mesh_n_entity_get(pmesh,
                                                      i_part,
                                                      PDM_MESH_ENTITY_EDGE);

    isos->n_vtx [i_part] = PDM_part_mesh_n_entity_get(pmesh,
                                                      i_part,
                                                      PDM_MESH_ENTITY_VTX);

    // Connectivities
    PDM_part_mesh_connectivity_get(pmesh,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                   &isos->cell_face    [i_part],
                                   &isos->cell_face_idx[i_part],
                                   PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_connectivity_get(pmesh,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                   &isos->face_edge    [i_part],
                                   &isos->face_edge_idx[i_part],
                                   PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_connectivity_get(pmesh,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                   &isos->face_vtx    [i_part],
                                   &isos->face_vtx_idx[i_part],
                                   PDM_OWNERSHIP_BAD_VALUE);

    int *edge_vtx_idx = NULL;
    PDM_part_mesh_connectivity_get(pmesh,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                   &isos->edge_vtx[i_part],
                                   &edge_vtx_idx,
                                   PDM_OWNERSHIP_BAD_VALUE);

    // Coordinates
    PDM_part_mesh_vtx_coord_get(pmesh,
                                i_part,
                                &isos->vtx_coord[i_part],
                                PDM_OWNERSHIP_BAD_VALUE);

    // Global IDs
    PDM_part_mesh_entity_ln_to_gn_get(pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_CELL,
                                      &isos->cell_gnum[i_part],
                                      PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_entity_ln_to_gn_get(pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      &isos->face_gnum[i_part],
                                      PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_entity_ln_to_gn_get(pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_EDGE,
                                      &isos->edge_gnum[i_part],
                                      PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_entity_ln_to_gn_get(pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_VTX,
                                      &isos->vtx_gnum[i_part],
                                      PDM_OWNERSHIP_BAD_VALUE);

    // Surfaces
    PDM_part_mesh_bound_concat_get(isos->pmesh,
                                   i_part,
                                   PDM_BOUND_TYPE_FACE,
                                   &isos->group_face_idx [i_part],
                                   &isos->group_face     [i_part],
                                   &isos->group_face_gnum[i_part],
                                   PDM_OWNERSHIP_BAD_VALUE);
  }
}


void
PDM_isosurface_mesh_nodal_set
(
 PDM_isosurface_t      *isos,
 PDM_part_mesh_nodal_t *pmn
)
{
  CHECK_IS_NOT_DIST(isos);
  PDM_ISOSURFACE_CHECK_ENTRY_MESH_COHERENCE(isos, -3);
  
  /**
   * TODO: add check on entry pmesh_nodal to be sure that all required data are in. 
   */

  isos->pmesh_nodal = pmn;
  isos->n_part      = PDM_part_mesh_nodal_n_part_get(pmn);
}


void
PDM_isosurface_redistribution_set
(
 PDM_isosurface_t        *isos,
 PDM_extract_part_kind_t  extract_kind,
 PDM_split_dual_t         part_method
)
{
  CHECK_IS_NOT_DIST(isos);
 
  isos->extract_kind = extract_kind;
  isos->part_method  = part_method;
}


void
PDM_isosurface_field_set
(
 PDM_isosurface_t *isos,
 int               id_isosurface,
 int               i_part,
 double           *field
)
{
  CHECK_IS_NOT_DIST(isos);
  CHECK_I_PART_SET(isos, i_part);

  PDM_ISOSURFACE_CHECK_ID(isos, id_isosurface);

  // > Check i_part is valid
  int n_part = isos->n_part;
  if (n_part==-1) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: mesh seems not to be defined.\n", isos->entry_mesh_type);
  }

  if (i_part>=n_part) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: trying to defined field for i_part >= n_part (%d >= %d).\n", i_part, n_part);
  }
 
  _isosurface_t *_iso = &isos->isosurfaces[id_isosurface];
  if (_iso->field==NULL) {
    PDM_malloc(_iso->field, isos->n_part, double *);
  }
  _iso->field[i_part] = field;
}


int
PDM_isosurface_local_parent_get
(
 PDM_isosurface_t     *isos,
 int                   id_isosurface,
 int                   i_part,
 PDM_mesh_entities_t   entity_type,
 int                 **entity_parent_idx,
 int                 **entity_parent,
 PDM_ownership_t       ownership
)
{
  // TODO: how to do in parallel for iso entities on entities between partition
  CHECK_IS_NOT_DIST(isos);

  PDM_ISOSURFACE_CHECK_ID      (isos, id_isosurface);
  PDM_ISOSURFACE_CHECK_COMPUTED(isos, id_isosurface);

  CHECK_I_PART_GET(isos, i_part);

  PDM_ISOSURFACE_CHECK_ENTITY_TYPE(entity_type);

  if (isos->extract_kind!=PDM_EXTRACT_PART_KIND_LOCAL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: extract_kind is not PDM_EXTRACT_PART_KIND_LOCAL.\n");
  }

  _isosurface_t *_iso = &isos->isosurfaces[id_isosurface];

  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    _iso->iso_owner_parent_lnum[entity_type][i_part] = ownership;
    _iso->iso_owner_parent_idx [entity_type][i_part] = ownership;
  }

  *entity_parent_idx = _iso->iso_entity_parent_idx [entity_type][i_part];
  *entity_parent     = _iso->iso_entity_parent_lnum[entity_type][i_part];

  return _iso->iso_n_entity[entity_type][i_part];
}


int
PDM_isosurface_parent_weight_get
(
  PDM_isosurface_t     *isos,
  int                   id_isosurface,
  int                   i_part,
  PDM_mesh_entities_t   entity_type,
  int                 **parent_idx,
  double              **parent_weight,
  PDM_ownership_t       ownership
)
{
  CHECK_IS_NOT_DIST(isos);

  PDM_ISOSURFACE_CHECK_ID      (isos, id_isosurface);
  PDM_ISOSURFACE_CHECK_COMPUTED(isos, id_isosurface);

  CHECK_I_PART_GET(isos, i_part);

  _isosurface_t *_iso = &isos->isosurfaces[id_isosurface];

  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    _iso->iso_owner_parent_idx [entity_type][i_part] = ownership;
    _iso->iso_owner_parent_wght[entity_type][i_part] = ownership;
  }

  *parent_idx    = _iso->iso_entity_parent_idx [entity_type][i_part];
  *parent_weight = _iso->iso_entity_parent_wght[entity_type][i_part];

  return _iso->iso_n_entity[entity_type][i_part];
}


int
PDM_isosurface_connectivity_get
(
 PDM_isosurface_t         *isos,
 int                       id_isosurface,
 int                       i_part,
 PDM_connectivity_type_t   connectivity_type,
 int                     **connect_idx,
 int                     **connect,
 PDM_ownership_t           ownership
)
{
  CHECK_IS_NOT_DIST(isos);

  PDM_ISOSURFACE_CHECK_ID      (isos, id_isosurface);
  PDM_ISOSURFACE_CHECK_COMPUTED(isos, id_isosurface);

  CHECK_I_PART_GET(isos, i_part);

  if (connectivity_type != PDM_CONNECTIVITY_TYPE_EDGE_VTX &&
      connectivity_type != PDM_CONNECTIVITY_TYPE_FACE_VTX) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: has no connectivity of type %d.\n", connectivity_type);
  }

  _isosurface_t *_iso = &isos->isosurfaces[id_isosurface];

  int n_entity = 0;

  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    _iso->iso_owner_connec[connectivity_type][i_part] = ownership;
  }

  *connect = _iso->iso_connec[connectivity_type][i_part];

  if (connectivity_type==PDM_CONNECTIVITY_TYPE_EDGE_VTX) {
    n_entity = _iso->iso_n_entity[PDM_MESH_ENTITY_EDGE][i_part];
  }
  else if (connectivity_type==PDM_CONNECTIVITY_TYPE_FACE_VTX) {
    n_entity     = _iso->iso_n_entity  [PDM_MESH_ENTITY_FACE][i_part];
    *connect_idx = _iso->iso_connec_idx[connectivity_type   ][i_part];
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: has no connectivity of type %d.\n", connectivity_type);
  }

  return n_entity;
}


int
PDM_isosurface_vtx_coord_get
(
 PDM_isosurface_t  *isos,
 int                id_isosurface,
 int                i_part,
 double           **vtx_coord,
 PDM_ownership_t    ownership
)
{
  CHECK_IS_NOT_DIST(isos);

  PDM_ISOSURFACE_CHECK_ID      (isos, id_isosurface);
  PDM_ISOSURFACE_CHECK_COMPUTED(isos, id_isosurface);

  CHECK_I_PART_GET(isos, i_part);

  _isosurface_t *_iso = &isos->isosurfaces[id_isosurface];

  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    _iso->iso_owner_vtx_coord[i_part] = ownership;
  }

  *vtx_coord = _iso->iso_vtx_coord[i_part];

  return _iso->iso_n_entity[PDM_MESH_ENTITY_VTX][i_part];
}


int
PDM_isosurface_ln_to_gn_get
(
 PDM_isosurface_t     *isos,
 int                   id_isosurface,
 int                   i_part,
 PDM_mesh_entities_t   entity_type,
 PDM_g_num_t         **ln_to_gn,
 PDM_ownership_t       ownership
)
{
  CHECK_IS_NOT_DIST(isos);

  PDM_ISOSURFACE_CHECK_ID      (isos, id_isosurface);
  PDM_ISOSURFACE_CHECK_COMPUTED(isos, id_isosurface);

  CHECK_I_PART_GET(isos, i_part);

  PDM_ISOSURFACE_CHECK_ENTITY_TYPE(entity_type);

  _isosurface_t *_iso = &isos->isosurfaces[id_isosurface];

  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    _iso->iso_owner_gnum[entity_type][i_part] = ownership;
  }

  *ln_to_gn = _iso->iso_entity_gnum[entity_type][i_part];

  return _iso->iso_n_entity[entity_type][i_part];
}


int
PDM_isosurface_group_get
(
  PDM_isosurface_t     *isos,
  int                   id_isosurface,
  int                   i_part,
  PDM_mesh_entities_t   entity_type,
  int                 **group_entity_idx,
  int                 **group_entity,
  PDM_g_num_t         **group_entity_gnum,
  PDM_ownership_t       ownership
)
{
  CHECK_IS_NOT_DIST(isos);

  PDM_ISOSURFACE_CHECK_ID      (isos, id_isosurface);
  PDM_ISOSURFACE_CHECK_COMPUTED(isos, id_isosurface);

  CHECK_I_PART_GET(isos, i_part);

  _isosurface_t *_iso = &isos->isosurfaces[id_isosurface];
  
  int n_group = 0;

  if (entity_type==PDM_MESH_ENTITY_EDGE) {
    if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
      _iso->iso_owner_edge_bnd[i_part] = ownership;
    }

     n_group           = _iso->iso_n_edge_group   ;
    *group_entity_idx  = _iso->iso_edge_group_idx [i_part];
    *group_entity      = _iso->iso_edge_group_lnum[i_part];
    *group_entity_gnum = _iso->iso_edge_group_gnum[i_part];
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: has no bounds for entity %d.\n",entity_type);
  }

  return n_group;
}


int
PDM_isosurface_isovalue_entity_idx_get
(
  PDM_isosurface_t     *isos,
  int                   id_isosurface,
  int                   i_part,
  PDM_mesh_entities_t   entity_type,
  int                 **isovalue_entity_idx,
  PDM_ownership_t       ownership
)
{
  CHECK_IS_NOT_DIST(isos);
  PDM_ISOSURFACE_CHECK_ID      (isos, id_isosurface);
  PDM_ISOSURFACE_CHECK_COMPUTED(isos, id_isosurface);
  CHECK_I_PART_GET(isos, i_part);

  PDM_ISOSURFACE_CHECK_ENTITY_TYPE(entity_type);

  _isosurface_t *_iso = &isos->isosurfaces[id_isosurface];

  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    _iso->iso_owner_isovalue_entity_idx[entity_type][i_part] = ownership;
  }

  *isovalue_entity_idx = _iso->isovalue_entity_idx[entity_type][i_part];

  return _iso->n_isovalues;
}


void
PDM_isosurface_n_part_out_set
(
  PDM_isosurface_t *isos,
  int               n_part_out
)
{
  CHECK_IS_NOT_DIST(isos);

  isos->iso_n_part = n_part_out;
}

#undef CHECK_IS_NOT_DIST

#ifdef  __cplusplus
}
#endif
