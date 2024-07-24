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
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/


static void
_check_is_not_dist
(
 PDM_isosurface_t *isos
)
{
  if (isos->is_dist_or_part==-1) {
    isos->is_dist_or_part=1;
  }
  else if (isos->is_dist_or_part==0) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t already set as distributed.\n");
  }
}


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
  _check_is_not_dist(isos);
  _check_entry_mesh_coherence(isos, -1);

  isos->n_part     = n_part;
  isos->iso_n_part = n_part;

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

  PDM_malloc(isos->n_group_face   , n_part, int          );
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

    isos->n_group_face   [i_part] = 0;
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
 int                     *connect_idx,
 int                     *connect
)
{
  /*
   * TODO: transform connectivity as pmesh or not ? 
   */

  _check_is_not_dist(isos);
  _check_entry_mesh_coherence(isos, -1);

  switch (connectivity_type) {
    case PDM_CONNECTIVITY_TYPE_CELL_FACE:
      isos->cell_face    [i_part] = connect;
      isos->cell_face_idx[i_part] = connect_idx;
      break;
    case PDM_CONNECTIVITY_TYPE_FACE_EDGE:
      isos->face_edge    [i_part] = connect;
      isos->face_edge_idx[i_part] = connect_idx;
      break;
    case PDM_CONNECTIVITY_TYPE_FACE_VTX:
      isos->face_vtx     [i_part] = connect;
      isos->face_vtx_idx [i_part] = connect_idx;
      break;
    case PDM_CONNECTIVITY_TYPE_EDGE_VTX:
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
 double           *vtx_coord
)
{
  _check_is_not_dist(isos);
  _check_entry_mesh_coherence(isos, -1);
  
  isos->vtx_coord[i_part] = vtx_coord;
}


void
PDM_isosurface_ln_to_gn_set
(
 PDM_isosurface_t    *isos,
 int                  i_part,
 PDM_mesh_entities_t  entity_type,
 int                  n_entity,
 PDM_g_num_t         *ln_to_gn
)
{
  _check_is_not_dist(isos);
  _check_entry_mesh_coherence(isos, -1);
  
  switch (entity_type) {
    case PDM_MESH_ENTITY_CELL:
      isos->n_cell   [i_part] = n_entity;
      isos->cell_gnum[i_part] = ln_to_gn;
      break;
    case PDM_MESH_ENTITY_FACE:
      isos->n_face   [i_part] = n_entity;
      isos->face_gnum[i_part] = ln_to_gn;
      break;
    case PDM_MESH_ENTITY_EDGE:
      isos->n_edge   [i_part] = n_entity;
      isos->edge_gnum[i_part] = ln_to_gn;
      break;
    case PDM_MESH_ENTITY_VTX:
      isos->n_vtx   [i_part]  = n_entity;
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
 int                  n_group,
 int                 *group_entity_idx,
 int                 *group_entity,
 PDM_g_num_t         *group_entity_ln_to_gn
)
{
  _check_is_not_dist(isos);
  _check_entry_mesh_coherence(isos, -1);

  switch (entity_type) {
    case PDM_MESH_ENTITY_FACE:
      isos->n_group_face     [i_part] = n_group;
      isos->  group_face_idx [i_part] = group_entity_idx;
      isos->  group_face     [i_part] = group_entity;
      isos->  group_face_gnum[i_part] = group_entity_ln_to_gn;
      break;
    // case PDM_MESH_ENTITY_EDGE:
    //   isos->n_group_edge    [i_part] = n_group;
    //   isos->  group_edge_idx[i_part] = group_entity_idx;
    //   isos->  group_edge    [i_part] = group_entity;
    //   break;
    // case PDM_MESH_ENTITY_VTX:
    //   isos->n_group_vtx    [i_part]  = n_group;
    //   isos->  group_vtx_idx[i_part]  = group_entity_idx;
    //   isos->  group_vtx    [i_part]  = group_entity;
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
  _check_is_not_dist(isos);
  _check_entry_mesh_coherence(isos, -2);
  
  /**
   * TODO: add check on entry pmesh to be sure that all required data are in. 
   */

  isos->pmesh = pmesh;
  // isos->n_part     = PDM_part_mesh_n_part_get(pmesh);
  // isos->iso_n_part = PDM_part_mesh_n_part_get(pmesh);

  /* Unpack part_mesh */
  isos->entry_mesh_type = -1; // héhé
  PDM_isosurface_n_part_set(isos, PDM_part_mesh_n_part_get(pmesh));

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
    isos->n_group_face[i_part] = PDM_part_mesh_n_bound_get(pmesh, PDM_BOUND_TYPE_FACE);

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
  _check_is_not_dist(isos);
  _check_entry_mesh_coherence(isos, -3);
  
  /**
   * TODO: add check on entry pmesh_nodal to be sure that all required data are in. 
   */

  isos->pmesh_nodal = pmn;
  isos->n_part     = PDM_part_mesh_nodal_n_part_get(pmn);
  isos->iso_n_part = PDM_part_mesh_nodal_n_part_get(pmn);
}


void
PDM_isosurface_redistribution_set
(
 PDM_isosurface_t        *isos,
 PDM_extract_part_kind_t  extract_kind,
 PDM_split_dual_t         part_method
)
{
  _check_is_not_dist(isos);
 
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
  _check_is_not_dist(isos);
  
  // > Check id isosurface is valid
  if (id_isosurface>isos->n_isosurface) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: trying to defined field for id_isosurface >= n_isosurface (%d >= %d).\n", id_isosurface, isos->n_isosurface);
  }

  // > Check i_part is valid
  int n_part = 0;
  if (isos->entry_mesh_type==-1) {
    n_part = isos->n_part;
    if (n_part==-1) {
      PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: entry_mesh_type = %d but n_part isn't defined.\n", isos->entry_mesh_type);
    }
  }
  else if (isos->entry_mesh_type==-2) {
    n_part = PDM_part_mesh_n_part_get(isos->pmesh);
  }
  else if (isos->entry_mesh_type==-3) {
    n_part = PDM_part_mesh_nodal_n_part_get(isos->pmesh_nodal);
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: Impossible to defined manually field without setting mesh first.\n", isos->entry_mesh_type);
  }

  if (i_part>=n_part) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: trying to defined field for i_part >= n_part (%d >= %d).\n", i_part, n_part);
  }
 
  isos->field[id_isosurface][i_part] = field;
}


void
PDM_isosurface_gradient_set
(
 PDM_isosurface_t *isos,
 int               id_isosurface,
 int               i_part,
 double           *gradient
)
{
  _check_is_not_dist(isos);
  
  // > Check id isosurface is valid
  if (id_isosurface>isos->n_isosurface) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: trying to defined field for id_isosurface >= n_isosurface (%d >= %d).\n", id_isosurface, isos->n_isosurface);
  }

  // > Check i_part is valid
  int n_part = 0;
  if (isos->entry_mesh_type==-1) {
    n_part = isos->n_part;
    if (n_part==-1) {
      PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: entry_mesh_type = %d but n_part isn't defined.\n", isos->entry_mesh_type);
    }
  }
  else if (isos->entry_mesh_type==-2) {
    n_part = PDM_part_mesh_n_part_get(isos->pmesh);
  }
  else if (isos->entry_mesh_type==-3) {
    n_part = PDM_part_mesh_nodal_n_part_get(isos->pmesh_nodal);
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: Impossible to defined manually gradient field without setting mesh first.\n", isos->entry_mesh_type);
  }

  if (i_part>=n_part) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: trying to defined gradient field for i_part >= n_part (%d >= %d).\n", i_part, n_part);
  }

  isos->gradient[id_isosurface][i_part] = gradient;
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
  _check_is_not_dist(isos);

  if (isos->extract_kind!=PDM_EXTRACT_PART_KIND_LOCAL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: extract_kind is not PDM_EXTRACT_PART_KIND_LOCAL.\n");
  }

  if (isos->iso_owner_lparent[id_isosurface]!=NULL) {
    int n_entity = 0;
    isos->iso_owner_lparent[id_isosurface][i_part][entity_type] = ownership;
    
    if (entity_type==PDM_MESH_ENTITY_VTX) {
      n_entity = isos->iso_n_vtx[id_isosurface][i_part];
      *entity_parent_idx = isos->iso_vtx_lparent_idx[id_isosurface][i_part];
      *entity_parent     = isos->iso_vtx_lparent    [id_isosurface][i_part];
    }
    else if (entity_type==PDM_MESH_ENTITY_EDGE) {
      n_entity = isos->iso_n_edge[id_isosurface][i_part];
      *entity_parent_idx = isos->iso_edge_lparent_idx[id_isosurface][i_part];
      *entity_parent     = isos->iso_edge_lparent    [id_isosurface][i_part];
    }
    else if (entity_type==PDM_MESH_ENTITY_FACE) {
      n_entity = isos->iso_n_face   [id_isosurface][i_part];
      *entity_parent_idx = isos->iso_face_lparent_idx[id_isosurface][i_part];
      *entity_parent     = isos->iso_face_lparent    [id_isosurface][i_part];
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: has no mesh entity of type %d.\n", entity_type);
    }
    return n_entity;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_local_parent_get: Isosurface with id %d seems not computed.\n", id_isosurface);
    return 0;
  }
}


int
PDM_isosurface_vtx_parent_weight_get
(
 PDM_isosurface_t  *isos,
 int                id_isosurface,
 int                i_part,
 double           **vtx_parent_weight,
 PDM_ownership_t    ownership
)
{
  _check_is_not_dist(isos);

  if (isos->iso_owner_vtx_parent_weight[id_isosurface]!=NULL) {
    isos->iso_owner_vtx_parent_weight[id_isosurface][i_part] = ownership;
    *vtx_parent_weight = isos->iso_vtx_parent_weight[id_isosurface][i_part];
    
    return isos->iso_n_vtx[id_isosurface][i_part];
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_vtx_parent_weight_get: Isosurface with id %d seems not computed.\n", id_isosurface);
    return 0;
  }
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
  _check_is_not_dist(isos);

  if (isos->iso_owner_connec[id_isosurface]!=NULL) {
    int n_entity = 0;
    isos->iso_owner_connec[id_isosurface][i_part][connectivity_type] = ownership;

    if (connectivity_type==PDM_CONNECTIVITY_TYPE_EDGE_VTX) {
      n_entity = isos->iso_n_edge  [id_isosurface][i_part];
      *connect = isos->iso_edge_vtx[id_isosurface][i_part];
    }
    else if (connectivity_type==PDM_CONNECTIVITY_TYPE_FACE_VTX) {
      n_entity     = isos->iso_n_face      [id_isosurface][i_part];
      *connect_idx = isos->iso_face_vtx_idx[id_isosurface][i_part];
      *connect     = isos->iso_face_vtx    [id_isosurface][i_part];
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: has no connectivity of type %d.\n", connectivity_type);
    }

    return n_entity;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_connectivity_get: Isosurface with id %d seems not computed.\n", id_isosurface);
    return 0;
  }
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
  _check_is_not_dist(isos);

  if (isos->iso_owner_vtx_coord[id_isosurface]!=NULL) {
    isos->iso_owner_vtx_coord[id_isosurface][i_part] = ownership;
    *vtx_coord = isos->iso_vtx_coord[id_isosurface][i_part];

    return isos->iso_n_vtx[id_isosurface][i_part];
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_vtx_coord_get: Isosurface with id %d seems not computed.\n", id_isosurface);
    return 0;
  }
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
  _check_is_not_dist(isos);

  if (isos->iso_owner_gnum[id_isosurface]!=NULL) {
    int n_entity = 0;
    isos->iso_owner_gnum[id_isosurface][i_part][entity_type] = ownership;
    if (entity_type==PDM_MESH_ENTITY_VTX ) {
      n_entity  = isos->iso_n_vtx   [id_isosurface][i_part];
      *ln_to_gn = isos->iso_vtx_gnum[id_isosurface][i_part];
    }
    else if (entity_type==PDM_MESH_ENTITY_EDGE) {
      n_entity  = isos->iso_n_edge   [id_isosurface][i_part];
      *ln_to_gn = isos->iso_edge_gnum[id_isosurface][i_part];
    }
    else if (entity_type==PDM_MESH_ENTITY_FACE) {
      n_entity  = isos->iso_n_face   [id_isosurface][i_part];
      *ln_to_gn = isos->iso_face_gnum[id_isosurface][i_part];
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: has no cell mesh entity.\n");
    }

    return n_entity;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_ln_to_gn_get: Isosurface with id %d seems not computed.\n", id_isosurface);
    return 0;
  }
}


int
PDM_isosurface_group_get
(
 PDM_isosurface_t     *isos,
 int                   id_isosurface,
 int                   i_part,
 PDM_mesh_entities_t   entity_type,
 int                  *n_group,
 int                 **group_entity_idx,
 int                 **group_entity,
 PDM_g_num_t         **group_entity_gnum,
 PDM_ownership_t       ownership
)
{
  _check_is_not_dist(isos);

  if (entity_type==PDM_MESH_ENTITY_EDGE) {
    isos->iso_owner_edge_bnd[id_isosurface][i_part] = ownership;
    *n_group           = isos->iso_n_edge_group   [id_isosurface];
    *group_entity_idx  = isos->iso_edge_group_idx [id_isosurface][i_part];
    *group_entity      = isos->iso_edge_group_lnum[id_isosurface][i_part];
    *group_entity_gnum = isos->iso_edge_group_gnum[id_isosurface][i_part];
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: has no bounds for entity %d.\n",entity_type);
  }

  return 0;
}


void
PDM_isosurface_enable_part_to_part
(
 PDM_isosurface_t     *isos,
 int                   id_isosurface,
 PDM_mesh_entities_t   entity_type
)
{
  _check_is_not_dist(isos);

  isos->compute_ptp[id_isosurface][entity_type] = 1;
}


// TODO: changer le nom "source to iso ptp" ?
void
PDM_isosurface_part_to_part_get
(
 PDM_isosurface_t     *isos,
 int                   id_isosurface,
 PDM_mesh_entities_t   entity_type,
 PDM_part_to_part_t  **ptp,
 PDM_ownership_t       ownership
)
{
  _check_is_not_dist(isos);
  
  if (isos->compute_ptp[id_isosurface][entity_type]!=1) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: part_to_part for entity %d of isosurface %d is not computed.\n", entity_type, id_isosurface);
  }
  isos->iso_owner_ptp[id_isosurface][entity_type] = ownership;

  if (entity_type==PDM_MESH_ENTITY_VTX) {
    *ptp = isos->iso_ptp_vtx[id_isosurface];
  }
  else if (entity_type==PDM_MESH_ENTITY_EDGE) {
    *ptp = isos->iso_ptp_edge[id_isosurface];
  }
  else if (entity_type==PDM_MESH_ENTITY_FACE) {
    *ptp = isos->iso_ptp_face[id_isosurface];
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: has no cell entity.\n");
  }

}


int
PDM_isosurface_isovalue_entity_idx_get
(
 PDM_isosurface_t     *isos,
 int                   id_isosurface,
 int                   i_part,
 PDM_mesh_entities_t   entity_type,
 int                 **isovalue_entity_idx
 )
{
  _check_is_not_dist(isos);

  if (entity_type==PDM_MESH_ENTITY_VTX) {
    *isovalue_entity_idx = isos->isovalue_vtx_idx[id_isosurface][i_part];
  }
  else if (entity_type==PDM_MESH_ENTITY_EDGE) {
    *isovalue_entity_idx = isos->isovalue_edge_idx[id_isosurface][i_part];
  }
  else if (entity_type==PDM_MESH_ENTITY_FACE) {
    *isovalue_entity_idx = isos->isovalue_face_idx[id_isosurface][i_part];
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: has no cell entity.\n");
  }

  return isos->n_isovalues[id_isosurface];
}


#ifdef  __cplusplus
}
#endif
