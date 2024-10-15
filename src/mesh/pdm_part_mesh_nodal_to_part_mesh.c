/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2024       ONERA

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
#include <math.h>
#include <assert.h>
#include <limits.h>

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_priv.h"

#include "pdm_error.h"
#include "pdm_logging.h"

#include "pdm_array.h"
#include "pdm_mesh_nodal.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_vtk.h"
#include "pdm_gnum.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_part_mesh_nodal_elmts_utils.h"

#include "pdm_part_mesh_nodal_to_part_mesh.h"
#include "pdm_part_mesh_nodal_to_part_mesh_priv.h"

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

#define CHECK_PMN_TO_PM(pmn_to_pm)                                       \
  if ((pmn_to_pm) == NULL) {                                             \
    PDM_error(__FILE__, __LINE__, 0,                                     \
              "Invalid PDM_part_mesh_nodal_to_part_mesh_t instance.\n"); \
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

// TODO -> move this to pdm.c ?
static const char *_entity_name[] = {
  "CELL",
  "FACE",
  "EDGE",
  "VTX"
};

/**
 * \brief Check coherence of requested data
 */
static void
_check_inputs
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm
)
{
  PDM_part_mesh_nodal_t *pmesh_nodal = pmn_to_pm->pmesh_nodal;

  PDM_geometry_kind_t geom_kind_parent = PDM_part_mesh_nodal_principal_geom_kind_get(pmesh_nodal);

  PDM_mesh_entities_t entity_type_parent = PDM_geometry_kind_to_entity_type(geom_kind_parent);

  PDM_bool_t build_entity[PDM_MESH_ENTITY_MAX];
  memset(build_entity, 0, PDM_MESH_ENTITY_MAX * sizeof(PDM_bool_t));

  build_entity[entity_type_parent]  = PDM_TRUE;
  build_entity[PDM_MESH_ENTITY_VTX] = PDM_TRUE;

  /* Switch ON entites involved in requested connectivities */
  for (PDM_connectivity_type_t connectivity_type = 0; connectivity_type < PDM_CONNECTIVITY_TYPE_MAX; connectivity_type++) {
    if (pmn_to_pm->build_connectivity[connectivity_type] == PDM_TRUE) {
      PDM_mesh_entities_t entity_type1, entity_type2;
      if (PDM_connectivity_type_to_entity_pair(connectivity_type, &entity_type1, &entity_type2) == 0) {
        // continue;
      }
      build_entity[entity_type1] = PDM_TRUE;
      build_entity[entity_type2] = PDM_TRUE;

      if (entity_type1 < entity_type_parent ||
          entity_type2 < entity_type_parent) {
        PDM_error(__FILE__, __LINE__, 0,
                  "Invalid requested connectivity type %d : highest dimension entity type is %d\n",
                  connectivity_type, entity_type_parent);
      }

      if (entity_type1 >= entity_type2) {
        PDM_error(__FILE__, __LINE__, 0,
                  "Invalid requested connectivity type %d : only strictly downward connectivities are supported at the moment\n",
                  connectivity_type);
      }
    }
  }

  /* Check gnums */
  for (PDM_mesh_entities_t entity_type = 0; entity_type < PDM_MESH_ENTITY_MAX; entity_type++) {
    if (pmn_to_pm->compute_g_nums[entity_type] == PDM_TRUE && build_entity[entity_type] == PDM_FALSE) {
      PDM_error(__FILE__, __LINE__, 0,
                "Gnums requested for entity type %d but this entity is not involved in the requested connectivities\n",
                 entity_type);
    }
  }

  /* Check groups and part_comm_graph */
  for (PDM_bound_type_t bound_type = PDM_BOUND_TYPE_FACE; bound_type < PDM_BOUND_TYPE_MAX; bound_type++) {

    PDM_mesh_entities_t entity_type = PDM_bound_type_to_entity_type(bound_type);

    if (pmn_to_pm->transfer_groups[bound_type] == PDM_TRUE && build_entity[entity_type] == PDM_FALSE) {
      PDM_error(__FILE__, __LINE__, 0,
                "Groups requested for entity type %d but this entity is not involved in the requested connectivities\n",
                 entity_type);
    }

    if (pmn_to_pm->build_part_comm_graph[bound_type] == PDM_TRUE && build_entity[entity_type] == PDM_FALSE) {
      PDM_error(__FILE__, __LINE__, 0,
                "Inter-partition comm graph requested for entity type %d but this entity is not involved in the requested connectivities\n",
                 entity_type);
    }
  }
}


/**
 * \brief Transfer highest-dimension entities
 *   - number of entities
 *   - global IDs (if requested)
 *   - entity_to_vtx connectivity (if requested)
 *   - link pmesh_nodal->pmesh (TODO)
 */
static void
_transfer_highest_dimension_entities
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm
)
{
  PDM_part_mesh_nodal_t *pmesh_nodal = pmn_to_pm->pmesh_nodal;
  PDM_part_mesh_t       *pmesh       = pmn_to_pm->pmesh;

  int n_part = pmn_to_pm->n_part;

  // Get highest-dimension geometry_kind and part_mesh_nodal_elmts
  PDM_geometry_kind_t geom_kind = PDM_part_mesh_nodal_principal_geom_kind_get(pmesh_nodal);

  PDM_part_mesh_nodal_elmts_t *pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmesh_nodal,
                                                                                    geom_kind);

  // Get entity type
  PDM_mesh_entities_t entity_type = PDM_geometry_kind_to_entity_type(geom_kind);

  // Get entity->vtx connectivity type
  PDM_connectivity_type_t connectivity_type = PDM_entity_pair_to_connectivity_type(entity_type, PDM_MESH_ENTITY_VTX);

  for (int i_part = 0; i_part < n_part; i_part++) {

    // Number of entities
    int n_elt_tot = PDM_part_mesh_nodal_elmts_n_elmts_get(pmne, i_part);

    PDM_part_mesh_n_entity_set(pmesh,
                               i_part,
                               entity_type,
                               n_elt_tot);


    /* Global IDs */
    if (pmn_to_pm->compute_g_nums[entity_type] == PDM_TRUE) {

      PDM_g_num_t *parent_ln_to_gn = NULL;
      PDM_malloc(parent_ln_to_gn, n_elt_tot, PDM_g_num_t);

      int  n_section   = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
      int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

      int i_parent = -1;

      for (int i_section = 0; i_section < n_section; i_section++) {

        int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne, sections_id[i_section], i_part);

        int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                                   sections_id[i_section],
                                                                   i_part,
                                                                   PDM_OWNERSHIP_BAD_VALUE);

        PDM_g_num_t *ln_to_gn = PDM_part_mesh_nodal_elmts_g_num_get(pmne,
                                                                    sections_id[i_section],
                                                                    i_part,
                                                                    PDM_OWNERSHIP_BAD_VALUE);

        if (ln_to_gn == NULL) {
          PDM_error(__FILE__, __LINE__, 0, "TODO -> use parent_entity_g_num instead? (what do we do with poly2d?)\n");
        }

        for (int i_elt = 0; i_elt < n_elt; i_elt++) {
          if (parent_num == NULL) {
            i_parent++;
          }
          else {
            i_parent = parent_num[i_elt];
          }

          parent_ln_to_gn[i_parent] = ln_to_gn[i_elt];
        }

      } // End loop on sections

      PDM_part_mesh_entity_ln_to_gn_set(pmesh,
                                        i_part,
                                        entity_type,
                                        parent_ln_to_gn,
                                        PDM_OWNERSHIP_KEEP);
    }


    /* Entity->Vtx connectivity */
    if (pmn_to_pm->build_connectivity[connectivity_type] == PDM_TRUE) {
      int *entity_to_vtx_idx = NULL;
      int *entity_to_vtx     = NULL;
      PDM_part_mesh_nodal_elmts_cell_vtx_connect_get(pmne,
                                                     i_part,
                                                     &entity_to_vtx_idx,
                                                     &entity_to_vtx);

      PDM_part_mesh_connectivity_set(pmesh,
                                     i_part,
                                     connectivity_type,
                                     entity_to_vtx,
                                     entity_to_vtx_idx,
                                     PDM_OWNERSHIP_KEEP);
    }

  } // End loop on parts

  if (pmn_to_pm->build_connectivity[connectivity_type] == PDM_TRUE) {
    pmn_to_pm->connectivity_done[connectivity_type] = PDM_TRUE;
  }

  if (pmn_to_pm->compute_g_nums[entity_type] == PDM_TRUE) {
    pmn_to_pm->g_nums_done[entity_type] = PDM_TRUE;
  }

}


/**
 * \brief Transfer vertices
 *   - coordinates
 *   - global IDs (if requested) (/!\ necessary if entity gnums are requested)
 */
static void
_transfer_vtx
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm
)
{
  PDM_part_mesh_nodal_t *pmesh_nodal = pmn_to_pm->pmesh_nodal;
  PDM_part_mesh_t       *pmesh       = pmn_to_pm->pmesh;

  int n_part = pmn_to_pm->n_part;

  for (int i_part = 0; i_part < n_part; i_part++) {

    int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmesh_nodal, i_part);
    double      *_vtx_coord    = PDM_part_mesh_nodal_vtx_coord_get(pmesh_nodal, i_part);
    PDM_g_num_t *_vtx_ln_to_gn = PDM_part_mesh_nodal_vtx_g_num_get(pmesh_nodal, i_part);
    double      *vtx_coord     = NULL;
    PDM_g_num_t *vtx_ln_to_gn  = NULL;

    if (pmesh_nodal->vtx[i_part]->owner == PDM_OWNERSHIP_KEEP &&
        pmn_to_pm->vtx_ownership_pmesh  == PDM_OWNERSHIP_KEEP) {
      // Copy vtx since both structures claim ownership
      PDM_malloc(vtx_coord, n_vtx * 3, double);
      memcpy(vtx_coord, _vtx_coord, sizeof(double) * n_vtx * 3);

      if (pmn_to_pm->compute_g_nums[PDM_MESH_ENTITY_VTX] == PDM_TRUE) {
        PDM_malloc(vtx_ln_to_gn, n_vtx, PDM_g_num_t);
        memcpy(vtx_ln_to_gn, _vtx_ln_to_gn, sizeof(PDM_g_num_t) * n_vtx);
      }
    }
    else {
      // Map pointers since only one structure claims ownership on vtx
      vtx_coord    = _vtx_coord;
      vtx_ln_to_gn = _vtx_ln_to_gn;
    }

    PDM_part_mesh_n_entity_set(pmesh,
                               i_part,
                               PDM_MESH_ENTITY_VTX,
                               n_vtx);

    PDM_part_mesh_vtx_coord_set(pmesh,
                                i_part,
                                vtx_coord,
                                pmn_to_pm->vtx_ownership_pmesh);

    if (pmn_to_pm->compute_g_nums[PDM_MESH_ENTITY_VTX] == PDM_TRUE) {
      PDM_part_mesh_entity_ln_to_gn_set(pmesh,
                                        i_part,
                                        PDM_MESH_ENTITY_VTX,
                                        vtx_ln_to_gn,
                                        pmn_to_pm->vtx_ownership_pmesh);
    }

  } // End loop on parts

  if (pmn_to_pm->compute_g_nums[PDM_MESH_ENTITY_VTX] == PDM_TRUE) {
    pmn_to_pm->g_nums_done[PDM_MESH_ENTITY_VTX] = PDM_TRUE;
  }
}


/**
 * \brief Generate downward connectivity
 */
static void
_generate_downward_connectivity
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm,
  PDM_connectivity_type_t             connectivity_type
)
{
  PDM_part_mesh_nodal_t *pmesh_nodal = pmn_to_pm->pmesh_nodal;
  PDM_part_mesh_t       *pmesh       = pmn_to_pm->pmesh;

  int n_part = pmn_to_pm->n_part;

  PDM_geometry_kind_t geom_kind_parent = PDM_part_mesh_nodal_principal_geom_kind_get(pmesh_nodal);

  PDM_mesh_entities_t entity_type1, entity_type2;
  PDM_connectivity_type_to_entity_pair(connectivity_type, &entity_type1, &entity_type2);

  PDM_geometry_kind_t geom_kind1 = PDM_entity_type_to_geometry_kind(entity_type1);
  PDM_geometry_kind_t geom_kind2 = PDM_entity_type_to_geometry_kind(entity_type2);

  PDM_part_mesh_nodal_elmts_t *pmne1 = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmesh_nodal, geom_kind1);
  PDM_part_mesh_nodal_elmts_t *pmne2 = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmesh_nodal, geom_kind2);


  int use_fake_pmne = 0;

  /* Special case (3D + face->edge) => we build a fake part_mesh_nodal_elmts */
  if (geom_kind_parent  == PDM_GEOMETRY_KIND_VOLUMIC &&
      connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_EDGE) {
    use_fake_pmne = 1;

    pmne1 = PDM_part_mesh_nodal_elmts_create(2, n_part, pmn_to_pm->comm);

    int id_section = PDM_part_mesh_nodal_elmts_add(pmne1, PDM_MESH_NODAL_POLY_2D);

    for (int i_part = 0; i_part < n_part; i_part++) {

      int n_face = PDM_part_mesh_n_entity_get(pmesh,
                                              i_part,
                                              PDM_MESH_ENTITY_FACE);

      int *face_vtx_idx = NULL;
      int *face_vtx     = NULL;
      PDM_part_mesh_connectivity_get(pmesh,
                                     i_part,
                                     PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                     &face_vtx,
                                     &face_vtx_idx,
                                     PDM_OWNERSHIP_BAD_VALUE);
      assert(face_vtx_idx != NULL);

      PDM_part_mesh_nodal_elmts_section_poly2d_set(pmne1,
                                                   id_section,
                                                   i_part,
                                                   n_face,
                                                   face_vtx_idx,
                                                   face_vtx,
                                                   NULL,
                                                   NULL,
                                                   PDM_OWNERSHIP_USER);
    }
  }


  /* Decompose pmne1 into lower-dimension entities */
  int **child_to_parent_idx    = NULL;
  int **child_to_parent        = NULL;
  int  *n_entity2              = NULL;
  int **entity2_to_vtx_idx     = NULL;
  int **entity2_to_vtx         = NULL;
  int **entity1_to_entity2_idx = NULL;
  int **entity1_to_entity2     = NULL;
  PDM_part_mesh_nodal_elmts_compute_child_parent(pmne1,
                                                 pmne2,
                                                 entity_type2,
                                                 PDM_TRUE,
                                                 &child_to_parent_idx,
                                                 &child_to_parent,
                                                 &n_entity2,
                                                 &entity2_to_vtx_idx,
                                                 &entity2_to_vtx,
                                                 &entity1_to_entity2_idx,
                                                 &entity1_to_entity2);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_free(child_to_parent_idx[i_part]);
    PDM_free(child_to_parent    [i_part]);
  }
  PDM_free(child_to_parent_idx);
  PDM_free(child_to_parent    );


  /* Second connectivity type : entity2->vtx (byproduct of decomposition) */
  PDM_connectivity_type_t connectivity_type2 = PDM_entity_pair_to_connectivity_type(entity_type2,
                                                                                    PDM_MESH_ENTITY_VTX);
  pmn_to_pm->connectivity_done[connectivity_type ] = PDM_TRUE;
  pmn_to_pm->connectivity_done[connectivity_type2] = PDM_TRUE;

  for (int i_part = 0; i_part < n_part; i_part++) {
    // Entity1->Entity2
    PDM_part_mesh_connectivity_set(pmesh,
                                   i_part,
                                   connectivity_type,
                                   entity1_to_entity2    [i_part],
                                   entity1_to_entity2_idx[i_part],
                                   PDM_OWNERSHIP_KEEP);

    // Entity2
    PDM_part_mesh_n_entity_set(pmesh,
                               i_part,
                               entity_type2,
                               n_entity2[i_part]);

    // Entity2->Vtx
    if (entity_type2 != PDM_MESH_ENTITY_VTX) {
      PDM_part_mesh_connectivity_set(pmesh,
                                     i_part,
                                     connectivity_type2,
                                     entity2_to_vtx    [i_part],
                                     entity2_to_vtx_idx[i_part],
                                     PDM_OWNERSHIP_KEEP);
    }
  }
  PDM_free(n_entity2             );
  PDM_free(entity2_to_vtx_idx    );
  PDM_free(entity2_to_vtx        );
  PDM_free(entity1_to_entity2_idx);
  PDM_free(entity1_to_entity2    );

  if (use_fake_pmne) {
    PDM_part_mesh_nodal_elmts_free(pmne1);
  }
}


/**
 * \brief Generate global IDs for lower-dimension entities using vtx global IDs
 */
static void
_generate_gnum
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm,
  PDM_mesh_entities_t                 entity_type
)
{
  PDM_part_mesh_t *pmesh = pmn_to_pm->pmesh;

  int n_part = pmn_to_pm->n_part;

  PDM_g_num_t **entity_parent_gnum = NULL;
  PDM_malloc(entity_parent_gnum, n_part, PDM_g_num_t *);

  /**
   *  First, let's make sure that the entity->vtx connectivity is consistent on all processes.
   *  To do so, the vertices of each entity are reordered so that
   *    - the 1st vertex of each entity is the one with the lowest global ID
   *    - the 2nd vertex is the one with the lowest global ID among the 1st vertex's direct neighbors
   *
   *  This procedure may flip the orientation of some entities, so we need to propagate the change
   *  in the direct downward connectivity (*->entity)
   */
  int  *n_entity              = NULL;
  int **entity_to_vtx_idx     = NULL;
  int **entity_to_vtx         = NULL;
  PDM_malloc(n_entity,          n_part, int  );
  PDM_malloc(entity_to_vtx_idx, n_part, int *);
  PDM_malloc(entity_to_vtx,     n_part, int *);

  // Get entity->vtx connectivity type
  PDM_connectivity_type_t connectivity_type_vtx  = PDM_entity_pair_to_connectivity_type(entity_type, PDM_MESH_ENTITY_VTX);

  // Get direct downward *->entity connectivity type
  PDM_connectivity_type_t connectivity_type_down = PDM_entity_pair_to_connectivity_type(entity_type-1, entity_type);

  // Get entity->vtx connectivity
  PDM_bool_t owner_entity_vtx_idx = PDM_FALSE;
  int max_n_vtx = 2;
  for (int i_part = 0; i_part < n_part; i_part ++) {
    n_entity[i_part] = PDM_part_mesh_n_entity_get(pmesh,
                                                  i_part,
                                                  entity_type);

    PDM_part_mesh_connectivity_get(pmesh,
                                   i_part,
                                   connectivity_type_vtx,
                                   &entity_to_vtx    [i_part],
                                   &entity_to_vtx_idx[i_part],
                                   PDM_OWNERSHIP_BAD_VALUE);
    if (entity_to_vtx_idx[i_part] == NULL) {
      if (entity_type != PDM_MESH_ENTITY_EDGE) {
        PDM_error(__FILE__, __LINE__, 0,
                  "We are supposed to have connectivity type %d, "
                  "yet n_entity = %d, entity_to_vtx_idx: %p, entity_to_vtx: %p (entity = %s)\n",
                  connectivity_type_vtx, n_entity[i_part], (void *) entity_to_vtx_idx[i_part], (void *) entity_to_vtx[i_part],
                  _entity_name[entity_type]);
      }
    }
    else {
      for (int i_entity = 0; i_entity < n_entity[i_part]; i_entity++) {
        int n_vtx = entity_to_vtx_idx[i_part][i_entity+1] - entity_to_vtx_idx[i_part][i_entity];
        max_n_vtx = PDM_MAX(max_n_vtx, n_vtx);
      }
    }
  }

  int *tmp_entity_to_vtx = NULL;
  PDM_malloc(tmp_entity_to_vtx, max_n_vtx, int);


  for (int i_part = 0; i_part < n_part; i_part++) {

    // Get vertex global IDs (mandatory)
    PDM_g_num_t *vtx_ln_to_gn = NULL;
    PDM_part_mesh_entity_ln_to_gn_get(pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_VTX,
                                      &vtx_ln_to_gn,
                                      PDM_OWNERSHIP_BAD_VALUE);
    assert(vtx_ln_to_gn != NULL);


    int *orientation = NULL;
    PDM_malloc(orientation,                n_entity[i_part],             int        );
    PDM_malloc(entity_parent_gnum[i_part], n_entity[i_part] * max_n_vtx, PDM_g_num_t);

    // Get direct downward connectivity
    int n_entity2 = PDM_part_mesh_n_entity_get(pmesh,
                                               i_part,
                                               entity_type-1);

    int *entity2_to_entity_idx = NULL;
    int *entity2_to_entity     = NULL;
    PDM_part_mesh_connectivity_get(pmesh,
                                   i_part,
                                   connectivity_type_down,
                                   &entity2_to_entity,
                                   &entity2_to_entity_idx,
                                   PDM_OWNERSHIP_BAD_VALUE);

    // Reorder each entity's vertices
    if (entity_type == PDM_MESH_ENTITY_EDGE) {

      for (int i_entity = 0; i_entity < n_entity[i_part]; i_entity++) {
        int         *ev = &entity_to_vtx     [i_part][2*i_entity];
        PDM_g_num_t *ep = &entity_parent_gnum[i_part][2*i_entity];

        int i_vtx0 = ev[0] - 1;
        int i_vtx1 = ev[1] - 1;

        PDM_g_num_t g_vtx0 = vtx_ln_to_gn[i_vtx0];
        PDM_g_num_t g_vtx1 = vtx_ln_to_gn[i_vtx1];

        if (g_vtx0 < g_vtx1) {
          orientation[i_entity] = 1;
          ep[0] = g_vtx0;
          ep[1] = g_vtx1;
        }
        else {
          orientation[i_entity] = -1;
          ev[0] = i_vtx1 + 1;
          ev[1] = i_vtx0 + 1;
          ep[0] = g_vtx1;
          ep[1] = g_vtx0;
        }

      } // End loop on edges

    } // End if edges

    else { // Faces
      assert(entity_type == PDM_MESH_ENTITY_FACE);
      for (int i_entity = 0; i_entity < n_entity[i_part]; i_entity++) {

        int n_vtx = entity_to_vtx_idx[i_part][i_entity+1] - entity_to_vtx_idx[i_part][i_entity];

        int         *ev = &entity_to_vtx     [i_part][entity_to_vtx_idx[i_part][i_entity]];
        PDM_g_num_t *ep = &entity_parent_gnum[i_part][max_n_vtx * i_entity];

        // pick 1st vtx (lowest global ID)
        PDM_g_num_t min_vtx;
#ifdef PDM_LONG_G_NUM
        min_vtx = LONG_MAX;
#else
        min_vtx = INT_MAX;
#endif
        int start = 0;
        for (int i = 0; i < n_vtx; i++) {
          int i_vtx = ev[i] - 1;
          if (vtx_ln_to_gn[i_vtx] < min_vtx) {
            start   = i;
            min_vtx = vtx_ln_to_gn[i_vtx];
          }
        }

        memcpy(tmp_entity_to_vtx,
               &entity_to_vtx[i_part][entity_to_vtx_idx[i_part][i_entity]],
               sizeof(int) * n_vtx);

        // pick 2nd vtx (neighbor of 1st vtx with lowest global ID)
        int i_prev = ev[(start + n_vtx - 1)%n_vtx] - 1;
        int i_next = ev[(start + 1)        %n_vtx] - 1;

        // reorder vertices of current entity
        if (vtx_ln_to_gn[i_next] < vtx_ln_to_gn[i_prev]) {
          orientation[i_entity] = 1;
          for (int i = 0; i < n_vtx; i++) {
            ev[i] = tmp_entity_to_vtx[(start+i)%n_vtx];
          }
        }
        else {
          orientation[i_entity] = -1;
          for (int i = 0; i < n_vtx; i++) {
            ev[i] = tmp_entity_to_vtx[(start+n_vtx-i)%n_vtx];
          }
        }

        // Store nuplet of vertex global IDs in new order
        for (int i = 0; i < n_vtx; i++) {
          ep[i] = vtx_ln_to_gn[ev[i] - 1];
        }
        for (int i = n_vtx; i < max_n_vtx; i++) {
          ep[i] = 0; // ¯\_(ツ)_/¯
        }

      } // End loop on entities
    } // End if faces

    // Update downward connectivity
    if (pmn_to_pm->build_connectivity[connectivity_type_down] == PDM_TRUE) {
      for (int i_entity2 = 0; i_entity2 < n_entity2; i_entity2++) {
        for (int idx_entity = entity2_to_entity_idx[i_entity2]; idx_entity < entity2_to_entity_idx[i_entity2+1]; idx_entity++) {
          int i_entity = PDM_ABS(entity2_to_entity[idx_entity]) - 1;
          entity2_to_entity[idx_entity] *= orientation[i_entity];
        }
      }
    }

    PDM_free(orientation);
  } // End loop on parts


  PDM_free(tmp_entity_to_vtx);


  // Generate entity global IDs
  PDM_gen_gnum_t *gen_gnum = PDM_gnum_create(3,         // unused
                                             n_part,
                                             PDM_FALSE, // unused
                                             1.,        // unused
                                             pmn_to_pm->comm,
                                             PDM_OWNERSHIP_USER);

  PDM_gnum_set_parents_nuplet(gen_gnum, max_n_vtx);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_gnum_set_from_parents(gen_gnum,
                              i_part,
                              n_entity          [i_part],
                              entity_parent_gnum[i_part]);
  }

  PDM_gnum_compute(gen_gnum);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_free(entity_parent_gnum[i_part]);

    PDM_g_num_t *entity_ln_to_gn = PDM_gnum_get(gen_gnum, i_part);

    PDM_part_mesh_entity_ln_to_gn_set(pmesh,
                                      i_part,
                                      entity_type,
                                      entity_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);
  }

  pmn_to_pm->g_nums_done[entity_type] = PDM_TRUE;

  // Free memory
  PDM_free(entity_parent_gnum);
  PDM_gnum_free(gen_gnum);

  if (owner_entity_vtx_idx == PDM_TRUE) {
    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_free(entity_to_vtx_idx[i_part]);
    }
  }

  PDM_free(n_entity         );
  PDM_free(entity_to_vtx_idx);
  PDM_free(entity_to_vtx    );
}


/**
 * \brief Transfer groups
 */
static void
_transfer_groups
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm,
  PDM_bound_type_t                    bound_type
)
{
  PDM_part_mesh_nodal_t *pmesh_nodal = pmn_to_pm->pmesh_nodal;
  PDM_part_mesh_t       *pmesh       = pmn_to_pm->pmesh;

  int n_part = pmn_to_pm->n_part;

  // Get geometry kind
  PDM_mesh_entities_t entity_type = PDM_bound_type_to_entity_type(bound_type);

  PDM_geometry_kind_t geom_kind = PDM_entity_type_to_geometry_kind(entity_type);

  // Number of groups
  int n_group = PDM_part_mesh_nodal_n_group_get(pmesh_nodal, geom_kind);

  PDM_part_mesh_n_bound_set(pmesh, bound_type, n_group);

  // Group->entity
  for (int i_group = 0; i_group < n_group; i_group++) {
    for (int i_part = 0; i_part < n_part; i_part++) {

      int          n_group_elmt   = 0;
      int         *group_elmt     = NULL;
      PDM_g_num_t *group_ln_to_gn = NULL;
      PDM_part_mesh_nodal_group_get(pmesh_nodal,
                                    geom_kind,
                                    i_part,
                                    i_group,
                                    &n_group_elmt,
                                    &group_elmt,
                                    &group_ln_to_gn,
                                    PDM_OWNERSHIP_BAD_VALUE);

      // Make deep copy
      // OK since pmne decomposition preserves local numbering of bound entities
      int         *copy_group_elmt     = NULL;
      PDM_g_num_t *copy_group_ln_to_gn = NULL;
      PDM_malloc(copy_group_elmt, n_group_elmt, int);
      memcpy(copy_group_elmt, group_elmt, sizeof(int) * n_group_elmt);
      if (pmn_to_pm->compute_g_nums[entity_type] == PDM_TRUE) {
        PDM_malloc(copy_group_ln_to_gn, n_group_elmt, PDM_g_num_t);
        memcpy(copy_group_ln_to_gn, group_ln_to_gn, sizeof(PDM_g_num_t) * n_group_elmt);
      }

      PDM_part_mesh_bound_set(pmesh,
                              i_part,
                              i_group,
                              bound_type,
                              n_group_elmt,
                              copy_group_elmt,
                              copy_group_ln_to_gn,
                              PDM_OWNERSHIP_KEEP);
    }
  }

  pmn_to_pm->groups_done[bound_type] = PDM_TRUE;
}


/**
 * \brief Generate intermediate entities
 *   - requested downward connectivities
 *   - global IDs (if requested)
 *   - groups (if requested)
 *   - link elmt->entity (if requested)
 *   - inter-partition communication graph (if requested) (TODO)
 */
static void
_generate_entities
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm
)
{
  /**
   * /!\ comment gérer le cas cell->face, (cell->edge,) face->edge, edge->vtx ??
   */

  for (PDM_connectivity_type_t connectivity_type = 0; connectivity_type < PDM_CONNECTIVITY_TYPE_MAX; connectivity_type++) {
    if (pmn_to_pm->build_connectivity[connectivity_type] == PDM_FALSE ||
        pmn_to_pm->connectivity_done [connectivity_type] == PDM_TRUE) {
      continue;
    }

    // Connectivity and keep link elt->entity
    _generate_downward_connectivity(pmn_to_pm, connectivity_type);

    // Global IDs of lower-dimension entity
    PDM_mesh_entities_t entity_type1, entity_type2;
    PDM_connectivity_type_to_entity_pair(connectivity_type, &entity_type1, &entity_type2);

    if (pmn_to_pm->compute_g_nums[entity_type2] == PDM_TRUE &&
        pmn_to_pm->g_nums_done   [entity_type2] == PDM_FALSE) {
      _generate_gnum(pmn_to_pm, entity_type2);
    }

    // Groups
    PDM_bound_type_t bound_type = PDM_entity_type_to_bound_type(entity_type2);
    if (pmn_to_pm->transfer_groups[bound_type] == PDM_TRUE &&
        pmn_to_pm->groups_done    [bound_type] == PDM_FALSE) {
      _transfer_groups(pmn_to_pm, bound_type);
    }
  }
}


/**
 * \brief Store link elt->entity in pmesh_nodal struct
 *        By construction, it is identical to parent_num (if it exists).
 */
static void
_store_link_pmn_to_pm
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm
)
{
  if (pmn_to_pm->keep_link_elmt_to_entity == PDM_FALSE) {
    return;
  }

  PDM_part_mesh_nodal_t *pmesh_nodal = pmn_to_pm->pmesh_nodal;
  int n_part = pmn_to_pm->n_part;

  // Get highest-dimension geometry_kind and part_mesh_nodal_elmts
  PDM_geometry_kind_t geom_kind_main = PDM_part_mesh_nodal_principal_geom_kind_get(pmesh_nodal);

  for (PDM_geometry_kind_t geom_kind = geom_kind_main; geom_kind < PDM_GEOMETRY_KIND_CORNER; geom_kind++) {

    PDM_part_mesh_nodal_elmts_t *pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmesh_nodal,
                                                                                      geom_kind);
    if (pmne == NULL) {
      continue;
    }

    int  n_section   = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
    int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

    for (int i_part = 0; i_part < n_part; i_part++) {

      int i_parent = -1;

      for (int i_section = 0; i_section < n_section; i_section++) {

        int *elt_to_entity = NULL;
        elt_to_entity = PDM_part_mesh_nodal_elmts_section_elt_to_entity_get(pmne,
                                                                            sections_id[i_section],
                                                                            i_part,
                                                                            PDM_OWNERSHIP_BAD_VALUE);
        if (elt_to_entity != NULL) {
          printf("Warning : elt_to_entity already exists => What should we do??? (keep untouched for now)\n");
          fflush(stdout);
          continue;
        }

        int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne, sections_id[i_section], i_part);

        int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                                   sections_id[i_section],
                                                                   i_part,
                                                                   PDM_OWNERSHIP_BAD_VALUE);

        PDM_malloc(elt_to_entity, n_elt, int);

        if (parent_num != NULL) {
          memcpy(elt_to_entity, parent_num, sizeof(int) * n_elt);
        }
        else {
          for (int i_elt = 0; i_elt < n_elt; i_elt++) {
            elt_to_entity[i_elt] = ++i_parent;
          }
        }

        PDM_part_mesh_nodal_elmts_section_elt_to_entity_set(pmne,
                                                            sections_id[i_section],
                                                            i_part,
                                                            elt_to_entity,
                                                            PDM_OWNERSHIP_KEEP);

      } // End loop on sections

    } // End loop on parts

  } // End loop on geom_kind
}

/**
 * \brief Get rid of unrequested data
 */
static void
_clean_up
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm
)
{
  PDM_part_mesh_t *pmesh = pmn_to_pm->pmesh;

  int n_part = pmn_to_pm->n_part;

  for (PDM_connectivity_type_t connectivity_type = 0; connectivity_type < PDM_CONNECTIVITY_TYPE_MAX; connectivity_type++) {
    if (pmn_to_pm->build_connectivity[connectivity_type] == PDM_FALSE) {

      for (int i_part = 0; i_part < n_part; i_part++) {
        int *entity1_to_entity2_idx = NULL;
        int *entity1_to_entity2     = NULL;
        PDM_part_mesh_connectivity_get(pmesh,
                                       i_part,
                                       connectivity_type,
                                       &entity1_to_entity2,
                                       &entity1_to_entity2_idx,
                                       PDM_OWNERSHIP_USER);
        PDM_free(entity1_to_entity2_idx);
        PDM_free(entity1_to_entity2);

        PDM_part_mesh_connectivity_set(pmesh,
                                       i_part,
                                       connectivity_type,
                                       NULL,
                                       NULL,
                                       PDM_OWNERSHIP_KEEP);
      }

    }
  }

}

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

PDM_part_mesh_nodal_to_part_mesh_t *
PDM_part_mesh_nodal_to_part_mesh_create
(
  PDM_part_mesh_nodal_t *pmesh_nodal,
  PDM_bool_t             keep_link_elmt_to_entity,
  PDM_ownership_t        vtx_ownership_pmesh
)
{
  if (pmesh_nodal == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Invalid PDM_part_mesh_nodal_t instance\n");
  }

  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm = NULL;

  PDM_malloc(pmn_to_pm, 1, PDM_part_mesh_nodal_to_part_mesh_t);

  memset(pmn_to_pm, 0, sizeof(PDM_part_mesh_nodal_to_part_mesh_t));

  pmn_to_pm->pmesh_nodal              = pmesh_nodal;
  pmn_to_pm->keep_link_elmt_to_entity = keep_link_elmt_to_entity;
  pmn_to_pm->vtx_ownership_pmesh      = vtx_ownership_pmesh;
  pmn_to_pm->owner_pmesh              = PDM_OWNERSHIP_BAD_VALUE;
  pmn_to_pm->comm                     = pmesh_nodal->comm;
  pmn_to_pm->n_part                   = pmesh_nodal->n_part;

  return pmn_to_pm;
}


void
PDM_part_mesh_nodal_to_part_mesh_connectivity_enable
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm,
  PDM_connectivity_type_t             connectivity_type
)
{
  CHECK_PMN_TO_PM(pmn_to_pm);
  // TODO: throw error/warning if requested connectivity is not handled

  pmn_to_pm->build_connectivity[connectivity_type] = PDM_TRUE;
}


void
PDM_part_mesh_nodal_to_part_mesh_g_nums_enable
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm,
  PDM_mesh_entities_t                 entity_type
)
{
  CHECK_PMN_TO_PM(pmn_to_pm);

  pmn_to_pm->compute_g_nums[entity_type] = PDM_TRUE;
}


void
PDM_part_mesh_nodal_to_part_mesh_groups_enable
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm,
  PDM_bound_type_t                    bound_type
)
{
  CHECK_PMN_TO_PM(pmn_to_pm);

  pmn_to_pm->transfer_groups[bound_type] = PDM_TRUE;
}


void
PDM_part_mesh_nodal_to_part_mesh_part_comm_graph_enable
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm,
  PDM_bound_type_t                    bound_type
)
{
  CHECK_PMN_TO_PM(pmn_to_pm);

  pmn_to_pm->build_part_comm_graph[bound_type] = PDM_TRUE;
}


void
PDM_part_mesh_nodal_to_part_mesh_compute
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm
)
{
  CHECK_PMN_TO_PM(pmn_to_pm);

  if (pmn_to_pm->pmesh != NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_mesh_nodal_to_part_mesh already computed\n");
  }

  PDM_part_mesh_nodal_t *pmesh_nodal = pmn_to_pm->pmesh_nodal;
  if (pmesh_nodal == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Invalid PDM_part_mesh_nodal_t instance\n");
  }

  /* Check coherence in requested connectivities, gnums, etc */
  _check_inputs(pmn_to_pm);

  pmn_to_pm->pmesh = PDM_part_mesh_create(pmn_to_pm->n_part, pmn_to_pm->comm);

  /* Entities of highest dimension */
  _transfer_highest_dimension_entities(pmn_to_pm);

  /* Vertices */
  _transfer_vtx(pmn_to_pm);

  /* Connectivities and lower-dimension entities */
  _generate_entities(pmn_to_pm);

  /* Store link pmesh_nodal -> pmesh */
  _store_link_pmn_to_pm(pmn_to_pm);

  /* Get rid of temporary connectivities */
  _clean_up(pmn_to_pm);
}


void
PDM_part_mesh_nodal_to_part_mesh_part_mesh_get
(
  PDM_part_mesh_nodal_to_part_mesh_t  *pmn_to_pm,
  PDM_part_mesh_t                    **pmesh,
  PDM_ownership_t                      ownership
)
{
  CHECK_PMN_TO_PM(pmn_to_pm);

  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    pmn_to_pm->owner_pmesh = ownership;
  }

  *pmesh = pmn_to_pm->pmesh;
}


void
PDM_part_mesh_nodal_to_part_mesh_free
(
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm
)
{
  CHECK_PMN_TO_PM(pmn_to_pm);

  if (pmn_to_pm->owner_pmesh == PDM_OWNERSHIP_KEEP) {
    PDM_part_mesh_free(pmn_to_pm->pmesh);
  }

  PDM_free(pmn_to_pm);
}


#undef CHECK_PMN_TO_PM

#ifdef  __cplusplus
}
#endif
