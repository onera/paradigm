#ifndef PDM_ISOSURFACE_PRIV_H
#define PDM_ISOSURFACE_PRIV_H

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
#include "pdm_part_mesh_nodal.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_part_mesh.h"
#include "pdm_dmesh.h"
#include "pdm_part_to_part.h"
#include "pdm_block_to_part.h"
#include "pdm_extract_part.h"
#include "pdm_isosurface.h"

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

#define PDM_ISOSURFACE_CHECK_ENTRY_MESH_COHERENCE(isos, entry_type) \
  if ((isos)->entry_mesh_type==0) { \
    (isos)->entry_mesh_type=(entry_type); \
  } \
  else if ((isos)->entry_mesh_type!=(entry_type)) { \
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t:entry_mesh_type already set to %d.\n", (isos)->entry_mesh_type); \
  }

#define PDM_ISOSURFACE_CHECK_ID(isos, id_isosurface) \
  if ((id_isosurface) >= (isos)->n_isosurface) { \
    PDM_error(__FILE__, __LINE__, 0, "Invalid id_isosurface %d (n_isosurface = %d).\n", (id_isosurface), (isos)->n_isosurface); \
  }

#define PDM_ISOSURFACE_CHECK_COMPUTED(isos, id_isosurface) \
  if ((isos)->is_computed[(id_isosurface)] == PDM_FALSE) { \
    PDM_error(__FILE__, __LINE__, 0, "Isosurface with id %d was not computed.\n", (id_isosurface)); \
  }

#define PDM_ISOSURFACE_CHECK_ENTITY_TYPE(entity_type) \
  if ((entity_type) != PDM_MESH_ENTITY_VTX  && \
      (entity_type) != PDM_MESH_ENTITY_EDGE && \
      (entity_type) != PDM_MESH_ENTITY_FACE) { \
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: has no mesh entity of type %d.\n", (entity_type)); \
  }

/*============================================================================
 * Type definitions
 *============================================================================*/


struct _pdm_isosurface_t {

  // =========
  // > General
  PDM_MPI_Comm comm;


  // ========================
  // > Entry mesh information

  // > Isosurface switch
  // TODO: is entry_mesh_type sign necessary ???
  int entry_is_part; // -1: invalid, 0: dist, 1: part
  int entry_mesh_type; //  0: undef, 1: dist_alamano, 2: dmesh, 3: dmesh_nodal, -1: part_alamano, -2: pmesh, -3: pmesh_nodal

  // > Mesh information
  int entry_mesh_dim;



  // ====================
  // > Isosurface options

  // > Isosurfaces
  int n_isosurface;

  // > Isosurface tolerance
  double ISOSURFACE_EPS;

  // > Isosurface type
  PDM_iso_surface_kind_t *kind;

  // > Isovalues
  int     *n_isovalues;
  double **isovalues;

  // > Equation args
  PDM_isosurface_field_function_t *field_function;
  double                         **eq_coeffs;
  int                             *use_gradient;

  // > Function args
  _pdm_isosurface_field_function_t *iso_func;

  // > Redistribution
  PDM_extract_part_t      **extrp;
  PDM_extract_part_kind_t   extract_kind;
  PDM_split_dual_t          part_method;

  // > Link with entry mesh
  PDM_bool_t *compute_ptp[PDM_MESH_ENTITY_MAX];
  
  // ========================
  // > Distributed entry data
  int dist_to_part_computed;
  PDM_block_to_part_t *btp_vtx;

  // > Mesh structs
  PDM_dmesh_t       *dmesh;
  PDM_dmesh_nodal_t *dmesh_nodal;

  // > Distribution
  PDM_g_num_t *distrib_cell;
  PDM_g_num_t *distrib_face;
  PDM_g_num_t *distrib_edge;
  PDM_g_num_t *distrib_vtx;

  // > Vertices
  double *dvtx_coord;

  // > Connectivities
  int         *dcell_face_idx;
  PDM_g_num_t *dcell_face;
  int         *dface_edge_idx;
  PDM_g_num_t *dface_edge;
  int         *dface_vtx_idx;
  PDM_g_num_t *dface_vtx;
  PDM_g_num_t *dedge_vtx;

  int we_have_edges;

  // > Boundaries
  int          n_dgroup_face;
  int         *dgroup_face_idx;
  PDM_g_num_t *dgroup_face;

  // > Field
  double **dfield;
  double **dgradient;

  double ***extract_field;
  PDM_part_mesh_nodal_t *extract_pmesh_nodal;

  // ========================
  // > Partitioned entry data

  // > Mesh structs
  PDM_part_mesh_t       *pmesh;
  PDM_part_mesh_nodal_t *pmesh_nodal;

  // > Partition
  int  n_part;
  int *n_cell;
  int *n_face;
  int *n_edge;
  int *n_vtx;
  PDM_g_num_t **cell_gnum;
  PDM_g_num_t **face_gnum;
  PDM_g_num_t **edge_gnum;
  PDM_g_num_t ** vtx_gnum;

  // > Vertices
  double **vtx_coord;

  // > Connectivities
  int **cell_face_idx;
  int **cell_face;
  int **face_edge_idx;
  int **face_edge;
  int **face_vtx_idx;
  int **face_vtx;
  int **edge_vtx;

  // > Nodal connectivities

  int          *extract_n_vtx;
  double      **extract_vtx_coord;
  PDM_g_num_t **extract_vtx_gnum; // from initial mesh
  int         **extract_vtx_lnum; // from initial mesh

  int          *extract_n_tri;
  int         **extract_tri_vtx;
  PDM_g_num_t **extract_tri_gnum; // from initial mesh
  int         **extract_tri_lnum; // from initial mesh
  int          *extract_tri_n_group; // from initial mesh
  int         **extract_tri_tag; // from initial mesh

  int          *extract_n_tetra;
  int         **extract_tetra_vtx;
  PDM_g_num_t **extract_tetra_gnum; // from initial mesh
  int         **extract_tetra_lnum; // from initial mesh

  // > Boundaries
  int          *n_group_face;
  int         **group_face_idx;
  int         **group_face;
  PDM_g_num_t **group_face_gnum;

  // > Field
  double ***field;
  double ***gradient;



  // =========================
  // > Partitioned output data

  int iso_mesh_dimension;

  // > Isosurface entities
  int            **iso_n_entity          [PDM_MESH_ENTITY_MAX];
  PDM_g_num_t   ***iso_entity_gnum       [PDM_MESH_ENTITY_MAX];
  int           ***iso_entity_parent_idx [PDM_MESH_ENTITY_MAX];
  int           ***iso_entity_parent_lnum[PDM_MESH_ENTITY_MAX];
  PDM_g_num_t   ***iso_entity_parent_gnum[PDM_MESH_ENTITY_MAX];
  int           ***isovalue_entity_idx   [PDM_MESH_ENTITY_MAX];
  int           ***iso_connec_idx        [PDM_CONNECTIVITY_TYPE_MAX];
  int           ***iso_connec            [PDM_CONNECTIVITY_TYPE_MAX];

  // > Vertices
  double         ***iso_vtx_coord;
  double         ***iso_vtx_parent_weight;

  // > Edges
  int             *iso_n_edge_group;
  int           ***iso_edge_group_idx;
  int           ***iso_edge_group_lnum;
  PDM_g_num_t   ***iso_edge_group_gnum;

  // > Part_to_part between iso entities and entry mesh entities
  PDM_part_to_part_t **iso_ptp[PDM_MESH_ENTITY_MAX];

  // > Owners
  PDM_ownership_t **iso_owner_vtx_coord;
  PDM_ownership_t **iso_owner_vtx_parent_weight;
  PDM_ownership_t **iso_owner_edge_bnd;
  PDM_ownership_t **iso_owner_gnum               [PDM_MESH_ENTITY_MAX];
  PDM_ownership_t **iso_owner_parent_lnum        [PDM_MESH_ENTITY_MAX];
  PDM_ownership_t **iso_owner_parent_idx         [PDM_MESH_ENTITY_MAX];
  PDM_ownership_t **iso_owner_isovalue_entity_idx[PDM_MESH_ENTITY_MAX];
  PDM_ownership_t  *iso_owner_ptp                [PDM_MESH_ENTITY_MAX];
  PDM_ownership_t **iso_owner_connec             [PDM_CONNECTIVITY_TYPE_MAX];


  // =========================
  // > Distributed output data

  int            *iso_dn_entity          [PDM_MESH_ENTITY_MAX];
  int           **iso_dentity_parent_idx [PDM_MESH_ENTITY_MAX];
  PDM_g_num_t   **iso_dentity_parent_gnum[PDM_MESH_ENTITY_MAX];
  int           **iso_dconnec_idx        [PDM_CONNECTIVITY_TYPE_MAX];
  PDM_g_num_t   **iso_dconnec            [PDM_CONNECTIVITY_TYPE_MAX];

  // > Vertices
  double         **iso_dvtx_coord;
  double         **iso_dvtx_parent_weight;

  // > Edges
  int            **iso_dedge_group_idx;
  PDM_g_num_t    **iso_dedge_group_gnum;

  // > Owners
  PDM_ownership_t  *iso_owner_dvtx_coord;
  PDM_ownership_t  *iso_owner_dvtx_parent_weight;
  PDM_ownership_t  *iso_owner_dconnec    [PDM_CONNECTIVITY_TYPE_MAX];
  PDM_ownership_t  *iso_owner_dparent_idx[PDM_MESH_ENTITY_MAX];
  PDM_ownership_t  *iso_owner_dparent    [PDM_MESH_ENTITY_MAX];
  PDM_ownership_t  *iso_owner_dedge_bnd;


  // ===============
  // > Internal data

  PDM_bool_t *is_computed;


  // ========
  // > Timers
  double                  timer1, timer2; // donner des noms de ce qu'on mesure

};


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/



#ifdef  __cplusplus
}
#endif

#endif // PDM_ISOSURFACE_PRIV_H
