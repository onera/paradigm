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

/*============================================================================
 * Type definitions
 *============================================================================*/


struct _pdm_isosurface_t {

  // ========================
  // > General
  PDM_MPI_Comm comm;


  // ========================
  // > Entry mesh information

  // > Isosurface switch
  // TODO: is entry_mesh_type sign necessary ???
  int is_dist_or_part; // -1: undef, 0: dist, 1: part
  int entry_mesh_type; //  0: undef, 1: dist_alamano, 2: dmesh, 3: dmesh_nodal, -1: part_alamano, -2: pmesh, -3: pmesh_nodal

  // > Mesh information
  int entry_mesh_dim;



  // ====================
  // > Isosurface options

  // > Isosurfaces
  int                     n_isosurface;

  // > Isosurface type
  int                     mesh_dimension;
  PDM_iso_surface_kind_t *kind;

  // > Isovalues
  int                    *n_isovalues;
  double                **isovalues;

  // > Equation args
  PDM_isosurface_field_function_t *field_function;
  double                         **eq_coeffs;
  int                             *use_gradient;

  // > Function args
  _pdm_isosurface_field_function_t *iso_func;

  // > Redistribution
  PDM_extract_part_kind_t extract_kind;
  PDM_split_dual_t        part_method;

  // > Link with entry mesh
  int **compute_ptp;
  
  // ========================
  // > Distributed entry data

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

  // > Boundaries
  int          n_dgroup_face;
  int         *dgroup_face_idx;
  PDM_g_num_t *dgroup_face;

  // int          n_dgroup_edge;
  // int         *dgroup_edge_idx;
  // PDM_g_num_t *dgroup_edge;
  
  // int          n_dgroup_vtx;
  // int         *dgroup_vtx_idx;
  // PDM_g_num_t *dgroup_vtx;

  // > Field
  double **dfield;
  double **dgradient;



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

  // > Boundaries
  int  *n_group_face;
  int **group_face_idx;
  int **group_face;

  // int  *n_group_edge;
  // int **group_edge_idx;
  // int **group_edge;
  
  // int  *n_group_vtx;
  // int **group_vtx_idx;
  // int **group_vtx;

  // > Field
  double ***field;
  double ***gradient;



  // =========================
  // > Partitioned output data

  int iso_mesh_dimension;
  int iso_n_part;

  // > Vertices
  int             **iso_n_vtx;
  double         ***iso_vtx_coord;
  PDM_g_num_t    ***iso_vtx_gnum;
  int            ***iso_vtx_lparent_idx;
  int            ***iso_vtx_lparent;
  double         ***iso_vtx_parent_weight;
  int            ***isovalue_vtx_idx;

  // > Edges
  int            **iso_n_edge;
  int           ***iso_edge_vtx;
  PDM_g_num_t   ***iso_edge_gnum;
  int           ***iso_edge_lparent_idx;
  int           ***iso_edge_lparent;
  int             *iso_n_edge_group;
  int           ***iso_edge_group_idx;
  int           ***iso_edge_group_lnum;
  PDM_g_num_t   ***iso_edge_group_gnum;
  int           ***isovalue_edge_idx;

  // > Faces
  int            **iso_n_face;
  int           ***iso_face_vtx_idx;
  int           ***iso_face_vtx;
  PDM_g_num_t   ***iso_face_gnum;
  int           ***iso_face_lparent_idx;
  int           ***iso_face_lparent;
  int           ***isovalue_face_idx;

  // > Part_to_part between iso entities and entry mesh entites
  PDM_part_to_part_t **iso_ptp_vtx;
  PDM_part_to_part_t **iso_ptp_edge;
  PDM_part_to_part_t **iso_ptp_face;

  // > Owners
  PDM_ownership_t  **iso_owner_vtx_coord;
  PDM_ownership_t  **iso_owner_vtx_parent_weight;
  PDM_ownership_t ***iso_owner_gnum;
  PDM_ownership_t ***iso_owner_connec;
  PDM_ownership_t ***iso_owner_lparent;
  PDM_ownership_t  **iso_owner_edge_bnd;
  PDM_ownership_t  **iso_owner_ptp;


  // ===============
  // > Internal data




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
