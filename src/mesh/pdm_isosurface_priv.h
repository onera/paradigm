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


// void (*_pdm_isosurface_field_function_t)
// (
//  const double  x,
//  const double  y,
//  const double  z,
//  double       *value
// );


struct _pdm_isosurface_t {

  // ========================
  // > General
  PDM_MPI_Comm comm;


  // ========================
  // > Isosurface API

  // > Isosurface type
  int                     mesh_dimension;
  PDM_Mesh_nodal_elt_t    elt_type;
  PDM_iso_surface_kind_t  kind;
  int                     is_dist_or_part; // -1: undef, 0: dist, 1: part

  // > Isovalues
  int                     n_isovalues;
  double                 *isovalues;

  // > Equation args
  double                 *eq_coeffs;
  int                     use_gradient;

  // // > Function args
  // _pdm_isosurface_field_function_t *iso_func;

  
  // ========================
  // > Distributed entry data

  // > 

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
  int          n_group_face;
  int         *dgroup_face_idx;
  PDM_g_num_t *dgroup_face;

  int          n_group_edge;
  int         *dgroup_edge_idx;
  PDM_g_num_t *dgroup_edge;
  
  int          n_group_vtx;
  int         *dgroup_vtx_idx;
  PDM_g_num_t *dgroup_vtx;

  // > Field
  double **dfield;
  double **dgradient;



  // ========================
  // > Partitioned entry data

  // > Distribution

  // > Vertices

  // > Connectivities

  // > Boundaries


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