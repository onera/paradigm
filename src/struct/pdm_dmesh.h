/*
 * \file
 */

#ifndef __PDM_DMESH_H__
#define __PDM_DMESH_H__

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2017       ONERA

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

typedef struct _pdm_dmesh_t PDM_dmesh_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Build a distributed mesh structure
 *
 * \param [in]   dn_cell             Number of distributed cells
 * \param [in]   dn_face             Number of distributed faces
 * \param [in]   dn_vtx              Number of distributed vertices
 * \param [in]   dn_bnd              Number of boundaries
 * \param [in]   n_join              Number of interfaces with other zones
 *
 * \return     Identifier
 */

PDM_dmesh_t*
PDM_dmesh_create
(
       PDM_ownership_t owner,
 const int             dn_cell,
 const int             dn_face,
 const int             dn_edge,
 const int             dn_vtx,
       PDM_MPI_Comm    comm
);

/**
 *
 * \brief Get the dimensions of the distributed mesh
 *
 * \param [in]    id                id of the requested dmesh
 * \param [out]   dn_cell            Number of distributed cells
 * \param [out]   dn_face            Number of distributed faces
 * \param [out]   dn_vtx             Number of distributed vertices
 * \param [out]   dn_bnd             Number of boundaries
 * \param [out]   n_join             Number of interfaces with other zones
 */

void
PDM_dmesh_dims_get
(
 PDM_dmesh_t *dmeshm,
 int         *dn_cell,
 int         *dn_face,
 int         *dn_edge,
 int         *dn_vtx
);

/**
 *
 * \brief Get the dimensions of the distributed mesh
 *
 * \param [in]    id                id of the dmesh requested
 * \param [in]    entity_type       Kind of entity
 * \param [in]    dn_cell           Number of distributed cells
 */
void
PDM_dmesh_dn_entity_set
(
 PDM_dmesh_t         *dmesh,
 PDM_mesh_entities_t  entity_type,
 int                  dn_entity
);

/**
 *
 * \brief Get the dimensions of the distributed mesh
 *
 * \param [in]    id                id of the dmesh requested
 * \param [in]    entity_type       Kind of entity
 */
int
PDM_dmesh_dn_entity_get
(
 PDM_dmesh_t         *dmesh,
 PDM_mesh_entities_t  entity_type
);

void
PDM_dmesh_vtx_coord_get
(
 PDM_dmesh_t      *dmesh,
 double          **dvtx_coord,
 PDM_ownership_t   ownership
);

void
PDM_dmesh_vtx_coord_set
(
 PDM_dmesh_t      *dmesh,
 double           *dvtx_coord,
 PDM_ownership_t   ownership
);

void
PDM_dmesh_connectivity_set
(
 PDM_dmesh_t              *dmesh,
 PDM_connectivity_type_t   connectivity_type,
 PDM_g_num_t              *connect,
 int                      *connect_idx,
 PDM_ownership_t           ownership
);

int
PDM_dmesh_connectivity_get
(
 PDM_dmesh_t              *dmesh,
 PDM_connectivity_type_t   connectivity_type,
 PDM_g_num_t             **connect,
 int                     **connect_idx,
 PDM_ownership_t           ownership
);

int
PDM_dmesh_bound_get
(
 PDM_dmesh_t       *dmesh,
 PDM_bound_type_t   bound_type,
 PDM_g_num_t      **connect,
 int              **connect_idx,
 PDM_ownership_t    ownership
);

int
PDM_dmesh_distrib_get
(
 PDM_dmesh_t              *dmesh,
 PDM_mesh_entities_t       entity_type,
 PDM_g_num_t             **distrib
);

/**
 *
 * \brief Free
 *
 * \param [in]   id           Identifier
 *
 */

void
PDM_dmesh_free
(
 PDM_dmesh_t        *dmesh
);


const double *
PDM_dmesh_global_extents_get
(
 PDM_dmesh_t         *dmesh
 );


void
PDM_dmesh_bound_set
(
 PDM_dmesh_t      *dmesh,
 PDM_bound_type_t  bound_type,
 int               n_bound,
 PDM_g_num_t      *connect,
 int              *connect_idx,
 PDM_ownership_t   ownership
);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DMESH_H__ */
