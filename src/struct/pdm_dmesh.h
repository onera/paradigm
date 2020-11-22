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

typedef enum {
  PDM_CONNECTIVITY_TYPE_CELL_ELMT   = 0,
  PDM_CONNECTIVITY_TYPE_CELL_CELL   = 1,
  PDM_CONNECTIVITY_TYPE_CELL_FACE   = 2,
  PDM_CONNECTIVITY_TYPE_CELL_EDGE   = 3,
  PDM_CONNECTIVITY_TYPE_CELL_VTX    = 4,
  PDM_CONNECTIVITY_TYPE_FACE_ELMT   = 5,
  PDM_CONNECTIVITY_TYPE_FACE_CELL   = 6,
  PDM_CONNECTIVITY_TYPE_FACE_FACE   = 7,
  PDM_CONNECTIVITY_TYPE_FACE_EDGE   = 8,
  PDM_CONNECTIVITY_TYPE_FACE_VTX    = 9,
  PDM_CONNECTIVITY_TYPE_EDGE_ELMT   = 10,
  PDM_CONNECTIVITY_TYPE_EDGE_CELL   = 11,
  PDM_CONNECTIVITY_TYPE_EDGE_FACE   = 12,
  PDM_CONNECTIVITY_TYPE_EDGE_EDGE   = 13,
  PDM_CONNECTIVITY_TYPE_EDGE_VTX    = 14,
  PDM_CONNECTIVITY_TYPE_VTX_ELMT    = 15,
  PDM_CONNECTIVITY_TYPE_VTX_CELL    = 16,
  PDM_CONNECTIVITY_TYPE_VTX_FACE    = 17,
  PDM_CONNECTIVITY_TYPE_VTX_EDGE    = 18,
  PDM_CONNECTIVITY_TYPE_VTX_VTX     = 19,
  PDM_CONNECTIVITY_TYPE_ELMT_CELL   = 20,
  PDM_CONNECTIVITY_TYPE_ELMT_FACE   = 21,
  PDM_CONNECTIVITY_TYPE_ELMT_EDGE   = 22,
  PDM_CONNECTIVITY_TYPE_ELMT_VTX    = 23,
  PDM_CONNECTIVITY_TYPE_MAX         = 24
} PDM_connectivity_type_t;

// Reflexion similaire pour les group


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
 const int             n_bnd,
 const int             n_join
);

/**
 *
 * \brief Set the arrays into the distributed mesh structure
 *
 * \param [in]   id                 id of the dmesh to be set
 * \param [in]   dvtx_coord          Coordinates of  vertices (size = 3 * dn_vtx)
 * \param [in]   dface_vtx_idx        Face-vertex connectivity index of
 *                                    faces (size = dn_face + 1)
 * \param [in]   dface_vtx           Face-vertex connectivity of faces
 *                                    (size = dface_vtx_idx[dn_face])
 * \param [in]   dface_cell          Face-cell connectivity of faces (size =
 *                                    2 * dn_face). If iface is a boundary face,
 *                                    dface_cell[2*iface + 1] = 0
 * \param [in]   dface_bound_idx      Index of faces list of each boundary
 *                                    (size = dn_bnd + 1)
 * \param [in]   dface_bound         Faces list of each boundary
 *                                    (size = dface_bound_idx[dn_bnd])
 * \param [in]   joins_glob_id       Global id of each join (size = n_join)
 * \param [in]   dface_join_idx       Index of faces list of each join
 *                                    (size = n_join + 1)
 * \param [in]   dface_join          Faces list of each join
 *                                    (size = dface_join_idx[n_join])
 */

void
PDM_dmesh_set
(
 PDM_dmesh_t        *dmeshm,
 const double       *dvtx_coord,
 const int          *dface_vtx_idx,
 const PDM_g_num_t  *dface_vtx,
 const PDM_g_num_t  *dface_cell,
 const int          *dface_bound_idx,
 const PDM_g_num_t  *dface_bound,
 const int          *joins_glob_id,
 const int          *dface_join_idx,
 const PDM_g_num_t  *dface_join
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
 int         *dn_vtx,
 int         *dn_bnd,
 int         *n_joins
);

/**
 *
 * \brief Get the data (arrays) of the distributed mesh
 *
 * \param [in]    id                 id of the requested dmesh
 * \param [out]   dvtx_coord          Coordinates of  vertices
 * \param [out]   dface_vtx_idx        Face-vertex connectivity indices
 * \param [out]   dface_vtx           Face-vertex connectivity
 * \param [out]   dface_cell          Face-cell connectivity of faces
 * \param [out]   dface_bound_idx      Indices of faces list of each boundary
 * \param [out]   dface_bound         Faces list of each boundary
 * \param [out]   joins_glob_id       Global Id of each join
 * \param [out]   dface_join_idx       Indices of faces list of each join
 * \param [out]   dface_join          Faces list of each join
 */

void
PDM_dmesh_data_get
(
 PDM_dmesh_t         *dmeshm,
 const double       **dvtx_coord,
 const int          **dface_vtx_idx,
 const PDM_g_num_t  **dface_vtx,
 const PDM_g_num_t  **dface_cell,
 const int          **dface_bound_idx,
 const PDM_g_num_t  **dface_bound,
 const int          **joins_glob_id,
 const int          **dface_join_idx,
 const PDM_g_num_t  **dface_join
);

void
PDM_dmesh_connectivity_get
(
 PDM_dmesh_t              *dmesh,
 PDM_connectivity_type_t   connectivity_type,
 PDM_g_num_t             **connect,
 int                     **connect_idx,
 PDM_ownership_t           ownership
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


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DMESH_H__ */
