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
 * \param [in]   dNCell             Number of distributed cells
 * \param [in]   dNFace             Number of distributed faces
 * \param [in]   dNVtx              Number of distributed vertices
 * \param [in]   dNBnd              Number of boundaries
 * \param [in]   dNJoin             Number of interfaces with other zones
 *
 * \return     Identifier
 */

int
PDM_dmesh_create
(
 const int          dNCell,
 const int          dNFace,
 const int          dNVtx,
 const int          dNBnd,
 const int          dNJoin
);

/**
 *
 * \brief Set the arrays into the distributed mesh structure
 *
 * \param [in]   id                 id of the dmesh to be set
 * \param [in]   dVtxCoord          Coordinates of  vertices (size = 3 * dNVtx)
 * \param [in]   dFaceVtxIdx        Face-vertex connectivity index of
 *                                    faces (size = dNFace + 1)
 * \param [in]   dFaceVtx           Face-vertex connectivity of faces
 *                                    (size = dFaceVtxIdx[dNFace])
 * \param [in]   dFaceCell          Face-cell connectivity of faces (size =
 *                                    2 * dNFace). If iface is a boundary face,
 *                                    dFaceCell[2*iface + 1] = 0
 * \param [in]   dFaceBoundIdx      Index of faces list of each boundary
 *                                    (size = dNBnd + 1)
 * \param [in]   dFaceBound         Faces list of each boundary
 *                                    (size = dfaceBoundIdx[dNBnd])
 * \param [in]   dJoinZoneOpp       Global Id of the opposed zone for each join
 *                                    (size = dNJoin)
 * \param [in]   dFaceJoinIdx       Index of faces list of each join
 *                                    (size = dNJoin + 1)
 * \param [in]   dFaceJoin          Faces list of each join
 *                                    (size = dfaceJoinIdx[dNJoin])
 * \param [in]   dFaceTag           Distributed face tag (size = dNFace)
 */

void
PDM_dmesh_set
(
 const int           id,
 const double       *dVtxCoord,
 const int          *dFaceVtxIdx,
 const PDM_g_num_t  *dFaceVtx,
 const PDM_g_num_t  *dFaceCell,
 const int          *dFaceBoundIdx,
 const PDM_g_num_t  *dFaceBound,
 const int          *dJoinZoneOpp,
 const int          *dFaceJoinIdx,
 const PDM_g_num_t  *dFaceJoin,
 const int          *dFaceTag
);

/**
 *
 * \brief Get the dimensions of the distributed mesh
 *
 * \param [in]    id                id of the requested dmesh
 * \param [out]   dNCell            Number of distributed cells
 * \param [out]   dNFace            Number of distributed faces
 * \param [out]   dNVtx             Number of distributed vertices
 * \param [out]   dNBnd             Number of boundaries
 * \param [out]   dNJoin            Number of interfaces with other zones
 */

void
PDM_dmesh_dims_get
(
 const int   id,
 int        *dNCell,
 int        *dNFace,
 int        *dNVtx,
 int        *dNBnd,
 int        *dNJoins
);

/**
 *
 * \brief Get the data (arrays) of the distributed mesh
 *
 * \param [in]    id                 id of the requested dmesh
 * \param [out]   dVtxCoord          Coordinates of  vertices
 * \param [out]   dFaceVtxIdx        Face-vertex connectivity indices
 * \param [out]   dFaceVtx           Face-vertex connectivity
 * \param [out]   dFaceCell          Face-cell connectivity of faces
 * \param [out]   dFaceBoundIdx      Indices of faces list of each boundary
 * \param [out]   dFaceBound         Faces list of each boundary
 * \param [out]   dJoinZoneOpp       Global Id of the opposed zone for each join
 * \param [out]   dFaceJoinIdx       Indices of faces list of each join
 * \param [out]   dFaceJoin          Faces list of each join
 * \param [out]   dFaceTag           Faces tags
 */

void
PDM_dmesh_data_get
(
 const int      id,
 double       **dVtxCoord,
 int          **dFaceVtxIdx,
 PDM_g_num_t  **dFaceVtx,
 PDM_g_num_t  **dFaceCell,
 int          **dFaceBoundIdx,
 PDM_g_num_t  **dFaceBound,
 int          **dJoinZoneOpp,
 int          **dFaceJoinIdx,
 PDM_g_num_t  **dFaceJoin,
 int          **dFaceTag
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
 const int id
);


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DMESH_H__ */
