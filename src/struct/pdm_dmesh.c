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

/*============================================================================
 * Interface structure to represent a distributed mesh
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_handles.h"
#include "pdm_printf.h"
#include "pdm_error.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_dmesh.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local structure definitions
 *============================================================================*/


/**
 * \struct _pdm_dmesh_t
 * \brief  Define a distributed mesh. Arrays are shared: this structure does not
 *         holds the data.
 *
 */

typedef struct
{
  int               dNCell;          /*!< Number of distributed cells         */
  int               dNFace;          /*!< Number of distributed faces         */
  int               dNVtx;           /*!< Number of distributed vertices      */
  int               dNBnd;           /*!< Number of boundaries                */
  int               dNJoin;          /*!< Number of interfaces with other zone*/
  const PDM_g_num_t *_dFaceCell;     /*!< Face-cell connectivity of distributed
                                        faces (size = 2 * dNFace)
                                        if iface is a boundary face,
                                        _dFaceCell[2*iface + 1] = 0           */
  const int         *_dFaceVtxIdx;   /*!< Face-vertex connectivity index of
                                        distributed faces (size = dNFace + 1) */
  const PDM_g_num_t *_dFaceVtx;      /*!< Face-vertex connectivity of
                                        distributed faces (size = dFaceVtxIdx[
                                        dNFace])                              */
  const double      *_dVtxCoord;     /*!< Coordinates of ditributed vertices
                                        (size = 3 * dNVtx)                    */
  const int         *_dFaceBoundIdx; /*!< Index of distributed faces list of
                                        each boundary (size = dNBnd + 1)
                                        or NULL                               */
  const PDM_g_num_t *_dFaceBound;    /*!< Distributed faces list of each
                                       boundary (size = dfaceBoundIdx[dNBnd])
                                        or NULL                               */
  const int         *_dJoinGIds;     /*!< Tuple JoinGId, JoinGIdDonnor for
                                        each join (size = 2*dNJoin) or NULL   */
  const int         *_dFaceJoinIdx;  /*!< Index of distributed faces list of
                                        each join (size = dNJoin + 1)
                                        or NULL                               */
  const PDM_g_num_t *_dFaceJoin;     /*!< Distributed faces list of each
                                       join (size = dfaceJoinIdx[dNJoin])
                                        or NULL                               */
} _pdm_dmesh_t;

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_dmeshes   = NULL;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Return ppart object from it identifier
 *
 * \param [in]   ppartId        ppart identifier
 *
 */

static _pdm_dmesh_t *
_get_from_id
(
 int  id
)
{

  _pdm_dmesh_t *dmesh = (_pdm_dmesh_t *) PDM_Handles_get (_dmeshes, id);

  if (dmesh == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_dmesh error : Bad identifier\n");
  }

  return dmesh;
}


/*=============================================================================
 * Public function definitions
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
)
{

  /*
   * Search a ppart free id
   */
  if (_dmeshes == NULL) {
    _dmeshes = PDM_Handles_create (4);
  }

  _pdm_dmesh_t *dmesh = (_pdm_dmesh_t *) malloc(sizeof(_pdm_dmesh_t));
  int id = PDM_Handles_store (_dmeshes, dmesh);

  dmesh->dNCell         = dNCell;
  dmesh->dNFace         = dNFace;
  dmesh->dNVtx          = dNVtx;
  dmesh->dNBnd          = dNBnd;
  dmesh->dNJoin         = dNJoin;
  dmesh->_dFaceCell     = NULL;
  dmesh->_dFaceVtxIdx   = NULL;
  dmesh->_dFaceVtx      = NULL;
  dmesh->_dVtxCoord     = NULL;
  dmesh->_dFaceBoundIdx = NULL;
  dmesh->_dFaceBound    = NULL;
  dmesh->_dJoinGIds     = NULL;
  dmesh->_dFaceJoinIdx  = NULL;
  dmesh->_dFaceJoin     = NULL;

  return id;

}

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
 * \param [in]   dJoinGIds          Tuple JoinGId, JoinGIdDonnor for
 *                                    each join (size = 2*dNJoin)
 * \param [in]   dFaceJoinIdx       Index of faces list of each join
 *                                    (size = dNJoin + 1)
 * \param [in]   dFaceJoin          Faces list of each join
 *                                    (size = dfaceJoinIdx[dNJoin])
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
 const int          *dJoinGIds,
 const int          *dFaceJoinIdx,
 const PDM_g_num_t  *dFaceJoin
)
{
  _pdm_dmesh_t *dmesh = _get_from_id (id);

  dmesh->_dVtxCoord     = dVtxCoord;
  dmesh->_dFaceVtxIdx   = dFaceVtxIdx;
  dmesh->_dFaceVtx      = dFaceVtx;
  dmesh->_dFaceCell     = dFaceCell;
  dmesh->_dFaceBoundIdx = dFaceBoundIdx;
  dmesh->_dFaceBound    = dFaceBound;
  dmesh->_dJoinGIds     = dJoinGIds;
  dmesh->_dFaceJoinIdx  = dFaceJoinIdx;
  dmesh->_dFaceJoin     = dFaceJoin;
}

/**
 *
 * \brief Get the dimensions of the distributed mesh
 *
 * \param [in]    id                id of the dmesh requested
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
 int        *dNJoin
)
{
  _pdm_dmesh_t *dmesh = _get_from_id (id);
  *dNCell = dmesh->dNCell;
  *dNFace = dmesh->dNFace;
  *dNVtx = dmesh->dNVtx;
  *dNBnd = dmesh->dNBnd;
  *dNJoin = dmesh->dNJoin;
}

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
 * \param [out]   dJoinGIds          Global Ids of the join and opposed join
 * \param [out]   dFaceJoinIdx       Indices of faces list of each join
 * \param [out]   dFaceJoin          Faces list of each join
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
 int          **dJoinGIds,
 int          **dFaceJoinIdx,
 PDM_g_num_t  **dFaceJoin
)
{
  _pdm_dmesh_t *dmesh = _get_from_id (id);

  *dVtxCoord     = dmesh->_dVtxCoord;
  *dFaceVtxIdx   = dmesh->_dFaceVtxIdx;
  *dFaceVtx      = dmesh->_dFaceVtx;
  *dFaceCell     = dmesh->_dFaceCell;
  *dFaceBoundIdx = dmesh->_dFaceBoundIdx;
  *dFaceBound    = dmesh->_dFaceBound;
  *dJoinGIds     = dmesh->_dJoinGIds;
  *dFaceJoinIdx  = dmesh->_dFaceJoinIdx;
  *dFaceJoin     = dmesh->_dFaceJoin;
}

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
)
{
  _pdm_dmesh_t *dmesh = _get_from_id (id);

  dmesh->dNCell   = 0;
  dmesh->dNFace   = 0;
  dmesh->dNVtx    = 0;
  dmesh->dNBnd    = 0;
  dmesh->dNJoin   = 0;
  dmesh->_dFaceCell     = NULL;
  dmesh->_dFaceVtxIdx   = NULL;
  dmesh->_dFaceVtx      = NULL;
  dmesh->_dVtxCoord     = NULL;
  dmesh->_dFaceBoundIdx = NULL;
  dmesh->_dFaceBound    = NULL;
  dmesh->_dJoinGIds     = NULL;
  dmesh->_dFaceJoinIdx  = NULL;
  dmesh->_dFaceJoin     = NULL;

  free (dmesh);

  PDM_Handles_handle_free (_dmeshes, id, PDM_FALSE);

  const int n_dmesh = PDM_Handles_n_get (_dmeshes);

  if (n_dmesh == 0) {
    _dmeshes = PDM_Handles_free (_dmeshes);
  }

}


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
