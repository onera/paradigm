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
 * TODO : write module description here
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_part.h"
#include "pdm_part_priv.h"
#include "pdm_handles.h"
#include "pdm_printf.h"
#include "pdm_error.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_multipart.h"

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
 * \struct _pdm_blockdata_t
 * \brief  This structure makes the link to the distributed data
 *         known by the process for a given block.
 *         Note that arrays are shared references : the structure does not hold
 *         the memory.
 *
 */

typedef struct
{
  int               dNCell;          /*!< Number of distributed cells         */
  int               dNFace;          /*!< Number of distributed faces         */
  int               dNVtx;           /*!< Number of distributed vertices      */
  int               nFaceGroup;      /*!< Number of boundaries                */
  const PDM_g_num_t *_dFaceCell;     /*!< Face-cell connectivity of distributed
                                        faces (size = 2 * dNFace, shared array)
                                        if iface is a boundary face,
                                        _dFaceCell[2*iface + 1] = 0           */
  const int         *_dFaceVtxIdx;   /*!< Face-vertex connectivity index of
                                        distributed faces (size = dNFace + 1,
                                        shared array)                         */
  const PDM_g_num_t *_dFaceVtx;      /*!< Face-vertex connectivity of
                                        distributed faces (size = dFaceVtxIdx[
                                        dNFace],shared array)                 */
  const double      *_dVtxCoord;     /*!< Coordinates of ditributed vertices
                                        (size = 3 * dNVtx, shared array)      */
  const int         *_dFaceGroupIdx; /*!< Index of distributed faces list of
                                        each boundary (size = nBound + 1)
                                        or NULL                               */
  const PDM_g_num_t *_dFaceGroup;    /*!< Distributed faces list of each
                                       boundary (size = dfaceBoundIdx[nBound])
                                        or NULL                               */
} _pdm_blockdata_t;

/**
 * \struct _pdm_multipart_t
 * \brief  This structure describe a multipart. In addition to splitting
 *         parameters, it stores the multiples blocks and part as well as
 *         the global numbering.
 *
 */

typedef struct  {

  int               n_block;          /*!< Number of blocks */
  int               n_part;           /*!< Number of partitions per proc in each block */
  PDM_bool_t        merge_blocks;     /*!< Merge before partitionning or not */
  PDM_part_split_t  split_method;     /*!< Partitioning method */
  PDM_MPI_Comm      comm;             /*!< MPI communicator */
  _pdm_blockdata_t  **meshBlocks;     /*!< Blocks (size = n_block) */
  _part_t           **meshParts;      /*!< Partitions built on this process (size = ?) */

} _pdm_multipart_t;

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_multiparts   = NULL;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Return multipart object from it identifier
 *
 * \param [in]   multipartId    multipart identifier
 *
 */

static _pdm_multipart_t *
_get_from_id
(
 int  id
)
{

  _pdm_multipart_t *multipart = (_pdm_multipart_t *) PDM_Handles_get (_multiparts, id);

  if (multipart == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_multipart error : Bad identifier\n");
  }

  return multipart;
}

static inline _pdm_blockdata_t*
_block_create
(
 void
)
{
  _pdm_blockdata_t *block = (_pdm_blockdata_t *) malloc(sizeof(_pdm_blockdata_t));
  block->dNCell = 0;
  block->dNFace = 0;
  block->dNVtx = 0;
  block->nFaceGroup = 0;
  block->_dFaceCell = NULL;
  block->_dFaceVtxIdx = NULL;
  block->_dFaceVtx = NULL;
  block->_dVtxCoord = NULL;
  block->_dFaceGroupIdx = NULL;
  block->_dFaceGroup = NULL;
  return block;
}
static void
_block_free
(
 _pdm_blockdata_t *block
)
{
  free(block);
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Build a multipart structure
 *
 * \param [in]   n_block      Number of blocks in the original mesh
 * \param [in]   n_part       Number of partition per proc in each block
 * \param [in]   merge_blocks Merge or not the blocks before splitting
 * \param [in]   split_method Choice of library used to split the mesh
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Identifier
 */

int
PDM_multipart_create
(
 const int              n_block,
 const int              n_part,
 const PDM_bool_t       merge_blocks,
 const PDM_part_split_t split_method,
 const PDM_MPI_Comm     comm
)
{

  /*
   * Search a ppart free id
   */

  if (_multiparts == NULL) {
    _multiparts = PDM_Handles_create (4);
  }

  _pdm_multipart_t *_multipart = (_pdm_multipart_t *) malloc(sizeof(_pdm_multipart_t));
  int id = PDM_Handles_store (_multiparts, _multipart);

  _multipart->n_block     = n_block;
  _multipart->n_part      = n_part;
  _multipart->merge_blocks= merge_blocks;
  _multipart->split_method= split_method;
  _multipart->comm        = comm;

  _multipart->meshBlocks = (_pdm_blockdata_t **) malloc(_multipart->n_block * sizeof(_pdm_blockdata_t *));
  for (int i = 0; i < _multipart->n_block; i++)
    _multipart->meshBlocks[i] = NULL;

  int totalPartNumber = (_multipart->n_block)*(_multipart->n_part);
  _multipart->meshParts = (_part_t **) malloc(totalPartNumber * sizeof(_part_t *));
  for (int i = 0; i < totalPartNumber; i++)
    _multipart->meshParts[i] = NULL;

  PDM_printf("Created from PDM_multipart_create. You requested a multipart with %d blocks \n", n_block);
  return id;

}

// void
// PROCF (pdm_multipart_create, PDM_MULTIPART_CREATE)
// (
//  const int          *n_block,
//  const PDM_MPI_Fint *fcomm,
//        int          *id
// )
// {
//   const PDM_MPI_Comm c_comm = PDM_MPI_Comm_f2c (*fcomm);

//   *id = PDM_multipart_create (*n_block, c_comm);
// }

/* TODO : copy doc of the function */

 void PDM_multipart_set_block
(
  const int                   id,
  const int                   i_block,
  const int                   dNCell,
  const int                   dNFace,
  const int                   dNVtx,
  const int                   nFaceGroup,
  const PDM_g_num_t          *dFaceCell,
  const int                  *dFaceVtxIdx,
  const PDM_g_num_t          *dFaceVtx,
  const double               *dVtxCoord,
  const int                  *dFaceGroupIdx,
  const PDM_g_num_t          *dFaceGroup
)
{
  PDM_printf("Set  block n° %d \n", i_block);
  _pdm_multipart_t *_multipart = _get_from_id (id);

  assert(_multipart->meshBlocks[i_block] == NULL);
  _multipart->meshBlocks[i_block] = _block_create();
  _pdm_blockdata_t *block =  _multipart->meshBlocks[i_block];

  block->dNCell = dNCell;
  block->dNFace = dNFace;
  block->dNVtx = dNVtx;
  block->nFaceGroup = nFaceGroup;

  block->_dFaceCell = dFaceCell;
  block->_dFaceVtxIdx = dFaceVtxIdx;
  block->_dFaceVtx = dFaceVtx;
  block->_dVtxCoord = dVtxCoord;
  block->_dFaceGroupIdx = dFaceGroupIdx;
  block->_dFaceGroup = dFaceGroup;
  PDM_printf("Block n° %d filled\n", i_block);
}

void
PDM_multipart_run_ppart
(
 const int id
)
{
  _pdm_multipart_t *_multipart = _get_from_id (id);

  if (_multipart->merge_blocks)
  {
    // 1. Generate global numerotation using all blocks
    // 2. Call the partitionner once on the global numbering
  }
  else
  {
    // 2. Loop over the blocks and call the partitionner
    for (int iblock = 0; iblock<_multipart->n_block; iblock++)
    {
      PDM_printf("You requested no merge : partitionning block %d/%d \n", iblock+1, _multipart->n_block);
      _pdm_blockdata_t *block =  _multipart->meshBlocks[iblock];
      int ppartId = 0;
      int have_dCellPart = 0;
      int *dCellPart = (int *) malloc(block->dNCell*sizeof(int));
      // We probably should create a subcomm to imply only the procs sharing this block
      PDM_part_create(&ppartId,
              _multipart->comm,
              _multipart->split_method,
              "PDM_PART_RENUM_CELL_NONE",
              "PDM_PART_RENUM_FACE_NONE",
              0,                          // nPropertyCell
              NULL,                       // renum_properties_cell
              0,                          // nPropertyFace
              NULL,                       // renum_properties_face
              _multipart->n_part,
              block->dNCell,
              block->dNFace,
              block->dNVtx,
              block->nFaceGroup,
              NULL,                       // dCellFaceIdx
              NULL,                       // dCellFace
              NULL,                       // dCellTag
              NULL,                       // dCellWeight
              have_dCellPart,
              dCellPart,                  // dCellPart
              block->_dFaceCell,
              block->_dFaceVtxIdx,
              block->_dFaceVtx,
              NULL,                       // dFaceTag
              block->_dVtxCoord,
              NULL,                       // dVtxTag
              block->_dFaceGroupIdx,
              block->_dFaceGroup);
      PDM_printf("New call to partitionner, ppardId is %d \n", ppartId);


      // Store partitions in array
      for (int ipart = 0; ipart < _multipart->n_part; ipart++) {

        assert(_multipart->meshParts[iblock+ipart] == NULL);
        _multipart->meshParts[iblock+ipart] = _part_create();
        _part_t *part =  _multipart->meshParts[iblock+ipart];

        int nProc;
        int nTPart;
        int sCellFace;
        int sFaceVtx;
        int sFaceGroup;

        PDM_part_part_dim_get(ppartId,
                             ipart,
                             &(part->nCell),
                             &(part->nFace),
                             &(part->nFacePartBound),
                             &(part->nVtx),
                             &nProc,
                             &nTPart,
                             &sCellFace,
                             &sFaceVtx,
                             &sFaceGroup,
                             &(part->nFaceGroup));

        PDM_part_part_val_get(ppartId,
                             ipart,
                             &(part->cellTag),
                             &(part->cellFaceIdx),
                             &(part->cellFace),
                             &(part->cellLNToGN),
                             &(part->faceTag),
                             &(part->faceCell),
                             &(part->faceVtxIdx),
                             &(part->faceVtx),
                             &(part->faceLNToGN),
                             &(part->facePartBoundProcIdx),
                             &(part->facePartBoundPartIdx),
                             &(part->facePartBound),
                             &(part->vtxTag),
                             &(part->vtx),
                             &(part->vtxLNToGN),
                             &(part->faceGroupIdx),
                             &(part->faceGroup),
                             &(part->faceGroupLNToGN));

      PDM_printf("ncell is %d \n", part->nCell);
      PDM_printf("cellLNToGN from multipart");
      for (int k=0; k<part->nCell; k++)
        PDM_printf(" %d ", part->cellLNToGN[k]);
      PDM_printf("\n");

      }
      free(dCellPart);
    }
  }
  // 3. rebuild the joins
}

void
PDM_multipart_part_dim_get
(
const   int  mpartId,
const   int  iblock,
const   int  ipart,
 int        *nCell,
 int        *nFace,
 int        *nFacePartBound,
 int        *nVtx,
 int        *nProc,
 int        *nTPart,
 int        *sCellFace,
 int        *sFaceVtx,
 int        *sFaceGroup,
 int        *nFaceGroup
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpartId);
  int numProcs;
  PDM_MPI_Comm_size(_multipart->comm, &numProcs);

  _part_t *meshPart = NULL;
  if (iblock < _multipart->n_block && ipart < _multipart->n_part)
    meshPart = _multipart->meshParts[iblock+ipart];

  if (meshPart == NULL) {
    PDM_printf("PDM_part_part_get error : unknown partition\n");
    exit(1);
  }

  *nCell           = meshPart->nCell;
  *nFace           = meshPart->nFace;
  *nFacePartBound  = meshPart->nFacePartBound;
  *nProc           = numProcs;
  // *nTPart          = _multipart->tNPart;
  *nTPart          = 0;
  *nVtx            = meshPart->nVtx;
  *sCellFace       = meshPart->cellFaceIdx[*nCell];
  *sFaceVtx        = meshPart->faceVtxIdx[*nFace];
  *sFaceGroup      = 0;
  // if (ppart->nFaceGroup > 0)
    // *sFaceGroup    = meshPart->faceGroupIdx[ppart->nFaceGroup];
    // *nFaceGroup    = ppart->nFaceGroup;
}

void
PDM_multipart_part_val_get
(
const int            mpartId,
const int            iblock,
const int            ipart,
      int          **cellTag,
      int          **cellFaceIdx,
      int          **cellFace,
      PDM_g_num_t  **cellLNToGN,
      int          **faceTag,
      int          **faceCell,
      int          **faceVtxIdx,
      int          **faceVtx,
      PDM_g_num_t  **faceLNToGN,
      int          **facePartBoundProcIdx,
      int          **facePartBoundPartIdx,
      int          **facePartBound,
      int          **vtxTag,
      double       **vtx,
      PDM_g_num_t  **vtxLNToGN,
      int          **faceGroupIdx,
      int          **faceGroup,
      PDM_g_num_t  **faceGroupLNToGN
)
{
   _pdm_multipart_t *_multipart = _get_from_id (mpartId);

  _part_t *meshPart = NULL;
  if (iblock < _multipart->n_block && ipart < _multipart->n_part)
    meshPart = _multipart->meshParts[iblock+ipart];

  if (meshPart == NULL) {
    PDM_printf("PDM_part_part_val_get error : unknown partition\n");
    exit(1);
  }
  *cellTag              = meshPart->cellTag;
  *cellFaceIdx          = meshPart->cellFaceIdx;
  *cellFace             = meshPart->cellFace;
  *cellLNToGN           = meshPart->cellLNToGN;
  *faceTag              = meshPart->faceTag;
  *faceCell             = meshPart->faceCell;
  *faceVtxIdx           = meshPart->faceVtxIdx;
  *faceVtx              = meshPart->faceVtx;
  *faceLNToGN           = meshPart->faceLNToGN;
  *facePartBoundProcIdx = meshPart->facePartBoundProcIdx;
  *facePartBoundPartIdx = meshPart->facePartBoundPartIdx;
  *facePartBound        = meshPart->facePartBound;
  *vtxTag               = meshPart->vtxTag;
  *vtx                  = meshPart->vtx;
  *vtxLNToGN            = meshPart->vtxLNToGN;
  *faceGroupIdx         = meshPart->faceGroupIdx;
  *faceGroup            = meshPart->faceGroup;
  *faceGroupLNToGN      = meshPart->faceGroupLNToGN;
}

/**
 *
 * \brief Free
 *
 * \param [in]   id           Identifier
 *
 */

void
PDM_multipart_free
(
 const int id
)
{
  _pdm_multipart_t *_multipart = _get_from_id (id);

  for (int i = 0; i < _multipart->n_block; i++) {
    if (_multipart->meshBlocks[i] != NULL)
      _block_free(_multipart->meshBlocks[i]);
    _multipart->meshBlocks[i] = NULL;
  }
  if (_multipart->meshBlocks != NULL)
    free(_multipart->meshBlocks);
  _multipart->meshBlocks = NULL;

  int totalPartNumber = (_multipart->n_block)*(_multipart->n_part);
  for (int i = 0; i < totalPartNumber; i++) {
    if (_multipart->meshParts[i] != NULL)
    _multipart->meshParts[i] = NULL;
  }
  if (_multipart->meshParts != NULL)
    free(_multipart->meshParts);
  _multipart->meshParts = NULL;

  //TODO : remplacer par le vrai ID
  PDM_part_free(0);
  free (_multipart);

  PDM_Handles_handle_free (_multiparts, id, PDM_FALSE);

  const int n_multipart = PDM_Handles_n_get (_multiparts);

  if (n_multipart == 0) {
    _multiparts = PDM_Handles_free (_multiparts);
  }
  PDM_printf("Cleaned from PDM_multipart_free\n");
}

void
PROCF (pdm_multipart_free, PDM_MULTIPART_FREE)
(
 const int *id
)
{
  PDM_multipart_free (*id);
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
