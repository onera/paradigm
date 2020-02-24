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
#include "pdm_handles.h"
#include "pdm_dmesh.h"
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
  int               *dmeshesIds;     /*!< Ids of distributed blocks (size = n_block)  */
  int               *partIds;         /*!< Ids of partitions built on each block of this
                                           process (size = n_block)                    */
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

  _multipart->dmeshesIds = (int *) malloc(_multipart->n_block * sizeof(int));
  for (int i = 0; i < _multipart->n_block; i++)
    _multipart->dmeshesIds[i] = -1;

  _multipart->partIds = (int *) malloc(_multipart->n_block * sizeof(int));
  for (int iblock = 0; iblock < _multipart->n_block; iblock++)
    _multipart->partIds[iblock] = -1;

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
void PDM_multipart_register_block
(
 const int        mpart_id,
 const int        block_id,
 const int        block_data_id
)
{
  PDM_printf("In multipart %d, set block n°%d using blockdata %d \n",
             mpart_id, block_id, block_data_id);

  _pdm_multipart_t *_multipart = _get_from_id (mpart_id);
  assert(block_id < _multipart->n_block);
  _multipart->dmeshesIds[block_id] = block_data_id;
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
      int blockId = _multipart->dmeshesIds[iblock];
      PDM_printf("block id for block %d is %d\n", iblock, blockId);
      int dNCell;
      int dNFace;
      int dNVtx;
      int dNBounds;
      double       *dVtxCoord;
      int          *dFaceVtxIdx;
      PDM_g_num_t  *dFaceVtx;
      PDM_g_num_t  *dFaceCell;
      int          *dFaceGroupIdx;
      PDM_g_num_t  *dFaceGroup;

      PDM_dmesh_dims_get(blockId, &dNCell, &dNFace, &dNVtx, &dNBounds);
      PDM_dmesh_data_get(blockId, &dVtxCoord, &dFaceVtxIdx, &dFaceVtx, &dFaceCell, &dFaceGroupIdx, &dFaceGroup);


      int ppartId = 0;
      int have_dCellPart = 0;
      int *dCellPart = (int *) malloc(dNCell*sizeof(int));
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
              dNCell,
              dNFace,
              dNVtx,
              dNBounds,
              NULL,                       // dCellFaceIdx
              NULL,                       // dCellFace
              NULL,                       // dCellTag
              NULL,                       // dCellWeight
              have_dCellPart,
              dCellPart,                  // dCellPart
              dFaceCell,
              dFaceVtxIdx,
              dFaceVtx,
              NULL,                       // dFaceTag
              dVtxCoord,
              NULL,                       // dVtxTag
              dFaceGroupIdx,
              dFaceGroup);
      PDM_printf("Partitionning done, ppardId is %d \n", ppartId);
      //Store the partition id for future access
      _multipart->partIds[iblock] = ppartId;

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

  assert(iblock < _multipart->n_block && ipart < _multipart->n_part);
  int ppartId = _multipart->partIds[iblock];

  PDM_part_part_dim_get(ppartId,
                        ipart,
                        nCell,
                        nFace,
                        nFacePartBound,
                        nVtx,
                        nProc,
                        nTPart,
                        sCellFace,
                        sFaceVtx,
                        sFaceGroup,
                        nFaceGroup);

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

  assert(iblock < _multipart->n_block && ipart < _multipart->n_part);
  int ppartId = _multipart->partIds[iblock];

  PDM_part_part_val_get(ppartId,
                        ipart,
                        cellTag,
                        cellFaceIdx,
                        cellFace,
                        cellLNToGN,
                        faceTag,
                        faceCell,
                        faceVtxIdx,
                        faceVtx,
                        faceLNToGN,
                        facePartBoundProcIdx,
                        facePartBoundPartIdx,
                        facePartBound,
                        vtxTag,
                        vtx,
                        vtxLNToGN,
                        faceGroupIdx,
                        faceGroup,
                        faceGroupLNToGN
                        );

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

  free(_multipart->dmeshesIds);

  for (int iblock = 0; iblock<_multipart->n_block; iblock++)
    PDM_part_free(_multipart->partIds[iblock]);
  free(_multipart->partIds);

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
