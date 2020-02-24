#ifndef __PDM_MULTIPART_H__
#define __PDM_MULTIPART_H__

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
);

// void
// PROCF (pdm_multipart_create, PDM_MULTIPART_CREATE)
// (
//  const int *n_block,
//  const PDM_MPI_Fint *fcomm,
//        int *id
// );

/**
 *
 * \brief Set a block in the multipart structure
 *
 * \param [in]   id           Identifier
 * \param [in]   i_block      Number of block to set
 * \param [in]   dFaceCell    Face to cell connectivity for the block
 * TODO LIST PARAMS
 *
 */
void PDM_multipart_register_block
(
 const int        mpart_id,
 const int        block_id,
 const int        block_data_id
);

/**
 *
 * \brief Call the partitionner (via PDM_part_create) on the multipart object
 *
 * \param [in]   id           Identifier
 *
 */

void
PDM_multipart_run_ppart
(
 const int id
);

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
);

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
);

void
PDM_multipart_part_color_get
(
const int            mpartId,
const int            iblock,
const int            ipart,
      int          **cellColor,
      int          **faceColor,
      int          **threadColor,
      int          **hyperPlaneColor
);

void
PDM_multipart_time_get
(
const int       mpartId,
const int       iblock,
      double  **elapsed,
      double  **cpu,
      double  **cpu_user,
      double  **cpu_sys
);


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
);

void
PROCF (pdm_multipart_free, PDM_MULTIPART_FREE)
(
 const int *id
);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MULTIPART_H__ */
