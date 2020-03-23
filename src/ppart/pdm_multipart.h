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
 * \param [in]   n_zone       Number of zones in the original mesh
 * \param [in]   n_part       Number of partition per proc in each zone
 * \param [in]   merge_blocks Merge or not the zones before splitting
 * \param [in]   split_method Choice of library used to split the mesh
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Identifier
 */

int
PDM_multipart_create
(
 const int              n_zone,
 const int             *n_part,
 const PDM_bool_t       merge_blocks,
 const PDM_part_split_t split_method,
 const PDM_MPI_Comm     comm
);


/**
 *
 * \brief Set a block in the multipart structure
 *
 * \param [in]   id           Identifier
 * \param [in]   i_block      Number of block to set
 * \param [in]   dface_cell    Face to cell connectivity for the block
 * TODO LIST PARAMS
 *
 */

void PDM_multipart_register_block
(
 const int        mpart_id,
 const int        zoneGId,
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
const   int  zoneGId,
const   int  i_part,
 int        *n_cell,
 int        *n_face,
 int        *n_face_part_bound,
 int        *n_vtx,
 int        *n_proc,
 int        *n_total_part,
 int        *scell_face,
 int        *sface_vtx,
 int        *sFaceBound,
 int        *n_faceBound,
 int        *sFaceJoin,
 int        *n_faceJoin
);

void
PDM_multipart_part_val_get
(
const int            mpartId,
const int            zoneGId,
const int            i_part,
      int          **cell_tag,
      int          **cell_face_idx,
      int          **cell_face,
      PDM_g_num_t  **cell_ln_to_gn,
      int          **face_tag,
      int          **face_cell,
      int          **face_vtx_idx,
      int          **face_vtx,
      PDM_g_num_t  **face_ln_to_gn,
      int          **face_part_bound_proc_idx,
      int          **face_part_bound_part_idx,
      int          **face_part_bound,
      int          **vtx_tag,
      double       **vtx,
      PDM_g_num_t  **vtx_ln_to_gn,
      int          **faceBoundIdx,
      int          **faceBound,
      PDM_g_num_t  **faceBoundLNToGN,
      int          **faceJoinIdx,
      int          **faceJoin,
      PDM_g_num_t  **faceJoinLNToGN
);

void
PDM_multipart_part_color_get
(
const int            mpartId,
const int            zoneGId,
const int            i_part,
      int          **cellColor,
      int          **faceColor,
      int          **threadColor,
      int          **hyperPlaneColor
);

void
PDM_multipart_time_get
(
const int       mpartId,
const int       zoneGId,
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


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MULTIPART_H__ */
