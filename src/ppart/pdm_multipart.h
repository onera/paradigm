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
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Identifier
 */

int
PDM_multipart_create
(
 const int          n_block,
 const PDM_MPI_Comm comm
);

void
PROCF (pdm_multipart_create, PDM_MULTIPART_CREATE)
(
 const int *n_block,
 const PDM_MPI_Fint *fcomm,
       int *id
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
