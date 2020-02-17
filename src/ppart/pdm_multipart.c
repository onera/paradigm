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
 * \brief  TODO write description
 *
 */

typedef struct  {

  int          n_block;     /*!< Number of partitions */
  PDM_MPI_Comm comm;        /*!< MPI communicator */
  PDM_g_num_t  n_g_elt;     /*!< Global number of elements */
  int          *n_elts;     /*!< Number of elements in partitions */
  PDM_g_num_t **g_nums;     /*!< Global numbering of elements */

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
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Identifier
 */

int
PDM_multipart_create
(
 const int          n_block,
 const PDM_MPI_Comm comm
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
  _multipart->comm        = comm;
  _multipart->n_g_elt     = -1;
  _multipart->g_nums      = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t * ) * n_block);
  _multipart->n_elts      = (int *) malloc (sizeof(int) * n_block);

  for (int i = 0; i < n_block; i++) {
    _multipart->g_nums[i] = NULL;
  }
  PDM_printf("Created from PDM_multipart_create\n");
  return id;

}

void
PROCF (pdm_multipart_create, PDM_MULTIPART_CREATE)
(
 const int          *n_block,
 const PDM_MPI_Fint *fcomm,
       int          *id
)
{
  const PDM_MPI_Comm c_comm = PDM_MPI_Comm_f2c (*fcomm);

  *id = PDM_multipart_create (*n_block, c_comm);
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
    free (_multipart->g_nums[i]);
  }


  free (_multipart->g_nums);
  free (_multipart->n_elts);

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
