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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_handles.cuh"
#include "pdm_error.h"
#include "pdm_cuda_error.cuh"

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

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/



/**
 * \brief Create handles storage
 *
 * \param [in] init_size      Initial size of storage array
 *
 * \return  New handles storage
 */

PDM_Handles_t *
PDM_Handles_create_GPU
(
 const int init_size
)
{
  PDM_Handles_t *nouv =  (PDM_Handles_t*)malloc (sizeof(PDM_Handles_t));

  nouv->s_array = init_size;
  nouv->n_handles = 0;
  nouv->array = (const void **)malloc(sizeof(void*) * init_size);
  nouv->idx = (int*)malloc(sizeof(int) * init_size);
  nouv->idx_inv = (int*)malloc(sizeof(int) * init_size);

  for (int i = 0; i < nouv->s_array; i++) {
    nouv->array[i] = NULL;
    nouv->idx[i] = -1;
    nouv->idx_inv[i] = -1;
  }

  return nouv;
}


/**
 * \brief Store a new handle pointer
 *
 * \param [in] handles      Current handles storage
 * \param [in] handle_ptr   Handle pointer
 *
 * \return  Handle index
 */

int
PDM_Handles_store_GPU
(
 PDM_Handles_t *handles,
 const void *handle_ptr
)
{
  if (handles->n_handles >= handles->s_array) {
    int p_s_array = handles->s_array;
    handles->s_array *= 2;
    handles->array   = (const void**)realloc(handles->array, sizeof(void*) * handles->s_array);
    handles->idx     = (int*)realloc(handles->idx, sizeof(int) * handles->s_array);
    handles->idx_inv = (int*)realloc(handles->idx_inv, sizeof(int) * handles->s_array);

    for (int i = p_s_array; i < handles->s_array; i++) {
      handles->array[i] = NULL;
      handles->idx[i] = -1;
      handles->idx_inv[i] = -1;
    }
  }

  int idx = 0;
  while (handles->array[idx] != NULL)
    idx++;

  handles->array[idx] = handle_ptr;
  handles->idx[handles->n_handles] = idx;
  handles->idx_inv[idx] = handles->n_handles;

  handles->n_handles += 1;

  return idx;
}


/**
 * \brief Get handle pointer
 *
 * \param [in] handles  Current handles storage
 * \param [in] handle_idx   Handle index
 *
 * \return  Handle pointer
 */

__device__
const void *
PDM_Handles_get_GPU
( 
 PDM_Handles_t *handles,
  const int handle_idx
)
{
  if (handle_idx >= handles->s_array) {
    PDM_error_GPU (__FILE__, __LINE__, 0, "PDM_Handle : Bad identifier\n");
  }

  return handles->array[handle_idx];
}


#ifdef __cplusplus
}
#endif /* __cplusplus */