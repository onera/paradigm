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

#include "pdm_handles.h"
#include "pdm_error.h"
#include "pdm_priv.h"

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
 * \struct _PDM_Handles_t
 * \brief  Helper struct to hold a collection of private internal structure
 *
 */

struct _PDM_Handles_t {

  int *idx;            /*!< list of Index */
  int *idx_inv;        /*!< Index -> indice (1 to n_handles) */
  const void **array;  /*!< Storage array */
  int s_array;         /*!< Size of array */
  int n_handles;       /*!< Number of stored handles */

};

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
PDM_Handles_create
(
 const int init_size
)
{
  PDM_Handles_t *new_handle;
  PDM_malloc(new_handle, 1, PDM_Handles_t);

  new_handle->s_array   = init_size;
  new_handle->n_handles = 0;
  PDM_malloc(new_handle->array  , init_size, const void*);
  PDM_malloc(new_handle->idx    , init_size,       int  );
  PDM_malloc(new_handle->idx_inv, init_size,       int  );

  for (int i = 0; i < new_handle->s_array; i++) {
    new_handle->array[i]   = NULL;
    new_handle->idx[i]     = -1;
    new_handle->idx_inv[i] = -1;
  }

  return new_handle;
}


/**
 * \brief Free handles storage
 *
 * \param [in] handles  Current handles storage
 *
 * \return               NULL
 */

PDM_Handles_t *
PDM_Handles_free
(
 PDM_Handles_t *handles
)
{
  if (handles != NULL) {
    PDM_free(handles->array);
    PDM_free(handles->idx);
    PDM_free(handles->idx_inv);
    PDM_free(handles);
  }
  return NULL;
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
PDM_Handles_store
(
 PDM_Handles_t *handles,
 const void *handle_ptr
)
{
  if (handles->n_handles >= handles->s_array) {
    int p_s_array = handles->s_array;
    handles->s_array *= 2;
    PDM_realloc(handles->array   ,handles->array   , handles->s_array, const void*);
    PDM_realloc(handles->idx     ,handles->idx     , handles->s_array,int  );
    PDM_realloc(handles->idx_inv ,handles->idx_inv , handles->s_array,int  );

    for (int i = p_s_array; i < handles->s_array; i++) {
      handles->array[i]   = NULL;
      handles->idx[i]     = -1;
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

const void *
PDM_Handles_get
(
 PDM_Handles_t *handles,
  const int handle_idx
)
{
  if (handle_idx >= handles->s_array) {
    PDM_error (__FILE__, __LINE__, 0, "PDM_Handle : Bad identifier\n");
  }

  return handles->array[handle_idx];
}


/**
 * \brief Get handles index
 *
 * \param [in] handles  Current handles storage
 *
 * \return  Handles index
 */

const int *
PDM_Handles_idx_get
(
 PDM_Handles_t *handles
)
{
  return handles->idx;
}


/**
 * \brief Free a handle
 *
 * \param [in] handles      Current handles storage
 * \param [in] handle_idx   Handle index
 * \param [in] st_free_data Free data or not
 *
 * \return  Handle pointer
 */

void
PDM_Handles_handle_free
(
 PDM_Handles_t *handles,
 const int handle_idx,
 const PDM_bool_t st_free_data
)
{
  if (handle_idx >= handles->s_array) {
    PDM_error (__FILE__, __LINE__, 0, "PDM_Handle : Bad identifier\n");
  }

  if (st_free_data) {
    if (handles->array[handle_idx] != NULL) {
      PDM_free(handles->array[handle_idx]);
    }
  }

  handles->array[handle_idx] = NULL;

  int ind = handles->idx_inv[handle_idx];

  for (int i = ind + 1; i < handles->n_handles; i++) {
    handles->idx[i-1] = handles->idx[i];
    handles->idx_inv[handles->idx[i-1]] = i-1;
  }

  handles->n_handles += -1;

}


/**
 * \brief Get number of stored handles
 *
 * \param [in] handles  Current handles storage
 *
 * \return  Number of stored handles
 */

int
PDM_Handles_n_get
(
 PDM_Handles_t *handles
)
{
  return handles->n_handles;
}
