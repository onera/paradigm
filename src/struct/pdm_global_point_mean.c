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

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_morton.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_handles.h"
#include "pdm_mpi.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_global_point_mean.h"

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
 * \struct _pdm_global_point_mean_t
 * \brief  Define a global point mean
 * 
 */

typedef struct  {

  int          n_part;            /*!< Number of partitions */
  PDM_MPI_Comm comm;              /*!< MPI communicator */
  int          *n_elts;           /*!< Number of elements in partitions */
  PDM_g_num_t **g_nums;           /*!< Global numbering of elements */ 
  PDM_part_to_block_t *ptb;       /*!< Part to block structure */
  PDM_block_to_part_t *btp;       /*!< Block to part structure */
  double      **local_field;      /*!< Local field */
  double      **local_weight;      /*!< Weight */
  double      **global_mean_field;/*!< Global mean field */
  
} _pdm_global_point_mean_t;


/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_gpms   = NULL;

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

static _pdm_global_point_mean_t *
_get_from_id
(
 int  id
)
{
  
  _pdm_global_point_mean_t *gpm = 
          (_pdm_global_point_mean_t *) PDM_Handles_get (_gpms, id);
    
  if (gpm == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_global_point_mean error : Bad identifier\n");
  }

  return gpm;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a structure that compute a global mean   
 *
 * \param [in]   n_part       Number of local partitions 
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Identifier    
 */

int
PDM_global_point_mean_create
(
 const int n_part,
 const PDM_MPI_Comm comm
)
{
  if (_gpms == NULL) {
    _gpms = PDM_Handles_create (4);
  }

  _pdm_global_point_mean_t *_gpm = 
          (_pdm_global_point_mean_t *) malloc(sizeof(_pdm_global_point_mean_t));
  int id = PDM_Handles_store (_gpms, _gpm);

  _gpm->n_part = n_part;
  _gpm->comm   = comm;
  _gpm->g_nums = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t * ) * n_part); 
  _gpm->n_elts = (int *) malloc (sizeof(int) * n_part);
  _gpm->ptb    = NULL;
  _gpm->btp    = NULL;
  
  for (int i = 0; i < n_part; i++) {
    _gpm->g_nums[i] = NULL;
  }
  
  _gpm->local_field = (double **) malloc (sizeof(double * ) * n_part); 
  _gpm->local_weight = (double **) malloc (sizeof(double * ) * n_part); 
  _gpm->global_mean_field = (double **) malloc (sizeof(double * ) * n_part);
  
  return id;
}

void
PROCF (pdm_global_point_mean_create, PDM_GLOBAL_POINT_MEAN_CREATE)
(
 const int *n_part,
 const PDM_MPI_Fint *fcomm,
       int *id
)
{
  const PDM_MPI_Comm c_comm = PDM_MPI_Comm_f2c (*fcomm);

  *id = PDM_global_point_mean_create (*n_part, c_comm);
}


/**
 *
 * \brief Set absolute number   
 *
 * \param [in]   id           Identifier
 * \param [in]   i_part       Current partition
 * \param [in]   n_point      Number of points in the partition
 * \param [in]   numabs       Absolute number of points
 *
 */

void
PDM_global_point_mean_set
(
 const int          id,
 const int          i_part,
 const int          n_point,
 const PDM_g_num_t *numabs
)
{
  _pdm_global_point_mean_t *_gpm = _get_from_id (id);
  
  _gpm->g_nums[i_part] = NULL;
  _gpm->n_elts[i_part]  = 0;
  
}

void
PROCF (pdm_global_point_mean_set, PDM_GLOBAL_POINT_MEAN_SET)
(
 const int         *id,
 const int         *i_part,
 const int         *n_point,
 const PDM_g_num_t *numabs
)
{
  PDM_global_point_mean_set (*id, *i_part, *n_point, numabs);
}


/**
 *
 * \brief Free a global point mean structure   
 *
 * \param [in]   id           Identifier
 *
 * \return     Identifier    
 */

int
PDM_global_point_mean_free
(
 const int          id
)
{
  _pdm_global_point_mean_t *_gpm = _get_from_id (id);
  
  if (_gpm->ptb != NULL) {
    _gpm->ptb = PDM_part_to_block_free (_gpm->ptb);
  }

  if (_gpm->btp != NULL) {
    _gpm->btp = PDM_block_to_part_free (_gpm->btp);    
  }
          
  free (_gpm->g_nums);
  free (_gpm->n_elts);
  free (_gpm->local_field);
  free (_gpm->local_weight);
  free (_gpm->global_mean_field);

}

void
PROCF (pdm_global_point_mean_free, PDM_GLOBAL_POINT_MEAN_FREE)
(
 const int         *id
)
{
  PDM_global_point_mean_free (*id);
}


/**
 *
 * \brief Set local field and it associated weight    
 *
 * \param [in]   id                    Identifier
 * \param [in]   i_part                Current partition
 * \param [in]   stride                Stride of the field
 * \param [in]   local_field           Local value of field
 * \param [in]   local_weight          Local weight used to compute the mean value
 * \param [in]   global_mean_field_ptr Pointer where global mean field
 *                                     will be stored after computing
 */

void
PDM_global_point_mean_field_set
(
 const int          id,
 const int          i_part,
 const int          stride, 
 const double      *local_field,
 const double      *local_weight,
 double            *global_mean_field_ptr
)
{
}

void
PROCF (pdm_global_point_mean_field_set, PDM_GLOBAL_POINT_MEAN_FIELD_SET)
(
 const int         *id,
 const int         *i_part,
 const int         *stride, 
 const double      *local_field,
 const double      *local_weight,
 double            *global_mean_field_ptr
)
{
  PDM_global_point_mean_field_set (*id, *i_part, *stride, 
                                   local_field, local_weight, global_mean_field_ptr);
}


/**
 *
 * \brief Compute the global average field   
 *
 * \param [in]   id           Identifier
 *
 */

void
PDM_global_point_mean_field_compute
(
 const int          id
)
{
}

void
PROCF (pdm_global_point_mean_field_compute, PDM_GLOBAL_POINT_MEAN_FIELD_COMPUTE)
(
 const int         *id
)
{
}



/**
 *
 * \brief Set local field and it associated weight    
 *
 * \param [in]   id                 Identifier
 * \param [in]   i_part             Current partition
 * \param [in]   global_mean_field  Global point mean field
 *
 */

void
PDM_global_point_mean_field_get
(
 const int          id,
 const int          i_part,
 const double      *global_mean_field
)
{
}

void
PROCF (pdm_global_point_mean_field_get, PDM_GLOBAL_POINT_MEAN_FIELD_GET)
(
 const int         *id,
 const int         *i_part,
 const double      *global_mean_field
)
{
}


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
