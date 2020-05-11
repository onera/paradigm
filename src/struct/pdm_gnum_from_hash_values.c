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
#include <stdint.h>
#include <stdbool.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_handles.h"
#include "pdm_binary_search.h"
#include "pdm_mpi.h"
#include "pdm_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_gnum_from_hash_values.h"

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
 * \struct _pdm_gnum_from_hv_t
 * \brief  Define a global numberring
 *
 */

typedef struct {
  int             n_part;          /*!< Number of partitions                     */
  PDM_bool_t      equilibrate;     /*!< Equilibrate the hash values distribution */
  PDM_MPI_Comm    comm;            /*!< MPI communicator                         */
  int             n_rank;          /*!< MPI communicator size                    */

  int            *n_elts;          /*!< Number of elements in partitions         */
  size_t        **part_hkeys;
  unsigned char **part_hdata;
  int           **part_hstri;

  // size_t         *blk_hkeys;
  unsigned char  *blk_hdata;
  int            *blk_hstri;

  PDM_g_num_t     n_g_elt;        /*!< Global number of elements                 */
  PDM_g_num_t   **g_nums;         /*!< Global numbering of elements              */

  PDM_g_num_t     *distribution;

} _pdm_gnum_from_hv_t;

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_gnums_from_hv   = NULL;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Return ppart object from it identifier
 *
 * \param [in]   ppart_id        ppart identifier
 *
 */

static _pdm_gnum_from_hv_t *
_get_from_id
(
 int  id
)
{

  _pdm_gnum_from_hv_t *gnum_from_hv = (_pdm_gnum_from_hv_t *) PDM_Handles_get (_gnums_from_hv, id);

  if (gnum_from_hv == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_gnum_from_hv error : Bad identifier\n");
  }

  return gnum_from_hv;
}

/**
 *
 * \brief Compute with equilibrate algorithm
 *
 * \param [in]   _gnum          Current _pdm_gnum_from_hv_t structure
 *
 */
static void
_compute_distribution_equilibrate
(
 _pdm_gnum_from_hv_t *_gnum
)
{
  printf("_gnum_from_hv_compute_equilibrate Not implemented \n");
  abort();

}

/**
 *
 * \brief Setup a naive distribution from min max of data
 */
static void
setup_distribution_from_min_max
(
 size_t       min_elt,
 size_t       max_elt,
 PDM_g_num_t* distribution,
 int          n_dist
)
{
  PDM_g_num_t nelmt = max_elt - min_elt + 1;

  assert(nelmt > 0);

  PDM_g_num_t quotient  = nelmt/n_dist;
  PDM_g_num_t remainder = nelmt%n_dist;

  printf(PDM_FMT_G_NUM"\n", quotient);
  printf(PDM_FMT_G_NUM"\n", remainder);

  distribution[0] = 0;
  for(int i = 1; i < n_dist+1; ++i) {
    distribution[i] = quotient;
    PDM_g_num_t i1 = i - 1;
    if(i1 < remainder){
      distribution[i] += 1;
    }
  }

  for(int i = 0; i < n_dist; ++i) {
    distribution[i+1] += distribution[i];
  }


  printf(" distribution:: ");
  for(int i = 0; i < n_dist+1; ++i) {
    printf(PDM_FMT_G_NUM" ", distribution[i]);
  }

}

/**
 *
 * \brief Compute with equilibrate algorithm
 *
 * \param [in]   _gnum          Current _pdm_gnum_from_hv_t structure
 *
 */
static void
_compute_distribution
(
 _pdm_gnum_from_hv_t *_gnum_from_hv
)
{

  size_t max_key_loc = 0;
  size_t min_key_loc = SIZE_MAX;

  size_t max_key = 0;
  size_t min_key = SIZE_MAX;

  for(int i_part = 0; i_part < _gnum_from_hv->n_part; i_part++){
    for(int ielt = 0; ielt < _gnum_from_hv->n_elts[i_part]; ++ielt){
      max_key_loc = PDM_MAX(max_key_loc, _gnum_from_hv->part_hkeys[i_part][ielt]);
      min_key_loc = PDM_MIN(min_key_loc, _gnum_from_hv->part_hkeys[i_part][ielt]);
    }
  }
  int ierr;
  ierr = PDM_MPI_Allreduce(&min_key_loc, &min_key, 1, PDM_MPI_UNSIGNED_LONG, PDM_MPI_MIN, _gnum_from_hv->comm);
  assert(ierr == 0);

  ierr = PDM_MPI_Allreduce(&max_key_loc, &max_key, 1, PDM_MPI_UNSIGNED_LONG, PDM_MPI_MAX, _gnum_from_hv->comm);
  assert(ierr == 0);

  printf(" max_key:: %lu \n", max_key);
  printf(" min_key:: %lu \n", min_key);

  /* Prepare distribution from min and max elements */
  setup_distribution_from_min_max(min_key, max_key, _gnum_from_hv->distribution, _gnum_from_hv->n_rank);



}

/**
 *
 * \brief Compute but without equilibrate
 *
 * \param [in]   _gnum          Current _pdm_gnum_from_hv_t structure
 *
 */
static void
_gnum_from_hv_compute
(
 _pdm_gnum_from_hv_t *_gnum_from_hv
)
{
  printf("_gnum_from_hv_compute \n");

  if(_gnum_from_hv->equilibrate) {
    _compute_distribution_equilibrate(_gnum_from_hv);
  } else {
    _compute_distribution(_gnum_from_hv);
  }

  /*
   * Remapping of partition data in block data according to the hash values distribution
   */



  /*
   * Generate global numbering from the block_data
   */
  // PDM_generate_global_id_from();

  /*
   * Reverse all_to_all exchange in order to remap global id on current partition
   */

}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Build a global numbering structure
 *
 * \param [in]   n_part       Number of local partitions
 * \param [in]   equilibrate  Use algorithm to equilibrate the block treatment (hash value is not a priori equi-reparti)
 * \param [in]   tolerance    Geometric tolerance (used if merge double points is activated)
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Identifier
 */

int
PDM_gnum_from_hash_values_create
(
 const int          n_part,
 const PDM_bool_t   equilibrate,
 const PDM_MPI_Comm comm
)
{
  int i_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  /*
   * Search a gnum_from_hash_values free id
   */
  if (_gnums_from_hv == NULL) {
    _gnums_from_hv = PDM_Handles_create (4);
  }

  _pdm_gnum_from_hv_t *_gnum_from_hv = (_pdm_gnum_from_hv_t *) malloc(sizeof(_pdm_gnum_from_hv_t));
  int id = PDM_Handles_store (_gnums_from_hv, _gnum_from_hv);

  _gnum_from_hv->n_part      = n_part;
  _gnum_from_hv->equilibrate = equilibrate;
  _gnum_from_hv->comm        = comm;
  _gnum_from_hv->n_rank      = n_rank;
  _gnum_from_hv->n_g_elt     = -1;
  _gnum_from_hv->g_nums      = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t * ) * n_part);

  _gnum_from_hv->n_elts      = (int            *) malloc (sizeof(int            ) * n_part);
  _gnum_from_hv->part_hkeys  = (size_t        **) malloc (sizeof(size_t        *) * n_part);
  _gnum_from_hv->part_hstri  = (int           **) malloc (sizeof(int           *) * n_part);
  _gnum_from_hv->part_hdata  = (unsigned char **) malloc (sizeof(unsigned char *) * n_part);

  for (int i = 0; i < n_part; i++) {
    _gnum_from_hv->g_nums[i]     = NULL;
    _gnum_from_hv->part_hkeys[i] = NULL;
    _gnum_from_hv->part_hstri[i] = NULL;
    _gnum_from_hv->part_hdata[i] = NULL;
  }

  _gnum_from_hv->distribution = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t ) * ( n_rank + 1 ));

  return id;

}

void
PROCF (pdm_gnum_from_hash_values_create, PDM_GNUM_FROM_HVALUES_CREATE)
(
 const int          *n_part,
 const int          *equilibrate,
 const PDM_MPI_Fint *fcomm,
       int          *id
)
{
  const PDM_MPI_Comm c_comm = PDM_MPI_Comm_f2c (*fcomm);

  *id = PDM_gnum_from_hash_values_create (*n_part, (PDM_bool_t) *equilibrate, c_comm);
}

/**
 *
 * \brief Set hash values for one partition
 *
 * \param [in]   id           Identifier
 * \param [in]   i_part       Current partition
 * \param [in]   n_elts       Number of elements
 * \param [in]   part_hkey    For each elements the hash value associated
 * \param [in]   part_strid   Stride between each data in part_hdata
 * \param [in]   part_hdata   Partition data which compute the hash value, we need it to setup in a block way
 *
 */

void
PDM_gnum_set_hash_values
(
 const int            id,
 const int            i_part,
 const int            n_elts,
 const size_t        *part_hkeys,
 const int           *part_hstri,
 const unsigned char *part_hdata
)
{
  _pdm_gnum_from_hv_t *_gnum_from_hv = _get_from_id (id);

  assert(_gnum_from_hv->part_hkeys != NULL);
  assert(_gnum_from_hv->part_hstri != NULL);
  assert(_gnum_from_hv->part_hdata != NULL);

  assert(_gnum_from_hv->part_hkeys[i_part] == NULL);
  assert(_gnum_from_hv->part_hstri[i_part] == NULL);
  assert(_gnum_from_hv->part_hdata[i_part] == NULL);

  _gnum_from_hv->n_elts[i_part]      = n_elts;
  _gnum_from_hv->part_hkeys[i_part]  = (size_t        *) part_hkeys;
  _gnum_from_hv->part_hstri[i_part]  = (int           *) part_hstri;
  _gnum_from_hv->part_hdata[i_part]  = (unsigned char *) part_hdata;

}

void
PROCF (pdm_gnum_set_hash_values, PDM_GNUM_SET_FROM_HASH_VALUES)
(
 const int           *id,
 const int           *i_part,
 const int           *n_elts,
 const size_t        *part_hkeys,
 const int           *part_hstri,
 const unsigned char *part_hdata
)
{
  PDM_gnum_set_hash_values (*id, *i_part, *n_elts, part_hkeys, part_hstri, part_hdata);
}


/**
 *
 * \brief Compute
 *
 * \param [in]   id           Identifier
 *
 */


void
PDM_gnum_from_hv_compute
(
 const int id
)
{
  _pdm_gnum_from_hv_t *_gnum_from_hv = _get_from_id (id);

  printf("PDM_gnum_from_hv_compute::oooooooooo \n");

  _gnum_from_hv_compute(_gnum_from_hv);

}

void
PROCF (PDM_gnum_from_hv_compute, PDM_GNUM_FROM_HV_COMPUTE)
(
 const int *id
)
{
  PDM_gnum_from_hv_compute (*id);
}


/**
 *
 * \brief Set from coordinates
 *
 * \param [in]   id           Identifier
 * \param [in]   i_part       Current partition
 * \param [in]   n_elts       Number of elements
 * \param [in]   coords       Coordinates (size = 3 * \ref n_elts)
 *
 */

PDM_g_num_t *
PDM_gnum_from_hv_get
(
 const int id,
 const int i_part
)
{
  _pdm_gnum_from_hv_t *_gnum_from_hv = _get_from_id (id);

  return _gnum_from_hv->g_nums[i_part];
}

void
PROCF (pdm_gnum_from_hv_get, PDM_GNUM_FROM_HV_GET)
(
 const int *id,
 const int *i_part,
 PDM_g_num_t *gnum
)
{
  _pdm_gnum_from_hv_t *_gnum_from_hv = _get_from_id (*id);

  const PDM_g_num_t *tmp = PDM_gnum_from_hv_get (*id, *i_part);
  for (int i = 0; i < _gnum_from_hv->n_elts[*i_part]; i++) {
    gnum[i] = tmp[i];
  }
}


/**
 *
 * \brief Free
 *
 * \param [in]   id           Identifier
 *
 */

void
PDM_gnum_from_hv_free
(
 const int id,
 const int partial
)
{
  _pdm_gnum_from_hv_t *_gnum_from_hv = _get_from_id (id);

  if (partial != 1) {
    for (int i = 0; i < _gnum_from_hv->n_part; i++) {
      free (_gnum_from_hv->g_nums[i]);
    }
  }

  free (_gnum_from_hv->g_nums);
  free (_gnum_from_hv->n_elts);
  free (_gnum_from_hv->part_hkeys);
  free (_gnum_from_hv->part_hstri);
  free (_gnum_from_hv->part_hdata);
  free (_gnum_from_hv->distribution);

  free (_gnum_from_hv);

  PDM_Handles_handle_free (_gnums_from_hv, id, PDM_FALSE);

  const int n_gnum_from_hv = PDM_Handles_n_get (_gnums_from_hv);

  if (n_gnum_from_hv == 0) {
    _gnums_from_hv = PDM_Handles_free (_gnums_from_hv);
  }

}

void
PROCF (pdm_gnum_from_hv_free, PDM_GNUM_FROM_HV_FREE)
(
 const int *id,
 const int *partial
)
{
  PDM_gnum_from_hv_free (*id, *partial);
}



void
PDM_generate_global_id_from
(
 const int              blk_size,
 const unsigned char   *blk_data,
 const int             *blk_stri,
 gnum_from_hv_compare   fcompare,
 gnum_from_hv_equal     fequal,
 PDM_g_num_t          **gnum
)
{
  printf(" TODO \n");
  abort();
  // int nBlock = blockPaths.size();
  // std::vector<int> orderName(nBlock);
  // std::iota(begin(orderName), end(orderName), 0);
  // std::sort(begin(orderName), end(orderName), [&](const int& i1, const int& i2){
  //   return blockPaths[i1] < blockPaths[i2];
  // });

  // // -------------------------------------------------------------------
  // // 2 - Give an local number for each element in blockPaths
  // std::vector<int> globalNameNum(nBlock);
  // int nextNameId =  0;
  // int nLocNameId =  0;
  // std::string lastName;
  // for(int i = 0; i < nBlock; i++){
  //   if(blockPaths[orderName[i]] == lastName){
  //     globalNameNum[orderName[i]] = nextNameId;
  //   } else {
  //     nextNameId++;
  //     nLocNameId++;
  //     globalNameNum[orderName[i]] = nextNameId;
  //     lastName = blockPaths[orderName[i]];
  //   }
  // }

  // // -------------------------------------------------------------------
  // // 3 - Setup global numbering by simply shift
  // int shiftG;
  // int ierr = MPI_Scan(&nLocNameId, &shiftG, 1, MPI_INT, MPI_SUM, comm);
  // assert(ierr == 0);
  // shiftG -= nLocNameId;

  // for(int i = 0; i < nBlock; i++){
  //   globalNameNum[i] += shiftG;
  // }
}



/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
