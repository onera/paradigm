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

//
typedef struct  {
  int* idx;
  int* arr;
} user_defined_sort_t;

// static int
// my_compare_edge
// (
// const void* a,
// const void* b,
//       void* ctxt)
// {
//   int i = *(const int *) a;
//   int j = *(const int *) b;

//   user_defined_sort_t* my_struct = (user_defined_sort_t *) ctxt;
//   int ni = 2
//   int nj = 2
//   int* arr_i = my_struct->arr[2*my_struct->idx[i]];
//   int* arr_j = my_struct->arr[2*my_struct->idx[j]];
// }


static int
my_compare
(
const void* a,
const void* b,
      void* ctxt)
{
  int i = *(const int *) a;
  int j = *(const int *) b;

  user_defined_sort_t* my_struct = (user_defined_sort_t*) ctxt;

  int ni = my_struct->idx[i+1] - my_struct->idx[i];
  int nj = my_struct->idx[j+1] - my_struct->idx[j];

  printf("my_compare:: %d %d - %d %d \n", i, j, ni, nj);

  int* arr_i = &my_struct->arr[my_struct->idx[i]];
  int* arr_j = &my_struct->arr[my_struct->idx[j]];

  if(ni < nj){
    return 1;
  } else if (ni == nj){
    for(int k = 0; k < ni; ++k){
      if(arr_i[k] < arr_j[k]) {
        return 1;
      }
    }
  }

  return 0;
}

static int
my_equal
(
const void* a,
const void* b,
      void* ctxt)
{
  int i = *(const int *) a;
  int j = *(const int *) b;


  user_defined_sort_t* my_struct = (user_defined_sort_t*) ctxt;

  int ni = my_struct->idx[i+1] - my_struct->idx[i];
  int nj = my_struct->idx[j+1] - my_struct->idx[j];

  printf("my_equal:: %d %d - %d %d - %d %d \n", i, j, ni, nj, my_struct->idx[i], my_struct->idx[j]);

  int* arr_i = &my_struct->arr[my_struct->idx[i]];
  int* arr_j = &my_struct->arr[my_struct->idx[j]];


  if(ni != nj){
    return 0;
  } else if (ni == nj){

    /* Dans notre cas on veut sort les entiers avant de les comparers */
    int* sort_arr_i = (int*) malloc( ni * sizeof(int));
    int* sort_arr_j = (int*) malloc( ni * sizeof(int));

    for(int k = 0; k < ni; ++k){
      sort_arr_i[k] = arr_i[k];
      sort_arr_j[k] = arr_j[k];
    }
    PDM_quick_sort_int(sort_arr_i, 0, ni-1);
    PDM_quick_sort_int(sort_arr_j, 0, ni-1);

    for(int k = 0; k < ni; ++k){
      printf(" \t sort_arr_i[%d] = %d | sort_arr_j[%d] = %d \n", k, sort_arr_i[k], k, sort_arr_j[k]);
      if(sort_arr_i[k] != sort_arr_j[k]) {
        free(sort_arr_i);
        free(sort_arr_j);
        return 0;
      }
    }

    free(sort_arr_i);
    free(sort_arr_j);
  }

  return 1;
}



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
  size_t          s_data;

  gnum_from_hv_compare fcompare;
  gnum_from_hv_equal   fequal;

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

  distribution[0] = min_elt;
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
  printf("\n");

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
  int* n_key_send  = (int *) malloc( _gnum_from_hv->n_rank * sizeof(int));
  int* n_key_recv  = (int *) malloc( _gnum_from_hv->n_rank * sizeof(int));
  int* n_data_send = (int *) malloc( _gnum_from_hv->n_rank * sizeof(int));
  int* n_data_recv = (int *) malloc( _gnum_from_hv->n_rank * sizeof(int));

  for(int i = 0; i < _gnum_from_hv->n_rank; ++i){
    n_key_send[i] = 0;
    n_data_send[i] = 0;
  }

  /*
   * Prepare send
   */
  for(int i_part = 0; i_part < _gnum_from_hv->n_part; ++i_part){
    for(int ielt = 0; ielt < _gnum_from_hv->n_elts[i_part]; ++ielt){

      PDM_g_num_t g_key = (PDM_g_num_t) _gnum_from_hv->part_hkeys[i_part][ielt];
      printf(" Search for :: %d\n", (int)g_key);
      int t_rank = PDM_binary_search_gap_long(g_key, _gnum_from_hv->distribution, _gnum_from_hv->n_rank+1);

      printf(" Found in t_rank :: %d\n", (int)t_rank);
      // n_data_send[t_rank] += _gnum_from_hv->s_data * _gnum_from_hv->part_hstri[i_part][ielt];
      n_data_send[t_rank] += _gnum_from_hv->part_hstri[i_part][ielt];
      n_key_send[t_rank]++;

    }
  }

  /*
   * Exchange
   */
  PDM_MPI_Alltoall(n_key_send , 1, PDM_MPI_INT, n_key_recv , 1, PDM_MPI_INT, _gnum_from_hv->comm);
  PDM_MPI_Alltoall(n_data_send, 1, PDM_MPI_INT, n_data_recv, 1, PDM_MPI_INT, _gnum_from_hv->comm);

  /*
   * Prepare utility array to setup the second exchange
   */
  int* i_key_send  = (int *) malloc( (_gnum_from_hv->n_rank+1) * sizeof(int));
  int* i_key_recv  = (int *) malloc( (_gnum_from_hv->n_rank+1) * sizeof(int));
  int* i_data_send = (int *) malloc( (_gnum_from_hv->n_rank+1) * sizeof(int));
  int* i_data_recv = (int *) malloc( (_gnum_from_hv->n_rank+1) * sizeof(int));

  i_key_send[0] = 0;
  i_key_recv[0] = 0;
  i_data_send[0] = 0;
  i_data_recv[0] = 0;
  for(int i = 0; i < _gnum_from_hv->n_rank; ++i){
    i_key_send[i+1] = i_key_send[i] + n_key_send[i];
    i_key_recv[i+1] = i_key_recv[i] + n_key_recv[i];
    i_data_send[i+1] = i_data_send[i] + n_data_send[i];
    i_data_recv[i+1] = i_data_recv[i] + n_data_recv[i];
    n_key_send[i]  = 0;
    n_data_send[i] = 0;
  }

  int s_send_keys = i_key_send[_gnum_from_hv->n_rank];
  int s_recv_keys = i_key_recv[_gnum_from_hv->n_rank];
  int s_send_data = i_data_send[_gnum_from_hv->n_rank] * _gnum_from_hv->s_data;
  int s_recv_data = i_data_recv[_gnum_from_hv->n_rank] * _gnum_from_hv->s_data;

  printf("s_send_keys::%d\n", s_send_keys);
  printf("s_recv_keys::%d\n", s_recv_keys);
  printf("s_send_data::%d\n", s_send_data);
  printf("s_recv_data::%d\n", s_recv_data);
  printf("i_data_send[_gnum_from_hv->n_rank]::%d\n", i_data_send[_gnum_from_hv->n_rank]);
  printf("i_data_recv[_gnum_from_hv->n_rank]::%d\n", i_data_recv[_gnum_from_hv->n_rank]);

  /*
   * Allocate
   */
  int *send_buffer_keys = (int *) malloc(sizeof(int) * s_send_keys );
  int *recv_buffer_keys = (int *) malloc(sizeof(int) * s_recv_keys );

  int *send_buffer_stri = (int *) malloc(sizeof(int) * ( s_send_keys + 1) );
  int *recv_buffer_stri = (int *) malloc(sizeof(int) * ( s_recv_keys + 1) );

  unsigned char *send_buffer_data = (unsigned char *) malloc(sizeof(unsigned char) * s_send_data);
  unsigned char *recv_buffer_data = (unsigned char *) malloc(sizeof(unsigned char) * s_recv_data);

  /*
   * Setup send_buffer
   */
  int s_data = (int) _gnum_from_hv->s_data;
  for(int i_part = 0; i_part < _gnum_from_hv->n_part; ++i_part){

    unsigned char* _part_data = (unsigned char *) _gnum_from_hv->part_hdata[i_part];

    int idx = 0;
    for(int ielt = 0; ielt < _gnum_from_hv->n_elts[i_part]; ++ielt){

      PDM_g_num_t g_key = (PDM_g_num_t) _gnum_from_hv->part_hkeys[i_part][ielt];
      printf(" Search for :: %d\n", (int)g_key);
      int t_rank = PDM_binary_search_gap_long(g_key, _gnum_from_hv->distribution, _gnum_from_hv->n_rank+1);

      printf(" Found in t_rank :: %d\n", (int)t_rank);

      /* Send key and stride */
      int idx_send = i_key_send[t_rank]+n_key_send[t_rank]++;
      send_buffer_keys[idx_send] = g_key;
      send_buffer_stri[idx_send] = _gnum_from_hv->part_hstri[i_part][ielt];

      /* Send data */
      int n_data = s_data * _gnum_from_hv->part_hstri[i_part][ielt];
      int shift  = n_data_send[t_rank];
      printf(" n_data = %d | shift = %d | idx = %d \n", n_data, shift, idx);
      for(int i_data = 0; i_data < n_data; ++i_data) {
        int shift_tot = ( i_data_send[t_rank] + shift)*s_data;
        send_buffer_data[shift_tot+i_data] = _part_data[idx*s_data+i_data];
      }
      idx += _gnum_from_hv->part_hstri[i_part][ielt];

      // n_data_send[t_rank] += s_data * _gnum_from_hv->part_hstri[i_part][ielt];
      n_data_send[t_rank] += _gnum_from_hv->part_hstri[i_part][ielt];

    }
  }
  // abort();

  int* send_buffer_data_int = (int*) send_buffer_data;
  printf("send_buffer_data_int:: ");
  for(int i = 0; i < s_send_data/s_data; ++i){
    printf("%d ", send_buffer_data_int[i]);
  }
  printf("\n");


  for(int i = 0; i < _gnum_from_hv->n_rank+1; ++i){
    i_data_send[i] = i_data_send[i] * s_data;
    i_data_recv[i] = i_data_recv[i] * s_data;
    n_data_send[i] = n_data_send[i] * s_data;
    n_data_recv[i] = n_data_recv[i] * s_data;
  }


  printf("i_key_send:: ");
  for(int i = 0; i < _gnum_from_hv->n_rank+1; ++i){
    printf("%d ", i_key_send[i]);
  }
  printf("\n");
  printf("n_key_send:: ");
  for(int i = 0; i < _gnum_from_hv->n_rank; ++i){
    printf("%d ", n_key_send[i]);
  }
  printf("\n");
  printf("i_key_recv:: ");
  for(int i = 0; i < _gnum_from_hv->n_rank+1; ++i){
    printf("%d ", i_key_recv[i]);
  }
  printf("\n");
  printf("n_key_recv:: ");
  for(int i = 0; i < _gnum_from_hv->n_rank; ++i){
    printf("%d ", n_key_recv[i]);
  }
  printf("\n");


  printf("i_data_send:: ");
  for(int i = 0; i < _gnum_from_hv->n_rank+1; ++i){
    printf("%d ", i_data_send[i]);
  }
  printf("\n");
  printf("n_data_send:: ");
  for(int i = 0; i < _gnum_from_hv->n_rank; ++i){
    printf("%d ", n_data_send[i]);
  }
  printf("\n");
  printf("i_data_recv:: ");
  for(int i = 0; i < _gnum_from_hv->n_rank+1; ++i){
    printf("%d ", i_data_recv[i]);
  }
  printf("\n");
  printf("n_data_recv:: ");
  for(int i = 0; i < _gnum_from_hv->n_rank; ++i){
    printf("%d ", n_data_recv[i]);
  }
  printf("\n");

  /*
   * Exchange
   */
  printf(" Exchnage 1 \n");
  PDM_MPI_Alltoallv(send_buffer_keys, n_key_send, i_key_send, PDM_MPI_INT,
                    recv_buffer_keys, n_key_recv, i_key_recv, PDM_MPI_INT, _gnum_from_hv->comm);

  printf(" Exchnage 2 \n");
  PDM_MPI_Alltoallv(send_buffer_stri, n_key_send, i_key_send, PDM_MPI_INT,
                    recv_buffer_stri, n_key_recv, i_key_recv, PDM_MPI_INT, _gnum_from_hv->comm);


  printf(" Exchnage 1 end \n");
  printf(" Exchnage 2 \n");
  PDM_MPI_Alltoallv(send_buffer_data, n_data_send, i_data_send, PDM_MPI_BYTE,
                    recv_buffer_data, n_data_recv, i_data_recv, PDM_MPI_BYTE, _gnum_from_hv->comm);

  printf(" Exchnage 2 end \n");

  if(0 == 0){
    printf("send_buffer_keys:: ");
    for(int i = 0; i < s_send_keys; ++i){
      printf("%d ", send_buffer_keys[i]);
    }
    printf("\n");
    printf("recv_buffer_keys:: ");
    for(int i = 0; i < s_recv_keys; ++i){
      printf("%d ", recv_buffer_keys[i]);
    }
    printf("\n");
    printf("send_buffer_stri:: ");
    for(int i = 0; i < s_send_keys; ++i){
      printf("%d ", send_buffer_stri[i]);
    }
    printf("\n");
    printf("recv_buffer_stri:: ");
    for(int i = 0; i < s_recv_keys; ++i){
      printf("%d ", recv_buffer_stri[i]);
    }
    printf("\n");
  }

  int* recv_buffer_data_int = (int*) recv_buffer_data;
  printf("recv_buffer_data_int:: ");
  for(int i = 0; i < s_recv_data/s_data; ++i){
    printf("%d ", recv_buffer_data_int[i]);
  }
  printf("\n");
  for(int i = 0; i < s_recv_data/s_data; ++i){
    printf("recv_buffer_data_int[%d] = %d \n ", i, recv_buffer_data_int[i]);
  }
  /*
   * Rebuild a total stride
   */
  int tmp1 = recv_buffer_stri[0];
  recv_buffer_stri[0] = 0;
  for(int i = 0; i < s_recv_keys; ++i){
    int tmp2 = recv_buffer_stri[i+1];
    recv_buffer_stri[i+1] = recv_buffer_stri[i] + tmp1;
    tmp1 = tmp2;
  }

  if(0 == 0){
    printf("recv_buffer_stri:: ");
    for(int i = 0; i < s_recv_keys+1; ++i){
      printf("%d ", recv_buffer_stri[i]);
    }
    printf("\n");
  }

  /*
   * Generate global numbering from the block_data
   */
  // PDM_generate_global_id_from();
  int* order = (int *) malloc( sizeof(int) * s_recv_keys);
  for(int i = 0; i < s_recv_keys; ++i){
    order[i] = i;
  }

  _gnum_from_hv->fcompare = my_compare;
  _gnum_from_hv->fequal   = my_equal;

  user_defined_sort_t* my_struct = (user_defined_sort_t *) malloc( sizeof(user_defined_sort_t) );
  my_struct->idx = recv_buffer_stri;
  my_struct->arr = (int *) recv_buffer_data;


  PDM_sort_long_s(order, s_recv_keys, _gnum_from_hv->fcompare, (void*) my_struct);


  /*
   * Panic verbose
   */
  if(0 == 0 ){
    printf("order = ");
    for(int i = 0; i < s_recv_keys; ++i){
      printf("%d ", order[i]);
    }
    printf("\n");


  }

  if(0 == 0 ){
    printf("order = ");
    for(int i = 0; i < s_recv_keys; ++i){
      printf("%d --> ", (int)order[i]);
      int j   = order[i];
      for(int k = my_struct->idx[j]; k < my_struct->idx[j+1]; ++k ){
        printf(" %d ", my_struct->arr[k]);
      }
      printf("\n");
    }
    printf("\n");
  }

  /*
   * Use operator == to have an global numbering
   */
  PDM_g_num_t* blk_ln_to_gn = (PDM_g_num_t*) malloc( sizeof(PDM_g_num_t*) * s_recv_keys);
  PDM_g_num_t next_id = 0;
  PDM_g_num_t n_id    = 0;
  PDM_g_num_t last_id = order[0];
  for(int i = 0; i < s_recv_keys; ++i){
    printf(" generate g_id :: %d \n", i);
    if(_gnum_from_hv->fequal(&order[i], &last_id, (void*) my_struct)){
      printf(" \t Cas 1 :: order[%d] = %d | next_id : %d\n", i, order[i], next_id);
      blk_ln_to_gn[order[i]] = next_id;
    } else {
      next_id++;
      n_id++;
      printf(" \t Cas 2 :: order[%d] = %d | next_id : %d\n", i, order[i], next_id);
      blk_ln_to_gn[order[i]] = next_id;
      last_id = order[i];
    }
  }

  printf("n_id   :: %d \n", n_id);
  printf("next_id:: %d \n", next_id);

  /*
   * Panic verbose
   */
  if(0 == 0 ){
    printf("blk_ln_to_gn = ");
    for(int i = 0; i < s_recv_keys; ++i){
      printf("%d --> ", (int)blk_ln_to_gn[i]);
      int j   = blk_ln_to_gn[i];
      for(int k = my_struct->idx[j]; k < my_struct->idx[j+1]; ++k ){
        printf(" %d ", my_struct->arr[k]);
      }
      printf("\n");
    }
    printf("\n");
  }

  /*
   * Reverse all_to_all exchange in order to remap global id on current partition
   */
  // MPI_Alltoallv(blk_gid , n_key_recv, i_key_recv, PDM__PDM_MPI_G_NUM,
  //               part_gid, n_key_send, i_key_send, PDM__PDM_MPI_G_NUM, _gnum_from_hv->comm);

  free(n_key_send);
  free(n_key_recv);
  free(n_data_send);
  free(n_data_recv);
  free(i_key_send);
  free(i_key_recv);
  free(i_data_send);
  free(i_data_recv);
  free(send_buffer_keys);
  free(recv_buffer_keys);
  free(send_buffer_data);
  free(recv_buffer_data);
  free(send_buffer_stri);
  free(recv_buffer_stri);
  free(blk_ln_to_gn);

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
 const size_t       s_data,
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

  _gnum_from_hv->s_data      = s_data;
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
 const size_t       *s_data,
 const PDM_MPI_Fint *fcomm,
       int          *id
)
{
  const PDM_MPI_Comm c_comm = PDM_MPI_Comm_f2c (*fcomm);

  *id = PDM_gnum_from_hash_values_create (*n_part, (PDM_bool_t) *equilibrate, *s_data, c_comm);
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
