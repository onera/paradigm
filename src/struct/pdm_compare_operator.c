/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_sort.h"
#include "pdm_quick_sort.h"
#include "pdm_compare_operator.h"


#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Compare operator for connectivities
 *
 */
int
PDM_operator_compare_string
(
const void* a,
const void* b,
      void* ctxt
)
{
  int i = *(const int *) a;
  int j = *(const int *) b;

  PDM_user_defined_sort* us = (PDM_user_defined_sort*) ctxt;

  // log_trace("PDM_operator_compare_string::  %d %d - idx -> %d %d \n",  i, j, us->idx[i], us->idx[j]);

  int ni = us->idx[i+1] - us->idx[i];
  int nj = us->idx[j+1] - us->idx[j];

  // log_trace("PDM_operator_compare_string:: %d %d - %d %d \n", i, j, ni, nj);

  char* arr_i = (char *) &us->arr[us->idx[i]*sizeof(char)];
  char* arr_j = (char *) &us->arr[us->idx[j]*sizeof(char)];

  char *carr_i = NULL;
  char *carr_j = NULL;
  PDM_malloc(carr_i, ni + 1, char);
  PDM_malloc(carr_j, nj + 1, char);
  for(int k = 0; k < ni; ++k){
    carr_i[k] = (char)arr_i[k];
  }
  for(int k = 0; k < nj; ++k){
    carr_j[k] = (char)arr_j[k];
  }
  carr_i[ni] = '\0';
  carr_j[nj] = '\0';
  int i_comp = strcmp(carr_i, carr_j);
  PDM_free(carr_i);
  PDM_free(carr_j);
  if(i_comp >= 0){
    return 0;
  } else {
    return 1;
  }

  // Implementation different
  // if(ni < nj){
  //   return 0;
  // } else {
  //   char* arr_i = (char *) &us->arr[us->idx[i]*sizeof(char)];
  //   char* arr_j = (char *) &us->arr[us->idx[j]*sizeof(char)];
  //   int i_comp = strncmp(arr_i, arr_j, ni);
  //   if(i_comp > 0) {
  //     return 0;
  //   } else {
  //     return 1;
  //   }
  // }
}

/**
 *
 * \brief Compare operator for connectivities
 *
 */
int
PDM_operator_compare_connectivity
(
const void* a,
const void* b,
      void* ctxt
)
{
  int i = *(const int *) a;
  int j = *(const int *) b;

  PDM_user_defined_sort* us = (PDM_user_defined_sort*) ctxt;

  // log_trace("PDM_operator_compare_connectivity:: %d %d \n", i, j);
  // log_trace("PDM_operator_compare_connectivity:: %d %d with key %lu %lu \n", i, j, us->key[i], us->key[j]);

  if( us->key[i] == us->key[j]){

    // log_trace("PDM_operator_compare_connectivity:: %d %d \n", i, j);
    int ni = us->idx[i+1] - us->idx[i];
    int nj = us->idx[j+1] - us->idx[j];

    if(ni == nj){

      int* arr_i = (int *) &us->arr[us->idx[i]*sizeof(int)];
      int* arr_j = (int *) &us->arr[us->idx[j]*sizeof(int)];

      /* Dans notre cas on veut sort les entiers avant de les comparers */
      int *sort_arr_i = NULL;
      int *sort_arr_j = NULL;
      PDM_malloc(sort_arr_i, ni ,int);
      PDM_malloc(sort_arr_j, ni ,int);

      for(int k = 0; k < ni; ++k){
        sort_arr_i[k] = arr_i[k];
        sort_arr_j[k] = arr_j[k];
        // log_trace("PDM_operator_compare_connectivity:: %d %d \n", i, j);
      }
      PDM_quick_sort_int(sort_arr_i, 0, ni-1);
      PDM_quick_sort_int(sort_arr_j, 0, ni-1);

      // log_trace("Comparison of ::");
      // for(int k = 0; k < ni; ++k){
      //   log_trace(" %d %d ", sort_arr_i[k], sort_arr_j[k]);
      // }
      // log_trace("\n");

      for(int k = 0; k < ni; ++k){
        if(sort_arr_i[k] < sort_arr_j[k]) {
          PDM_free(sort_arr_i);
          PDM_free(sort_arr_j);
          return 1;
        } else if( sort_arr_i[k] > sort_arr_j[k] ) {
          PDM_free(sort_arr_i);
          PDM_free(sort_arr_j);
          return 0;
        }
      }
      PDM_free(sort_arr_i);
      PDM_free(sort_arr_j);
    } else {
      return ni < nj;
    }
  } else {
    return us->key[i] < us->key[j];
  }
  return 0;
}


int
PDM_compare_unsigned_ordered_nuplets_int
(
  const int size1,
  const int nuplet1[],
  const int size2,
  const int nuplet2[]
)
{
  if (size1 != size2) {
    return 0;
  }

  if (size1 == 2) {
    // Handle pairs separately
    if (nuplet2[0] == nuplet1[0] && nuplet2[1] == nuplet1[1]) {
      return 1;
    }
    else if (nuplet2[0] == nuplet1[1] && nuplet2[1] == nuplet1[0]) {
      return -1;
    }
    else {
      return 0;
    }
  }

  assert(size1 >= 3);

  for (int i = 0; i < size1; i++) {
    if (nuplet2[i] == nuplet1[0]) { // we found a common element ("pivot")
      int direct  = 1;
      int reverse = 1;
      for (int j = 1; j < size1; j++) {
        // Starting from pivot, traverse nuplet 2 in both directions
        // and compare to nuplet 1 element-wise
        if (nuplet2[(i+j)%size1] != nuplet1[j]) {
          direct = 0;
        }
        if (nuplet2[(i+size1-j)%size1] != nuplet1[j]) {
          reverse = 0;
        }
      }

      if (direct && !reverse) {
        return 1;
      } else if (reverse) {
        return -1;
      } else {
        return 0;
      }
    }
  }

  return 0;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
