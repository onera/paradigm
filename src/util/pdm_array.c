/*============================================================================
 * Some utils for manupulating arrays
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_array.h"
#include "pdm_config.h"
#include "pdm_priv.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Definition des types locaux
 *============================================================================*/

/*============================================================================
 * Definition des fonctions locales
 *============================================================================*/

/*============================================================================
 * Definition des fonctions publiques
 *============================================================================*/

/* Utils functions creating new arrays */


/*
 * Array copy
 */
int* PDM_array_copy_int(const int *array, const int size) {
  assert (size >= 0);
  int *array_cp = NULL;
  PDM_malloc(array_cp, size, int);
  assert (array != NULL);
  for (int i = 0; i < size; i++) array_cp[i] = array[i];
  return array_cp;
}
PDM_g_num_t* PDM_array_copy_gnum(const PDM_g_num_t *array, const int size) {
  assert (size >= 0);
  PDM_g_num_t *array_cp = NULL;
  PDM_malloc(array_cp, size, PDM_g_num_t);
  assert (array != NULL);
  for (int i = 0; i < size; i++) array_cp[i] = array[i];
  return array_cp;
}

/*
 * Allocate a new array of size size and fill it with 0
*/
int* PDM_array_zeros_int(const int size) {
  assert (size >= 0);
  int *array;
  PDM_malloc(array, size, int);
  assert (array != NULL);
  for (int i = 0; i < size; i++) array[i] = 0;
  return array;
}
PDM_g_num_t* PDM_array_zeros_gnum(const int size) {
  assert (size >= 0);
  PDM_g_num_t *array = NULL;
  PDM_malloc(array, size, PDM_g_num_t);
  assert (array != NULL);
  for (int i = 0; i < size; i++) array[i] = 0;
  return array;
}
double* PDM_array_zeros_double(const int size) {
  assert (size >= 0);
  double *array = NULL;
  PDM_malloc(array, size, double);
  assert (array != NULL);
  for (int i = 0; i < size; i++) array[i] = 0.;
  return array;
}

/*
 * Allocate a new array of size size and fill it with value
*/
int* PDM_array_const_int(const int size, const int value) {
  assert (size >= 0);
  int *array;
  PDM_malloc(array, size, int);
  assert (array != NULL);
  for (int i = 0; i < size; i++) array[i] = value;
  return array;
}
PDM_g_num_t* PDM_array_const_gnum(const int size, const PDM_g_num_t value) {
  assert (size >= 0);
  PDM_g_num_t *array;
  PDM_malloc(array, size, PDM_g_num_t);
  assert (array != NULL);
  for (int i = 0; i < size; i++) array[i] = value;
  return array;
}
double* PDM_array_const_double(const int size, const double value) {
  assert (size >= 0);
  double *array;
  PDM_malloc(array, size, double);
  assert (array != NULL);
  for (int i = 0; i < size; i++) array[i] = value;
  return array;
}

/*
 * Allocate a new array (size=size) and fill it with identity
*/
int* PDM_array_new_range_int(const int size) {
  int *array;
  PDM_malloc(array, size, int);
  for (int i = 0; i < size; ++i) {
    array[i] = i;
  }
  return array;
}


/*
 * Allocate a new int array starting at start, incremented by step and finishing at end max
*/
int* PDM_array_new_arange_int(const int start, const int end, const int step) {
  int size = (end-start)/step+1;
  int *array;
  PDM_malloc(array, size, int);
  for (int i = 0; i < size; ++i) {
    array[i] = start+i*step;
    
int* PDM_array_new_range_with_step_int(const int size, const int step) {
  int* array = NULL;;
  PDM_malloc(array, size, int);
  for (int i = 0; i < size; ++i) {
    array[i] = i*step;
  }
  return array;
}

/*
 * Allocate a new index array (size=size+1) and fill it from
 * the size array size_array (which must be of size size)
*/
int* PDM_array_new_idx_from_sizes_int(const int *size_array, const int size) {
  int *idx_array;
  PDM_malloc(idx_array, size+1, int);
  idx_array[0] = 0;
  for (int i = 0; i < size; i++) idx_array[i+1] = idx_array[i] + size_array[i];
  return idx_array;
}
PDM_g_num_t* PDM_array_new_idx_from_sizes_gnum(const int *size_array, const int size) {
  PDM_g_num_t *idx_array;
  PDM_malloc(idx_array, size+1, PDM_g_num_t);
  idx_array[0] = 0;
  for (int i = 0; i < size; i++) idx_array[i+1] = idx_array[i] + size_array[i];
  return idx_array;
}

int* PDM_array_new_idx_from_const_stride_int(const int stride, const int size) {
  int *idx_array;
  PDM_malloc(idx_array, size+1, int);
  idx_array[0] = 0;
  for (int i = 0; i < size; i++) idx_array[i+1] = idx_array[i] + stride;
  return idx_array;
}


/*
 * Allocate a new size array (size=size) and fill it from
 * the index array (which must be of size `size+1`)
*/
int* PDM_array_new_size_from_idx_int(const int *idx_array, const int size) {
  int *size_array = NULL;
  PDM_malloc(size_array, size, int);
  for (int i = 0; i < size; i++) size_array[i] = idx_array[i+1] - idx_array[i];
  return size_array;
}


/* Utils functions compararing arrays */

/*
 * Return 1 if the two arrays are equal, 0 otherwise
*/
int PDM_array_are_equal_int(const int *array1, const int *array2, const int size) {
  for (int i = 0; i < size; i++) {
    if (array1[i] != array2[i])
      return 0;
  }
  return 1;
}
int PDM_array_are_equal_gnum(const PDM_g_num_t *array1, const PDM_g_num_t *array2, const int size) {
  for (int i = 0; i < size; i++) {
    if (array1[i] != array2[i])
      return 0;
  }
  return 1;
}

/* Utils functions for reducing arrays */
int PDM_array_max_int(const int *array, const int size) {
  int max = -INT_MAX;
  for (int i = 0; i < size; i++) {
    if (array[i] > max)
      max = array[i];
  }
  return max;
}
int PDM_array_min_int(const int *array, const int size) {
  int min = INT_MAX;
  for (int i = 0; i < size; i++) {
    if (array[i] < min)
      min = array[i];
  }
  return min;
}

PDM_g_num_t PDM_array_max_gnum(const PDM_g_num_t *array, const int size) {
  PDM_g_num_t max = -PDM_G_NUM_MAX;
  for (int i = 0; i < size; i++) {
    if (array[i] > max)
      max = array[i];
  }
  return max;
}
PDM_g_num_t PDM_array_min_gnum(const PDM_g_num_t *array, const int size) {
  PDM_g_num_t min = PDM_G_NUM_MAX;
  for (int i = 0; i < size; i++) {
    if (array[i] < min)
      min = array[i];
  }
  return min;
}


/* Utils functions modifying arrays*/


/*
 * Fill the given array of size size with the value value
*/
void PDM_array_reset_int(int *array, const int size, const int value) {
    for (int i = 0; i < size; i++) array[i] = value;
}
void PDM_array_reset_gnum(PDM_g_num_t *array, const int size, const PDM_g_num_t value) {
    for (int i = 0; i < size; i++) array[i] = value;
}

/*
 * Fill the index array idx_array (preallocated at size=size+1) from the size array
 * size_array (of size size)
*/
void PDM_array_idx_from_sizes_int(const int *size_array, const int size, int *idx_array) {
  idx_array[0] = 0;
  for (int i = 0; i < size; i++) idx_array[i+1] = idx_array[i] + size_array[i];
}
void PDM_array_idx_from_sizes_gnum(const int *size_array, const int size, PDM_g_num_t *idx_array) {
  idx_array[0] = 0;
  for (int i = 0; i < size; i++) idx_array[i+1] = idx_array[i] + size_array[i];
}

/*
 * Accumulate the values of the input array array of size size, ie
 * array[0] remains unchanged, and array[i] = sum_{j=0...i} array[j]
*/
void PDM_array_accumulate_int(int *array, const int size) {
  for (int i = 1; i < size; i++) array[i] = array[i] + array[i-1];
}
void PDM_array_accumulate_gnum(PDM_g_num_t *array, const int size) {
  for (int i = 1; i < size; i++) array[i] = array[i] + array[i-1];
}


/* Arrays algorithms */


/*
 * Count the number of occurence of each color in a color_array.
 * Fill n_per_col, which must be pre allocated at size n_col
 * Colors are expected to start at 0.
*/
void PDM_array_count_per_col_int
(
 const int  n_col,
 const int  n_elem,
 const int *elem_col,
 int       *n_per_col
)
{
  PDM_array_reset_int(n_per_col, n_col, 0);
  for (int i_elem = 0; i_elem < n_elem; i_elem++)
    n_per_col[elem_col[i_elem]]++;
  int sum_elem = 0;
  for (int i_col = 0; i_col < n_col; i_col++)
    sum_elem += n_per_col[i_col];
  assert(sum_elem == n_elem);
}

/*
 * From an array of colors, compute the index order to be used to
 * read the array in increasing color order.
 * Fill ordered_idx, which must be pre allocated at size n_col+1, and
 * ordered, which must be preallocated at size n_elem.
 * Colors are expected to start at 0.
*/
void PDM_array_repart_per_col_int
(
 const int  n_col,
 const int  n_elem,
 const int *elem_col,
 int       *ordered_idx,
 int       *ordered
)
{
  int *count;
  PDM_malloc(count, n_col, int);
  PDM_array_count_per_col_int(n_col, n_elem, elem_col, count);

  ordered_idx[0] = 0;
  for (int i_col=0; i_col < n_col; i_col++) {
    ordered_idx[i_col+1] = ordered_idx[i_col] + count[i_col];
    count[i_col] = 0;
  }
  assert(ordered_idx[n_col] == n_elem);
  for (int i_elem = 0; i_elem < n_elem; i_elem++) {
    int col = elem_col[i_elem];
    ordered[ordered_idx[col] + count[col]++] = i_elem;
  }
  PDM_free(count);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
