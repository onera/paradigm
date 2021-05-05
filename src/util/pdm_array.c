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

#include "pdm_array.h"
#include "pdm_config.h"

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

/*
 * Create a new array int array of size size and fill it with 0
*/
inline int* PDM_array_zeros_int(const int size) {
  int *array = (int *) malloc(size * sizeof(int));
  assert (array != NULL);
  for (int i = 0; i < size; i++)
    array[i] = 0;
  return array;
}

inline int* PDM_array_const_int(const int size, const int value) {
  int *array = (int *) malloc(size * sizeof(int));
  assert (array != NULL);
  for (int i = 0; i < size; i++)
    array[i] = value;
  return array;
}

inline void PDM_array_reset_int(int *array, const int size, const int value) {
    for (int i = 0; i < size; i++)
        array[i] = value;
}

inline PDM_g_num_t* PDM_array_const_gnum(const int size, const PDM_g_num_t value) {
  PDM_g_num_t *array = (PDM_g_num_t *) malloc(size * sizeof(PDM_g_num_t));
  assert (array != NULL);
  for (int i = 0; i < size; i++)
    array[i] = value;
  return array;
}

inline void PDM_array_reset_gnum(PDM_g_num_t *array, const int size, const PDM_g_num_t value) {
    for (int i = 0; i < size; i++)
        array[i] = value;
}

inline int* PDM_array_new_idx_from_sizes_int(const int *size_array, const int size) {
  int *idx_array = (int *) malloc((size+1) * sizeof(int));
  idx_array[0] = 0;
  for (int i = 0; i < size; i++)
    idx_array[i+1] = idx_array[i] + size_array[i];
  return idx_array;
}

inline void PDM_array_idx_from_sizes_int(const int *size_array, const int size, int *idx_array) {
  idx_array[0] = 0;
  for (int i = 0; i < size; i++)
    idx_array[i+1] = idx_array[i] + size_array[i];
}

inline PDM_g_num_t* PDM_array_new_idx_from_sizes_gnum(const int *size_array, const int size) {
  PDM_g_num_t *idx_array = (PDM_g_num_t *) malloc((size+1) * sizeof(PDM_g_num_t));
  idx_array[0] = 0;
  for (int i = 0; i < size; i++)
    idx_array[i+1] = idx_array[i] + size_array[i];
  return idx_array;
}

inline void PDM_array_idx_from_sizes_gnum(const int *size_array, const int size, PDM_g_num_t *idx_array) {
  idx_array[0] = 0;
  for (int i = 0; i < size; i++)
    idx_array[i+1] = idx_array[i] + size_array[i];
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
