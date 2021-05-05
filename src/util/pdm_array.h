#ifndef __PDM_ARRAY_H__
#define __PDM_ARRAY_H__

/*----------------------------------------------------------------------------*/
#include "pdm.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Definition des types
 *============================================================================*/

/*============================================================================
 * Interfaces des fonctions publiques
 *============================================================================*/

int* PDM_array_zeros_int(const int size);
int* PDM_array_const_int(const int size, const int value);
void PDM_array_reset_int(int *array, const int size, const int value);

PDM_g_num_t* PDM_array_const_gnum(const int size, const PDM_g_num_t value);
void PDM_array_reset_gnum(PDM_g_num_t *array, const int size, const PDM_g_num_t value);

int* PDM_array_new_idx_from_sizes_int(const int *size_array, const int size);
void PDM_array_idx_from_sizes_int(const int *size_array, const int size, int *idx_array);
PDM_g_num_t* PDM_array_new_idx_from_sizes_gnum(const int *size_array, const int size);
void PDM_array_idx_from_sizes_gnum(const int *size_array, const int size, PDM_g_num_t *idx_array);


void PDM_array_count_per_col_int(const int n_col, const int n_elem, const int *elem_col, int *n_per_col);
void PDM_array_repart_per_col_int
(
 const int   n_col,
 const int   n_elem,
 const int *elem_col,
 int       *ordered_idx,
 int       *ordered
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FICHIER_SEQ_H__ */
