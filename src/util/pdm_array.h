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

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FICHIER_SEQ_H__ */
