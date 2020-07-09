#ifndef __PDM_TIMER_CUH__
#define __PDM_TIMER_CUH__

/*----------------------------------------------------------------------------*/

#include "pdm_timer.h"

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

/*----------------------------------------------------------------------------
 * Suspend la mesure du temps ecoule et incremente le temps ecoule
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

__device__ void PDM_timer_hang_on_GPU(PDM_timer_t *timer);

/*----------------------------------------------------------------------------
 * Reprend la mesure du temps ecoule
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

__device__ void PDM_timer_resume_GPU(PDM_timer_t *timer);

/*----------------------------------------------------------------------------
 * Retourne le temps elaps en secondes
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

__device__ double PDM_timer_elapsed_GPU(PDM_timer_t *timer);



void PDM_timer_free_GPU(PDM_timer_t *timer);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FICHIER_SEQ_CUH__ */
