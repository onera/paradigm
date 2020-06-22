#ifndef __PDM_TIMER_CUH__
#define __PDM_TIMER_CUH__

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Definition des types
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure de mesure des temps d'execution
 *----------------------------------------------------------------------------*/

typedef struct _pdm_timer_t PDM_timer_t;

/*============================================================================
 * Interfaces des fonctions publiques
 *============================================================================*/

void PDM_timer_free_GPU(PDM_timer_t *timer);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FICHIER_SEQ_CUH__ */
