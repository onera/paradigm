/*
 * \file
 */

#ifndef __PDM_TIMER_H__
#define __PDM_TIMER_H__

/*----------------------------------------------------------------------------*/
#include "pdm_mpi.h"


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

PDM_timer_t*
PDM_timer_create
(
  PDM_MPI_Comm comm
);


void
PDM_timer_start
(
        PDM_timer_t *timer,
  const char        *name,
        int          force_synchro
);


void
PDM_timer_end
(
        PDM_timer_t *timer,
  const char        *name,
        int          force_synchro
);

void
PDM_timer_gather
(
  PDM_timer_t *timer
);

void
PDM_timer_gather_dump
(
  PDM_timer_t *timer,
  char        *filename
);

char* PDM_timer_get_report_string(PDM_timer_t *timer, int mode);

void
PDM_timer_dump
(
        PDM_timer_t *timer
);

void
PDM_timer_free
(
        PDM_timer_t *timer
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FICHIER_SEQ_H__ */
