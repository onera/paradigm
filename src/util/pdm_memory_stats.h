#ifndef __PDM_MEMORY_STATS_H__
#define __PDM_MEMORY_STATS_H__

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _pdm_memory_stats_t PDM_memory_stats_t;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

PDM_memory_stats_t*
PDM_memory_stats_create
(
 int          n_memory_snapshot,
 PDM_MPI_Comm comm
);


void
PDM_memory_stats_add
(
 PDM_memory_stats_t *ms,
 int                 i_snapshot,
 const char         *name
);

void
PDM_memory_stats_log
(
 PDM_memory_stats_t* ms
);


void
PDM_memory_stats_free
(
 PDM_memory_stats_t* ms
);



#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MEMORY_STATS_H__ */
