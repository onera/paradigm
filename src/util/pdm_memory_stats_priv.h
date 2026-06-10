#ifndef __PDM_MEMORY_STATS_PRIV_H__
#define __PDM_MEMORY_STATS_PRIV_H__

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
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

/**
 * \struct _pdm_extract_part_t
 * \brief  Define a partition mesh. Arrays are shared
 *
 */

struct _pdm_memory_stats_t
{
  PDM_MPI_Comm           comm;
  int                    n_memory_snapshot;

  char**                 snapshot_name;
  long*                  curr_real_mem;
  long*                  peak_real_mem;
  long*                  curr_virt_mem;
  long*                  peak_virt_mem;

};


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MEMORY_STATS_PRIV_H__ */
