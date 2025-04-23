/*
 * \file
 */

#ifndef __PDM_EXCHANGE_HELPER_H__
#define __PDM_EXCHANGE_HELPER_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stddef.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef  __cplusplus
extern "C" {
#endif

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \struct PDM_exchange_helper_t
 * \brief  Helper struct to manage asynchronous and persistent communication
 *
 */

typedef struct _pdm_exchange_helper_t PDM_exchange_helper_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

PDM_exchange_helper_t *
PDM_exchange_helper_create
(
 const PDM_MPI_Comm    comm,
       int             n_request_init
);

// Async + persistent


#ifdef  __cplusplus
}
#endif

#endif  /* __PDM_EXCHANGE_HELPER_H__ */
