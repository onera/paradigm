#ifndef __PDM_EXCHANGE_HELPER_PRIV_H__
#define __PDM_EXCHANGE_HELPER_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

typedef enum {
  EXCHANGE_HELPER_STATUS_FREE,
  EXCHANGE_HELPER_STATUS_READY,
  EXCHANGE_HELPER_STATUS_ONGOING,
} _exch_helper_status_t;

/**
 * \struct _pdm_partgnum1_partgnum2_t
 *
 * \brief  Data transfer from partitions to blocks
 *
 */

struct _pdm_exchange_helper_t {

  PDM_MPI_Comm        comm;                         /*!< MPI communicator */
  int                 n_request;

  _exch_helper_status_t  *requests_status;
  int                    *n_sub_requests;
  PDM_MPI_Request       **sub_requests;

};


/*=============================================================================
 * Static global variables
 *============================================================================*/


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_EXCHANGE_HELPER_PRIV_H__ */
