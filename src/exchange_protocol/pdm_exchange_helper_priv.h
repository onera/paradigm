#ifndef __PDM_EXCHANGE_HELPER_PRIV_H__
#define __PDM_EXCHANGE_HELPER_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

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
  EXCHANGE_HELPER_STATUS_ONGOING
} _exch_helper_status_t;

/**
 * \struct _pdm_exchange_helper_t
 *
 * \brief  Helper for data_exchange
 *
 */

struct _pdm_exchange_helper_t {

  PDM_MPI_Comm             comm;               /*!< MPI communicator */
  int                      n_request;
  int                      topo_kind;

  long                     max_tag;
  int                      seed_tag;
  int                      next_tag;

  _exch_helper_status_t   *requests_status;
  int                     *is_persistent;
  int                     *n_sub_requests;
  PDM_MPI_Request        **sub_requests;

  void                   **send_buffer;
  void                   **recv_buffer;

  /* Very specific to RMA, allowing RMA + API Persistent */
  int                    **recv_n;
  int                    **recv_idx;

  PDM_MPI_Win             *win_send;
  PDM_MPI_Win             *win_recv;
  PDM_MPI_Group           *group_send;
  PDM_MPI_Group           *group_recv;
  int                    **target_disp;

  PDM_mpi_comm_kind_t     *k_comm;
  PDM_stride_t            *t_stride;
  size_t                  *s_data;
  int                     *cst_stride;
  PDM_MPI_Datatype        *mpi_type;


  /* High-User helper to keep pointer */
  int                   ***p_send_stride;
  void                  ***p_send_data;
  int                   ***p_recv_stride;
  void                  ***p_recv_data;

  int                    **d_send_stride;
  void                   **d_send_data;
  int                    **d_recv_stride;
  void                   **d_recv_data;


};


/*=============================================================================
 * Static global variables
 *============================================================================*/


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_EXCHANGE_HELPER_PRIV_H__ */
