#ifndef __PDM_PART_GRAPH_COMM_PRIV_H__
#define __PDM_PART_GRAPH_COMM_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_timer.h"
#include "pdm_error.h"
#include "pdm_printf.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation csback to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _PDM_ol_t
 * \brief  Overlay type
 *
 * _PDM_ol_t defines a overlaying structure
 *
 */
struct _pdm_part_comm_graph_t {

  PDM_MPI_Comm comm;
  int          n_part;
  int          n_g_part;

  int         *n_entity_graph;
  int        **pentity_graph;

  int         *send_idx;
  int         *recv_idx;
  int         *send_n;
  int         *recv_n;

  int        **part_to_send_buffer;
  int        **part_to_recv_buffer;

  int        **bound_owner;

  /* Communication variability */
  int          n_active_rank_send;
  int          n_active_rank_recv;
  int         *active_rank_send;
  int         *active_rank_recv;

  /* Asynchonous */

};


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_GRAPH_COMM_PRIV_H__ */
