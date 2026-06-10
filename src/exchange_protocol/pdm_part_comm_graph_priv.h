#ifndef __PDM_PART_GRAPH_COMM_PRIV_H__
#define __PDM_PART_GRAPH_COMM_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_exchange_helper.h"
#include "pdm_timer.h"
#include "pdm_error.h"
#include "pdm_printf.h"

#ifdef __cplusplus
extern "C" {
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

struct _pdm_part_comm_graph_t {

  PDM_MPI_Comm      comm;
  int               n_part;
  int               n_g_part;
  int               is_signed;

  int              *n_entity_graph;
  int             **pentity_graph;
  int               nuplet_size;
  int             **pentity_nuplet;

  PDM_ownership_t   owner_graph;
  PDM_ownership_t   owner_nuplet;

  int              *send_idx;
  int              *recv_idx;
  int              *send_n;
  int              *recv_n;

  int              *active_send_idx;
  int              *active_recv_idx;
  int              *active_send_n;
  int              *active_recv_n;

  int             **part_to_send_buffer;
  int             **part_to_recv_buffer;

  int             **bound_owner;

  /* Communication variability */
  int               n_active_rank_send;
  int               n_active_rank_recv;
  int              *active_rank_send;
  int              *active_rank_recv;

  /* Asynchronous */
  PDM_exchange_helper_t *exch_h;

};


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_GRAPH_COMM_PRIV_H__ */
