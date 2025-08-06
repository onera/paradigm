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



// Faire le exch classique pour gerer les variantes directement cachés : P2P/RMA/Collective
void
PDM_exchange_helper_exch
(
  PDM_exchange_helper_t *exch_helper,
  PDM_mpi_comm_kind_t    k_comm,
  int                    cst_stride,
  size_t                 s_data,
  int                   *send_idx,
  int                   *send_n,
  void                  *send_buffer,
  int                   *recv_idx,
  int                   *recv_n,
  void                  *recv_buffer
);

int
PDM_exchange_helper_iexch
(
  PDM_exchange_helper_t *exch_helper,
  PDM_mpi_comm_kind_t    k_comm,
  size_t                 s_data,
  int                    cst_stride,
  int                   *send_idx,
  int                   *send_n,
  void                  *send_buffer,
  int                   *recv_idx,
  int                   *recv_n,
  void                  *recv_buffer
);

int
PDM_exchange_helper_exch_init
(
  PDM_exchange_helper_t *exch_helper,
  PDM_mpi_comm_kind_t    k_comm,
  size_t                 s_data,
  int                    cst_stride,
  int                   *send_idx,
  int                   *send_n,
  void                  *send_buffer,
  int                   *recv_idx,
  int                   *recv_n,
  void                  *recv_buffer
);


/*
 * Idéal pour cwipi si on veut separé les send / recv
 * Doit couvrir les ONE-SIDED
 *
 */
// Persistent
// Pour les RMA : La fonction INIT prepare la window, on doit retourné le buffer !!!!
int
PDM_exchange_helper_exch_one_way_pack_init
(
  PDM_exchange_helper_t *exch_helper,
  PDM_mpi_comm_kind_t    kcomm,
  int                    cst_stride,
  size_t                 s_data,
  int                    tag,
  int                    n_active_rank,
  int                   *active_rank,
  int                   *buffer_idx,
  void                  *buffer,
  PDM_ownership_t        ownership
);

// API pour get le buffer interne envoie / reception ...

//
// Persistent
// Pour les RMA : La fonction INIT prepare la window, on doit retourné le buffer !!!!
int
PDM_exchange_helper_exch_one_way_init
(
  PDM_exchange_helper_t    *exch_helper,
  PDM_exchange_direction_t  direction,
  int                       cst_stride,
  size_t                    s_data,
  int                       n_active_rank,
  int                      *active_rank,
  int                       tag,
  void                     *buffer,
  PDM_ownership_t           ownership
);


void
PDM_exchange_helper_exch_start
(
  PDM_exchange_helper_t *exch_helper,
  int                    request
);


void
PDM_exchange_helper_exch_wait
(
  PDM_exchange_helper_t *exch_helper,
  int                    request
);

//
// Attention :
//   - MPI_Request_free(&request) necessaire a la fin des persistente
//
//  MPI_Wait(&request, &status) (ou MPI_Test(&request, &flag, &status)) :
//   Ces appels sont utilisés pour compléter la communication.
//   Quand MPI_Wait se termine ou MPI_Test indique que l'opération est finie (flag == true), la communication est considérée comme terminée.
//   La requête repasse alors dans un état inactif. Elle n'est pas encore libérée de la mémoire, mais elle n'est plus associée à une communication en cours et peut être redémarrée par un autre MPI_Start.
// En interne de paradim il faut check qu'on repasse pas la requete a NULL car elle pourrait restart après ....
// RMA


// Pour les RMA : La fonction INIT prepare la window, on doit retourné le buffer !!!!


void
PDM_exchange_helper_exch_free
(
  PDM_exchange_helper_t *exch_helper,
  int                    request_id
);

void
PDM_exchange_helper_free
(
  PDM_exchange_helper_t *exch_helper
);

#ifdef  __cplusplus
}
#endif

#endif  /* __PDM_EXCHANGE_HELPER_H__ */
