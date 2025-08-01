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

// MPI_Topo_test --> Permet de determiner si Alltoall ou neigtbor_alltoall automatiquement
//

// Pour ne pas avoir a syncrhoniser les rang pour l'utilisateur on peut faire :
//   PDM_exchange_helper_exch_iinit (donc non bloquant + synchrone)
//   PDM_exchange_helper_exch_warm_up (table de hash sur les tag enregistrés + check + création des windows / RMA / Persistent )
//   Attention en GPU c'est plus subtile car on pourra pas faire l'alloc dans paradigm ?

// Async + persistent

//
// Interface synchone à faire --> Toujours pratique
//   --> Nous permet de facilement delegué avec différents type d'échange
//   Il faut réfléchir a la stride variable --> Besoin d'échanger avant pour les tailles
//   En persistant pas trop de sens
// Si le comm est defini par un dist_graph_comm, on doit recuper n_send_rank, recv_rank, n_recv_rank, recv_rank via le comm
//   PDM_MPI_Dist_graph_create_adjacent
// Cas tordu : Le comm est un graph, mais on veut faire du p2p
// Persistent
int
PDM_exchange_helper_exch_init2
(
  PDM_exchange_helper_t *exch_helper,
  PDM_mpi_comm_kind_t    k_comm,
  int                    cst_stride,
  size_t                 s_data,
  int                    tag,
  int                   *send_buffer_idx,
  void                  *send_buffer,
  int                   *recv_buffer_idx,
  void                  *recv_buffer,
  PDM_ownership_t        ownership
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
  PDM_exchange_helper_t  *exch_helper,
  PDM_mpi_comm_kind_t     kcomm,
  int                     cst_stride,
  size_t                  s_data,
  int                     tag,
  void                  **buffer, // A voir
  PDM_ownership_t         ownership
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
