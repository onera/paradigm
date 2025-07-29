/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_exchange_helper.h"
#include "pdm.h"
#include "pdm_array.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include "pdm_mpi.h"
#include "pdm_exchange_helper_priv.h"
#include "pdm_priv.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/


/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/


static
int
_find_available_request
(
  PDM_exchange_helper_t *exch_helper
)
{

  for(int i = 0; i < exch_helper->n_request; ++i) {
    if(exch_helper->requests_status[i] == EXCHANGE_HELPER_STATUS_FREE) {
      return i;
    }
  }

  // Realloc
  int n_new_request = PDM_MAX(exch_helper->n_request * 2, 1);
  PDM_realloc(exch_helper->requests_status, exch_helper->requests_status, n_new_request, _exch_helper_status_t  );
  PDM_realloc(exch_helper->sub_requests   , exch_helper->sub_requests   , n_new_request, PDM_MPI_Request       *);
  PDM_realloc(exch_helper->n_sub_requests , exch_helper->n_sub_requests , n_new_request, int                    );
  PDM_realloc(exch_helper->send_buffer    , exch_helper->send_buffer    , n_new_request, void                  *);
  PDM_realloc(exch_helper->recv_buffer    , exch_helper->recv_buffer    , n_new_request, void                  *);

  PDM_realloc(exch_helper->t_stride       , exch_helper->t_stride       , n_new_request, PDM_stride_t           );
  PDM_realloc(exch_helper->s_data         , exch_helper->s_data         , n_new_request, size_t                 );
  PDM_realloc(exch_helper->cst_stride     , exch_helper->cst_stride     , n_new_request, int                    );
  PDM_realloc(exch_helper->p_send_stride  , exch_helper->p_send_stride  , n_new_request, int                  **);
  PDM_realloc(exch_helper->p_send_data    , exch_helper->p_send_data    , n_new_request, void                 **);
  PDM_realloc(exch_helper->p_recv_stride  , exch_helper->p_recv_stride  , n_new_request, int                  **);
  PDM_realloc(exch_helper->p_recv_data    , exch_helper->p_recv_data    , n_new_request, void                 **);

  PDM_realloc(exch_helper->d_send_stride  , exch_helper->d_send_stride  , n_new_request, int                   *);
  PDM_realloc(exch_helper->d_send_data    , exch_helper->d_send_data    , n_new_request, void                  *);
  PDM_realloc(exch_helper->d_recv_stride  , exch_helper->d_recv_stride  , n_new_request, int                   *);
  PDM_realloc(exch_helper->d_recv_data    , exch_helper->d_recv_data    , n_new_request, void                  *);

  for(int i = exch_helper->n_request; i < n_new_request; ++i) {
    exch_helper->requests_status[i] = EXCHANGE_HELPER_STATUS_FREE;
    exch_helper->n_sub_requests [i] = 0;
    exch_helper->sub_requests   [i] = NULL;

    exch_helper->send_buffer    [i] = NULL;
    exch_helper->recv_buffer    [i] = NULL;

    exch_helper->s_data         [i] = 0;
    exch_helper->cst_stride     [i] = 0;
    exch_helper->p_send_stride  [i] = NULL;
    exch_helper->p_send_data    [i] = NULL;
    exch_helper->p_recv_stride  [i] = NULL;
    exch_helper->p_recv_data    [i] = NULL;
    exch_helper->d_send_stride  [i] = NULL;
    exch_helper->d_send_data    [i] = NULL;
    exch_helper->d_recv_stride  [i] = NULL;
    exch_helper->d_recv_data    [i] = NULL;
  }

  exch_helper->n_request = n_new_request;
  return exch_helper->n_request-1;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_exchange_helper_t *
PDM_exchange_helper_create
(
  const PDM_MPI_Comm    comm,
        int             n_request_init
)
{
  PDM_exchange_helper_t *exch_helper = NULL;
  PDM_malloc(exch_helper, 1, PDM_exchange_helper_t);

  exch_helper->comm      = comm;
  exch_helper->n_request = n_request_init;

  PDM_malloc(exch_helper->requests_status, exch_helper->n_request, _exch_helper_status_t  );
  PDM_malloc(exch_helper->sub_requests   , exch_helper->n_request, PDM_MPI_Request       *);
  PDM_malloc(exch_helper->n_sub_requests , exch_helper->n_request, int                    );
  PDM_malloc(exch_helper->send_buffer    , exch_helper->n_request, void                  *);
  PDM_malloc(exch_helper->recv_buffer    , exch_helper->n_request, void                  *);

  PDM_malloc(exch_helper->t_stride       , exch_helper->n_request, PDM_stride_t           );
  PDM_malloc(exch_helper->s_data         , exch_helper->n_request, size_t                 );
  PDM_malloc(exch_helper->cst_stride     , exch_helper->n_request, int                    );
  PDM_malloc(exch_helper->p_send_stride  , exch_helper->n_request, int                  **);
  PDM_malloc(exch_helper->p_send_data    , exch_helper->n_request, void                 **);
  PDM_malloc(exch_helper->p_recv_stride  , exch_helper->n_request, int                  **);
  PDM_malloc(exch_helper->p_recv_data    , exch_helper->n_request, void                 **);

  PDM_malloc(exch_helper->d_send_stride  , exch_helper->n_request, int                   *);
  PDM_malloc(exch_helper->d_send_data    , exch_helper->n_request, void                  *);
  PDM_malloc(exch_helper->d_recv_stride  , exch_helper->n_request, int                   *);
  PDM_malloc(exch_helper->d_recv_data    , exch_helper->n_request, void                  *);

  for(int i = 0; i < exch_helper->n_request; ++i) {
    exch_helper->requests_status[i] = EXCHANGE_HELPER_STATUS_FREE;
    exch_helper->n_sub_requests [i] = 0;
    exch_helper->sub_requests   [i] = NULL;
    exch_helper->send_buffer    [i] = NULL;
    exch_helper->recv_buffer    [i] = NULL;

    exch_helper->s_data         [i] = 0;
    exch_helper->cst_stride     [i] = 0;
    exch_helper->p_send_stride  [i] = NULL;
    exch_helper->p_send_data    [i] = NULL;
    exch_helper->p_recv_stride  [i] = NULL;
    exch_helper->p_recv_data    [i] = NULL;
    exch_helper->d_send_stride  [i] = NULL;
    exch_helper->d_send_data    [i] = NULL;
    exch_helper->d_recv_stride  [i] = NULL;
    exch_helper->d_recv_data    [i] = NULL;
  }

  return exch_helper;
}

void
PDM_exchange_helper_exch
(
  PDM_exchange_helper_t *exch_helper,
  PDM_mpi_comm_kind_t    kcomm,
  int                    cst_stride,
  size_t                 s_data,
  int                   *send_idx,
  int                   *send_n,
  void                  *send_buffer,
  int                   *recv_idx,
  int                   *recv_n,
  void                  *recv_buffer
)
{
  int s_data_tot = s_data * cst_stride;

  PDM_MPI_Datatype mpi_type;
  PDM_MPI_Type_create_contiguous(s_data_tot, PDM_MPI_BYTE, &mpi_type);
  PDM_MPI_Type_commit(&mpi_type);

  if(kcomm == PDM_MPI_COMM_KIND_COLLECTIVE) {
    PDM_MPI_Alltoallv(send_buffer,
                      send_n,
                      send_idx,
                      mpi_type,
                      recv_buffer,
                      recv_n,
                      recv_idx,
                      mpi_type,
                      exch_helper->comm);
  } else if(kcomm == PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE) {
    PDM_MPI_Neighbor_alltoallv(send_buffer,
                               send_n,
                               send_idx,
                               mpi_type,
                               recv_buffer,
                               recv_n,
                               recv_idx,
                               mpi_type,
                               exch_helper->comm);
  } else if(kcomm == PDM_MPI_COMM_KIND_P2P) {
    PDM_MPI_Alltoallv_p2p(send_buffer,
                          send_n,
                          send_idx,
                          mpi_type,
                          recv_buffer,
                          recv_n,
                          recv_idx,
                          mpi_type,
                          exch_helper->comm);
  } else {
    PDM_error(__FILE__, __LINE__, 0,
              "Error PDM_exchange_helper_exch not yet implemented with kcomm = %i\n", kcomm);
  }

  PDM_MPI_Type_free(&mpi_type);
}

int
PDM_exchange_helper_exch_init
(
  PDM_exchange_helper_t *exch_helper,
  PDM_mpi_comm_kind_t    kcomm,
  size_t                 s_data,
  int                    cst_stride,
  int                   *send_idx,
  int                   *send_n,
  void                  *send_buffer,
  int                   *recv_idx,
  int                   *recv_n,
  void                  *recv_buffer
)
{
  // Ownership ou pas ??
  // Vu qu'on gère des request deja, et qu'on souhaite allegé les structures internes, je dirai que oui,
  // Le ptp construit les buffer, et on les gardes en interne ici plutot que dans le ptp ?

  int request_id = _find_available_request(exch_helper);

  int s_data_tot = s_data * cst_stride;

  PDM_MPI_Datatype mpi_type;
  PDM_MPI_Type_create_contiguous(s_data_tot, PDM_MPI_BYTE, &mpi_type);
  PDM_MPI_Type_commit(&mpi_type);

  if(kcomm == PDM_MPI_COMM_KIND_COLLECTIVE) {
    exch_helper->n_sub_requests[request_id] = 1;
    PDM_malloc(exch_helper->sub_requests[request_id], exch_helper->n_sub_requests[request_id], PDM_MPI_Request);
    PDM_MPI_Alltoallv_init(send_buffer,
                           send_n,
                           send_idx,
                           mpi_type,
                           recv_buffer,
                           recv_n,
                           recv_idx,
                           mpi_type,
                           exch_helper->comm,
                           &exch_helper->sub_requests[request_id][0]);
  } else if(kcomm == PDM_MPI_COMM_KIND_P2P) {
    PDM_MPI_Alltoallv_p2p_init(send_buffer,
                               send_n,
                               send_idx,
                               mpi_type,
                               recv_buffer,
                               recv_n,
                               recv_idx,
                               mpi_type,
                               exch_helper->comm,
                               &exch_helper->n_sub_requests[request_id],
                               &exch_helper->sub_requests  [request_id]);
  } else {
    PDM_error(__FILE__, __LINE__, 0,
              "Error PDM_exchange_helper_exch not yet implemented with kcomm = %i\n", kcomm);
  }
  PDM_MPI_Type_free(&mpi_type);

  exch_helper->requests_status[request_id] = EXCHANGE_HELPER_STATUS_READY;

  return request_id;
}



void
PDM_exchange_helper_exch_start
(
  PDM_exchange_helper_t *exch_helper,
  int                    request_id
)
{
  if(exch_helper->requests_status[request_id] != EXCHANGE_HELPER_STATUS_READY) {
    PDM_error(__FILE__, __LINE__, 0,
              "Error PDM_exchange_helper_exch_start with status = %i for request_id = %i, you should initialize exch with PDM_exchange_helper_exch_init or PDM_exchange_helper_iexch\n", exch_helper->requests_status[request_id], request_id);
  }

  exch_helper->requests_status[request_id] = EXCHANGE_HELPER_STATUS_ONGOING;
  PDM_MPI_Startall(exch_helper->n_sub_requests[request_id],
                   exch_helper->sub_requests  [request_id]);
}


void
PDM_exchange_helper_exch_wait
(
  PDM_exchange_helper_t *exch_helper,
  int                    request_id
)
{
  if(exch_helper->requests_status[request_id] != EXCHANGE_HELPER_STATUS_ONGOING) {
    PDM_error(__FILE__, __LINE__, 0,
              "Error PDM_exchange_helper_exch_wait with status = %i for request_id = %i, you should initialize exch with PDM_exchange_helper_exch_init or PDM_exchange_helper_iexch\n", exch_helper->requests_status[request_id], request_id);
  }

  for(int i = 0; i < exch_helper->n_sub_requests[request_id]; ++i) {
    PDM_MPI_Wait(&exch_helper->sub_requests[request_id][i]);
  }
  exch_helper->requests_status[request_id] = EXCHANGE_HELPER_STATUS_READY;
}

void
PDM_exchange_helper_exch_free
(
  PDM_exchange_helper_t *exch_helper,
  int                    request_id
)
{
  // Ca le remet a null automatiquemnt
  for(int i = 0; i < exch_helper->n_sub_requests[request_id]; ++i) {
    PDM_MPI_Request_free(&exch_helper->sub_requests[request_id][i]);
  }
  PDM_free(exch_helper->sub_requests[request_id]);
  exch_helper->requests_status[request_id] = EXCHANGE_HELPER_STATUS_FREE;
  exch_helper->n_sub_requests [request_id] = 0;
}


int
PDM_exchange_helper_exch_one_way_init
(
  PDM_exchange_helper_t  *exch_helper,
  PDM_mpi_comm_kind_t     kcomm,
  int                     cst_stride,
  size_t                  s_data,
  int                     tag,
  void                  **buffer,
  PDM_ownership_t         ownership
)
{
  PDM_UNUSED(exch_helper);
  PDM_UNUSED(kcomm);
  PDM_UNUSED(cst_stride);
  PDM_UNUSED(s_data);
  PDM_UNUSED(tag);
  PDM_UNUSED(buffer);
  PDM_UNUSED(ownership);

  return -1;
}



void
PDM_exchange_helper_free
(
  PDM_exchange_helper_t *exch_helper
)
{
  for(int i_req = 0; i_req < exch_helper->n_request; ++i_req) {
    for(int i = 0; i < exch_helper->n_sub_requests[i_req]; ++i) {
      PDM_MPI_Request_free(&exch_helper->sub_requests[i_req][i]);
    }
  }
  PDM_free(exch_helper->requests_status);
  PDM_free(exch_helper->n_sub_requests );
  PDM_free(exch_helper->sub_requests   );
  PDM_free(exch_helper->send_buffer    );
  PDM_free(exch_helper->recv_buffer    );

  PDM_free(exch_helper->t_stride       );
  PDM_free(exch_helper->s_data         );
  PDM_free(exch_helper->cst_stride     );
  PDM_free(exch_helper->p_send_stride  );
  PDM_free(exch_helper->p_send_data    );
  PDM_free(exch_helper->p_recv_stride  );
  PDM_free(exch_helper->p_recv_data    );
  PDM_free(exch_helper->d_send_stride  );
  PDM_free(exch_helper->d_send_data    );
  PDM_free(exch_helper->d_recv_stride  );
  PDM_free(exch_helper->d_recv_data    );


  PDM_free(exch_helper);
}


#ifdef __cplusplus
}
#endif
