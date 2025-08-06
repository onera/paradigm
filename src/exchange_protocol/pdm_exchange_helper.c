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
_get_next_tag
(
  PDM_exchange_helper_t *exch_helper,
  PDM_mpi_comm_kind_t    k_comm
)
{
  int tag = -10000;
  if (k_comm == PDM_MPI_COMM_KIND_P2P) {
    tag  = exch_helper->seed_tag;
    tag += (exch_helper->next_tag++);
    tag %= exch_helper->max_tag;
  }
  return tag;
}

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

  PDM_realloc(exch_helper->win_send       , exch_helper->win_send       , n_new_request, PDM_MPI_Win            );
  PDM_realloc(exch_helper->win_recv       , exch_helper->win_recv       , n_new_request, PDM_MPI_Win            );
  PDM_realloc(exch_helper->group_send     , exch_helper->group_send     , n_new_request, PDM_MPI_Group          );
  PDM_realloc(exch_helper->group_recv     , exch_helper->group_recv     , n_new_request, PDM_MPI_Group          );
  PDM_realloc(exch_helper->target_disp    , exch_helper->target_disp    , n_new_request, int                   *);

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

    exch_helper->win_send       [i] = PDM_MPI_WIN_NULL;
    exch_helper->win_recv       [i] = PDM_MPI_WIN_NULL;
    exch_helper->group_send     [i] = PDM_MPI_GROUP_NULL;
    exch_helper->group_recv     [i] = PDM_MPI_GROUP_NULL;
    exch_helper->target_disp    [i] = NULL;

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

  PDM_MPI_Comm_dup(comm, &exch_helper->comm);
  exch_helper->n_request = n_request_init;

  PDM_malloc(exch_helper->requests_status, exch_helper->n_request, _exch_helper_status_t  );
  PDM_malloc(exch_helper->sub_requests   , exch_helper->n_request, PDM_MPI_Request       *);
  PDM_malloc(exch_helper->n_sub_requests , exch_helper->n_request, int                    );
  PDM_malloc(exch_helper->send_buffer    , exch_helper->n_request, void                  *);
  PDM_malloc(exch_helper->recv_buffer    , exch_helper->n_request, void                  *);

  // RMA
  PDM_malloc(exch_helper->win_send       , exch_helper->n_request, PDM_MPI_Win            );
  PDM_malloc(exch_helper->win_recv       , exch_helper->n_request, PDM_MPI_Win            );
  PDM_malloc(exch_helper->group_send     , exch_helper->n_request, PDM_MPI_Group          );
  PDM_malloc(exch_helper->group_recv     , exch_helper->n_request, PDM_MPI_Group          );
  PDM_malloc(exch_helper->target_disp    , exch_helper->n_request, int                   *);

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
    exch_helper->win_send       [i] = PDM_MPI_WIN_NULL;
    exch_helper->win_recv       [i] = PDM_MPI_WIN_NULL;
    exch_helper->group_send     [i] = PDM_MPI_GROUP_NULL;
    exch_helper->group_recv     [i] = PDM_MPI_GROUP_NULL;
    exch_helper->target_disp    [i] = NULL;

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

  // Create tag for P2P
  void  *max_tag_tmp;
  int    flag = 0;

  // Mandatory to call with PDM_MPI_COMM_WORLD becuase only this one keep attributes (openMPI implemntation for exemple)
  PDM_MPI_Comm_get_attr_tag_ub(PDM_MPI_COMM_WORLD, &max_tag_tmp, &flag);
  exch_helper->max_tag  = (long) (*((int *) max_tag_tmp));
  exch_helper->seed_tag = 1;
  exch_helper->next_tag = 1;

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

  int topo_kind = PDM_MPI_UNDEFINED;
  PDM_MPI_Topo_test(exch_helper->comm, &topo_kind);

  if(kcomm == PDM_MPI_COMM_KIND_COLLECTIVE && topo_kind == PDM_MPI_UNDEFINED) {
    PDM_MPI_Alltoallv(send_buffer,
                      send_n,
                      send_idx,
                      mpi_type,
                      recv_buffer,
                      recv_n,
                      recv_idx,
                      mpi_type,
                      exch_helper->comm);
  } else if(kcomm == PDM_MPI_COMM_KIND_COLLECTIVE && topo_kind == PDM_MPI_DIST_GRAPH) {
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
)
{
  int tag        = _get_next_tag          (exch_helper, k_comm);
  int request_id = _find_available_request(exch_helper);

  int s_data_tot = s_data * cst_stride;

  int n_rank;
  PDM_MPI_Comm_size(exch_helper->comm, &n_rank);

  PDM_MPI_Datatype mpi_type;
  PDM_MPI_Type_create_contiguous(s_data_tot, PDM_MPI_BYTE, &mpi_type);
  PDM_MPI_Type_commit(&mpi_type);

  if(k_comm == PDM_MPI_COMM_KIND_COLLECTIVE) {
    exch_helper->n_sub_requests[request_id] = 1;
    PDM_malloc(exch_helper->sub_requests[request_id], exch_helper->n_sub_requests[request_id], PDM_MPI_Request);
    PDM_MPI_Ialltoallv(send_buffer,
                       send_n,
                       send_idx,
                       mpi_type,
                       recv_buffer,
                       recv_n,
                       recv_idx,
                       mpi_type,
                       exch_helper->comm,
                       &exch_helper->sub_requests[request_id][0]);
  } else if (k_comm == PDM_MPI_COMM_KIND_P2P) {
    PDM_MPI_Ialltoallv_p2p(send_buffer,
                           send_n,
                           send_idx,
                           mpi_type,
                           0,    // n_send_rank
                           NULL, // send_rank
                           recv_buffer,
                           recv_n,
                           recv_idx,
                           mpi_type,
                           0,    // n_recv_rank
                           NULL, // recv_rank
                           tag,
                           exch_helper->comm,
                           &exch_helper->n_sub_requests[request_id],
                           &exch_helper->sub_requests  [request_id]);
  } else if(k_comm == PDM_MPI_COMM_KIND_WIN_RMA) {

    /*
     * Creates a memory window (`win_send`) on `send_buffer`.
     */
    PDM_MPI_Win_create(send_buffer, s_data_tot * send_idx[n_rank], s_data_tot, exch_helper->comm, &exch_helper->win_send[request_id]);

    /*
     * Initialization of MPI synchronization groups for the Active Target model.
     * This step is necessary to identify communication partners.
     * It is typically done once at the beginning to be reused later.
     */
    PDM_MPI_Group world_group;
    PDM_MPI_Comm_group(exch_helper->comm, &world_group);

    /*
     * Creates a process group (`group_recv`).
     * This group contains all processes from whom the current process will receive data.
     */
    int *tmp_rank_id = NULL;
    PDM_malloc(tmp_rank_id, n_rank, int);

    int n_recv = 0;
    for(int i = 0; i < n_rank; ++i) {
      if(recv_n[i] > 0) {
        tmp_rank_id[n_recv++] = i;
      }
    }

    PDM_MPI_Group_incl(world_group, n_recv, tmp_rank_id, &exch_helper->group_recv[request_id]);
    PDM_MPI_Group_free(&world_group);

    PDM_malloc(exch_helper->target_disp[request_id], n_rank, int);

    /*
     * Uses an MPI_Alltoall collective communication to exchange offsets.
     * Each process sends its own offset (`send_idx`) and receives the offsets from all others.
     * These offsets (stored in `exch_helper->target_disp`) are crucial
     * for the current process to know where to read/write in the remote windows.
     */
    PDM_MPI_Alltoall(send_idx                            , 1, PDM_MPI_INT,
                     exch_helper->target_disp[request_id], 1, PDM_MPI_INT, exch_helper->comm);

    if(0 == 1) {
      PDM_log_trace_array_int(exch_helper->target_disp[request_id], n_rank, "exch_helper->target_disp[request_id]");
      PDM_log_trace_array_int(send_idx, n_rank+1, "send_idx ::");
      PDM_log_trace_array_int(recv_idx, n_rank+1, "recv_idx ::");
      PDM_log_trace_array_int(tmp_rank_id, n_recv, "tmp_rank_id ::");
    }
    PDM_free(tmp_rank_id);

    /*
     * This is the beginning of an MPI Active Target RMA communication epoch.
     * These two calls, `Win_post` and `Win_start`, are typically called back-to-back,
     * but they have very different roles. They can be thought of as a handshake
     * between the partners defined by the groups.
     *
     * 1. Declares that this process's window (`win_send`) is available.
     *    The current process "posts" its window, indicating that the processes in `group_recv`
     *    are now allowed to access its memory. This enables the group members to start
     *    their RMA operations (Rget/Rput) to/from this window.
     *
     */
    PDM_MPI_Win_post (exch_helper->group_recv[request_id], 0, exch_helper->win_send[request_id]);

    /*
     * 2. Starts an access epoch to the partners' windows.
     *    The current process "starts" an access epoch, declaring that it will initiate
     *    RMA operations to the windows of the processes in `group_recv`.
     *    This call is a prerequisite before launching any `MPI_Rget` or `MPI_Rput` functions.
     */
    PDM_MPI_Win_start(exch_helper->group_recv[request_id], 0, exch_helper->win_send[request_id]);

    /*
     * Once the `post` and `start` epochs are established, data transfers can be launched.
     *
     * Launches the asynchronous Ialltoallv_p2p_rma operations.
     * These operations will use the local `win_send` window, the remote offsets
     * (`target_disp`), and the local `recv_buffer` to read data from partners
     * (with MPI_Rget) or write data to them (with MPI_Rput).
     * These calls are non-blocking and generate requests that must be completed later
     * with a `PDM_MPI_Waitall` (or equivalent).
     */
    PDM_MPI_Ialltoallv_p2p_rma(exch_helper->win_send   [request_id],
                               exch_helper->target_disp[request_id],
                               recv_buffer,
                               recv_n,
                               recv_idx,
                               mpi_type,
                               exch_helper->comm,
                               &exch_helper->n_sub_requests[request_id],
                               &exch_helper->sub_requests  [request_id]);
  } else {
    PDM_error(__FILE__, __LINE__, 0,
              "Error PDM_exchange_helper_iexch not yet implemented with kcomm = %i\n", k_comm);
  }

  PDM_MPI_Type_free(&mpi_type);

  exch_helper->requests_status[request_id] = EXCHANGE_HELPER_STATUS_ONGOING;

  return request_id;
}

// A réfléchir pour appeler plus facilement les IAlltoall
// int
// PDM_exchange_helper_exch_init
// (
//   PDM_exchange_helper_t *exch_helper,
//   PDM_mpi_comm_kind_t    k_comm,
//   size_t                 s_data,
//   int                    cst_stride,
//   int                    n_send_rank,
//   int                   *send_rank,
//   int                   *send_idx,
//   int                   *send_n,
//   void                  *send_buffer,
//   int                    n_recv_rank,
//   int                   *recv_rank,
//   int                   *recv_idx,
//   int                   *recv_n,
//   void                  *recv_buffer
// );
// Si le send_rank ou le recv_rank est NULL --> Alltoall classique
// Sinon on peut faire du ISend/Irecv plus fin

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
)
{
  int tag        = _get_next_tag          (exch_helper, k_comm);
  int request_id = _find_available_request(exch_helper);

  int s_data_tot = s_data * cst_stride;

  PDM_MPI_Datatype mpi_type;
  PDM_MPI_Type_create_contiguous(s_data_tot, PDM_MPI_BYTE, &mpi_type);
  PDM_MPI_Type_commit(&mpi_type);

  if(k_comm == PDM_MPI_COMM_KIND_COLLECTIVE) {
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
  } else if(k_comm == PDM_MPI_COMM_KIND_P2P) {
    PDM_MPI_Alltoallv_p2p_init(send_buffer,
                               send_n,
                               send_idx,
                               mpi_type,
                               recv_buffer,
                               recv_n,
                               recv_idx,
                               mpi_type,
                               tag,
                               exch_helper->comm,
                               &exch_helper->n_sub_requests[request_id],
                               &exch_helper->sub_requests  [request_id]);
  } else {
    PDM_error(__FILE__, __LINE__, 0,
              "Error PDM_exchange_helper_exch not yet implemented with kcomm = %i\n", k_comm);
  }
  PDM_MPI_Type_free(&mpi_type);

  exch_helper->requests_status[request_id] = EXCHANGE_HELPER_STATUS_READY;

  return request_id;
}

int
PDM_mpi_comm_kind_is_persistent
(
  PDM_mpi_comm_kind_t    k_comm
)
{
  return 0;
}

// Si comm graph -> Neighbor automatique ? Voire dans le part_migrate
// RMA --> init, start, free (Permet d'économiser la création des groupes et du target_disp)
// Etat des lieux des k_comm dans PDM --> Pour moi pas cohérent
// Echange ONE-Way à faire ( send / recv )
// Gestion des actives rank pour le p2p
// Test exch_helper


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
  printf("request_id = %i - n_sub_requests = %i \n", request_id, exch_helper->n_sub_requests[request_id]);

  for(int i = 0; i < exch_helper->n_sub_requests[request_id]; ++i) {
    PDM_MPI_Wait(&exch_helper->sub_requests[request_id][i]);
  }

  /*
   * Completes the window's access epoch.
   * This call signals to the target processes (those with whom this process communicated)
   * that this process has finished all its RMA operations (Rget/Rput) to their windows.
   * This allows target processes, which are waiting with PDM_MPI_Win_wait, to proceed.
   */
  PDM_MPI_Win_complete(exch_helper->win_send[request_id]);

  /* Completes the window's exposure epoch.
   * This call waits for all source processes (those who accessed this process's memory)
   * to have finished their work and signaled their completion via their PDM_MPI_Win_complete.
   * It ensures that all incoming data "pushed" into this process's memory is
   * now visible and ready to be used by the local process.
   */
  PDM_MPI_Win_wait(exch_helper->win_send[request_id]);

  // RM Group and target_dsip


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


  PDM_MPI_Win_free  (&exch_helper->win_send  [request_id]);
  PDM_MPI_Win_free  (&exch_helper->win_recv  [request_id]);
  PDM_MPI_Group_free(&exch_helper->group_send[request_id]);
  PDM_MPI_Group_free(&exch_helper->group_recv[request_id]);
  PDM_free(exch_helper->target_disp[request_id]);


}


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
)
{
  PDM_UNUSED(exch_helper);
  PDM_UNUSED(direction);
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

  PDM_free(exch_helper->win_send       );
  PDM_free(exch_helper->win_recv       );
  PDM_free(exch_helper->group_send     );
  PDM_free(exch_helper->group_recv     );
  PDM_free(exch_helper->target_disp    );

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

  PDM_MPI_Comm_free(&exch_helper->comm);

  PDM_free(exch_helper);
}


#ifdef __cplusplus
}
#endif
