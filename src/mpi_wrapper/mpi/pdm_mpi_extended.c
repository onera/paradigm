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
#include "pdm_mpi.h"
#include "pdm_mpi_priv.h"
#include "pdm_mem_tool.h"

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

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

int
PDM_MPI_Sends_init
(
  const void              *sendbuf,
        int               *sendcounts,
        int               *sdispls,
        PDM_MPI_Datatype   datatype,
        int                n_active_send,
        int               *active_send,
        int                tag,
        PDM_MPI_Comm       comm,
        PDM_MPI_Request  **out_requests
)
{
  PDM_MPI_Request *requests = NULL;
  PDM_malloc(requests, n_active_send, PDM_MPI_Request);

  int size_send_type;
  MPI_Type_size(datatype, &size_send_type);

  int code = MPI_SUCCESS;
  for (int i = 0; i < n_active_send; i++) {
    void *buf = (void *) ((unsigned char*) sendbuf + sdispls[i] * size_send_type);
    int t_rank = active_send[i];
    code = MPI_Send_init(buf,
                         sendcounts[i],
                         datatype,
                         t_rank,
                         tag,
                         comm,
                         &requests[i]);
    if (code != MPI_SUCCESS) {
      break;
    }
  }
  *out_requests = requests;

  return code;
}



int
PDM_MPI_Isends
(
  const void              *sendbuf,
        int               *sendcounts,
        int               *sdispls,
        PDM_MPI_Datatype   datatype,
        int                n_active_send,
        int               *active_send,
        int                tag,
        PDM_MPI_Comm       comm,
        PDM_MPI_Request  **out_requests
)
{
  PDM_MPI_Request *requests = NULL;
  PDM_malloc(requests, n_active_send, PDM_MPI_Request);

  int size_send_type;
  MPI_Type_size(datatype, &size_send_type);

  int code = MPI_SUCCESS;
  for (int i = 0; i < n_active_send; i++) {
    void *buf = (void *) ((unsigned char*) sendbuf + sdispls[i] * size_send_type);
    int t_rank = active_send[i];
    code = MPI_Isend(buf,
                     sendcounts[i],
                     datatype,
                     t_rank,
                     tag,
                     comm,
                     &requests[i]);
    if (code != MPI_SUCCESS) {
      break;
    }
  }
  *out_requests = requests;

  return code;
}


int
PDM_MPI_Recvs_init
(
        void              *recvbuf,
        int               *recvcounts,
        int               *rdispls,
        PDM_MPI_Datatype   datatype,
        int                n_active_recv,
        int               *active_recv,
        int                tag,
        PDM_MPI_Comm       comm,
        PDM_MPI_Request  **out_requests
)
{
  PDM_MPI_Request *requests = NULL;
  PDM_malloc(requests, n_active_recv, PDM_MPI_Request);

  int size_recv_type;
  MPI_Type_size(datatype, &size_recv_type);

  int code = MPI_SUCCESS;
  for (int i = 0; i < n_active_recv; i++) {
    void *buf = (void *) ((unsigned char*) recvbuf + rdispls[i] * size_recv_type);
    int t_rank = active_recv[i];
    code = MPI_Recv_init(buf,
                         recvcounts[i],
                         datatype,
                         t_rank,
                         tag,
                         comm,
                         &requests[i]);
    if (code != MPI_SUCCESS) {
      break;
    }
  }
  *out_requests = requests;
  return code;
}


int
PDM_MPI_Irecvs
(
  const void              *recvbuf,
        int               *recvcounts,
        int               *sdispls,
        PDM_MPI_Datatype   datatype,
        int                n_active_recv,
        int               *active_recv,
        int                tag,
        PDM_MPI_Comm       comm,
        PDM_MPI_Request  **out_requests
)
{
  PDM_MPI_Request *requests = NULL;
  PDM_malloc(requests, n_active_recv, PDM_MPI_Request);

  int size_recv_type;
  MPI_Type_size(datatype, &size_recv_type);

  int code = MPI_SUCCESS;
  for (int i = 0; i < n_active_recv; i++) {
    void *buf = (void *) ((unsigned char*) recvbuf + sdispls[i] * size_recv_type);
    int t_rank = active_recv[i];
    code = MPI_Irecv(buf,
                     recvcounts[i],
                     datatype,
                     t_rank,
                     tag,
                     comm,
                     &requests[i]);
    if (code != MPI_SUCCESS) {
      break;
    }
  }
  *out_requests = requests;

  return code;
}

void
PDM_MPI_Partofactiverank
(
  int          *sendcounts,
  int          *recvcounts,
  PDM_MPI_Comm  comm,
  double       *part_active_rank
)
{
  int size;
  MPI_Comm_size(comm, &size);

  int rank;
  MPI_Comm_rank(comm, &rank);

  int n_active_rank = 0;
  for (int i = 0; i < size; i++) {
    if ((sendcounts[i] > 0) || (recvcounts[i] > 0)) {
      n_active_rank++;
    }
  }

  double _part_active_rank = (double) n_active_rank / (double) size;

  PDM_MPI_Allreduce (&_part_active_rank,
                     part_active_rank,
                     1,
                     PDM_MPI_DOUBLE,
                     PDM_MPI_MAX,
                     comm);

}

int
PDM_MPI_Ialltoallv_select_p2p
(
  void              *sendbuf,
  int               *sendcounts,
  int               *sdispls,
  PDM_MPI_Datatype   sendtype,
  int                n_send_rank,
  int               *send_rank,
  void              *recvbuf,
  int               *recvcounts,
  int               *rdispls,
  PDM_MPI_Datatype   recvtype,
  int                n_recv_rank,
  int               *recv_rank,
  int                tag,
  PDM_MPI_Comm       comm,
  int               *n_send_recv_request,
  PDM_MPI_Request  **out_requests
)
{

  PDM_MPI_Request *requests = NULL;
  int n_request = n_recv_rank + n_send_rank;
  PDM_malloc(requests, n_request, PDM_MPI_Request);

  int size_send_type;
  int size_recv_type;
  MPI_Type_size(sendtype, &size_send_type);
  MPI_Type_size(recvtype, &size_recv_type);

  int code = MPI_SUCCESS;
  n_request = 0;
  for (int i = 0; i < n_recv_rank; i++) {
    void *buf = (void *) ((unsigned char*) recvbuf + rdispls[i] * size_recv_type);
    int t_rank = recv_rank[i];
    code = MPI_Irecv(buf,
                     recvcounts[i],
                     recvtype,
                     t_rank,
                     tag,
                     comm,
                     &requests[n_request]);
    n_request++;
    if (code != MPI_SUCCESS) {
      break;
    }
  }

  for (int i = 0; i < n_send_rank; i++) {
    void *buf = (void *) ((unsigned char*) sendbuf + sdispls[i] * size_send_type);
    int t_rank = send_rank[i];
    code = MPI_Isend(buf,
                     sendcounts[i],
                     sendtype,
                     t_rank,
                     tag,
                     comm,
                     &requests[n_request]);
    n_request++;
    if (code != MPI_SUCCESS) {
      break;
    }
  }

  if (code != MPI_SUCCESS) {
    return code;
  }

  *n_send_recv_request = n_request;
  *out_requests        = requests;

  return code;
}

int
PDM_MPI_Ialltoallv_p2p
(
  void              *sendbuf,
  int               *sendcounts,
  int               *sdispls,
  PDM_MPI_Datatype   sendtype,
  int                n_send_rank,
  int               *send_rank,
  void              *recvbuf,
  int               *recvcounts,
  int               *rdispls,
  PDM_MPI_Datatype   recvtype,
  int                n_recv_rank,
  int               *recv_rank,
  int                tag,
  PDM_MPI_Comm       comm,
  int               *n_send_recv_request,
  PDM_MPI_Request  **out_requests
)
{

  int code = MPI_SUCCESS;
  int n_rank;
  MPI_Comm_size(comm, &n_rank);

  // Short-cut if send_rank and recv_rank is specified
  if(send_rank != NULL && recv_rank != NULL) {
    return PDM_MPI_Ialltoallv_select_p2p(sendbuf,
                                         sendcounts,
                                         sdispls,
                                         sendtype,
                                         n_send_rank,
                                         send_rank,
                                         recvbuf,
                                         recvcounts,
                                         rdispls,
                                         recvtype,
                                         n_recv_rank,
                                         recv_rank,
                                         tag,
                                         comm,
                                         n_send_recv_request,
                                         out_requests);
  }

  // Count number of request
  int n_request = 0;
  for (int i = 0; i < n_rank; i++) {
    if (recvcounts[i] != 0) {
      n_request++;
    }
    if (sendcounts[i] != 0) {
      n_request++;
    }
  }

  PDM_MPI_Request *requests = NULL;
  PDM_malloc(requests, n_request, PDM_MPI_Request);

  int size_send_type;
  int size_recv_type;
  MPI_Type_size(sendtype, &size_send_type);
  MPI_Type_size(recvtype, &size_recv_type);

  n_request = 0;
  for (int i = 0; i < n_rank; i++) {
    if (recvcounts[i] != 0) {
      void *buf = (void *) ((unsigned char*) recvbuf + rdispls[i] * size_recv_type);
      code = MPI_Irecv(buf,
                       recvcounts[i],
                       recvtype,
                       i,
                       tag,
                       comm,
                       &requests[n_request]);
      n_request++;
      if (code != MPI_SUCCESS) {
        break;
      }
    }

    if (sendcounts[i] != 0) {
      void *buf = (void *) ((unsigned char*) sendbuf + sdispls[i] * size_send_type);
      code = MPI_Isend(buf,
                       sendcounts[i],
                       sendtype,
                       i,
                       tag,
                       comm,
                       &requests[n_request]);
      n_request++;
      if (code != MPI_SUCCESS) {
        break;
      }
    }
  }

  if (code != MPI_SUCCESS) {
    return code;
  }

  *n_send_recv_request = n_request;
  *out_requests        = requests;

  return code;
}



int
PDM_MPI_Ialltoallv_p2p_rma
(
  PDM_MPI_Win        send_win,
  int               *target_disp,
  void              *recvbuf,
  int               *recvcounts,
  int               *rdispls,
  PDM_MPI_Datatype   recvtype,
  PDM_MPI_Comm       comm,
  int               *n_send_recv_request,
  PDM_MPI_Request  **out_requests
)
{
  int code = MPI_SUCCESS;

  int n_rank;
  MPI_Comm_size(comm, &n_rank);

  // Count number of request
  int n_request = 0;
  for (int i = 0; i < n_rank; i++) {
    if (recvcounts[i] != 0) {
      n_request++;
    }
  }

  PDM_MPI_Request *requests = NULL;
  PDM_malloc(requests, n_request, PDM_MPI_Request);

  int size_recv_type;
  MPI_Type_size(recvtype, &size_recv_type);

  n_request = 0;
  for (int i = 0; i < n_rank; i++) {
    if (recvcounts[i] != 0) {
      void *buf = (void *) ((unsigned char*) recvbuf + rdispls[i] * size_recv_type);
      code = MPI_Rget(buf,
                      recvcounts[i],
                      recvtype,
                      i,
                      target_disp[i],
                      recvcounts[i],
                      recvtype,
                      send_win,
                      &requests[n_request]);

      n_request++;
      if (code != MPI_SUCCESS) {
        break;
      }
    }
  }

  if (code != MPI_SUCCESS) {
    return code;
  }

  *n_send_recv_request = n_request;
  *out_requests        = requests;

  return code;
}



int
PDM_MPI_Alltoallv_p2p_init
(
  void              *sendbuf,
  int               *sendcounts,
  int               *sdispls,
  PDM_MPI_Datatype   sendtype,
  void              *recvbuf,
  int               *recvcounts,
  int               *rdispls,
  PDM_MPI_Datatype   recvtype,
  int                tag,
  PDM_MPI_Comm       comm,
  int               *n_send_recv_request,
  PDM_MPI_Request  **out_requests
)
{
  int code = MPI_SUCCESS;

  int n_rank;
  MPI_Comm_size(comm, &n_rank);

  // Count number of request
  int n_request = 0;
  for (int i = 0; i < n_rank; i++) {
    if (recvcounts[i] != 0) {
      n_request++;
    }
    if (sendcounts[i] != 0) {
      n_request++;
    }
  }

  PDM_MPI_Request *requests = NULL;
  PDM_malloc(requests, n_request, PDM_MPI_Request);

  int size_send_type;
  int size_recv_type;
  MPI_Type_size(sendtype, &size_send_type);
  MPI_Type_size(recvtype, &size_recv_type);

  n_request = 0;
  for (int i = 0; i < n_rank; i++) {
    if (recvcounts[i] != 0) {
      void *buf = (void *) ((unsigned char*) recvbuf + rdispls[i] * size_recv_type);
      code = MPI_Recv_init(buf,
                           recvcounts[i],
                           recvtype,
                           i,
                           tag,
                           comm,
                           &requests[n_request]);
      n_request++;
      if (code != MPI_SUCCESS) {
        break;
      }
    }
  }

  for (int i = 0; i < n_rank; i++) {
    if (sendcounts[i] != 0) {
      void *buf = (void *) ((unsigned char*) sendbuf + sdispls[i] * size_send_type);
      code = MPI_Send_init(buf,
                           sendcounts[i],
                           sendtype,
                           i,
                           tag,
                           comm,
                           &requests[n_request]);
      n_request++;
      if (code != MPI_SUCCESS) {
        break;
      }
    }
  }

  if (code != MPI_SUCCESS) {
    return code;
  }

  *n_send_recv_request = n_request;
  *out_requests        = requests;

  return code;
}


int
PDM_MPI_Alltoallv_p2p
(
  void             *sendbuf,
  int              *sendcounts,
  int              *sdispls,
  PDM_MPI_Datatype  sendtype,
  void             *recvbuf,
  int              *recvcounts,
  int              *rdispls,
  PDM_MPI_Datatype  recvtype,
  PDM_MPI_Comm      comm
)
{
  int code = MPI_SUCCESS;

  int size;
  MPI_Comm_size(comm, &size);

  MPI_Request *request_r;
  MPI_Request *request_s;
  PDM_malloc(request_r, size, MPI_Request);
  PDM_malloc(request_s, size, MPI_Request);

  int n_request_r = 0;
  int n_request_s = 0;

  int size_send_type;
  MPI_Type_size(sendtype, &size_send_type);

  int size_recv_type;
  MPI_Type_size(recvtype, &size_recv_type);

  for (int i = 0; i < size; i++) {
    if (recvcounts[i] != 0) {
      void *buf = (void *) ((unsigned char*) recvbuf + rdispls[i] * size_recv_type);
      code = MPI_Irecv(buf,
                       recvcounts[i],
                       recvtype,
                       i,
                       0,
                       comm,
                       request_r + n_request_r);
      n_request_r++;
      if (code != MPI_SUCCESS) {
        break;
      }
    }

    if (sendcounts[i] != 0) {
      void *buf = (void *) ((unsigned char*) sendbuf + sdispls[i] * size_send_type);
      code = MPI_Isend(buf,
                       sendcounts[i],
                       sendtype,
                       i,
                       0,
                       comm,
                       request_s + n_request_s);
      n_request_s++;
      if (code != MPI_SUCCESS) {
        break;
      }
    }
  }

  if (code != MPI_SUCCESS) {
    return code;
  }
  for (int i = 0; i < n_request_r; i++) {
    code = MPI_Wait(request_r + i, MPI_STATUS_IGNORE);
    if (code != MPI_SUCCESS) {
      break;
    }
  }

  if (code != MPI_SUCCESS) {
    return code;
  }
  for (int i = 0; i < n_request_s; i++) {
    code = MPI_Wait(request_s + i, MPI_STATUS_IGNORE);
    if (code != MPI_SUCCESS) {
      break;
    }
  }

  PDM_free(request_r);
  PDM_free(request_s);

  return code;
}


int
PDM_MPI_Alltoallv_p2p_l
(
  void             *sendbuf,
  int              *sendcounts,
  size_t           *sdispls,
  PDM_MPI_Datatype  sendtype,
  void             *recvbuf,
  int              *recvcounts,
  size_t           *rdispls,
  PDM_MPI_Datatype  recvtype,
  PDM_MPI_Comm      comm
)
{
  int code = MPI_SUCCESS;

  int size;
  MPI_Comm_size(comm, &size);

  MPI_Request *request_r;
  MPI_Request *request_s;
  PDM_malloc(request_r, size, MPI_Request);
  PDM_malloc(request_s, size, MPI_Request);

  int n_request_r = 0;
  int n_request_s = 0;

  int size_send_type;
  MPI_Type_size(sendtype, &size_send_type);

  int size_recv_type;
  MPI_Type_size(recvtype, &size_recv_type);

  for (int i = 0; i < size; i++) {
    if (recvcounts[i] != 0) {
      void *buf = (void *) ((unsigned char*) recvbuf + rdispls[i] * size_recv_type);
      code = MPI_Irecv(buf,
                       recvcounts[i],
                       recvtype,
                       i,
                       0,
                       comm,
                       request_r + n_request_r);
      n_request_r++;
      if (code != MPI_SUCCESS) {
        break;
      }
    }

    if (sendcounts[i] != 0) {
      void *buf = (void *) ((unsigned char*) sendbuf + sdispls[i] * size_send_type);
      code = MPI_Isend(buf,
                       sendcounts[i],
                       sendtype,
                       i,
                       0,
                       comm,
                       request_s + n_request_s);
      n_request_s++;
      if (code != MPI_SUCCESS) {
        break;
      }
    }
  }

  if (code != MPI_SUCCESS) {
    return code;
  }
  for (int i = 0; i < n_request_r; i++) {
    code = MPI_Wait(request_r + i, MPI_STATUS_IGNORE);
    if (code != MPI_SUCCESS) {
      break;
    }
  }

  if (code != MPI_SUCCESS) {
    return code;
  }
  for (int i = 0; i < n_request_s; i++) {
    code = MPI_Wait(request_s + i, MPI_STATUS_IGNORE);
    if (code != MPI_SUCCESS) {
      break;
    }
  }

  PDM_free(request_r);
  PDM_free(request_s);

  return code;
}

int
PDM_MPI_Alltoallv_l
(
  void             *sendbuf,
  int              *sendcounts,
  size_t           *sdispls,
  PDM_MPI_Datatype  sendtype,
  void             *recvbuf,
  int              *recvcounts,
  size_t           *rdispls,
  PDM_MPI_Datatype  recvtype,
  PDM_MPI_Comm      comm
)
{
  int code = MPI_SUCCESS;

  int size;
  MPI_Comm_size(comm, &size);

  INT_MAX;
  int coeff = 4;
  int large = 0;
  if ((sdispls[size-1] > (size_t) (INT_MAX/coeff)) || (rdispls[size-1] > (size_t) (INT_MAX/coeff))) {
    large = 1;
  }

  int s_large = 0;
  MPI_Allreduce (&large, &s_large, 1, MPI_INT,  MPI_SUM,  comm);
  large = s_large;

  if (!large) {

    int *_sdispls;
    int *_rdispls;
    PDM_malloc(_sdispls, size, int);
    PDM_malloc(_rdispls, size, int);

    for (int i = 0; i < size; i++) {
      _sdispls[i] = (int) sdispls[i];
      _rdispls[i] = (int) rdispls[i];
    }

    MPI_Alltoallv(sendbuf,
                  sendcounts,
                  _sdispls,
                  sendtype,
                  recvbuf,
                  recvcounts,
                  _rdispls,
                  recvtype,
                  comm);

    PDM_free(_sdispls);
    PDM_free(_rdispls);
  } else {

    MPI_Request *request_r;
    MPI_Request *request_s;
    PDM_malloc(request_r, size, MPI_Request);
    PDM_malloc(request_s, size, MPI_Request);

    int size_send_type;
    MPI_Type_size(sendtype, &size_send_type);

    int size_recv_type;
    MPI_Type_size(recvtype, &size_recv_type);

    for (int i = 0; i < size; i++) {
      if (recvcounts[i] != 0) {
        void *buf = (void *) ((unsigned char*) recvbuf + rdispls[i] * size_recv_type);
        code = MPI_Irecv(buf,
                         recvcounts[i],
                         recvtype,
                         i,
                         0,
                         comm,
                         request_r + i);
        if (code != MPI_SUCCESS) {
          break;
        }
      }
      if (sendcounts[i] != 0) {
        void *buf = (void *) ((unsigned char*) sendbuf + sdispls[i] * size_send_type);
        code = MPI_Isend(buf,
                         sendcounts[i],
                         sendtype,
                         i,
                         0,
                         comm,
                         request_s + i);
        if (code != MPI_SUCCESS) {
          break;
        }
      }
    }

    if (code != MPI_SUCCESS) {
      return code;
    }

    for (int i = 0; i < size; i++) {
      if (recvcounts[i] != 0) {
        code = MPI_Wait(request_r + i, MPI_STATUS_IGNORE);
      }
      if (code != MPI_SUCCESS) {
        break;
      }
    }

    if (code != MPI_SUCCESS) {
      return code;
    }

    for (int i = 0; i < size; i++) {
      if (sendcounts[i] != 0) {
        code = MPI_Wait(request_s + i, MPI_STATUS_IGNORE);
      }
      if (code != MPI_SUCCESS) {
        break;
      }
    }

    PDM_free(request_r);
    PDM_free(request_s);
  }

  return code;
}

#ifdef __cplusplus
}
#endif
