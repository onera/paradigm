/*
 * \file
 */

#ifndef __PDM_MPI_EXTENDED_H__
#define __PDM_MPI_EXTENDED_H__

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief Initializes multiple non-blocking send operations to a selected set of processes.
 *
 * This function is a wrapper for MPI_Send_init, designed to simplify the
 * initialization of multiple non-blocking sends to a predefined list of
 * destination ranks. It dynamically allocates an array of MPI requests and
 * initializes a separate send request for each active destination.
 *
 * This function does not perform any communication; it only prepares the
 * requests. The communication must be started later using functions like
 * MPI_Start, and completed with functions like MPI_Waitall or MPI_Testall.
 *
 * \param[in]  sendbuf        Send buffer
 * \param[in]  sendcounts     Number of elements to send to each active destination process (size = \p n_tgt_rank)
 * \param[in]  sdispls        Displacement (relative to \p sendbuf) to the data relative to each active destination process (size = \p n_tgt_rank)
 * \param[in]  datatype       Type of the data elements in \p sendbuf
 * \param[in]  n_tgt_rank     Number of active destination processes
 * \param[in]  tgt_rank       Ranks in \p comm of active destination processes (size = \p n_send_rank)
 * \param[in]  tag            Message tag to be used for the send operations
 * \param[in]  comm           MPI communicator
 * \param[out] out_requests   MPI requests
 *
 * \return An MPI error code. PDM_MPI_SUCCESS on success, otherwise an error code.
 *
 * \pre The arrays sendcounts, sdispls, and tgt_rank must be allocated and correctly sized (at least n_tgt_rank).
 * \post A non-null pointer to an array of MPI requests is returned via out_requests, which must be freed by the user after communication is complete.
 *
 * \see PDM_MPI_Recvs_init for the corresponding receive operation.
 * \see MPI_Send_init, MPI_Start, MPI_Waitall
 */
int
PDM_MPI_Sends_init
(
  const void              *sendbuf,
        int               *sendcounts,
        int               *sdispls,
        PDM_MPI_Datatype   datatype,
        int                n_tgt_rank,
        int               *tgt_rank,
        int                tag,
        PDM_MPI_Comm       comm,
        PDM_MPI_Request  **out_requests
);

/**
 * \brief Initializes and starts multiple non-blocking send operations to selected processes.
 *
 * This function is a convenience wrapper that combines the initialization of multiple non-blocking send requests and their immediate start. It acts as  single call to perform a series of non-blocking sends to a predefined list of destination ranks.
 *
 * The function allocates an array of MPI requests and initiates a separate `MPI_Isend` for each active destination, but unlike `PDM_MPI_Sends_init`, it does not create persistent requests.
 *
 * \param[in]  sendbuf        Send buffer
 * \param[in]  sendcounts     Number of elements to send to each active destination process (size = \p n_tgt_rank)
 * \param[in]  sdispls        Displacement (relative to \p sendbuf) to the data relative to each active destination process (size = \p n_tgt_rank)
 * \param[in]  datatype       Type of the data elements in \p sendbuf
 * \param[in]  n_tgt_rank     Number of active destination processes
 * \param[in]  tgt_rank       Ranks in \p comm of active destination processes (size = \p n_send_rank)
 * \param[in]  tag            Message tag to be used for the send operations
 * \param[in]  comm           MPI communicator
 * \param[out] out_requests   MPI requests
 *
 * \return An MPI error code. MPI_SUCCESS on success, otherwise an error code.
 *
 * \pre The arrays sendcounts, sdispls, and tgt_rank must be allocated and correctly sized (at least n_tgt_rank).
 * \post A non-null pointer to an array of MPI requests is returned via out_requests, which must be freed by the user after communication is completed (e.g., with MPI_Waitall).
 *
 * \see PDM_MPI_Sends_init for persistent requests.
 * \see PDM_MPI_Irecvs for the corresponding receive operation.
 * \see MPI_Isend, MPI_Waitall, MPI_Testall
 */
int
PDM_MPI_Isends
(
  const void              *sendbuf,
        int               *sendcounts,
        int               *sdispls,
        PDM_MPI_Datatype   datatype,
        int                n_tgt_rank,
        int               *tgt_rank,
        int                tag,
        PDM_MPI_Comm       comm,
        PDM_MPI_Request  **out_requests
);


/**
 * \brief Initializes and starts multiple non-blocking receive operations from selected processes.
 *
 * This function is a convenience wrapper that combines the initialization of multiple non-blocking receive requests and their immediate start. It acts as a single call to perform a series of non-blocking receives from a predefined list of source ranks.
 *
 * The function allocates an array of MPI requests and initiates a separate `MPI_Irecv` for each active source.
 *
 * \param[out] recvbuf        Receive buffer
 * \param[in]  recvcounts     Number of elements to receive from each active source process (size = \p n_src_rank)
 * \param[in]  rdispls        Displacement (relative to \p recvbuf) to the data relative to each active source process (size = \p n_src_rank)
 * \param[in]  datatype       Type of the data elements in \p recvbuf
 * \param[in]  n_src_rank     Number of active source processes
 * \param[in]  src_rank       Ranks in \p comm of active source processes (size = \p n_src_rank)
 * \param[in]  tag            Message tag to be used for the receive operations
 * \param[in]  comm           MPI communicator
 * \param[out] out_requests   MPI requests
 *
 * \return An MPI error code. MPI_SUCCESS on success, otherwise an error code.
 *
 * \pre The arrays recvcounts, sdispls, and src_rank must be allocated and correctly sized (at least n_src_rank).
 * \post A non-null pointer to an array of MPI requests is returned via out_requests, which must be freed by the user after communication is completed (e.g., with MPI_Waitall).
 *
 * \see PDM_MPI_Recvs_init for persistent requests.
 * \see PDM_MPI_Isends for the corresponding send operation.
 * \see MPI_Irecv, MPI_Waitall, MPI_Testall
 */
int
PDM_MPI_Irecvs
(
  const void              *recvbuf,
        int               *recvcounts,
        int               *rdispls,
        PDM_MPI_Datatype   datatype,
        int                n_src_rank,
        int               *src_rank,
        int                tag,
        PDM_MPI_Comm       comm,
        PDM_MPI_Request  **out_requests
);


/**
 * \brief Initializes multiple non-blocking receive operations from a selected set of processes.
 *
 * This function is a wrapper for MPI_Recv_init, designed to simplify the
 * initialization of multiple non-blocking receives from a predefined list of
 * source ranks. It dynamically allocates an array of MPI requests and
 * initializes a separate receive request for each active source.
 *
 * This function does not perform any communication; it only prepares the
 * requests. The communication must be started later using functions like
 * MPI_Start, and completed with functions like MPI_Waitall or MPI_Testall.
 *
 * \param[in]  recvbuf        Receive buffer
 * \param[in]  recvcounts     Number of elements to receive from each active source process (size = \p n_src_rank)
 * \param[in]  rdispls        Displacement (relative to \p recvbuf) to the data relative to each active source process (size = \p n_src_rank)
 * \param[in]  datatype       Type of the data elements in \p recvbuf
 * \param[in]  n_src_rank     Number of active source processes
 * \param[in]  src_rank       Ranks in \p comm of active source processes (size = \p n_src_rank)
 * \param[in]  tag            Message tag to be used for the receive operations
 * \param[in]  comm           MPI communicator
 * \param[out] out_requests   MPI requests

 *
 * \return An MPI error code. MPI_SUCCESS on success, otherwise an error code.
 *
 * \pre The arrays recvcounts, rdispls, and src_rank must be allocated and correctly sized (at least n_src_rank).
 * \post A non-null pointer to an array of MPI requests is returned via out_requests, which must be freed by the user after communication is complete.
 *
 * \see PDM_MPI_Sends_init for the corresponding send operation.
 * \see MPI_Recv_init, MPI_Start, MPI_Waitall
 */
int
PDM_MPI_Recvs_init
(
  void              *recvbuf,
  int               *recvcounts,
  int               *rdispls,
  PDM_MPI_Datatype   datatype,
  int                n_src_rank,
  int               *src_rank,
  int                tag,
  PDM_MPI_Comm       comm,
  PDM_MPI_Request  **out_requests
);

/**
 * \brief Calculates the global proportion of "active" processes based on communication counts.
 *
 * An active process is defined as any rank in the communicator for which the provided sendcounts or recvcounts arrays indicate a non-zero transfer.
 * The function determines this count locally by iterating over the provided arrays. It then calculates the local proportion (n_active_rank / total_size) and uses PDM_MPI_Allreduce (with the MAX operation) to ensure that the final, correct proportion (which should be identical on all ranks) is distributed back to all processes.
 *
 * This function is typically used in performance metrics or to optimize collective
 * operations based on the communication sparsity defined by the counts arrays (e.g., in an Alltoallv).
 *
 * @note It is assumed that the sendcounts and recvcounts arrays are globally consistent
 * (i.e., they hold the full communication matrix counts) across all participating ranks.
 *
 * \param [in]  sendcounts        Number of elements to send to each rank in \p comm
 * \param [in]  recvcounts        Number of elements to receive from each rank in \p comm
 * \param [in]  comm              MPI communicator
 * \param [out] part_active_rank  Global proportion of active ranks (between 0 and 1)
 */
void
PDM_MPI_Partofactiverank
(
  int          *sendcounts,
  int          *recvcounts,
  PDM_MPI_Comm  comm,
  double       *part_active_rank
);

/**
 * \brief Performs a non-blocking all-to-all communication for a selected, sparse set of ranks using P2P messages.
 *
 * This function emulates MPI_Ialltoallv for a sparse communication pattern by posting non-blocking
 * P2P sends (MPI_Issend) and receives (MPI_Irecv) to a user-defined list of ranks. It is particularly
 * useful when the communication partners are known explicitly and represent a small subset of the
 * communicator's ranks. The function returns an array of requests that must be waited on
 * (e.g., with PDM_MPI_Waitall) to ensure the completion of the communication.
 *
 * \param[in]  sendbuf              Send buffer
 * \param[in]  sendcounts           Number of elements to send to each destination process (size = \p n_send_rank)
 * \param[in]  sdispls              Displacement (relative to \p sendbuf) to the data relative to each destination process (size = \p n_send_rank)
 * \param[in]  sendtype             Type of the data elements in \p sendbuf
 * \param[in]  n_send_rank          Number of destination processes
 * \param[in]  send_rank            Ranks in \p comm of destination processes (size = \p n_send_rank)
 * \param[out] recvbuf              Receive buffer
 * \param[in]  recvcounts           Number of elements to receive from each source process (size = \p n_recv_rank)
 * \param[in]  rdispls              Displacement (relative to \p recvbuf) to the data relative to each source process (size = \p n_recv_rank)
 * \param[in]  recvtype             Type of the data elements in \p recvbuf
 * \param[in]  n_recv_rank          Number of source processes
 * \param[in]  recv_rank            Ranks in \p comm of source processes (size = \p n_recv_rank)
 * \param[in]  tag                  The message tag for P2P communication
 * \param[in]  comm                 MPI communicator
 * \param[out] n_send_recv_request  Total number of MPI requests
 * \param[out] out_requests         MPI requests
 *
 * \return PDM_SUCCESS or an error code from the underlying MPI calls.
 */
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
);


/**
 * \brief Emulates MPI_Ialltoallv using non-blocking P2P messages.
 *
 * This function provides a flexible implementation of a non-blocking all-to-all communication.
 * It dynamically posts non-blocking P2P sends and receives based on non-zero entries in the
 * sendcounts and recvcounts arrays. If explicit rank lists are provided (send_rank/recv_rank),
 * it delegates to the more specialized PDM_MPI_Ialltoallv_select_p2p function. Otherwise, it
 * assumes a dense communication pattern and iterates through all ranks of the communicator.
 * This function is useful for scenarios where a native collective may not be optimal
 * (e.g., for sparse communication patterns or debugging).
 *
 * \param[in]  sendbuf              Send buffer
 * \param[in]  sendcounts           Number of elements to send to each destination process (size = \p n_send_rank if provided, else size of \p comm)
 * \param[in]  sdispls              Displacement (relative to \p sendbuf) to the data relative to each destination process (size = \p n_send_rank if provided, else size of \p comm)
 * \param[in]  sendtype             Type of the data elements in \p sendbuf
 * \param[in]  n_send_rank          (Optional) Number of destination processes (only used if \p send_rank is not NULL)
 * \param[in]  send_rank            (Optional) Ranks in \p comm of destination processes (size = \p n_send_rank if provided, else size of \p comm)
 *                                  If NULL, a dense check on \p sendcounts is performed.
 * \param[out] recvbuf              Receive buffer
 * \param[in]  recvcounts           Number of elements to receive from each source process (size = \p n_recv_rank if provided, else size of \p comm)
 * \param[in]  rdispls              Displacement (relative to \p recvbuf) to the data relative to each source process (size = \p n_recv_rank if provided, else size of \p comm)
 * \param[in]  recvtype             Type of the data elements in \p recvbuf
 * \param[in]  n_recv_rank          (Optional) Number of source processes (only used if \p recv_rank is not NULL)
 * \param[in]  recv_rank            (Optional) Ranks in \p comm of source processes (size = \p n_recv_rank if provided, else size of \p comm)
 *                                  If NULL, a dense check on \p recvcounts is performed.
 * \param[in]  tag                  The message tag for P2P communication
 * \param[in]  comm                 MPI communicator
 * \param[out] n_send_recv_request  Total number of MPI requests
 * \param[out] out_requests         MPI requests
 *
 * \return PDM_SUCCESS or an error code from the underlying MPI calls.
 */
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
);

/**
 * \brief Emulates MPI_Ialltoallv using non-blocking One-Sided (RMA) communication via MPI_Rget.
 *
 * This function performs an **asynchronous, non-blocking** all-to-all communication based on the
 * Remote Memory Access (RMA) paradigm, specifically by using **MPI_Rget** operations.
 *
 * For every rank from which the local process expects to receive data (recvcounts > 0),
 * an `MPI_Rget` is posted to pull data from the remote process's exposed memory window
 * (send_win) into the local receive buffer (recvbuf). The operation is non-blocking,
 * and the returned requests must be completed using `PDM_MPI_Wait` or equivalent synchronization
 * before the received data can be accessed and before the RMA epoch ends.
 *
 * \param[in]  send_win             MPI Window handle (\ref PDM_MPI_Win) on the *target* process from which data will be retrieved
 * \param[in]  target_disp          Displacement in the remote target window (\p send_win) to the start of the data retrieved from each rank
 * \param[out] recvbuf              Local receive buffer where the data will be placed (origin buffer of the Rget)
 * \param[in]  recvcounts           Number of elements to receive from each rank
 * \param[in]  rdispls              Displacement in the local (\p recvbuf) to the start of each received message
 * \param[in]  recvtype             Type of the data elements in \p recvbuf and used for the target window access
 * \param[in]  comm                 MPI communicator
 * \param[out] n_send_recv_request  Total number of allocated MPI_Rget requests (equal to the number of non-zero entries in \p recvcounts)
 * \param[out] out_requests         Non-blocking MPI_Rget requests
 *
 * \return PDM_SUCCESS or an error code from the underlying MPI calls.
 */
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
);



/**
 * \brief Performs a standard blocking All-to-all communication with variable counts and large (size_t) displacements.
 *
 * This function serves as the **Large-Offset (L) wrapper** around the native MPI_Alltoallv (or equivalent MPI_Ialltoallv/Wait)
 * but uses size_t for the displacement arrays (sdispls, rdispls) to support extremely large data transfers
 * that exceed the limitations of standard integer displacement arrays. The underlying implementation will typically
 * rely on a specialized MPI function (like MPI_Ialltoallv with custom datatype offsets or MPIX_Alltoallv_ll) if available,
 * or fall back to an internal large-offset P2P emulation. The call is blocking.
 *
 * \param[in]  sendbuf     Send buffer
 * \param[in]  sendcounts  Number of elements to send to each rank
 * \param[in]  sdispls     Displacement in \p sendbuf to the start of each message (*large* offset)
 * \param[in]  sendtype    The datatype of send buffer elements
 * \param[out] recvbuf     Receive buffer
 * \param[in]  recvcounts  Number of elements to receive from each rank
 * \param[in]  rdispls     Displacement in \p recvbuf to the start of each message (*large* offset)
 * \param[in]  recvtype    The datatype of receive buffer elements
 * \param[in]  comm        MPI communicator
 *
 * \return PDM_SUCCESS or an error code from the underlying MPI calls.
 */
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
);

/**
 * \brief Emulates MPI_Alltoallv using blocking P2P messages.
 *
 * This function performs a **blocking** all-to-all communication by manually posting point-to-point (P2P) sends and receives.
 * The function returns only after all data has been safely sent and received. It is typically used for debugging or to
 * bypass potential performance issues with native collective implementations in sparse communication scenarios.
 *
 * \param[in]  sendbuf     Send buffer
 * \param[in]  sendcounts  Number of elements to send to each rank
 * \param[in]  sdispls     Displacement in \p sendbuf to the start of each message
 * \param[in]  sendtype    The datatype of send buffer elements
 * \param[out] recvbuf     Receive buffer
 * \param[in]  recvcounts  Number of elements to receive from each rank
 * \param[in]  rdispls     Displacement in \p recvbuf to the start of each message
 * \param[in]  recvtype    The datatype of receive buffer elements
 * \param[in]  comm        MPI  communicator
 *
 * \return PDM_SUCCESS or an error code from the underlying MPI calls.
 */
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
);


/**
 * \brief Emulates MPI_Alltoallv for very large data transfers using blocking P2P messages and size_t offsets.
 *
 * This is the Large-Offset (L) version of the blocking P2P all-to-all communication. It uses size_t for
 * displacement arrays (sdispls and rdispls) to support extremely large buffers or complex memory layouts
 * where 32-bit integer offsets are insufficient. The operation is blocking, returning only after all transfers are complete.
 *
 * \param[in]  sendbuf     Send buffer
 * \param[in]  sendcounts  Number of elements to send to each rank
 * \param[in]  sdispls     Displacement in \p sendbuf to the start of each message (*large* offset)
 * \param[in]  sendtype    The datatype of send buffer elements
 * \param[out] recvbuf     Receive buffer
 * \param[in]  recvcounts  Number of elements to receive from each rank
 * \param[in]  rdispls     Displacement in \p recvbuf to the start of each message (*large* offset)
 * \param[in]  recvtype    The datatype of receive buffer elements
 * \param[in]  comm        MPI communicator
 *
 * \return PDM_SUCCESS or an error code from the underlying MPI calls.
 */
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
);

/**
 * \brief Initializes a persistent P2P communication for AlltoallV exchanges.
 *
 * This function sets up the structure for a repeated all-to-all communication using non-blocking P2P persistent requests
 * (based on MPI_Send_init/MPI_Recv_init). The actual data transfer must be initiated later using PDM_MPI_Start/PDM_MPI_Startall
 * and completed with PDM_MPI_Wait/PDM_MPI_Waitall. This is highly efficient when the communication pattern (counts and displacements)
 * remains constant across multiple steps.
 *
 * \param[in]  sendbuf              Send buffer
 * \param[in]  sendcounts           Number of elements to send to each rank
 * \param[in]  sdispls              Displacement in \p sendbuf to the start of each message
 * \param[in]  sendtype             The datatype of send buffer elements
 * \param[out] recvbuf              Receive buffer
 * \param[in]  recvcounts           Number of elements to receive from each rank
 * \param[in]  rdispls              Displacement in \p recvbuf to the start of each message
 * \param[in]  recvtype             The datatype of receive buffer elements
 * \param[in]  tag                  The message tag for P2P persistent communication
 * \param[in]  comm                 MPI communicator
 * \param[out] n_send_recv_request  Total number of allocated persistent requests
 * \param[out] requests             MPI requests for persistent communications
 *
 * \return PDM_SUCCESS or an error code from the underlying MPI calls.
 */
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
  PDM_MPI_Request  **requests
);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MPI_EXTENDED_H__ */
