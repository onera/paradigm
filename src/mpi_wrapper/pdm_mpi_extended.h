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
 * \param[in]  sendbuf Pointer to the start of the send buffer.
 * \param[in]  sendcounts Array of size n_active_send, where sendcounts[i] is
 *             the number of elements to send to the i-th active
 *             destination.
 * \param[in]  sdispls Array of size n_active_send, where sdispls[i] is the
 *             displacement from sendbuf to the starting element for
 *             the i-th active destination.
 * \param[in]  datatype The type of the data elements in the send buffer.
 * \param[in]  n_active_send The number of active destination processes.
 * \param[in]  active_send Array of size n_active_send, containing the ranks
 *             of the destination processes.
 * \param[in]  tag The message tag to be used for the send operations.
 * \param[in]  comm The MPI communicator to be used.
 * \param[out] out_requests Pointer to a pointer that will be set to the
 * dynamically allocated array of MPI requests.
 *
 * \return An MPI error code. PDM_MPI_SUCCESS on success, otherwise an error code.
 *
 * \pre The arrays sendcounts, sdispls, and active_send must be allocated
 * and correctly sized (at least n_active_send).
 * \post A non-null pointer to an array of MPI requests is returned via out_requests,
 * which must be freed by the user after communication is complete.
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
        int                n_active_send,
        int               *active_send,
        int                tag,
        PDM_MPI_Comm       comm,
        PDM_MPI_Request  **out_requests
);

/**
 * \brief Initializes and starts multiple non-blocking send operations to selected processes.
 *
 * This function is a convenience wrapper that combines the initialization of
 * multiple non-blocking send requests and their immediate start. It acts as
 * a single call to perform a series of non-blocking sends to a predefined
 * list of destination ranks.
 *
 * The function allocates an array of MPI requests and initiates a separate
 * `MPI_Isend` for each active destination, but unlike `PDM_MPI_Sends_init`,
 * it does not create persistent requests.
 *
 * \param[in] sendbuf Pointer to the start of the send buffer.
 * \param[in] sendcounts Array of size n_active_send, where sendcounts[i] is
 * the number of elements to send to the i-th active destination.
 * \param[in] sdispls Array of size n_active_send, where sdispls[i] is the
 * displacement from sendbuf to the starting element for the i-th active
 * destination.
 * \param[in] datatype The type of the data elements in the send buffer.
 * \param[in] n_active_send The number of active destination processes.
 * \param[in] active_send Array of size n_active_send, containing the ranks
 * of the destination processes.
 * \param[in] tag The message tag to be used for the send operations.
 * \param[in] comm The MPI communicator to be used.
 * \param[out] out_requests Pointer to a pointer that will be set to the
 * dynamically allocated array of non-blocking requests.
 *
 * \return An MPI error code. MPI_SUCCESS on success, otherwise an error code.
 *
 * \pre The arrays sendcounts, sdispls, and active_send must be allocated
 * and correctly sized (at least n_active_send).
 * \post A non-null pointer to an array of MPI requests is returned via
 * out_requests, which must be freed by the user after communication is
 * completed (e.g., with MPI_Waitall).
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
        int                n_active_send,
        int               *active_send,
        int                tag,
        PDM_MPI_Comm       comm,
        PDM_MPI_Request  **out_requests
);


/**
 * \brief Initializes and starts multiple non-blocking receive operations from selected processes.
 *
 * This function is a convenience wrapper that combines the initialization of
 * multiple non-blocking receive requests and their immediate start. It acts as
 * a single call to perform a series of non-blocking receives from a predefined
 * list of source ranks.
 *
 * The function allocates an array of MPI requests and initiates a separate
 * `MPI_Irecv` for each active source.
 *
 * \param[out] recvbuf Pointer to the start of the receive buffer.
 * \param[in] recvcounts Array of size n_active_recv, where recvcounts[i] is
 * the number of elements to receive from the i-th active source.
 * \param[in] sdispls Array of size n_active_recv, where sdispls[i] is the
 * displacement from recvbuf to the starting element for the i-th active source.
 * \param[in] datatype The type of the data elements in the receive buffer.
 * \param[in] n_active_recv The number of active source processes.
 * \param[in] active_recv Array of size n_active_recv, containing the ranks
 * of the source processes.
 * \param[in] tag The message tag to be used for the receive operations.
 * \param[in] comm The MPI communicator to be used.
 * \param[out] out_requests Pointer to a pointer that will be set to the
 * dynamically allocated array of non-blocking requests.
 *
 * \return An MPI error code. MPI_SUCCESS on success, otherwise an error code.
 *
 * \pre The arrays recvcounts, sdispls, and active_recv must be allocated
 * and correctly sized (at least n_active_recv).
 * \post A non-null pointer to an array of MPI requests is returned via
 * out_requests, which must be freed by the user after communication is
 * completed (e.g., with MPI_Waitall).
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
        int               *sdispls,
        PDM_MPI_Datatype   datatype,
        int                n_active_recv,
        int               *active_recv,
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
 * \param[in]  recvbuf Pointer to the start of the receive buffer.
 * \param[in]  recvcounts Array of size n_active_recv, where recvcounts[i] is
 *             the number of elements to receive from the i-th active
 *             source.
 * \param[in]  rdispls Array of size n_active_recv, where rdispls[i] is the
 *             displacement from recvbuf to the starting element for
 *             the i-th active source.
 * \param[in]  datatype The type of the data elements in the receive buffer.
 * \param[in]  n_active_recv The number of active source processes.
 * \param[in]  active_recv Array of size n_active_recv, containing the ranks
 *             of the source processes.
 * \param[in]  tag The message tag to be used for the receive operations.
 * \param[in]  comm The MPI communicator to be used.
 * \param[out] out_requests Pointer to a pointer that will be set to the
 *             dynamically allocated array of MPI requests.
 *
 * \return An MPI error code. MPI_SUCCESS on success, otherwise an error code.
 *
 * \pre The arrays recvcounts, rdispls, and active_recv must be allocated
 *      and correctly sized (at least n_active_recv).
 * \post A non-null pointer to an array of MPI requests is returned via out_requests,
 *       which must be freed by the user after communication is complete.
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
  int                n_active_recv,
  int               *active_recv,
  int                tag,
  PDM_MPI_Comm       comm,
  PDM_MPI_Request  **out_requests
);

void
PDM_MPI_Partofactiverank
(
  int          *sendcounts,
  int          *recvcounts,
  PDM_MPI_Comm  comm,
  double       *part_active_rank
);

/**
 * @brief Performs a non-blocking all-to-all communication for a selected, sparse set of ranks using P2P messages.
 *
 * This function emulates MPI_Ialltoallv for a sparse communication pattern by posting non-blocking
 * P2P sends (MPI_Issend) and receives (MPI_Irecv) to a user-defined list of ranks. It is particularly
 * useful when the communication partners are known explicitly and represent a small subset of the
 * communicator's ranks. The function returns an array of requests that must be waited on
 * (e.g., with PDM_MPI_Waitall) to ensure the completion of the communication.
 *
 * \param sendbuf       The send buffer.
 * \param sendcounts    An array of integers specifying the number of elements to send to each rank in @a send_rank.
 * \param sdispls       An array of integers specifying the displacement in @a sendbuf for each message.
 * \param sendtype      The datatype of send buffer elements.
 * \param n_send_rank   The number of ranks to send to.
 * \param send_rank     An array of integers containing the ranks to send to.
 * \param recvbuf       The receive buffer.
 * \param recvcounts    An array of integers specifying the number of elements to receive from each rank in @a recv_rank.
 * \param rdispls       An array of integers specifying the displacement in @a recvbuf for each message.
 * \param recvtype      The datatype of receive buffer elements.
 * \param n_recv_rank   The number of ranks to receive from.
 * \param recv_rank     An array of integers containing the ranks to receive from.
 * \param tag           The message tag for P2P communication.
 * \param comm          The communicator.
 * \param n_send_recv_request A pointer to an integer that will be filled with the total number of requests.
 * \param out_requests  A pointer to a PDM_MPI_Request array that will be allocated and filled with the requests.
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
 * @brief Emulates MPI_Ialltoallv using non-blocking P2P messages.
 *
 * This function provides a flexible implementation of a non-blocking all-to-all communication.
 * It dynamically posts non-blocking P2P sends and receives based on non-zero entries in the
 * sendcounts and recvcounts arrays. If explicit rank lists are provided (@a send_rank/@a recv_rank),
 * it delegates to the more specialized PDM_MPI_Ialltoallv_select_p2p function. Otherwise, it
 * assumes a dense communication pattern and iterates through all ranks of the communicator.
 * This function is useful for scenarios where a native collective may not be optimal
 * (e.g., for sparse communication patterns or debugging).
 *
 * \param sendbuf       The send buffer.
 * \param sendcounts    An array of integers specifying the number of elements to send to each rank.
 * \param sdispls       An array of integers specifying the displacement in @a sendbuf for each message.
 * \param sendtype      The datatype of send buffer elements.
 * \param n_send_rank   (Optional) The number of ranks to send to. Used only if @a send_rank is not NULL.
 * \param send_rank     (Optional) An array of integers containing the ranks to send to. If NULL, a dense check on sendcounts is performed.
 * \param recvbuf       The receive buffer.
 * \param recvcounts    An array of integers specifying the number of elements to receive from each rank.
 * \param rdispls       An array of integers specifying the displacement in @a recvbuf for each message.
 * \param recvtype      The datatype of receive buffer elements.
 * \param n_recv_rank   (Optional) The number of ranks to receive from. Used only if @a recv_rank is not NULL.
 * \param recv_rank     (Optional) An array of integers containing the ranks to receive from. If NULL, a dense check on recvcounts is performed.
 * \param tag           The message tag for P2P communication.
 * \param comm          The communicator.
 * \param n_send_recv_request A pointer to an integer that will be filled with the total number of requests.
 * \param out_requests  A pointer to a PDM_MPI_Request array that will be allocated and filled with the requests.
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
 * @brief Emulates MPI_Ialltoallv using non-blocking One-Sided (RMA) communication via MPI_Rget.
 *
 * This function performs an **asynchronous, non-blocking** all-to-all communication based on the
 * Remote Memory Access (RMA) paradigm, specifically by using **MPI_Rget** operations.
 *
 * For every rank from which the local process expects to receive data (@a recvcounts > 0),
 * an `MPI_Rget` is posted to pull data from the remote process's exposed memory window
 * (@a send_win) into the local receive buffer (@a recvbuf). The operation is non-blocking,
 * and the returned requests must be completed using `PDM_MPI_Wait` or equivalent synchronization
 * before the received data can be accessed and before the RMA epoch ends.
 *
 * \param send_win       The MPI Window handle (PDM_MPI_Win) on the **target** process from which data will be retrieved.
 * \param target_disp    An array of integers specifying the **displacement** (offset) in the remote target window (@a send_win) for the start of the data retrieval from each rank.
 * \param recvbuf        The local receive buffer where the data will be placed (origin buffer of the Rget).
 * \param recvcounts     An array of integers specifying the number of elements to **receive** (retrieve) from each rank.
 * \param rdispls        An array of integers specifying the **displacement** in the local @a recvbuf for each received message.
 * \param recvtype       The datatype of receive buffer elements and the type used for the target window access.
 * \param comm           The communicator.
 * \param n_send_recv_request A pointer to an integer that will be filled with the total number of allocated `MPI_Rget` requests (equal to the number of non-zero entries in @a recvcounts).
 * \param out_requests   A pointer to a PDM_MPI_Request array that will be allocated and filled with the non-blocking `MPI_Rget` requests.
 * @return PDM_SUCCESS or an error code from the underlying MPI calls.
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
 * @brief Emulates MPI_Ialltoallv for very large data transfers using non-blocking P2P messages and size_t offsets.
 *
 * This is the Large-Offset (L) version of the non-blocking P2P all-to-all communication, designed to handle cases
 * where the buffer displacements exceed the limits of a standard 32-bit integer (e.g., for data structures larger than 2GB).
 * It uses @a size_t for displacement arrays (@a sdispls and @a rdispls). Communication still relies on dynamic non-blocking
 * P2P sends and receives.
 *
 * \param sendbuf        The send buffer.
 * \param sendcounts     An array of integers specifying the number of elements to send to each rank.
 * \param sdispls        An array of size_t specifying the **displacement** in @a sendbuf for each message (Large offset).
 * \param sendtype       The datatype of send buffer elements.
 * \param recvbuf        The receive buffer.
 * \param recvcounts     An array of integers specifying the number of elements to receive from each rank.
 * \param rdispls        An array of size_t specifying the **displacement** in @a recvbuf for each message (Large offset).
 * \param recvtype       The datatype of receive buffer elements.
 * \param comm           The communicator.
 * \param request_s      A pointer to a PDM_MPI_Request array that will be allocated and filled with the send requests.
 * \param request_r      A pointer to a PDM_MPI_Request array that will be allocated and filled with the receive requests.
 * \param n_request_s    A pointer to an integer that will be filled with the total number of send requests.
 * \param n_request_r    A pointer to an integer that will be filled with the total number of receive requests.
 * \return PDM_SUCCESS or an error code from the underlying MPI calls.
 */
int
PDM_MPI_Ialltoallv_p2p_l
(
  void              *sendbuf,
  int               *sendcounts,
  size_t            *sdispls,
  PDM_MPI_Datatype   sendtype,
  void              *recvbuf,
  int               *recvcounts,
  size_t            *rdispls,
  PDM_MPI_Datatype   recvtype,
  PDM_MPI_Comm       comm,
  PDM_MPI_Request  **request_s,
  PDM_MPI_Request  **request_r,
  int               *n_request_s,
  int               *n_request_r
);


/**
 * @brief Performs a standard blocking All-to-all communication with variable counts and large (size_t) displacements.
 *
 * This function serves as the **Large-Offset (L) wrapper** around the native MPI_Alltoallv (or equivalent MPI_Ialltoallv/Wait)
 * but uses @a size_t for the displacement arrays (@a sdispls, @a rdispls) to support extremely large data transfers
 * that exceed the limitations of standard integer displacement arrays. The underlying implementation will typically
 * rely on a specialized MPI function (like MPI_Ialltoallv with custom datatype offsets or MPIX_Alltoallv_ll) if available,
 * or fall back to an internal large-offset P2P emulation. The call is blocking.
 *
 * \param sendbuf        The send buffer.
 * \param sendcounts     An array of integers specifying the number of elements to send to each rank.
 * \param sdispls        An array of size_t specifying the **displacement** in @a sendbuf for each message (Large offset).
 * \param sendtype       The datatype of send buffer elements.
 * \param recvbuf        The receive buffer.
 * \param recvcounts     An array of integers specifying the number of elements to receive from each rank.
 * \param rdispls        An array of size_t specifying the **displacement** in @a recvbuf for each message (Large offset).
 * \param recvtype       The datatype of receive buffer elements.
 * \param comm           The communicator.
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
 * @brief Emulates MPI_Alltoallv using blocking P2P messages.
 *
 * This function performs a **blocking** all-to-all communication by manually posting point-to-point (P2P) sends and receives.
 * The function returns only after all data has been safely sent and received. It is typically used for debugging or to
 * bypass potential performance issues with native collective implementations in sparse communication scenarios.
 *
 * \param sendbuf        The send buffer.
 * \param sendcounts     An array of integers specifying the number of elements to send to each rank.
 * \param sdispls        An array of integers specifying the displacement in @a sendbuf for each message.
 * \param sendtype       The datatype of send buffer elements.
 * \param recvbuf        The receive buffer.
 * \param recvcounts     An array of integers specifying the number of elements to receive from each rank.
 * \param rdispls        An array of integers specifying the displacement in @a recvbuf for each message.
 * \param recvtype       The datatype of receive buffer elements.
 * \param comm           The communicator.
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
 * @brief Emulates MPI_Alltoallv for very large data transfers using blocking P2P messages and size_t offsets.
 *
 * This is the Large-Offset (L) version of the blocking P2P all-to-all communication. It uses @a size_t for
 * displacement arrays (@a sdispls and @a rdispls) to support extremely large buffers or complex memory layouts
 * where 32-bit integer offsets are insufficient. The operation is blocking, returning only after all transfers are complete.
 *
 * \param sendbuf        The send buffer.
 * \param sendcounts     An array of integers specifying the number of elements to send to each rank.
 * \param sdispls        An array of size_t specifying the **displacement** in @a sendbuf for each message (Large offset).
 * \param sendtype       The datatype of send buffer elements.
 * \param recvbuf        The receive buffer.
 * \param recvcounts     An array of integers specifying the number of elements to receive from each rank.
 * \param rdispls        An array of size_t specifying the **displacement** in @a recvbuf for each message (Large offset).
 * \param recvtype       The datatype of receive buffer elements.
 * \param comm           The communicator.
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
 * @brief Initializes a persistent P2P communication for AlltoallV exchanges.
 *
 * This function sets up the structure for a repeated all-to-all communication using non-blocking P2P persistent requests
 * (based on MPI_Send_init/MPI_Recv_init). The actual data transfer must be initiated later using PDM_MPI_Start/PDM_MPI_Startall
 * and completed with PDM_MPI_Wait/PDM_MPI_Waitall. This is highly efficient when the communication pattern (counts and displacements)
 * remains constant across multiple steps.
 *
 * \param sendbuf        The send buffer.
 * \param sendcounts     An array of integers specifying the number of elements to send to each rank.
 * \param sdispls        An array of integers specifying the displacement in @a sendbuf for each message.
 * \param sendtype       The datatype of send buffer elements.
 * \param recvbuf        The receive buffer.
 * \param recvcounts     An array of integers specifying the number of elements to receive from each rank.
 * \param rdispls        An array of integers specifying the displacement in @a recvbuf for each message.
 * \param recvtype       The datatype of receive buffer elements.
 * \param tag            The message tag for P2P persistent communication.
 * \param comm           The communicator.
 * \param n_send_recv_request A pointer to an integer that will be filled with the total number of allocated persistent requests.
 * \param requests       A pointer to a PDM_MPI_Request array that will be allocated and filled with the persistent requests.
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
