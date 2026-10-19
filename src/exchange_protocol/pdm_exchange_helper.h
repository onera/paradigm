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

/**
 * \brief Creates and initializes an exchange helper object.
 *
 * This function allocates memory for a new PDM_exchange_helper_t object
 * and initializes its internal state. It is the constructor for the exchange
 * helper and must be called before using any other functions in this module.
 *
 * \param[in] comm           The MPI communicator to be used for all exchange operations.
 * \param[in] n_request_init The initial number of non-blocking requests that
 *                           can be stored. This value can be a hint; the helper may dynamically resize.
 *
 * \return A pointer to the newly created PDM_exchange_helper_t object on success,
 * or a null pointer on failure.
 *
 * \see PDM_exchange_helper_free to destroy the object.
 */
PDM_exchange_helper_t *
PDM_exchange_helper_create
(
  const PDM_MPI_Comm    comm
);


/**
 * \brief Performs a blocking exchange of data between processes.
 *
 * This function handles both sending and receiving data in a single,
 * blocking call. It is typically used for `all-to-all` or `point-to-point`
 * exchanges where the sending and receiving data sizes and indices are
 * known beforehand.
 *
 * \param[in,out] exch_helper The initialized exchange helper object.
 * \param[in]     k_comm      The communication kind (e.g., all-to-all, point-to-point).
 * \param[in]     mpi_type    An predefinid mpi type see \ref PDM_MPI_Datatype
 * \param[in]     send_idx    Array of destination ranks for sending data.
 * \param[in]     send_n      Array of data counts to be sent to each destination.
 * \param[in]     send_buffer Pointer to the data buffer for sending.
 * \param[in]     recv_idx    Array of source ranks for receiving data.
 * \param[in]     recv_n      Array of data counts to be received from each source.
 * \param[out]    recv_buffer Pointer to the data buffer for receiving.
 */
void
PDM_exchange_helper_mpi_type_exch
(
  PDM_exchange_helper_t *exch_helper,
  PDM_mpi_comm_kind_t    k_comm,
  PDM_MPI_Datatype       mpi_type,
  int                   *send_idx,
  int                   *send_n,
  void                  *send_buffer,
  int                   *recv_idx,
  int                   *recv_n,
  void                  *recv_buffer
);

/**
 * \brief Performs a blocking exchange of data between processes.
 *
 * This function handles both sending and receiving data in a single,
 * blocking call. It is typically used for `all-to-all` or `point-to-point`
 * exchanges where the sending and receiving data sizes and indices are
 * known beforehand.
 *
 * \param[in,out] exch_helper The initialized exchange helper object.
 * \param[in]     k_comm The communication kind (e.g., all-to-all, point-to-point).
 * \param[in]     cst_stride The stride for non-contiguous send/receive buffers.
 * \param[in]     s_data The size of a single data element in bytes.
 * \param[in]     send_idx Array of destination ranks for sending data.
 * \param[in]     send_n Array of data counts to be sent to each destination.
 * \param[in]     send_buffer Pointer to the data buffer for sending.
 * \param[in]     recv_idx Array of source ranks for receiving data.
 * \param[in]     recv_n Array of data counts to be received from each source.
 * \param[out]    recv_buffer Pointer to the data buffer for receiving.
 */
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

/**
 * \brief Initializes and starts a non-blocking exchange of data.
 *
 * This function initiates a non-blocking exchange and immediately returns.
 * The operation must be completed later by calling a wait function. This
 * allows other work to be performed concurrently with communication.
 *
 * \param[in,out] exch_helper The initialized exchange helper object.
 * \param[in]     k_comm The communication kind (e.g., all-to-all, point-to-point).
 * \param[in]     s_data The size of a single data element in bytes.
 * \param[in]     cst_stride The stride for non-contiguous send/receive buffers.
 * \param[in]     send_idx Array of destination ranks for sending data.
 * \param[in]     send_n Array of data counts to be sent to each destination.
 * \param[in]     send_buffer Pointer to the data buffer for sending.
 * \param[in]     recv_idx Array of source ranks for receiving data.
 * \param[in]     recv_n Array of data counts to be received from each source.
 * \param[out]    recv_buffer Pointer to the data buffer for receiving.
 *
 * \return An ID for the non-blocking request, to be used with wait function.
 *
 * \see PDM_exchange_helper_exch_wait to complete the operation.
 */
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

/**
 * \brief Prepares a persistent, non-blocking exchange request.
 *
 * This function initializes a persistent communication request for a
 * two-way exchange. The communication itself is not started. This is useful
 * when the same communication pattern (e.g., same ranks and counts) will be
 * repeated multiple times.
 *
 * \param[in,out] exch_helper The initialized exchange helper object.
 * \param[in] k_comm The communication kind (e.g., all-to-all, point-to-point).
 * \param[in] s_data The size of a single data element in bytes.
 * \param[in] cst_stride The stride for non-contiguous send/receive buffers.
 * \param[in] send_idx Array of destination ranks for sending data.
 * \param[in] send_n Array of data counts to be sent to each destination.
 * \param[in] send_buffer Pointer to the data buffer for sending.
 * \param[in] recv_idx Array of source ranks for receiving data.
 * \param[in] recv_n Array of data counts to be received from each source.
 * \param[out] recv_buffer Pointer to the data buffer for receiving.
 *
 * \return An ID for the persistent request, to be used with start/wait/free functions.
 *
 * \see PDM_exchange_helper_exch_start and PDM_exchange_helper_exch_wait.
 */
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


/**
 * \brief Starts a non-blocking one-way communication operation.
 *
 * This function initiates a one-way exchange of data in a non-blocking manner
 * and returns immediately with a request ID. The actual communication happens
 * asynchronously in the background. The operation can be a set of send or
 * receive operations, as specified by the direction parameter.
 *
 * The caller must use a `wait` function with the returned request ID to ensure
 * the communication has completed before reusing the data buffer.
 *
 * \param[in,out] exch_helper The initialized exchange helper object.
 * \param[in]     direction The direction of communication, either sending or receiving.
 * \param[in]     s_data The size of a single data element in bytes.
 * \param[in]     cst_stride The stride for non-contiguous data within the buffer.
 * \param[in]     n_active_rank The number of ranks involved in the communication.
 * \param[in]     active_rank Array of size n_active_rank, containing the ranks
 *                of the destination (for sends) or source (for receives) processes.
 * \param[in]     send_or_recv_idx Array of size n_active_rank, where each element
 *                is the displacement from the buffer to the starting element for the corresponding active rank.
 * \param[in]     send_or_recv_n Array of size n_active_rank, where each element
 *                is the number of elements to send or receive for the corresponding active rank.
 * \param[in]     tag The MPI message tag to be used for the exchange.
 * \param[in,out] buffer Pointer to the data buffer for the send or receive operation.
 *
 * \return The integer ID of the non-blocking request, which must be passed to
 *         a corresponding wait function to complete the communication. \see PDM_exchange_helper_exch_wait
 */
int
PDM_exchange_helper_iexch_one_way
(
  PDM_exchange_helper_t    *exch_helper,
  PDM_exchange_direction_t  direction,
  size_t                    s_data,
  int                       cst_stride,
  int                       n_active_rank,
  int                      *active_rank,
  int                      *send_or_recv_idx,
  int                      *send_or_recv_n,
  int                       tag,
  void                     *buffer
);


/**
 * \brief Prepares a persistent one-way communication request.
 *
 * This function is similar to `PDM_exchange_helper_exch_init` but is
 * specialized for one-way communication (either all sends or all receives).
 * This is useful for asymmetric communication patterns.
 *
 * \param[in,out] exch_helper The initialized exchange helper object.
 * \param[in]     direction The direction of communication (send or receive).
 * \param[in]     s_data The size of a single data element in bytes.
 * \param[in]     cst_stride The stride for non-contiguous send/receive buffers.
 * \param[in]     n_active_rank The number of active ranks for this operation.
 * \param[in]     active_rank Array of active ranks.
 * \param[in]     send_or_recv_idx Array of offsets for the data buffer.
 * \param[in]     send_or_recv_n Array of counts for each active rank.
 * \param[in]     tag The message tag to be used.
 * \param[in,out] buffer Pointer to the data buffer.
 *
 * \return An ID for the persistent request.
 *
 * \see PDM_exchange_helper_exch_start and PDM_exchange_helper_exch_wait.
 */
int
PDM_exchange_helper_exch_one_way_init
(
  PDM_exchange_helper_t    *exch_helper,
  PDM_exchange_direction_t  direction,
  size_t                    s_data,
  int                       cst_stride,
  int                       n_active_rank,
  int                      *active_rank,
  int                      *send_or_recv_idx,
  int                      *send_or_recv_n,
  int                       tag,
  void                     *buffer
);

/**
 * \brief Starts a previously initialized persistent communication request.
 *
 * This function activates a persistent request created with `..._init` functions.
 * The communication is initiated in a non-blocking manner.
 *
 * \param[in,out] exch_helper The initialized exchange helper object.
 * \param[in]     request The ID of the persistent request to be started.
 *
 * \see PDM_exchange_helper_exch_init to create the request.
 */
void
PDM_exchange_helper_exch_start
(
  PDM_exchange_helper_t *exch_helper,
  int                    request
);

/**
 * \brief Waits for a non-blocking communication request to complete.
 *
 * This function blocks until the non-blocking communication specified by
 * the request ID has finished. It is used to complete requests from
 * `PDM_exchange_helper_iexch` or `PDM_exchange_helper_exch_start`.
 *
 * \param[in,out] exch_helper The initialized exchange helper object.
 * \param[in]     request The ID of the non-blocking request to wait for.
 *
 * \see PDM_exchange_helper_iexch to initiate a non-blocking exchange.
 */
void
PDM_exchange_helper_exch_wait
(
  PDM_exchange_helper_t *exch_helper,
  int                    request
);

/**
 * \brief Checks the completion status of a multi-segment exchange request.
 *
 * This function queries the status of a specific exchange request which may
 * consist of multiple underlying MPI sub-requests. It mimics the behavior of
 * \ref MPI_Testall: it returns 1 only if all sub-communications associated
 * with the \p request_id have completed.
 *
 * If the status is not \ref EXCHANGE_HELPER_STATUS_ONGOING, it returns 1 immediately.
 * Otherwise, it iterates through all sub-requests to progress the communication
 * and verify their status.
 *
 * \param[in,out] exch_helper \ref PDM_exchange_helper_t structure managing the requests.
 * \param[in]     request_id  The ID of the composite request to check.
 *
 * \return An integer acting as a boolean flag:
 * - 1 if all sub-requests are completed (or if the request was already finished).
 * - 0 if at least one sub-request is still pending.
 *
 * \note Calling this function helps progress the MPI communication engine for
 * all pending sub-requests within the specified exchange.
 */
int
PDM_exchange_helper_exch_test
(
  PDM_exchange_helper_t *exch_helper,
  int                    request_id
);

/**
 * \brief Frees a persistent communication request.
 *
 * This function releases all resources associated with a specific
 * persistent request. It should be called after a request is no longer needed.
 *
 * \param[in,out] exch_helper The initialized exchange helper object.
 * \param[in]     request_id The ID of the persistent request to free.
 *
 * \pre The communication associated with request_id must be completed before freeing it.
 *
 * \see PDM_exchange_helper_exch_init to create the request.
 */
void
PDM_exchange_helper_exch_free
(
  PDM_exchange_helper_t *exch_helper,
  int                    request_id
);

/**
 * \brief Frees all memory associated with an exchange helper object.
 *
 * This is the destructor for the exchange helper. It frees all internal
 * memory, including any non-freed persistent requests. The object should
 * not be used after this call.
 *
 * \param[in,out] exch_helper A pointer to the exchange helper object to be freed.
 *
 */
void
PDM_exchange_helper_free
(
  PDM_exchange_helper_t *exch_helper
);

#ifdef  __cplusplus
}
#endif

#endif  /* __PDM_EXCHANGE_HELPER_H__ */
