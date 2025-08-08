/*
 * \file
 */

#ifndef __PDM_PART_COMM_GRAPH_H__
#define __PDM_PART_COMM_GRAPH_H__

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

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _pdm_part_comm_graph_t PDM_part_comm_graph_t;


/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Build a \ref PDM_part_comm_graph_t instance
 * \param [in]   n_part                 Number of partition on current process
 * \param [in]   pn_entity_graph        Number of bound (size = \p n_part)
 * \param [in]   pentity_graph          Graph comm identifier (size = 4 * \p pn_entity_graph[i_part]) :
                                            For each entity :
                                              - entity local number (1-based)
                                              - Connected process   (0-based)
                                              - Connected partition on the connected process (1-based)
                                              - Connected entity local number in the connected partition (1-based)
 * \param [in]   ownership              Ownership for \p pentity_graph
 * \param [in]   comm                   MPI communicator

 * \return   Initialized \ref PDM_part_comm_graph_t instance
 */
PDM_part_comm_graph_t*
PDM_part_comm_graph_create
(
  int               n_part,
  int              *pn_entity_graph,
  int             **pentity_graph,
  PDM_ownership_t   ownership,
  PDM_MPI_Comm      comm
);


/**
 *
 * \brief Build a \ref PDM_part_comm_graph_t instance using additional information represented as a n-uplet
 * \param [in]   n_part                 Number of partition on current process
 * \param [in]   pn_entity_graph        Number of bound (size = \p n_part)
 * \param [in]   pentity_graph          Graph comm identifier (size = 4 * \p pn_entity_graph[i_part]) :
                                            For each entity :
                                              - entity local number (1-based)
                                              - Connected process   (0-based)
                                              - Connected partition on the connected process (1-based)
                                              - Connected entity local number in the connected partition (1-based)
 * \param [in]   ownership_graph        Ownership for \p pentity_graph
 * \param [in]   nuplet_size            N-uplet size
 * \param [in]   pentity_nuplet         Additional nuplets (size = \p nuplet_size * \p pn_entity_graph[i_part])
 * \param [in]   ownership_nuplet       Ownership for \p pentity_nuplet
 * \param [in]   is_signed              Use signed nuplets
 * \param [in]   comm                   MPI communicator
 *
 * \return   Initialized \ref PDM_part_comm_graph_t instance
 */
PDM_part_comm_graph_t*
PDM_part_comm_graph_with_nuplet_create
(
  int               n_part,
  int              *pn_entity_graph,
  int             **pentity_graph,
  PDM_ownership_t   ownership_graph,
  int               nuplet_size,
  int             **pentity_nuplet,
  PDM_ownership_t   ownership_nuplet,
  PDM_bool_t        is_signed,
  PDM_MPI_Comm      comm
);


/**
 *
 * \brief Exchange data between graph comm with synchronous blocking exchange
 * \param [in]   pcg                 \ref PDM_part_comm_graph_t structure
 * \param [in]   s_data              Data size
 * \param [in]   t_stride            Kind of stride (see \ref PDM_stride_t)
 * \param [in]   cst_stride          Constant stride
 * \param [in]   send_entity_stride  Stride of send data (following pentity_graph)
 * \param [in]   send_entity_data    Send data           (following pentity_graph)
 * \param [out]  recv_entity_stride  Stride of recv data (following pentity_graph)
 * \param [out]  recv_entity_data    Recv data           (following pentity_graph)
 *
 */
void
PDM_part_comm_graph_exch
(
  PDM_part_comm_graph_t   *pcg,
  size_t                   s_data,
  PDM_stride_t             t_stride,
  int                      cst_stride,
  int                    **send_entity_stride,
  void                   **send_entity_data,
  int                   ***recv_entity_stride,
  void                  ***recv_entity_data
);


/**
 *
 * \brief Start exchange data between graph comm with asyncrhonous exchange,
 *        after this call you need to use PDM_part_comm_graph_exch_start and wait.
 *        When you finish, you need to free the persistent exchange with PDM_part_comm_graph_exch_free
 * \param [in]   pcg                 \ref PDM_part_comm_graph_t structure
 * \param [in]   kcomm               Kind of MPI communication
 * \param [in]   s_data              Data size
 * \param [in]   t_stride            Kind of stride (see \ref PDM_stride_t)
 * \param [in]   cst_stride          Constant stride
 * \param [in]   send_entity_stride  Stride of send data (following pentity_graph)
 * \param [in]   send_entity_data    Send data           (following pentity_graph)
 * \param [out]  recv_entity_stride  Stride of recv data (following pentity_graph)
 * \param [out]  recv_entity_data    Recv data           (following pentity_graph)
 *
 * \return Request id
 *
 */
int
PDM_part_comm_graph_iexch
(
 PDM_part_comm_graph_t   *pcg,
 PDM_mpi_comm_kind_t      kcomm,
 size_t                   s_data,
 PDM_stride_t             t_stride,
 int                      cst_stride,
 int                    **send_entity_stride,
 void                   **send_entity_data,
 int                   ***recv_entity_stride,
 void                  ***recv_entity_data
);

/**
 *
 * \brief Prepare exchange data between graph comm with persistent exchange,
 *        after this call you need to use PDM_part_comm_graph_exch_start and wait.
 *        When you finish, you need to free the persistent exchange with PDM_part_comm_graph_exch_free
 * \param [in]   pcg                 \ref PDM_part_comm_graph_t structure
 * \param [in]   kcomm               Kind of MPI communication
 * \param [in]   s_data              Data size
 * \param [in]   t_stride            Kind of stride (see \ref PDM_stride_t)
 * \param [in]   cst_stride          Constant stride
 * \param [in]   send_entity_stride  Stride of send data (following pentity_graph)
 * \param [in]   send_entity_data    Send data           (following pentity_graph)
 * \param [out]  recv_entity_stride  Stride of recv data (following pentity_graph)
 * \param [out]  recv_entity_data    Recv data           (following pentity_graph)
 *
 * \return Request id
 *
 */
int
PDM_part_comm_graph_exch_init
(
 PDM_part_comm_graph_t   *pcg,
 PDM_mpi_comm_kind_t      kcomm,
 size_t                   s_data,
 PDM_stride_t             t_stride,
 int                      cst_stride,
 int                    **send_entity_stride,
 void                   **send_entity_data,
 int                   ***recv_entity_stride,
 void                  ***recv_entity_data
);

/**
 *
 * \brief Start exchange ( initalize by PDM_part_comm_graph_exch_init )
 * \param [in]   pcg                 \ref PDM_part_comm_graph_t structure
 * \param [in]   request_id          Request id
 *
 */
void
PDM_part_comm_graph_exch_start
(
  PDM_part_comm_graph_t   *pcg,
  int                      request_id
);


/**
 *
 * \brief Wait exchange ( initalize by PDM_part_comm_graph_exch_init and launch by PDM_part_comm_graph_exch_start )
 * \param [in]   pcg                 \ref PDM_part_comm_graph_t structure
 * \param [in]   request_id          Request id
 *
 */
void
PDM_part_comm_graph_exch_wait
(
  PDM_part_comm_graph_t   *pcg,
  int                      request_id
);


/**
 *
 * \brief Free persistent exchange ( initalize by PDM_part_comm_graph_exch_init )
 * \param [in]   pcg                 \ref PDM_part_comm_graph_t structure
 * \param [in]   request_id          Request id
 *
 */
void
PDM_part_comm_graph_exch_free
(
  PDM_part_comm_graph_t   *pcg,
  int                      request_id
);


/**
 *
 * \brief Get internal indirection to fill send buffer throw MPI from user data layout
 * \param [in]   pcg                 \ref PDM_part_comm_graph_t structure
 * \param [in]   part_to_send_buffer Indirection table to fill directly send buffer
 *
 */
void
PDM_part_comm_graph_part_to_send_buffer_get
(
  PDM_part_comm_graph_t   *pcg,
  int                   ***out_part_to_send_buffer
);

/**
 *
 * \brief Get internal indirection to fill recv buffer throw MPI from user data layout
 * \param [in]   pcg                 \ref PDM_part_comm_graph_t structure
 * \param [in]   part_to_recv_buffer Indirection table to fill directly recv buffer
 *
 */
void
PDM_part_comm_graph_part_to_recv_buffer_get
(
  PDM_part_comm_graph_t   *pcg,
  int                   ***part_to_recv_buffer
);

int
PDM_part_comm_graph_exch_one_way_raw_init
(
 PDM_part_comm_graph_t      *pcg,
 PDM_exchange_direction_t    direction,
 size_t                      s_data,
 int                         cst_stride,
 int                        *raw_buffer,
 int                         tag
);

void
PDM_part_comm_graph_exch_one_way_raw_start
(
 PDM_part_comm_graph_t      *pcg,
 int                         request_id
);

void
PDM_part_comm_graph_exch_one_way_raw_wait
(
  PDM_part_comm_graph_t   *pcg,
  int                      request_id
);

void
PDM_part_comm_graph_exch_one_way_raw_free
(
 PDM_part_comm_graph_t      *pcg,
 int                         request_id
);

/**
 *
 * \brief Get the owner array computed inside the structure, useful to manage reduction of array for example
 * \param [in]   pcg           \ref PDM_part_comm_graph_t structure
 * \param [in]   i_part        Id of current partition
 *
 * \return   Array of size pentity_graph[i_part] that contains 0 if not owner and 1 if owner. Ownership is determined by the lowest rank that holds the entity
 */
const int*
PDM_part_comm_graph_owner_get
(
  PDM_part_comm_graph_t *pcg,
  int                    i_part
);


/**
 *
 * \brief Reorder internally the graph with the table \p old_to_new.
 *        This method is useful when we want to change the local order of entities and update the exchange protocol.
 *        This method changes the internal data for future exchanges.
 * \param [in]   pcg            \ref PDM_part_comm_graph_t structure
 * \param [in]   old_to_new     Permutation id old to new (0-based)
 */
void
PDM_part_comm_graph_reorder
(
  PDM_part_comm_graph_t  *pcg,
  int                   **old_to_new
);


/**
 *
 * \brief Gather local and distant data through part_comm_graph communicator
 *
 * \param [in]  pcg             \ref PDM_part_comm_graph_t structure
 * \param [in]  size_data       Data size
 * \param [in]  t_stride        Kind of stride (see \ref PDM_stride_t)
 * \param [in]  n_entity        n_entity for data (should be > pn_entity_bound)
 * \param [in]  data_stride     Index of data to gather (size[i_part] = n_entity[i_part])
 * \param [in]  data            Data to gather (size[i_part] = cumsum(data_stride[i_part]))
 * \param [out] out_data_stride Gathered data index (size[i_part] = n_entity[i_part]+1)
 * \param [out] out_data        Gathered data (size[i_part] = cumsum(out_data_stride[i_part])))
 *
 */
void
PDM_part_comm_graph_gather_strided_data
(
  PDM_part_comm_graph_t   *pcg,
  const size_t             size_data,
  PDM_stride_t             t_stride,
  int                     *n_entity,
  int                    **data_stride,
  void                   **data,
  int                   ***out_data_stride,
  void                  ***out_data
);


/**
 *
 * \brief Free \ref PDM_part_comm_graph_t structure
 *
 * \param pcg               \ref PDM_part_comm_graph_t structure
 *
 */
void
PDM_part_comm_graph_free
(
  PDM_part_comm_graph_t* pcg
);

/**
 *
 * \brief Get entity graph
 *
 * \param [in]  pcg           Pointer to \ref PDM_part_comm_graph_t instance
 * \param [in]  i_part        Partition identifier
 * \param [out] entity_graph  Entity graph (size = 4 * n_entity_graph)
 * \param [in]  ownership     Ownership for \p entity_graph
 *
 * \return Number of entities in graph in current partition
 */
int
PDM_part_comm_graph_entity_graph_get
(
  PDM_part_comm_graph_t  *pcg,
  int                     i_part,
  int                   **entity_graph,
  PDM_ownership_t         ownership
);

/**
 *
 * \brief Inplace reduce value on current graph. Allow synchronisation.
 *        Only PDM_MPI_DOUBLE and PDM_MPI_INT are allowed
 *
 * \param [in]    pcg            Pointer to \ref PDM_part_comm_graph_t instance
 * \param [in]    datatype       Mpi datatype (PDM_MPI_DOUBLE/PDM_MPI_INT)
 * \param [in]    stride         Constant stride
 * \param [in]    op             Reduction operation kind (SUM/MIN/MAX)
 * \param [inout] pdata          Buffer of data to synchronise (size = n_entity)
 */
void
PDM_part_comm_graph_all_reduce
(
  PDM_part_comm_graph_t   *pcg,
  PDM_MPI_Datatype         datatype,
  int                      stride,
  PDM_MPI_Op               op,
  unsigned char          **pdata
);

/**
 *
 * \brief Get entity nuplets
 *
 * \param [in]  pcg            Pointer to \ref PDM_part_comm_graph_t instance
 * \param [in]  i_part         Partition identifier
 * \param [out] entity_nuplet  Entity nuplets (size = nuplet_size * n_entity_graph)
 * \param [in]  ownership      Ownership for \p entity_nuplet
 *
 * \return Size of nuplet
 */
int
PDM_part_comm_graph_entity_nuplet_get
(
  PDM_part_comm_graph_t  *pcg,
  int                     i_part,
  int                   **entity_nuplet,
  PDM_ownership_t         ownership
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /*  __PDM_PART_COMM_GRAPH_H__ */
