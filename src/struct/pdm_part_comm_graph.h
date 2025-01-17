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
 * \brief Reorder internally all comm graph with the table \p old_to_new.
 *        This method is useful when who want to change the local order of entity and update exchange protocol
 *        This method change the internal data for future exchange and update with the new value of \p pentity_graph and \p old_to_new
 * \param [in]   pcg            \ref PDM_part_comm_graph_t structure
 * \param [in]   pentity_graph  Comm graph identifier (size = 4 * \p pn_entity_graph[i_part])
 * \param [in]   old_to_new     Permutation id old to new (0-based)
 */
void
PDM_part_comm_graph_reorder
(
  PDM_part_comm_graph_t  *pcg,
  int                   **pentity_graph,
  int                   **old_to_new
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
 * \brief Compute link between entity2 from entity1 link. Useful in order to deduce graph of faces with graph of vertices, for example. High level API
 *
 * \param [in]  ptpgc_entity1        \ref PDM_part_comm_graph_t structure for entity1
 * \param [in]  pn_entity1           Number of entity1 (size = n_part)
 * \param [in]  pn_entity2           Number of entity1 (size = n_part)
 * \param [in]  entity2_entity1_idx  Connectivity index (size = \p pn_entity2 + 1)
 * \param [in]  entity2_entity1      Connectivity array (size = \p entity2_entity1_idx[\p pn_entity2] )
 * \param [out] ptpgc_entity1        \ref PDM_part_comm_graph_t structure for entity2
 *
 */
void
PDM_part_comm_graph_entity1_to_part_comm_graph_entity2
(
  PDM_part_comm_graph_t   *ptpgc_entity1,
  int                     *pn_entity1,
  int                     *pn_entity2,
  int                    **entity2_entity1_idx,
  int                    **entity2_entity1,
  PDM_part_comm_graph_t  **ptpgc_entity2
);


/**
 *
 * \brief Compute link between entity2 from entity1 link. Useful in order to deduce graph of faces with graph of vertices, for example. Low level API
 *
 * \param [in]  comm                 MPI communicator
 * \param [in]  n_part               Number of partition on current process
 * \param [in]  pn_entity1_graph     Number of bound (size = \p n_part)
 * \param [in]  pentity1_graph       Graph comm identifier (size = 4 * \p pn_entity1_graph[i_part]) :
 * \param [in]  pn_entity1           Number of entity1 (size = \p n_part)
 * \param [in]  pn_entity2           Number of entity1 (size = \p n_part)
 * \param [in]  entity2_entity1_idx  Connectivity index (size = \p pn_entity2 + 1 )
 * \param [in]  entity2_entity1      Connectivity array (size = \p entity2_entity1_idx[\p pn_entity2] )
 * \param [out] pn_entity2_graph     Number of bound (size = \p n_part)
 * \param [out] pentity2_graph       Graph comm identifier (size = 4 * \p pn_entity2_graph[i_part]) :
 *
 */
void
PDM_part_comm_graph_entity1_to_entity2
(
  PDM_MPI_Comm             comm,
  int                      n_part,
  int                     *pn_entity1_graph,
  int                    **pentity1_graph,
  int                     *pn_entity1,
  int                     *pn_entity2,
  int                    **entity2_entity1_idx,
  int                    **entity2_entity1,
  int                    **out_pn_entity2_graph,
  int                   ***out_pentity2_graph
);


/**
 *
 * \brief Get entity graph
 *
 * \param [in]  pcg           Pointer to \ref PDM_part_comm_graph_t instance
 * \param [in]  i_part        Partition identifier
 * \param [out] entity_graph  Entity graph (size = 4 * n_entity_graph)
 *
 * \return Number of entities in graph in current partition
 */
int
PDM_part_comm_graph_entity_graph_get
(
  PDM_part_comm_graph_t  *pcg,
  int                     i_part,
  int                   **entity_graph
);


/**
 *
 * \brief Get entity nuplets
 *
 * \param [in]  pcg            Pointer to \ref PDM_part_comm_graph_t instance
 * \param [in]  i_part         Partition identifier
 * \param [out] entity_nuplet  Entity nuplets (size = nuplet_size * n_entity_graph)
 *
 * \return Size of nuplet
 */
int
PDM_part_comm_graph_entity_nuplet_get
(
  PDM_part_comm_graph_t  *pcg,
  int                     i_part,
  int                   **entity_nuplet
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /*  __PDM_PART_COMM_GRAPH_H__ */
