/*
 * \file
 */

#ifndef __PDM_PART_COMM_GRAPH_H__
#define __PDM_PART_COMM_GRAPH_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


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
 * \param [in]   pn_entity_graph        Number of bound (size = n_part)
 * \param [in]   pentity_graph          Graph comm identifier (size = 4* pn_entity_graph[i_part]) :
                                            For each entity :
                                              - entity local number
                                              - Connected process
                                              - Connected partition on the connected process
                                              - Connected entity local number in the connected partition
 * \return   Initialized \ref PDM_part_comm_graph_t instance
 */
PDM_part_comm_graph_t*
PDM_part_comm_graph_create
(
  int            n_part,
  int           *pn_entity_graph, // Pas besoin, juste besoin de la taille en fait ?
  int          **pentity_graph,
  PDM_MPI_Comm   comm
);



/**
 *
 * \brief Exchange data between graph comm with synchronous blocking exchange
 * \param [in]   ptpgc               \ref PDM_part_comm_graph_t structure
 * \param [in]   s_data              Data size
 * \param [in]   t_stride            Kind of stride (see \ref PDM_stride_t )
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
 PDM_part_comm_graph_t   *ptpgc,
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
 * \brief Get the owner array compute inside the structure, usefull to manage reduction of array for exemple
 * \param [in]   ptpgc         \ref PDM_part_comm_graph_t structure
 * \param [in]   i_part        Id of current partition
 *
 * \return   Array of size pentity_graph[i_part] that contains 0 if not owner and 1 if owner. Ownership is determine by the lowest rank that hold the entity
 */
const int*
PDM_part_comm_graph_owner_get
(
 PDM_part_comm_graph_t *ptpgc,
 int                    i_part
);


/**
 *
 * \brief Reorder internaly all graph comm with the table old_to_new.
 *        This method is usefull when who want to change the local order of entity and update exchange protocol
 *        This method change the internal data for futur exchange and update with the new value of pentity_graph and old_to_new
 * \param [in]   ptpgc          \ref PDM_part_comm_graph_t structure
 * \param [in]   pentity_graph  Graph comm identifier (size = 4* pn_entity_graph[i_part]) :
 * \param [in]   old_to_new     Permutation id old to new
 */
void
PDM_part_comm_graph_reorder
(
  PDM_part_comm_graph_t  *ptpgc,
  int                   **pentity_graph,
  int                   **old_to_new
);


/**
 *
 * \brief Free \ref PDM_part_comm_graph_t structure
 *
 * \param ptpgc               \ref PDM_part_comm_graph_t structure
 *
 */
void
PDM_part_comm_graph_free
(
 PDM_part_comm_graph_t* ptpgc
);


void
PDM_part_comm_graph_entity1_to_entity2
(
  PDM_part_comm_graph_t  *ptpgc_entity1,
  int                    *pn_entity1,
  int                    *pn_entity2,
  int                   **entity2_entity1_idx,
  int                   **entity2_entity1
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /*  __PDM_PART_COMM_GRAPH_H__ */
