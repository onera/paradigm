/*
 * \file
 */

#ifndef __PDM_PART_COMM_GRAPH_ALGORITHM_H__
#define __PDM_PART_COMM_GRAPH_ALGORITHM_H__

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

/*============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Compute link between entity2 from entity1 link. Useful in order to deduce graph of faces with graph of vertices, for example. High level API
 *
 * \param [in]  pcg_entity1          \ref PDM_part_comm_graph_t structure for entity1
 * \param [in]  pn_entity1           Number of entity1 (size = n_part)
 * \param [in]  pn_entity2           Number of entity1 (size = n_part)
 * \param [in]  entity2_entity1_idx  Connectivity index (size = \p pn_entity2 + 1)
 * \param [in]  entity2_entity1      Connectivity array (size = \p entity2_entity1_idx[\p pn_entity2] )
 * \param [out] pcg_entity1          \ref PDM_part_comm_graph_t structure for entity2
 *
 */
void
PDM_part_comm_graph_entity1_to_part_comm_graph_entity2
(
  PDM_part_comm_graph_t   *pcg_entity1,
  int                     *pn_entity1,
  int                     *pn_entity2,
  int                    **entity2_entity1_idx,
  int                    **entity2_entity1,
  PDM_part_comm_graph_t  **pcg_entity2
);


/**
 *
 * \brief Compute link between entity2 from entity1 link. Useful in order to deduce graph of faces with graph of vertices, for example. Low level API
 *
 * \param [in]  comm                 MPI communicator
 * \param [in]  n_part               Number of partition on current process
 * \param [in]  pn_entity1_graph     Number of bound (size = \p n_part)
 * \param [in]  pentity1_graph       Graph comm identifier (size = 4 * \p pn_entity1_graph[i_part]) :
 * \param [in]  nuplet_size          Nuplet size
 * \param [in]  pentity1_nuplet      Additional nuplets (NULL or size = \p nuplet_size * \p pn_entity_graph[i_part])
 * \param [in]  pn_entity1           Number of entity1 (size = \p n_part)
 * \param [in]  pn_entity2           Number of entity1 (size = \p n_part)
 * \param [in]  entity2_entity1_idx  Connectivity index (size = \p pn_entity2 + 1 )
 * \param [in]  entity2_entity1      Connectivity array (size = \p entity2_entity1_idx[\p pn_entity2] )
 * \param [out] out_pn_entity2_graph Number of bound (size = \p n_part)
 * \param [out] out_pentity2_graph   Graph comm identifier (size = 4 * \p pn_entity2_graph[i_part]) :
 * \param [out] out_pentity2_nuplet  Nuplets for entity2 (NULL or size = 4 * \p pn_entity2_graph[i_part]) :
 *
 */
void
PDM_part_comm_graph_entity1_to_entity2
(
  PDM_MPI_Comm             comm,
  int                      n_part,
  int                     *pn_entity1_graph,
  int                    **pentity1_graph,
  int                      nuplet_size,
  int                    **pentity1_nuplet,
  int                     *pn_entity1,
  int                     *pn_entity2,
  int                    **entity2_entity1_idx,
  int                    **entity2_entity1,
  int                    **out_pn_entity2_graph,
  int                   ***out_pentity2_graph,
  int                   ***out_pentity2_nuplet
);

/**
 * \brief Get selected of entities2 from list selected entities1
 *
 * \param [in]  n_selected_entity1    Number of selected entities1
 * \param [in]  selected_entity1      Local IDs of selected entities1 (1-based) (if NULL, assume all IDs from 1 to n_selected_entity1)
 * \param [in]  entity1_entity2_idx   Index for entity1->entity2 connectivity
 * \param [in]  entity1_entity2       Entity1->entity2 connectivity
 * \param [in]  pcg_entity2           \ref PDM_part_comm_graph_t instance for entities2
 * \param [out] out_n_entity2         Number of selected entities2
 * \param [out] out_selected_entity2  Local IDs of selected entities1 (1-based) (if NULL, assume all IDs from 1 to n_selected_entity1)
 *
 */
void
PDM_part_comm_graph_selected_entity1_to_selected_entity2
(
  int                     *n_selected_entity1,
  int                    **selected_entity1,
  int                    **entity1_entity2_idx,
  int                    **entity1_entity2,
  PDM_part_comm_graph_t   *pcg_entity2,
  int                    **out_n_entity2,
  int                   ***out_selected_entity2
);

/**
 * \brief Create an unique part_comm_graph from multiple ones.
 *
 * \param [in] comm  MPI communicator
 * \param [in] n_pcg Number of PDM_part_comm_graph_t objects to concatenate
 * \param [in] pcgs  PDM_part_comm_graph_t objects to concatenate
 *
 * \return Initialized concatenated \ref PDM_part_comm_graph_t instance
 *
 */
PDM_part_comm_graph_t *
PDM_part_comm_graph_concatenate
(
  PDM_MPI_Comm            comm,
  int                     n_pcg,
  PDM_part_comm_graph_t **pcgs
);

/**
 * \brief Split part_comm_graph into n_color part_comm_graph objects
 *
 * \param [in   ] pcg          Initial PDM_part_comm_graph_t object
 * \param [in   ] n_color      Number of colors
 * \param [in   ] entity_color List of n_part arrays of size n_entity_graph[i_part] each,
 *                             indicating the ID (0-based) of the output pcg
 * \param [inout] split_pcgs   Splitted part_comm_graph objects (size=n_color, must be allocated first)
 *
 */
void
PDM_part_comm_graph_split
(
  PDM_part_comm_graph_t  *pcg,
  const int               n_color,
  const int             **entity_color,
  PDM_part_comm_graph_t **split_pcgs
);

/**
 * \brief Remove entries in a part_comm_graph object
 *
 * \param [in] pcg   Initial PDM_part_comm_graph_t object
 * \param [in] flag  List of n_part arrays of size n_entity_graph[i_part] each,
 *                   indicating if entity is selected or not
 * \param [in] both  If True, keep entry if both sides are flagged True;
 *                   otherwise, keep entry if at least one side is flagged True.
 *
 * \return \ref PDM_part_comm_graph_t instance
 *
 */
PDM_part_comm_graph_t*
PDM_part_comm_graph_filter
(
  PDM_part_comm_graph_t   *pcg,
  const int              **flag,
  const int                both
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /*  __PDM_PART_COMM_GRAPH_ALGORITHM_H__ */
