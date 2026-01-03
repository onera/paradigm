/*
 * \file
 */

#ifndef __PDM_PART_GRAPH_DUAL_H__
#define __PDM_PART_GRAPH_DUAL_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stddef.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_part_comm_graph.h"

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
 * \brief Assemble a distributed dual graph in global numbering for partitioning
 * \param [in]    comm                   MPI communicator
 * \param [in]    n_part                 Number of partition on current process
 * \param [in]    n_node                 Number of nodes in each partition (size = \p n_part)
 * \param [in]    n_arc                  Number of arcs in each partition (size = \p n_part)
 * \param [in]    select_node            Filter for nodes to include in the graph (size = \p n_part, \p select_node[i_part] size = \p n_node[i_part]). If NULL, all nodes are used.
 * \param [in]    node_arc_idx           Index for node to arc adjacency (size = \p n_part, \p node_arc_idx[i_part] size = \p n_node[i_part] + 1)
 * \param [in]    node_arc               Adjacency from node to arc (size = \p n_part, \p node_arc[i_part] size = \p node_arc_idx[i_part][n_node])
 * \param [in]    arc_node_idx           Index for arc to node adjacency (size = \p n_part, \p arc_node_idx[i_part] size = \p n_arc[i_part] + 1)
 * \param [in]    arc_node               Adjacency from arc to node (size = \p n_part, \p arc_node[i_part] size = \p arc_node_idx[i_part][n_arc])
 * \param [in]    node_weight            Weights for nodes (multi-constraint support) (size = \p n_part, \p node_weight[i_part] size = n_con * \p n_node[i_part])
 * \param [in]    arc_weight             Weights for arcs (edge weights) (size = \p n_part, \p arc_weight[i_part] size = \p n_arc[i_part])
 * \param [in]    pcg_node               Communication graph for nodes. Required if \p pcg_arc is NULL.
 * \param [in]    pcg_arc                Communication graph for arcs. Required if \p pcg_node is NULL.
 * \param [out]   out_gnode_node_idx     Output CSR index for the global dual graph (size = total_n_node + 1)
 * \param [out]   out_gnode_node         Output CSR adjacency (global numbering)
 * \param [out]   out_garc_weight        Output edge weights for the dual graph
 * \param [out]   out_distrib_node       Global distribution of nodes across ranks (size = n_rank + 1)
 *
 * \details This function builds a dual graph where nodes are connected if they share an arc.
 * It resolves interfaces using \p pcg_node or \p pcg_arc to exchange global identifiers.
 * The resulting graph is sorted, unique, and formatted for partitioners like ParMETIS or PT-SCOTCH.
 */
void
PDM_part_assembly_dual_graph
(
  PDM_MPI_Comm            comm,
  int                     n_part,
  int                    *n_node,
  int                    *n_arc,
  int                   **select_node,
  int                   **node_arc_idx,
  int                   **node_arc,
  int                   **arc_node_idx,
  int                   **arc_node,
  int                   **node_weight,
  int                   **arc_weight,
  PDM_part_comm_graph_t  *pcg_node,
  PDM_part_comm_graph_t  *pcg_arc,
  int                    *out_n_tot_node,
  PDM_g_num_t           **out_gnode_node_idx,
  PDM_g_num_t           **out_gnode_node,
  int                   **out_garc_weight,
  PDM_g_num_t           **out_distrib_node
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /*  __PDM_PART_GRAPH_DUAL_H__ */
