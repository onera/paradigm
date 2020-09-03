#ifndef __PDM_PARA_GRAPH_DUAL_H__
#define __PDM_PARA_GRAPH_DUAL_H__

#include <stdio.h>
#include "pdm.h"
#include "pdm_mpi.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#if !defined (__hpux) && !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/
/**
 * \enum PDM_part_split_t
 * \brief Split method
 *
 */

typedef enum {
  PDM_SPLIT_DUAL_WITH_PARMETIS = 1,
  PDM_SPLIT_DUAL_WITH_PTSCOTCH = 2
} PDM_split_dual_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

void
PDM_para_graph_compress_connectivity
(
 PDM_g_num_t *dual_graph,
 PDM_g_num_t *dual_graph_idx,
 int         *dual_graph_n,
 int          dn_elt
);

/**
 *
 * \brief Compute the dual graph in parallel for a face cell connectivity
 *
 */
void
PDM_para_graph_dual_from_arc2node
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *graph_node_distrib,
 const PDM_g_num_t     *graph_arc_distrib,
 const PDM_g_num_t     *darc_to_node,
       PDM_g_num_t    **dual_graph_idx,
       PDM_g_num_t    **dual_graph,
 const int              compute_dnode_to_arc,
       int            **dnode_to_arc_idx,
       PDM_g_num_t    **dnode_to_arc
);


/**
 *
 * \brief Compute the dual graph in parallel for a cell face connectivity
 */
void
PDM_para_graph_dual_from_node2arc
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *graph_node_distrib,
 const PDM_g_num_t     *graph_arc_distrib,
 const int             *dnode_arc_idx,
 const PDM_g_num_t     *dnode_arc,
       PDM_g_num_t    **dual_graph_idx,
       PDM_g_num_t    **dual_graph

);

/**
 *
 * \brief Call the chosen graph partitioner to split the dual graph
*/
void
PDM_para_graph_split
(
 const PDM_split_dual_t  split_method,
 const PDM_g_num_t      *graph_node_distrib,
 const PDM_g_num_t      *dual_graph_idx,
 const PDM_g_num_t      *dual_graph,
 const int              *node_weight,
 const int              *arc_weight,
 const int               n_part,
 const double           *part_fraction,
 int                    *node_part_id,
 const PDM_MPI_Comm      comm
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_PARA_GRAPH_DUAL_H__ */
