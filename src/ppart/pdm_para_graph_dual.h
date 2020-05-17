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

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

void
PDM_compress_connectivity
(
 PDM_g_num_t *dual_graph,
 int         *dual_graph_idx,
 int         *dual_graph_n,
 int          dn_elt
);

/**
 *
 * \brief Compute the dual graph in parallel for a face cell connectivity
 *
 */
// void
// PDM_para_graph_dual_from_face_cell
// (
//  const PDM_MPI_Comm     comm,
//        PDM_g_num_t     *cell_distribution,
//        PDM_g_num_t     *face_distribution,
//        PDM_g_num_t     *dface_cell,
//        PDM_g_num_t    **dual_graph,
//        int            **dual_graph_idx
// );

/**
 *
 * \brief Compute the dual graph in parallel for a face cell connectivity
 *
 */
void
PDM_para_graph_dual_from_face_cell
(
 const PDM_MPI_Comm     comm,
 const int              dn_cell,
 const int              dn_face,
       PDM_g_num_t     *dface_cell,
       int            **dual_graph_idx,
       PDM_g_num_t    **dual_graph
);


/**
 *
 * \brief Compute the dual graph in parallel for a cell face connectivity
 */
void
PDM_para_graph_dual_from_cell_face
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *cell_distribution,
 const PDM_g_num_t     *face_distribution,
 const PDM_g_num_t     *dcell_face,
       int            **dual_graph_idx,
       PDM_g_num_t    **dual_graph

);


void
PDM_split_graph
(
 const PDM_MPI_Comm  comm,
 int                *dual_graph_idx,
 PDM_g_num_t        *dual_graph,
 int                *delmt_weight,
 int                *cell_part,
 int                 dn_elmt,
 int                 n_part
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_PARA_GRAPH_DUAL_H__ */
