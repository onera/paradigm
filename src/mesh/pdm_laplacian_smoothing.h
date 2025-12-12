#ifndef __PDM_LAPLACIAN_SMOOTHING_H__
#define __PDM_LAPLACIAN_SMOOTHING_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_part_comm_graph.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Types definition
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/**
 *
 * \brief Compute Inverse Distance Weighting-like edge weights
 *
 * \param [in]  n_part            Number of partition on current process
 * \param [in]  p_vtx_coord       Vertex coordinates (size = \p n_part, for each part size = 3 * \p p_n_vtx [i_part])
 * \param [in]  p_n_edge          Number of edges (size = \p n_part)
 * \param [in]  p_edge_vtx        Edge→vertex connectivity (size = \p n_part, for each part size = 2 * \p p_n_edge [i_part])
 * \param [in]  exponent          Edge length exponent
 * \param [out] out_p_edge_weight Edge weight (size = \p n_part, for each part size = \p p_n_edge [i_part])
 */
 void PDM_laplacian_smoothing_idw_weights_compute
(
  int       n_part,
  double  **p_vtx_coord,
  int      *p_n_edge,
  int     **p_edge_vtx,
  int       exponent,
  double ***out_p_edge_weight
);

/**
 *
 * \brief Compute Beltrami (cotangent) edge weights
 *
 * \param [in]  n_part            Number of partition on current process
 * \param [in]  p_n_vtx           Number of vertices (size = \p n_part)
 * \param [in]  p_vtx_coord       Vertex coordinates (size = \p n_part, for each part size = 3 * \p p_n_vtx [i_part])
 * \param [in]  p_n_elt           Number of elements (size = \p n_part)
 * \param [in]  p_elt_vtx_idx     Index of element→vertex connectivity (size = \p n_part, for each part size = \p p_n_elt [i_part]+1)
 * \param [in]  p_elt_vtx         Element→vertex connectivity (size = \p n_part, for each part size = \p p_elt_vtx_idx [i_part][\p p_n_elt [i_part]])
 * \param [in]  p_n_edge          Number of edges (size = \p n_part)
 * \param [in]  p_edge_vtx        Edge→vertex connectivity (size = \p n_part, for each part size = 2 * \p p_n_edge [i_part])
 * \param [out] out_p_edge_weight Edge weight (size = \p n_part, for each part size = \p p_n_edge [i_part])
 */
void PDM_laplacian_smoothing_beltrami_weights_compute
(
  int       n_part,
  int      *p_n_vtx,
  double  **p_vtx_coord,
  int      *p_n_elt,
  int     **p_elt_vtx_idx,
  int     **p_elt_vtx,
  int      *p_n_edge,
  int     **p_edge_vtx,
  double ***out_p_edge_weight
);

/**
 *
 * \brief Apply Laplacian smoothing to strided fields (interlaced).
 *        If no weight is provided, default unity weights are used.
 *
 * \param [in]     comm             MPI communicator
 * \param [in]     n_part           Number of partition on current process
 * \param [in]     p_n_vtx          Number of vertices (size = \p n_part)
 * \param [in]     p_n_vtx_frozen   Number of frozen vertices (size = \p n_part) or NULL
 * \param [in]     p_vtx_frozen     Local ID of frozen vertices (size = \p n_part, for each part size = \p p_n_vtx_frozen [i_part]) or NULL
 * \param [in]     pcg_vtx          Pointer to \ref PDM_part_comm_graph_t instance for vertices
 * \param [in]     p_n_edge         Number of edges (size = \p n_part)
 * \param [in]     p_edge_vtx       Edge→vertex connectivity (size = \p n_part, for each part size = 2 * \p p_n_edge [i_part])
 * \param [in]     p_edge_weight    Edge weight (size = \p n_part, for each part size = \p p_n_edge [i_part]) or NULL
 * \param [in]     pcg_edge         Pointer to \ref PDM_part_comm_graph_t instance for edges or NULL
 * \param [in]     damping          Damping constant (between 0. and 1.)
 * \param [in]     n_iter           Number of smoothing iterations
 * \param [in]     tol              Relative tolerance for convergence (ignored if negative)
 * \param [in]     stride           Field stride (interlaced values)
 * \param [in/out] p_vtx_field      Fields (size = \p n_part, for each part size = \p p_n_vtx [i_part])
 */
void
PDM_laplacian_smoothing_fields
(
  const PDM_MPI_Comm            comm,
        int                     n_part,
        int                    *p_n_vtx,
        int                    *p_n_vtx_frozen,
        int                   **p_vtx_frozen,
        PDM_part_comm_graph_t  *pcg_vtx,
        int                    *p_n_edge,
        int                   **p_edge_vtx,
        double                **p_edge_weight,
        PDM_part_comm_graph_t  *pcg_edge,
        double                  damping,
        int                     n_iter,
        double                  tol,
        int                     stride,
        double                **p_vtx_field
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_LAPLACIAN_SMOOTHING__ */
