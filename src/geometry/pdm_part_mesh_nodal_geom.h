#ifndef __PDM_PART_MESH_NODAL_GEOM_H__
#define __PDM_PART_MESH_NODAL_GEOM_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_part_mesh_nodal.h"


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


/**
 * \brief Compute *local* measure of dual volumes for simplex (tria/tetra only)
 *
 * \param [in]  dim             Number of partitions
 * \param [in]  n_part          Number of partitions
 * \param [in]  n_elt           Number of elements for each partitions
 * \param [in]  n_vtx           Vertex coordinates (size = 3 * n_vtx)
 * \param [in]  elt_vtx         Element vertices connectivities array (stride is 3 ou 4 implicitly)
 * \param [in]  vtx_coord       Vertex coordinates (size = 3 * n_vtx)
 * \param [out] out_vtx_volume  Entity center computed by mean with vertices coordinnates (size = n_vtx)
 */
void
PDM_compute_dual_volume_simplex
(
  int                 dim,
  int                 n_part,
  int                *n_elt,
  int                *n_vtx,
  int               **elt_vtx,
  double            **vtx_coord,
  double           ***out_vtx_volume
);


/**
 *
 * \brief Compute dual volumes
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  synchronize  Enable synchronization at partition boundaries
 * \param [in]  dual_vol     For each part, dual volume for each vertex
 *
 */
void
PDM_part_mesh_nodal_dual_volume_compute
(
  PDM_part_mesh_nodal_t   *pmn,
  PDM_bool_t               synchronize,
  double                ***dual_vol
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_MESH_NODAL_GEOM_H__ */
