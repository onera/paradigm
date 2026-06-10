#ifndef __PDM_HO_BASIS_H__
#define __PDM_HO_BASIS_H__

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_mesh_nodal.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/
/*----------------------------------------------------------------------------
 *
 * Callback to define the basis functions of an high order
 * element
 *
 * parameters:
 *   order             <-- element order
 *   n_nodes           <-- number of nodes of the element
 *   n_pts             <-- number of points
 *   uvw               <-- Parametric coordinates of points
 *   projected_uvw     --> Interpolation weights associated to uvw coordinates
 *
 *----------------------------------------------------------------------------*/

typedef
void (*PDM_ho_basis_fct_t)
(
  const int     entities_dim,
  const int     order,
  const int     n_nodes,
  const int     n_pts,
  const double *uvw,
  double       *weights
);

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Evaluate high-order basis functions
 *
 *
 * \param [in]  type      Element type structure
 * \param [in]  order     Element order
 * \param [in]  n_nodes   Number of nodes
 * \param [in]  n_pts     Number of points
 * \param [in]  uvw       Parametric coordinates of the points (size = elt_dim * \ref n_pts)
 * \param [out] weights   Weights (size = \ref n_pts * \ref n_nodes)
 *
 */
void
PDM_ho_basis
(
 const PDM_Mesh_nodal_elt_t  type,
 const int                   order,
 const int                   n_nodes,
 const int                   n_pts,
 const double               *uvw,
       double               *weights
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_HO_BASIS_H__ */
