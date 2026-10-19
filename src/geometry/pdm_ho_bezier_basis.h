#ifndef __PDM_HO_BEZIER_BASIS_H__
#define __PDM_HO_BEZIER_BASIS_H__

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_mesh_nodal.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/**
 *
 * \brief Evaluate high-order basis Bézier functions
 *
 *
 * \param [in]  type      Element type structure
 * \param [in]  order     Element order
 * \param [in]  n_pts     Number of points
 * \param [in]  uvw       Parametric coordinates of the points (size = elt_dim * \ref n_pts)
 * \param [out] weights   Weights (size = \ref n_pts * \ref n_nodes)
 *
 */

void
PDM_ho_bezier_basis
(
 const PDM_Mesh_nodal_elt_t  type,
 const int                   order,
 const int                   n_pts,
 const double               *uvw,
 double                     *weights
);

/**
 *
 * \brief Evaluate high-order basis Bézier functions derivatives
 *
 *
 * \param [in]  type      Element type structure
 * \param [in]  order     Element order
 * \param [in]  n_pts     Number of points
 * \param [in]  uvw       Parametric coordinates of the points (size = elt_dim * \ref n_pts)
 * \param [out] dw_du   Weights derivatives with respect to u
 * \param [out] dw_dv   Weights derivatives with respect to v
 * \param [out] dw_dw   Weights derivatives with respect to w
 *
 */

void
PDM_ho_bezier_basis_derivative
(
 const PDM_Mesh_nodal_elt_t  type,
 const int                   order,
 const int                   n_pts,
 const double               *uvw,
 double            *restrict dw_du,
 double            *restrict dw_dv,
 double            *restrict dw_dw
);



#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_HO_BEZIER_BASIS_H__ */
