/*
 * \file
 */

#ifndef __PDM_LINEAR_PROGRAMMING_H__
#define __PDM_LINEAR_PROGRAMMING_H__

/*============================================================================
 * Search octrees and quadtrees of boxes.
 *============================================================================*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


/*============================================================================
 * Type
 *============================================================================*/

typedef enum {

  PDM_LP_FEASIBLE   = 0,
  PDM_LP_UNFEASIBLE = 1,
  PDM_LP_UNBOUNDED  = 2

} PDM_lp_status_t;

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Solve the d-dimensional linear optimization problem
 *          maximize c.x
 *          subject to constraints ai.x <= bi
 *
 * \param [in]     dim   Dimension
 * \param [in]     n     Number of inequality constraints
 * \param [in]     a     a in ax <= b (size = \p n * \p dim)
 * \param [in]     b     b in ax <= b (size = \p n)
 * \param [in]     c     Constant in the objective function (size = \p dim)
 * \param [inout]  x     Initial point - Optimum (size = \p dim)
 *
 * \return Problem status
 */

PDM_lp_status_t
PDM_lp_solve_nd
(
 const int  dim,
 const int  n,
 double    *a,
 double    *b,
 double    *c,
 double    *x
 );

/**
 *
 * \brief Determine if the current box intersects a given volume
 *
 * \param [in]   n_plane          Number of planes in the current volume
 * \param [in]   plane_origin     Coordinates of a point on each plane
 * \param [in]   plane_normal     Normal vector of each plane
 * \param [in]   box_extents      Extents of the box (x_min, y_min, z_min, x_max, y_max, z_max)
 *
 */

int
PDM_lp_intersect_volume_box
(
 const int  n_plane,
 double    *plane_origin,
 double    *plane_normal,
 double    *box_extents
);


/**
 *
 * \brief Classify target points with respect to the convex hull of given source points.
 *
 * \warning Coordinates must always be defined in dimension 3, even if \p dim is lower
 *
 * \param [in]    dim         Spatial dimension (<= 3)
 * \param [in]    n_src       Number of source points
 * \param [in]    src_coord   Coordinates of source points (size = 3 * \p n_src)
 * \param [in]    n_tgt       Number of target points
 * \param [in]    tgt_coord   Coordinates of target points (size = 3 * \p n_tgt)
 * \param [inout] tgt_status  Status of each target point (0 = outside, 1 = inside) (size = \p n_tgt)
 *
 */
void
PDM_lp_pts_inside_convex_hull
(
  const int     dim,
  const int     n_src,
  const double *src_coord,
  const int     n_tgt,
  const double *tgt_coord,
        int    *tgt_status
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_LINEAR_PROGRAMMING_H__ */
