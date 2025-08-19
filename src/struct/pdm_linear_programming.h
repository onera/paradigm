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
 *                                 l <= x <= u
 *
 * \note The matrix \p a is defined in row-major (C) order, i.e. a_{i,j} = a[dim*i+j]
 *
 * \param [in]     dim   Dimension
 * \param [in]     n     Number of inequality constraints
 * \param [in]     a     a in ax <= b (size = \p n * \p dim)
 * \param [in]     b     b in ax <= b (size = \p n)
 * \param [in]     l     Lower bounds l <= x (size = \p dim)
 * \param [in]     u     Upper bounds x <= u (size = \p dim)
 * \param [in]     c     Constant in the objective function (size = \p dim)
 * \param [inout]  x     Initial point - Optimum (size = \p dim)
 *
 * \return Problem status
 */
PDM_lp_status_t
PDM_lp_solve_nd
(
  const int     dim,
  const int     n,
  const double *a,
  const double *b,
  const double *l,
  const double *u,
  const double *c,
        double *x
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_LINEAR_PROGRAMMING_H__ */
