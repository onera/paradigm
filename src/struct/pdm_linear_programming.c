/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_error.h"
#include "pdm_linear_programming.h"
#include "pdm_logging.h"
#include "pdm_priv.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Static global variables
 *============================================================================*/

static const double LP_EPS  = 1.e-15;
static const int    verbose = 0;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/**
 * Solve the 1-dimensional linear programming problem:
 *  maximize c*x
 *  subject to constraints ai*x <= bi
 *                         l <= x <= u
 */
static PDM_lp_status_t
_lp_solve_1d
(
  const int     n,
  const double *a,
  const double *b,
  const double  l,
  const double  u,
  const double  c,
        double *x
)
{
  double low  = l;
  double high = u;

  if (verbose) {
    log_trace(">> _lp_solve_1d\n");
    PDM_log_trace_array_double(a, n, "  a : ");
    PDM_log_trace_array_double(b, n, "  b : ");
    log_trace("  l,u = %20.16f / %20.16f\n", l, u);
  }


  for (int i = 0; i < n; i++) {
    if (a[i] < 0) {
      low  = PDM_MAX(low,  b[i]/a[i]);
    }
    else if (a[i] > 0) {
      high = PDM_MIN(high, b[i]/a[i]);
    }
  }

  if (high < low) {
    return PDM_LP_UNFEASIBLE;
  }
  else {
    if (c < 0) {
      *x = low;
    }
    else {
      *x = high;
    }
    return PDM_LP_FEASIBLE;
  }
}


/** 
 * Solve the dim-dimensional linear programming problem:
 *  maximize c.x
 *  subject to constraints ai.x <= bi
 *                         l <= x <= u
 * 
 * Algorithm from R. Seidel, "Linear programming and convex hulls made easy", 1990
 */
static PDM_lp_status_t
_lp_solve_nd
(
  const int     dim,
  const int     n,
  const double *a,
  const double *b,
  const double *l,
  const double *u,
  const double *c,
        double *x
)
{
  if (dim == 1) {
    return _lp_solve_1d(n,
                        a,
                        b,
                        l[0],
                        u[0],
                        c[0],
                        x);
  }

  if (verbose) log_trace("\n\n>> _lp_solve_nd, dim = %d\n", dim);


  // Compute optimum solution x for the constraints given by l and u
  for (int i = 0; i < dim; i++) {
    if (c[i] < 0) {
      x[i] = l[i];
    }
    else {
      x[i] = u[i];
    }
  }

  if (verbose) PDM_log_trace_array_double(x, dim, "initial solution : ");

  // Set for sub-problem
  double sub_a[(n+2)*(dim-1)];
  double sub_b[(n+2)];
  double sub_l[dim-1];
  double sub_u[dim-1];
  double sub_c[dim-1];
  double sub_x[dim-1];

  // int used[n];
  // memset(used, 0, n*sizeof(int));
 
  for (int i_constraint = 0; i_constraint < n; i_constraint++) {

    // int i = -1;
    // do {
    //   i = rand() % n; // pick random constraint
    // } while (used[i]);
    // used[i] = 1;
    int i = i_constraint;

    // i-th constraint
    const double *ai = a + dim*i;
    double ax_b = -b[i];
    for (int j = 0; j < dim; j++) {
      ax_b += ai[j] * x[j];
    }

    if (verbose) {
      log_trace("constraint #%d\n", i);
      PDM_log_trace_array_double(x,  dim, "  x  : ");
      PDM_log_trace_array_double(ai, dim, "  ai : ");
      log_trace("  bi = %f\n", b[i]);
      log_trace("  ai.x - bi = %f\n", ax_b);
    }

    // case: i-th constraint is satisfied
    if (ax_b <= 0) {
      continue;
    }

    // case: i-th constraint not satisfied
    // find the largest component in absolute value
    int    k   = -1;
    double aik = 0.;
    for (int j = 0; j < dim; j++) {
      double aij = PDM_ABS(ai[j]);
      if (aij > aik) {
        k   = j;
        aik = aij;
      }
    }

    if (verbose) {
      if (k < 0) {
        log_trace("  k = %d\n", k);
      }
      else {
        log_trace("  k = %d, aik = %f\n", k, ai[k]);
      }
    }

    if (k < 0) {
      // null constraint
      if (b[i] < 0) {
        return PDM_LP_UNFEASIBLE;
      }
      else {
        continue;
      }
    }

    // project previous constraints to lower dimension by eliminating k-th decision variable
    double inv_aik = 1. / ai[k];

    int jj;
    for (int m = 0; m < i; m++) {
      const double *am     = a     +  dim   *m;
      double       *sub_am = sub_a + (dim-1)*m;
      jj = 0;
      for (int j = 0; j < dim; j++) {
        if (j != k) {
          sub_am[jj] = am[j] - ai[j]*am[k]*inv_aik;
          jj++;
        }
      }
      sub_b[m] = b[m] - b[i]*am[k]*inv_aik;
    }

    // incorporate the lk <= xk <= uk constraints into sub_a/b
    jj = 0;
    double *sub_al = sub_a + (dim-1)* i;
    double *sub_au = sub_a + (dim-1)*(i+1);
    for (int j = 0; j < dim; j++) {
      if (j != k) {
        sub_al[jj] =  ai[j]*inv_aik;
        sub_au[jj] = -ai[j]*inv_aik;
        jj++;
      }
    }
    sub_b[i  ] =        b[i]*inv_aik - l[k];
    sub_b[i+1] = u[k] - b[i]*inv_aik;


    // lower/upper bounds for remaining decision variables
    jj = 0;
    for (int j = 0; j < dim; j++) {
      if (j != k) {
        sub_l[jj] = l[j];
        sub_u[jj] = u[j];
        jj++;
      }
    }

    // project c to lower dimension
    jj = 0;
    double mag_sub_c = 0.;
    for (int j = 0; j < dim; j++) {
      if (j != k) {
        sub_c[jj] = c[j] - c[k]*ai[j]*inv_aik;
        mag_sub_c += sub_c[jj]*sub_c[jj];
        sub_x[jj] = x[j];
        jj++;
      }
    }

    if (mag_sub_c < LP_EPS) {
      if (verbose) log_trace("c is null !\n");
      // what do we do?
    }

    // solve sub-problem
    PDM_lp_status_t stat = _lp_solve_nd(dim-1,
                                        i+2,
                                        sub_a,
                                        sub_b,
                                        sub_l,
                                        sub_u,
                                        sub_c,
                                        sub_x);

    if (stat == PDM_LP_UNFEASIBLE) {
      return PDM_LP_UNFEASIBLE;
    }

    // lift solution back into dimension dim
    jj = 0;
    double aix = 0.;
    for (int j = 0; j < dim; j++) {
      if (j == k) continue;
      x[j] = sub_x[jj];
      aix += ai[j] * x[j];
      jj++;
    }

    x[k] = (b[i] - aix) * inv_aik;
    if (verbose) PDM_log_trace_array_double(x, dim, "  new x : ");

  } // End loop on constraints

  return PDM_LP_FEASIBLE;
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

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
)
{
  // Solve
  return _lp_solve_nd(dim,
                      n,
                      a,
                      b,
                      l,
                      u,
                      c,
                      x);
}



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
)
{
  // objective function
  double c[3] = {1., 1., 1.};

  // constraints
  double a[n_plane*3];
  double b[n_plane];
  a[ 0] = -1.; a[ 1] =  0.; a[ 2] =  0.;
  a[ 3] =  0.; a[ 4] = -1.; a[ 5] =  0.;
  a[ 6] =  0.; a[ 7] =  0.; a[ 8] = -1.;
  a[ 9] =  1.; a[10] =  0.; a[11] =  0.;
  a[12] =  0.; a[13] =  1.; a[14] =  0.;
  a[15] =  0.; a[16] =  0.; a[17] =  1.;
  for (int i = 0; i < 3; i++) {
    b[i]   = -box_extents[i];
    b[i+3] = box_extents[i+3];
  }
  for (int iplane = 0; iplane < n_plane; iplane++) {
    for (int i = 0; i < 3; i++) {
      a[3*iplane + i] = -plane_normal[3*iplane + i];
    }
    b[iplane] = -PDM_DOT_PRODUCT(plane_normal + 3*iplane, plane_origin + 3*iplane);
  }

  // upper/lower bounds
  double l[3] = {box_extents[0], box_extents[1], box_extents[2]};
  double u[3] = {box_extents[3], box_extents[4], box_extents[5]};


  // Solve LP problem
  double x[3];
  PDM_lp_status_t stat =  PDM_lp_solve_nd(3, n_plane, a, b, l, u, c, x);

  return (stat != PDM_LP_UNFEASIBLE);
}



void
PDM_lp_pts_inside_convex_hull
(
  const int     dim,
  const int     n_src,
  const double *src_coord,
  const int     n_tgt,
  const double *tgt_coord,
        int    *tgt_status
)
{
  /**
   * Based upon https://www.cs.mcgill.ca/~fukuda/soft/polyfaq/node22.html
   * 
   * We seek for an hyperplane separating the tgt point from the src point cloud, 
   * i.e. find (h_0, ..., h_dim) such that 
   *    dot((h_0, ..., h_{dim-1}), tgt_coord) - h_dim >= 0, and
   *    dot((h_0, ..., h_{dim-1}), s        ) - h_dim <= 0 for all s in src_coord
   * 
   * If such hyperplane does not exist then tgt_coord is inside the convex hull.
   */
  if (dim > 3) {
    PDM_error(__FILE__, __LINE__, 0, "Invalid dim %d (must be <= 3)\n", dim);
  }

  // Setup LP arrays
  double lp_constraint[(dim+1) * (n_src+1)];
  double lp_rhs       [          (n_src+1)];
  double lp_lower     [ dim+1             ];
  double lp_upper     [ dim+1             ];
  double lp_objective [ dim+1             ];
  double lp_solution  [ dim+1             ];

  // Source points axis-aligned bounding box (AABB)
  double src_aabb[2*dim];
  for (int j = 0; j < dim; j++) {
    src_aabb[    j] =  HUGE_VAL;
    src_aabb[dim+j] = -HUGE_VAL;
  }

  for (int i = 0; i < n_src; i++) {
    for (int j = 0; j < dim; j++) {
      src_aabb[    j] = PDM_MIN(src_aabb[    j], src_coord[3*i+j]);
      src_aabb[dim+j] = PDM_MIN(src_aabb[dim+j], src_coord[3*i+j]);
    }
  }

  // Target-point-independent constraints
  for (int i = 0; i < n_src; i++) {
    for (int j = 0; j < dim; j++) {
      lp_constraint[(dim+1)*i+j] = src_coord[3*i+j];
    }
    lp_constraint[(dim+1)*i+dim] = -1;
    lp_rhs[i] = 0;
  }

  // Upper and lower bound constraints
  for (int j = 0; j < dim+1; j++) {
    lp_lower[j] = -1e9;
    lp_upper[j] =  1e9;
  }

  for (int i_tgt = 0; i_tgt < n_tgt; i_tgt++) {

    // Quick point-in-AABB test
    for (int j = 0; j < dim; j++) {
      if (tgt_coord[3*i_tgt+j] < src_aabb[j] || tgt_coord[3*i_tgt+j] > src_aabb[dim+j]) {
        // the current target point is outside the source AABB, hence it is outside the convex hull
        tgt_status[i_tgt] = 0;
        continue;
      }
    }

    // Objective function
    for (int j = 0; j < dim; j++) {
      lp_objective[j] = tgt_coord[3*i_tgt+j];
    }
    lp_objective[dim] = -1;

    // Constraint = opposite sign w.r.t. src points
    for (int j = 0; j <= dim; j++) {
      lp_constraint[(dim+1)*n_src+j] = -lp_objective[j];
    }
    lp_rhs[n_src] = -1;

    // Solve LP problem
    PDM_lp_status_t stat = PDM_lp_solve_nd(dim+1,
                                              n_src+1,
                                              lp_constraint,
                                              lp_rhs,
                                              lp_lower,
                                              lp_upper,
                                              lp_objective,
                                              lp_solution);

    tgt_status[i_tgt] = (stat == PDM_LP_UNFEASIBLE);
  } // End loop on target points
}