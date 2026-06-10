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
