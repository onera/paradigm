/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_convex.h"
#include "pdm_error.h"
#include "pdm_linear_programming.h"
#include "pdm_logging.h"
#include "pdm_priv.h"


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

void
PDM_points_inside_convex_hull
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
    PDM_error("Invalid dim %d (must be <= 3)", dim);
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
      src_aabb[dim+j] = PDM_MAX(src_aabb[dim+j], src_coord[3*i+j]);
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



int
PDM_intersect_convex_volume_box
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
  PDM_lp_status_t stat = PDM_lp_solve_nd(3, n_plane, a, b, l, u, c, x);

  return (stat != PDM_LP_UNFEASIBLE);
}
