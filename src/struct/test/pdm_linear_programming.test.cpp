
#include <stdlib.h>
#include <vector>
#include "doctest/doctest.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_mem_tool.h"
#include "pdm_linear_programming.h"
#include "pdm_logging.h"
#include "pdm_priv.h"

static double eps = 1e-6;
static double big = 1e6;


TEST_CASE("[pdm_linear_programming] - PDM_lp_solve_nd 2d") {

  /**
   * Maximize 15x + 10y
   * subject to constraints:
   *   x     >= 0
   *       y >= 0
   *   x + y <= 6
   *   x - y <= 4
   */

  const int dim = 2; // spatial dimension
  const int n   = 2; // number of constraints

  // Objective function
  double c[dim] = {15, 10};

  // Constraints
  double a[dim*n] = {
     0.25, 1,
     1.25, 0.5
  };

  double b[n] = {
    65,
    90
  };

  // Lower bounds
  double l[dim] = {
    0,
    0
  };

  // Upper bounds
  double u[dim] = {
    big,
    big
  };

  // Solve
  double x[dim];
  PDM_lp_status_t stat = PDM_lp_solve_nd(dim, n, a, b, l, u, c, x);

  // printf("[pdm_linear_programming] - PDM_lp_solve_nd 2d\n");
  // printf("  stat = %d\n", stat);
  // printf("  sol  = %20.16f %20.16f\n", x[0], x[1]);

  // Check
  double expected_sol[dim] = {
    51.11111111111111,
    52.22222222222222
  };
  CHECK(stat == PDM_LP_FEASIBLE);
  CHECK_EQ_C_ARRAY_FLOAT(x, expected_sol, dim, eps);
}



TEST_CASE("[pdm_linear_programming] - PDM_lp_solve_nd 3d") {

  /**
   * Maximize 20x + 10y + 15z
   * subject to constraints:
   *    x           >=  0
   *         y      >=  0
   *              z >=  0
   *   3x + 2y + 5z <= 55
   *   2x +  y +  z <= 26
   *    x +  y + 3z <= 30
   *   5x + 2y + 4z <= 57
   */

  const int dim = 3; // spatial dimension
  const int n   = 4; // number of constraints

  // Objective function
  double c[dim] = {20, 10, 15};

  // Constraints
  double a[dim*n] = {
    3,  2,  5,
    2,  1,  1,
    1,  1,  3,
    5,  2,  4
  };

  double b[n] = {
    55,
    26,
    30,
    57
  };

  
  // Lower bounds
  double l[dim] = {
    0,
    0,
    0
  };
  
  // Upper bounds
  double u[dim] = {
    big,
    big,
    big
  };


  // Solve
  double x[dim];
  PDM_lp_status_t stat = PDM_lp_solve_nd(dim, n, a, b, l, u, c, x);

  // printf("[pdm_linear_programming] - PDM_lp_solve_nd 3d\n");
  // printf("  stat = %d\n", stat);
  // printf("  sol  = %f %f %f\n", x[0], x[1], x[2]);

  // Check
  double expected_sol[dim] = {1.8, 20.8, 1.6};
  CHECK(stat == PDM_LP_FEASIBLE);
  CHECK_EQ_C_ARRAY_FLOAT(x, expected_sol, dim, eps);
}



TEST_CASE("[pdm_linear_programming] - PDM_lp_pts_inside_convex_hull 2d") {

  const int dim   = 2;
  const int n_src = 8;
  const int n_tgt = 3;

  double src_coord[n_src*3] = {
     0.0,  0.0, 0,
    12.0,  2.0, 0,
     4.0,  4.0, 0,
    14.0,  8.0, 0,
    10.0,  6.0, 0,
     4.0, 10.0, 0,
     8.0,  2.0, 0,
     0.0,  6.0, 0
  };


  double tgt_coord[n_tgt*3] = {
     2.0,  2.0, 0, // inside convex hull
    18.0,  2.0, 0, // outside AABB
    10.0, 10.0, 0  // inside  AABB but outside convex hull
  };

  int tgt_status[n_tgt];
  PDM_lp_pts_inside_convex_hull(dim,
                                n_src,
                                src_coord,
                                n_tgt,
                                tgt_coord,
                                tgt_status);

  int expected_tgt_status[n_tgt] = {1, 0, 0};
  CHECK_EQ_C_ARRAY(tgt_status, expected_tgt_status, n_tgt);
}



TEST_CASE("[pdm_linear_programming] - PDM_lp_pts_inside_convex_hull 3d") {

  const int dim   = 3;
  const int n_src = 12;
  const int n_tgt = 3;

  double src_coord[n_src*3] = {
     2.0, 4.0, -1.0,
     1.0, 2.0,  2.0,
    -1.0, 0.0,  1.0,
     1.0, 1.0,  1.0,
     0.0, 1.0,  0.0,
    -1.0, 3.0,  0.0,
    -1.0, 1.0,  2.0,
     3.0, 1.0, -1.0,
    -1.0, 1.0, -1.0,
     2.0, 2.0,  0.0,
     1.0, 2.0, -1.0,
     1.0, 1.0,  0.0
  };


  double tgt_coord[n_tgt*3] = {
    0.0, 1.0, 0.0, // inside convex hull
    4.0, 1.0, 1.0, // outside AABB
    1.0, 3.0, 1.0  // inside  AABB but outside convex hull
  };

  int tgt_status[n_tgt];
  PDM_lp_pts_inside_convex_hull(dim,
                                n_src,
                                src_coord,
                                n_tgt,
                                tgt_coord,
                                tgt_status);

  int expected_tgt_status[n_tgt] = {1, 0, 0};
  CHECK_EQ_C_ARRAY(tgt_status, expected_tgt_status, n_tgt);
}