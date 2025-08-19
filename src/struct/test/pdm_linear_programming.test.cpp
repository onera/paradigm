
#include <stdlib.h>
#include <vector>
#include "doctest/doctest.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_linear_programming.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include "pdm_priv.h"

static double eps = 1e-6;
static double big = 1e6;


TEST_CASE("[pdm_linear_programming] - PDM_lp_solve_nd 2d") {

  const int dim   = 2;
  const int n_max = 3;
  PDM_lp_status_t expected_stat;

  int    n = 0;
  double a         [dim*n_max];
  double b         [    n_max];
  double c         [dim      ];
  double expected_x[dim      ];

  SUBCASE("Feasible") {
    /**
     * Maximize 15*x0 + 10*x1
     * subject to constraints:
     * 0.25*x0 +     x1 <= 65
     * 1.25*x0 - 0.5*x1 <= 90
     *      x0          >= 0
     *               x1 >= 0
     */

    // Objective function
    c[0] = 15; c[1] = 10;

    // Constraints
    n = 2;
    a[0] = 0.25; a[1] = 1.0; b[0] = 65;
    a[2] = 1.25; a[3] = 0.5; b[1] = 90;

    // Expected result
    expected_stat = PDM_LP_FEASIBLE;
    expected_x[0] = 51.11111111111111; expected_x[1] = 52.22222222222222;
  }

  SUBCASE("Unfeasible") {
    /**
     * Maximize x0 + 2*x1
     * subject to constraints:
     *  3*x0 + 2*x1 <=  6
     *  2*x0 + 5*x1 <=  10
     * -4*x0 + 3*x1 <= -10
     *    x0        >=  0
     *           x1 >=  0
     */

    // Objective function
    c[0] = 1; c[1] = 2;

    // Constraints
    n = 3;
    a[0] =  3; a[1] = 2; b[0] =  6;
    a[2] =  2; a[3] = 5; b[1] =  10;
    a[4] = -4; a[5] = 3; b[2] = -10;

    // Expected result
    expected_stat = PDM_LP_UNFEASIBLE;
  }

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
  CHECK(stat == expected_stat);
  if (stat == PDM_LP_FEASIBLE) {
    CHECK_EQ_C_ARRAY_FLOAT(x, expected_x, dim, eps);
  }
}



TEST_CASE("[pdm_linear_programming] - PDM_lp_solve_nd 3d") {

  const int dim   = 3;
  const int n_max = 4;
  PDM_lp_status_t expected_stat;

  int    n = 0;
  double a         [dim*n_max];
  double b         [    n_max];
  double c         [dim      ];
  double expected_x[dim      ];

  SUBCASE("Feasible") {
    /**
     * Maximize 20*x0 + 10*x1 + 15*x2
     * subject to constraints:
     *   3*x0 + 2*x1 + 5*x2 <= 55
     *   2*x0 +   x1 +   x2 <= 26
     *     x0 +   x1 + 3*x2 <= 30
     *   5*x0 + 2*x1 + 4*x2 <= 57
     *     x0               >=  0
     *            x1        >=  0
     *                   x2 >=  0
     */

    // Objective function
    c[0] = 20; c[1] = 10; c[2] = 15;

    // Constraints
    n = 4;
    a[ 0] = 3; a[ 1] = 2; a[ 2] = 5; b[0] = 55;
    a[ 3] = 2; a[ 4] = 1; a[ 5] = 1; b[1] = 26;
    a[ 6] = 1; a[ 7] = 1; a[ 8] = 3; b[2] = 30;
    a[ 9] = 5; a[10] = 2; a[11] = 4; b[3] = 57;

    // Expected result
    expected_stat = PDM_LP_FEASIBLE;
    expected_x[0] = 1.8; expected_x[1] = 20.8; expected_x[2] = 1.6;
  }

  SUBCASE("Unfeasible") {
    /**
     * Maximize x0 + 2*x1 + 3*x2
     * subject to constraints:
     *     x0 + 2*x1 + 4*x2 <=  4
     *   6*x0 + 3*x1 +   x2 <=  3
     *  -3*x0 - 2*x1 + 3*x2 <= -6
     *     x0               >=  0
     *            x1        >=  0
     *                   x2 >=  0
     */

    // Objective function
    c[0] = 1; c[1] = 2; c[2] = 3;

    // Constraints
    n = 3;
    a[0] =  1; a[1] =  2; a[2] = 4; b[0] =  4;
    a[3] =  6; a[4] =  3; a[5] = 1; b[1] =  3;
    a[6] = -3; a[7] = -2; a[8] = 3; b[2] = -6;

    // Expected result
    expected_stat = PDM_LP_UNFEASIBLE;
  }

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
  // printf("  sol  = %20.16f %20.16f %20.16f\n", x[0], x[1], x[2]);

  // Check
  CHECK(stat == expected_stat);
  if (stat == PDM_LP_FEASIBLE) {
    CHECK_EQ_C_ARRAY_FLOAT(x, expected_x, dim, eps);
  }
}
