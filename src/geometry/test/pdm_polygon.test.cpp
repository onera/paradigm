#include "doctest/extensions/doctest_mpi.h"
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include "doctest/doctest.h"
#include "pdm.h"
#include "pdm_mem_tool.h"
#include "pdm_mpi.h"
#include "pdm_polygon.h"


MPI_TEST_CASE("PDM_polygon_evaluate_position", 1) {

  double x[3];
  const int n_pts=5;
  double closestPoint[3], expected_closestPoint[3];
  double minDist2, expected_minDist2;
  PDM_polygon_status_t location, expected_location;
  double tol = 1e-16;

  double rand_fact1 = 0.1 * (double) rand() / (double) RAND_MAX;
  double rand_fact2 = 0.1 * (double) rand() / (double) RAND_MAX;

  double pts[3*n_pts] = { 0.0, 0.0, 0.0,
                          1.0, 0.0, 0.0,
                          0.5, 1.0, 0.0,
                          2.0 + rand_fact1, 2.0 + rand_fact2, 0.0,
                         -0.5, 1.0, 0.0};

  SUBCASE("point inside case 1"){
    x[0] = 0.0; x[1] = 1.0; x[2] = 1.0;
    expected_location = PDM_POLYGON_INSIDE;
    expected_minDist2 = 1.0;
    expected_closestPoint[0] = 0.0; expected_closestPoint[1] = 1.0;  expected_closestPoint[2] = 0.0;
  }
  SUBCASE("point inside case 2"){
    x[0] = 2.0 + rand_fact1; x[1] = 2.0 + rand_fact2; x[2] = 0.0;
    expected_location = PDM_POLYGON_INSIDE;
    expected_minDist2 = 0.0;
    expected_closestPoint[0] = x[0]; expected_closestPoint[1] = x[1];  expected_closestPoint[2] = x[2];
  }
  SUBCASE("point outside case 1"){
    x[0] = 2.0; x[1] = -1.0; x[2] = 0.0;
    expected_location = PDM_POLYGON_OUTSIDE;
    expected_minDist2 = 2.0;
    expected_closestPoint[0] = 1.0; expected_closestPoint[1] = 0.0;  expected_closestPoint[2] = 0.0;
  }
  SUBCASE("point outside case 2"){
    x[0] = 1.0; x[1] = 0.625; x[2] = 1.0;
    expected_location = PDM_POLYGON_OUTSIDE;
    expected_minDist2 = 1.078125;
    expected_closestPoint[0] = 0.75; expected_closestPoint[1] = 0.5;  expected_closestPoint[2] = 0.0;
  }

  location = PDM_polygon_evaluate_position(x,
                                           n_pts,
                                           pts,
                                           closestPoint,
                                           &minDist2);

  CHECK(location == expected_location);
  CHECK(fabs(minDist2 - expected_minDist2) < tol);

  for (int i=0; i<3; i++){
    CHECK(fabs(closestPoint[i] - expected_closestPoint[i]) < tol);    
  }

}
