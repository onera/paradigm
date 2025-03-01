#include <stdlib.h>
#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm_mem_tool.h"
#include "pdm_mpi.h"
#include "pdm_priv.h"
#include "pdm_predicate.h"


TEST_CASE("PDM_predicate") {

  PDM_predicate_exactinit();

  double a[2] = {-3.28,  2.17};
  double b[2] = {-4.08, -3.63};
  double c[2] = { 2.58, -2.65};
  double d[2] = { 1.04,  0.61};
  double e[2] = {-2.14, -5.49};

  double sens_plus_d = PDM_predicate_incircle (a, b, c, d);
  double sens_plus_e = PDM_predicate_incircle (a, b, c, e);

  printf("sens (+) : d = %f, e = %f\n", sens_plus_d, sens_plus_e);
  CHECK(PDM_SIGN(sens_plus_d) ==  1);
  CHECK(PDM_SIGN(sens_plus_e) == -1);

  double sens_minus_d = PDM_predicate_incircle (a, c, b, d);
  double sens_minus_e = PDM_predicate_incircle (a, c, b, e);

  printf("sens (-) : d = %f, e = %f\n", sens_minus_d, sens_minus_e);
  CHECK(PDM_SIGN(sens_minus_d) == -1);
  CHECK(PDM_SIGN(sens_minus_e) ==  1);

  double p[2] = {5.97200000000000041922, 8.94400000000000083844};
  double u[2] = {7.00000000000000000000, 8.00000000000000000000};
  double v[2] = {5.76999999999999957367, 9.81000000000000049738};
  double w[2] = {8.73000000000000042633, 9.28999999999999914735};


  double sens_plus_u  = PDM_predicate_incircle (u, v, w, p);
  double sens_minus_u = PDM_predicate_incircle (v, u, w, p);

  printf("(+) : %g\n", sens_plus_u );
  printf("(-) : %g\n", sens_minus_u);
  CHECK(PDM_SIGN(sens_plus_u ) == -1);
  CHECK(PDM_SIGN(sens_minus_u) ==  1);


  // CHECK_EQ_C_ARRAY_FLOAT(vector_out,exp_out,12,EPS);

}
