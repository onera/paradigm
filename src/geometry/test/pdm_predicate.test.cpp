#include <stdlib.h>
#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm_mem_tool.h"
#include "pdm_mpi.h"
#include "pdm_priv.h"
#include "pdm_predicate.h"


TEST_CASE("PDM_predicate - in_circle") {

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
}

TEST_CASE("PDM_predicate - insphere") {

  PDM_predicate_exactinit();

  double a[3] = {-1.46, -2.60,  0.37};
  double b[3] = { 2.49, -1.07, -1.28};
  double c[3] = {-1.20,  2.34,  1.44};
  double d[3] = {-0.61, -0.47,  2.90};
  double e[3] = { 1.00,  2.00,  1.00};
  double f[3] = {-2.00, -3.00,  1.00};

  double volume = PDM_predicate_orient3d (a, b, c, d);
  printf("volume = %f\n", volume);

  CHECK(volume == doctest::Approx(46.765126).epsilon(0.01));

  double insphere_e = PDM_predicate_insphere (a, b, c, d, e);
  double insphere_f = PDM_predicate_insphere (a, b, c, d, f);

  printf("insphere? : e = %f, f = %f\n", insphere_e, insphere_f);
  CHECK(insphere_e == doctest::Approx(139.000714 ).epsilon(0.01));
  CHECK(insphere_f == doctest::Approx(-232.157103).epsilon(0.01));

}
