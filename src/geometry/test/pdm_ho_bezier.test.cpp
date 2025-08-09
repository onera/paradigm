#include <cstdio>
#include <math.h>
#include "doctest/doctest.h"
#include "pdm_ho_bezier.h"

TEST_CASE("PDM_ho_bezier_triangle_location") {

  double point_coord[3] = {
    -1.7516283337977998e+03, -2.4366417578898277e+03, -1.6904612902513327e+03
  };

  double node_coord[3*10] = {
    -8.0296495277832003e+02, -1.3890805515241000e+03, -9.6136786496847003e+02,
    -7.7952673129348398e+02, -1.3930219640567957e+03, -9.9049326381078606e+02,
    -7.5524779978206868e+02, -1.3962348401410800e+03, -1.0191639157580662e+03,
    -7.3017588396461997e+02, -1.3984962345714000e+03, -1.0470787613708001e+03,
    -7.7634577363059907e+02, -1.4168863752472719e+03, -9.6070418259120015e+02,
    -7.5232512212535721e+02, -1.4203639364413416e+03, -9.8951284758796032e+02,
    -7.2749705788034714e+02, -1.4229139556528405e+03, -1.0177074408379543e+03,
    -7.4877630148282299e+02, -1.4439180852325642e+03, -9.5946534676684917e+02,
    -7.2426425149072077e+02, -1.4467704408092804e+03, -9.8784935035719172e+02,
    -7.2039530384498005e+02, -1.4698580869662001e+03, -9.5755153831604002e+02
  };

  double proj_coord[3];
  double uvw[3];
  PDM_ho_bezier_triangle_location(3,
                                  10,
                                  node_coord,
                                  point_coord,
                                  proj_coord,
                                  uvw);

  // printf("proj_coord = %12.5e / %12.5e / %12.5e \n", proj_coord[0], proj_coord[1], proj_coord[2]);
  // printf("proj_coord = %12.5e / %12.5e / %12.5e \n", uvw[0], uvw[1], uvw[2]);

  double expexted_proj_coord[3] = {-7.42836e+02, -1.42693e+03, -9.91425e+02};
  double expected_uvw       [3] = { 1.67385e-01,  1.81742e-01,  6.50874e-01};

  CHECK(proj_coord[0] == doctest::Approx(expexted_proj_coord[0]).epsilon(0.01));
  CHECK(proj_coord[1] == doctest::Approx(expexted_proj_coord[1]).epsilon(0.01));
  CHECK(proj_coord[2] == doctest::Approx(expexted_proj_coord[2]).epsilon(0.01));

  CHECK(uvw[0] == doctest::Approx(expected_uvw[0]).epsilon(0.01));
  CHECK(uvw[1] == doctest::Approx(expected_uvw[1]).epsilon(0.01));
  CHECK(uvw[2] == doctest::Approx(expected_uvw[2]).epsilon(0.01));

}
