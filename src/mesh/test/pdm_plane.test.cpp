#include "doctest/extensions/doctest_mpi.h"

#include <limits.h>
#include <float.h>
#include <math.h>
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_priv.h"
#include "pdm_plane.h"
#include "pdm_logging.h"

static inline double
_rand(void) {
  return (double) rand() / (double) RAND_MAX;
}


MPI_TEST_CASE("PDM_plane_get_cartesian_plane", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  // Generate point cloud
  int n_part = 2;
  int n_pts[2] = {13, 28};

  srand(0);

  double **coord = NULL;
  PDM_malloc(coord, n_part, double *);
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_malloc(coord[i_part], n_pts[i_part] * 3, double);
    for (int i = 0; i < n_pts[i_part]; i++) {
      coord[i_part][3*i  ] = _rand();
      coord[i_part][3*i+1] = _rand();
      coord[i_part][3*i+2] = _rand();
    }
  }


  int expected_plane = -1;

  SUBCASE("XY") {
    expected_plane = 0;
    for (int i_part = 0; i_part < n_part; i_part++) {
      for (int i = 0; i < n_pts[i_part]; i++) {
        coord[i_part][3*i+2] = 42.0;
      }
    }
  }

  SUBCASE("YZ") {
    expected_plane = 1;
    for (int i_part = 0; i_part < n_part; i_part++) {
      for (int i = 0; i < n_pts[i_part]; i++) {
        coord[i_part][3*i  ] = 3.14;
      }
    }
  }

  SUBCASE("ZX") {
    expected_plane = 2;
    for (int i_part = 0; i_part < n_part; i_part++) {
      for (int i = 0; i < n_pts[i_part]; i++) {
        coord[i_part][3*i+1] = 1.23;
      }
    }
  }

  SUBCASE("Not cartesian plane") {
    expected_plane = -1;
  }


  double tolerance = 1e-14;

  int plane = PDM_plane_get_cartesian_plane(pdm_comm,
                                            n_part,
                                            n_pts,
                                            coord,
                                            tolerance);

  CHECK(plane == expected_plane);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_free(coord[i_part]);
  }
  PDM_free(coord);
}
