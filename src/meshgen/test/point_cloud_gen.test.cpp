#include <ext/alloc_traits.h>
#include <vector>
#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_point_cloud_gen.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_doctest.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mpi.h"


MPI_TEST_CASE("[PDM_point_cloud_gen] - cartesian - 1p", 1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  int    n_vtx_seg = 5;
  double length    = 1.;
  double h         = 0.5*length/(double) (n_vtx_seg - 1);

  const double xmin = -0.5*length;
  const double ymin = -0.5*length;
  const double zmin = -0.5*length;

  int          n_pts;
  double      *pts_coord;
  PDM_g_num_t *pts_g_num;
  PDM_point_cloud_gen_cartesian (pdm_comm,
                                 n_vtx_seg,
                                 n_vtx_seg,
                                 n_vtx_seg,
                                 xmin + h,
                                 ymin + h,
                                 zmin + h,
                                 xmin + length - h,
                                 ymin + length - h,
                                 zmin + length - h,
                                 &n_pts,
                                 &pts_coord,
                                 &pts_g_num);

  free(pts_coord);
  free(pts_g_num);

  MPI_CHECK(0, n_pts == 125);

}


MPI_TEST_CASE("[PDM_point_cloud_gen] - cartesian - 2p", 2) {
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  int    n_vtx_seg = 5;
  double length    = 1.;
  double h         = 0.5*length/(double) (n_vtx_seg - 1);

  const double xmin = -0.5*length;
  const double ymin = -0.5*length;
  const double zmin = -0.5*length;

  int          n_pts;
  double      *pts_coord;
  PDM_g_num_t *pts_g_num;
  PDM_point_cloud_gen_cartesian (pdm_comm,
                                 n_vtx_seg,
                                 n_vtx_seg,
                                 n_vtx_seg,
                                 xmin + h,
                                 ymin + h,
                                 zmin + h,
                                 xmin + length - h,
                                 ymin + length - h,
                                 zmin + length - h,
                                 &n_pts,
                                 &pts_coord,
                                 &pts_g_num);

  free(pts_coord);
  free(pts_g_num);

  MPI_CHECK(0, n_pts == 63);
  MPI_CHECK(1, n_pts == 62);

}

MPI_TEST_CASE("[PDM_point_cloud_gen] - random - 1p", 1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  PDM_g_num_t   n_g_pts_clouds = 10;
  double        length         = 1.;
  double radius                = length;

  int          n_pts_clouds;
  double      *pts_coord;
  PDM_g_num_t *pts_g_num;
  PDM_point_cloud_gen_random (pdm_comm,
                              0, // seed
                              0, // geometric_g_num
                              n_g_pts_clouds,
                              -radius, -radius, -radius,
                              radius, radius, radius,
                              &n_pts_clouds,
                              &pts_coord,
                              &pts_g_num);

  free(pts_coord);
  free(pts_g_num);

  MPI_CHECK(0, n_pts_clouds == 10);

}


MPI_TEST_CASE("[PDM_point_cloud_gen] - random - 2p", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  PDM_g_num_t   n_g_pts_clouds = 10;
  double        length         = 1.;
  double radius                = length;

  int          n_pts_clouds;
  double      *pts_coord;
  PDM_g_num_t *pts_g_num;
  PDM_point_cloud_gen_random (pdm_comm,
                              0, // seed
                              0, // geometric_g_num
                              n_g_pts_clouds,
                              -radius, -radius, -radius,
                              radius, radius, radius,
                              &n_pts_clouds,
                              &pts_coord,
                              &pts_g_num);

  free(pts_coord);
  free(pts_g_num);

  MPI_CHECK(0, n_pts_clouds == 5);
  MPI_CHECK(1, n_pts_clouds == 5);

}


MPI_TEST_CASE("[PDM_dpoint_cloud_gen_cartesian] - cartesian - 1p", 1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  int          n_vtx_seg = 5;
  double      *dpts_coord  = NULL;
  PDM_g_num_t *distrib_pts = NULL;
  PDM_dpoint_cloud_gen_cartesian(pdm_comm,
                                 n_vtx_seg,
                                 n_vtx_seg,
                                 n_vtx_seg,
                                 0., 0., 0.,
                                 1., 1., 1.,
                                 &dpts_coord,
                                 &distrib_pts);

  PDM_g_num_t distrib_expected[2] = {0, 125};

  CHECK_EQ_C_ARRAY( distrib_pts, distrib_expected, n_rank+1);

  free(dpts_coord);
  free(distrib_pts);

}


MPI_TEST_CASE("[PDM_dpoint_cloud_gen_cartesian] - cartesian - 2p", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  int          n_vtx_seg = 5;
  double      *dpts_coord  = NULL;
  PDM_g_num_t *distrib_pts = NULL;
  PDM_dpoint_cloud_gen_cartesian(pdm_comm,
                                 n_vtx_seg,
                                 n_vtx_seg,
                                 n_vtx_seg,
                                 0., 0., 0.,
                                 1., 1., 1.,
                                 &dpts_coord,
                                 &distrib_pts);

  PDM_g_num_t distrib_expected[3] = {0, 63, 125};

  CHECK_EQ_C_ARRAY( distrib_pts, distrib_expected, n_rank+1);

  free(dpts_coord);
  free(distrib_pts);

}

MPI_TEST_CASE("[PDM_dpoint_cloud_gen_random] - random - 1p", 1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  PDM_g_num_t n_g_pts   = 10;
  double      *dpts_coord  = NULL;
  PDM_g_num_t *distrib_pts = NULL;
  PDM_dpoint_cloud_gen_random(pdm_comm,
                              0, // seed
                              n_g_pts,
                              0., 0., 0.,
                              1., 1., 1.,
                              &dpts_coord,
                              &distrib_pts);

  PDM_g_num_t distrib_expected[2] = {0, 10};

  CHECK_EQ_C_ARRAY( distrib_pts, distrib_expected, n_rank+1);

  free(dpts_coord);
  free(distrib_pts);

}


MPI_TEST_CASE("[PDM_dpoint_cloud_gen_random] - random - 2p", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  PDM_g_num_t n_g_pts   = 10;
  double      *dpts_coord  = NULL;
  PDM_g_num_t *distrib_pts = NULL;
  PDM_dpoint_cloud_gen_random(pdm_comm,
                              0, // seed
                              n_g_pts,
                              0., 0., 0.,
                              1., 1., 1.,
                              &dpts_coord,
                              &distrib_pts);

  PDM_g_num_t distrib_expected[3] = {0, 5, 10};

  CHECK_EQ_C_ARRAY( distrib_pts, distrib_expected, n_rank+1);

  free(dpts_coord);
  free(distrib_pts);

}
