#include <memory>
#include <vector>
#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_closest_points.h"
#include "pdm_logging.h"
#include "pdm_mpi.h"

MPI_TEST_CASE("[pdm_closest_points] - 1p",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  std::vector<double>      tgt_coords = {0.5, 0.5, 0};
  std::vector<PDM_g_num_t> tgt_gnum   = {1};
  std::vector<double>      src_coords = {0.0, 0.0, 0};
  std::vector<PDM_g_num_t> src_gnum   = {2};

  int n_tgt_pts = tgt_gnum.size();
  int n_src_pts = src_gnum.size();

  int n_closest = 1;

  PDM_closest_point_t* cls = PDM_closest_points_create(pdm_comm,
                                                       n_closest,
                                                       PDM_OWNERSHIP_KEEP);

  PDM_closest_points_n_part_cloud_set(cls, 1, 1);

  PDM_closest_points_tgt_cloud_set(cls,
                                   0,
                                   n_tgt_pts,
                                   tgt_coords.data(),
                                   tgt_gnum  .data());

  PDM_closest_points_src_cloud_set(cls,
                                   0,
                                   n_src_pts,
                                   src_coords.data(),
                                   src_gnum  .data());

  PDM_closest_points_compute(cls);

  PDM_g_num_t *closest_src_gnum     = NULL;
  double      *closest_src_distance = NULL;
  PDM_closest_points_get(cls,
                         0,
                         &closest_src_gnum,
                         &closest_src_distance);
  if(0 == 1) {
    PDM_log_trace_array_long  (closest_src_gnum    , n_tgt_pts * n_closest, "closest_src_gnum     ::");
    PDM_log_trace_array_double(closest_src_distance, n_tgt_pts * n_closest, "closest_src_distance ::");
  }

  PDM_g_num_t closest_src_gnum_expected_p0[1] = {2  };
  double      closest_src_dist_expected_p0[1] = {0.5};

  MPI_CHECK_EQ_C_ARRAY(0, closest_src_gnum, closest_src_gnum_expected_p0, n_tgt_pts * n_closest);

  for (int i = 0; i < n_tgt_pts * n_closest; ++i) {
    CHECK(closest_src_distance[i] == doctest::Approx(closest_src_dist_expected_p0[i]).epsilon(0.01));
  }

  PDM_closest_points_free(cls);
}



MPI_TEST_CASE("[pdm_closest_points] - 2p",2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  std::vector<std::vector<double>>      tgt_coords = {{0.5, 0.5, 0}, {1., 1., 0}};
  // std::vector<std::vector<PDM_g_num_t>> tgt_gnum   = {{1}          , {4}}; // Partial block is not managed yet
  std::vector<std::vector<PDM_g_num_t>> tgt_gnum   = {{1}            , {2}};
  std::vector<std::vector<double>>      src_coords = {{0.25, 0.25, 0}, {0.75, 1., 0}};
  std::vector<std::vector<PDM_g_num_t>> src_gnum   = {{2}          , {1}};

  int n_tgt_pts = tgt_gnum[i_rank].size();
  int n_src_pts = src_gnum[i_rank].size();

  int n_closest = 1;

  PDM_closest_point_t* cls = PDM_closest_points_create(pdm_comm,
                                                       n_closest,
                                                       PDM_OWNERSHIP_KEEP);

  PDM_closest_points_n_part_cloud_set(cls, 1, 1);

  PDM_closest_points_tgt_cloud_set(cls,
                                   0,
                                   n_tgt_pts,
                                   tgt_coords[i_rank].data(),
                                   tgt_gnum  [i_rank].data());

  PDM_closest_points_src_cloud_set(cls,
                                   0,
                                   n_src_pts,
                                   src_coords[i_rank].data(),
                                   src_gnum  [i_rank].data());

  PDM_closest_points_compute(cls);

  PDM_g_num_t *closest_src_gnum     = NULL;
  double      *closest_src_distance = NULL;
  PDM_closest_points_get(cls,
                         0,
                         &closest_src_gnum,
                         &closest_src_distance);

  if(0 == 1) {
    PDM_log_trace_array_long  (closest_src_gnum    , n_tgt_pts * n_closest, "closest_src_gnum     ::");
    PDM_log_trace_array_double(closest_src_distance, n_tgt_pts * n_closest, "closest_src_distance ::");
  }

  PDM_g_num_t closest_src_gnum_expected_p0[1] = {2    };
  double      closest_src_dist_expected_p0[1] = {0.125};

  PDM_g_num_t closest_src_gnum_expected_p1[1] = {1     };
  double      closest_src_dist_expected_p1[1] = {0.0625};

  MPI_CHECK_EQ_C_ARRAY(0, closest_src_gnum, closest_src_gnum_expected_p0, n_tgt_pts * n_closest);
  MPI_CHECK_EQ_C_ARRAY(1, closest_src_gnum, closest_src_gnum_expected_p1, n_tgt_pts * n_closest);

  for (int i = 0; i < n_tgt_pts * n_closest; ++i) {
    MPI_CHECK(0, closest_src_distance[i] == doctest::Approx(closest_src_dist_expected_p0[i]).epsilon(0.01));
    MPI_CHECK(1, closest_src_distance[i] == doctest::Approx(closest_src_dist_expected_p1[i]).epsilon(0.01));
  }

  PDM_closest_points_free(cls);
}
