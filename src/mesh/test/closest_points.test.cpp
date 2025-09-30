#include <memory>
#include <vector>
#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_closest_points.h"
#include "pdm_logging.h"
#include "pdm_mpi.h"
#include "pdm_mem_tool.h"
#include "pdm_priv.h"

MPI_TEST_CASE("[pdm_closest_points] - 1p",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);

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
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);

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



MPI_TEST_CASE("[pdm_closest_points] - 1p - small", 1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int n_src = 27;
  PDM_g_num_t src_gnum[27] = {1,21,17,18,19,20,22,15,23,24,25,26,16,14,2,7,3,4,5,6,8,13,9,10,11,12,27};
  double src_coords[3*27] = {
    0.166667, 0.166667, 0.166667,
    0.833333, 0.166667, 0.833333,
    0.500000, 0.833333, 0.500000,
    0.833333, 0.833333, 0.500000,
    0.166667, 0.166667, 0.833333,
    0.500000, 0.166667, 0.833333,
    0.166667, 0.500000, 0.833333,
    0.833333, 0.500000, 0.500000,
    0.500000, 0.500000, 0.833333,
    0.833333, 0.500000, 0.833333,
    0.166667, 0.833333, 0.833333,
    0.500000, 0.833333, 0.833333,
    0.166667, 0.833333, 0.500000,
    0.500000, 0.500000, 0.500000,
    0.500000, 0.166667, 0.166667,
    0.166667, 0.833333, 0.166667,
    0.833333, 0.166667, 0.166667,
    0.166667, 0.500000, 0.166667,
    0.500000, 0.500000, 0.166667,
    0.833333, 0.500000, 0.166667,
    0.500000, 0.833333, 0.166667,
    0.166667, 0.500000, 0.500000,
    0.833333, 0.833333, 0.166667,
    0.166667, 0.166667, 0.500000,
    0.500000, 0.166667, 0.500000,
    0.833333, 0.166667, 0.500000,
    0.833333, 0.833333, 0.833333
  };

  int n_tgt = 56;
  PDM_g_num_t tgt_gnum[56] = {36,41,40,39,38,37,35,34,33,32,31,42,43,29,50,55,54,53,52,51,49,44,48,47,46,45,30,28,1,7,12,11,10,9,8,6,14,5,4,3,2,13,27,21,26,25,24,23,22,20,15,19,18,17,16,56};
  double tgt_coords[3*56] = {
    1.225000, 1.125000, 0.625000,
    1.475000, 1.375000, 0.625000,
    1.225000, 1.375000, 0.625000,
    0.975000, 1.375000, 0.625000,
    1.725000, 1.125000, 0.625000,
    1.475000, 1.125000, 0.625000,
    0.975000, 1.125000, 0.625000,
    1.725000, 0.875000, 0.625000,
    1.475000, 0.875000, 0.625000,
    1.225000, 0.875000, 0.625000,
    1.725000, 0.625000, 0.625000,
    1.725000, 1.375000, 0.625000,
    1.225000, 0.625000, 0.875000,
    1.225000, 0.625000, 0.625000,
    1.225000, 1.125000, 0.875000,
    1.475000, 1.375000, 0.875000,
    1.225000, 1.375000, 0.875000,
    0.975000, 1.375000, 0.875000,
    1.725000, 1.125000, 0.875000,
    1.475000, 1.125000, 0.875000,
    0.975000, 1.125000, 0.875000,
    1.475000, 0.625000, 0.875000,
    1.725000, 0.875000, 0.875000,
    1.475000, 0.875000, 0.875000,
    1.225000, 0.875000, 0.875000,
    1.725000, 0.625000, 0.875000,
    1.475000, 0.625000, 0.625000,
    1.725000, 1.375000, 0.375000,
    1.225000, 0.625000, 0.125000,
    0.975000, 1.125000, 0.125000,
    1.225000, 1.375000, 0.125000,
    0.975000, 1.375000, 0.125000,
    1.725000, 1.125000, 0.125000,
    1.475000, 1.125000, 0.125000,
    1.225000, 1.125000, 0.125000,
    1.725000, 0.875000, 0.125000,
    1.725000, 1.375000, 0.125000,
    1.475000, 0.875000, 0.125000,
    1.225000, 0.875000, 0.125000,
    1.725000, 0.625000, 0.125000,
    1.475000, 0.625000, 0.125000,
    1.475000, 1.375000, 0.125000,
    1.475000, 1.375000, 0.375000,
    0.975000, 1.125000, 0.375000,
    1.225000, 1.375000, 0.375000,
    0.975000, 1.375000, 0.375000,
    1.725000, 1.125000, 0.375000,
    1.475000, 1.125000, 0.375000,
    1.225000, 1.125000, 0.375000,
    1.725000, 0.875000, 0.375000,
    1.225000, 0.625000, 0.375000,
    1.475000, 0.875000, 0.375000,
    1.225000, 0.875000, 0.375000,
    1.725000, 0.625000, 0.375000,
    1.475000, 0.625000, 0.375000,
    1.725000, 1.375000, 0.875000
  };

  //Brute force expected results
  PDM_g_num_t *expected_closest_gnum = NULL;
  double      *expected_closest_dist = NULL;
  PDM_malloc(expected_closest_gnum, n_tgt, PDM_g_num_t);
  PDM_malloc(expected_closest_dist, n_tgt, double     );
  for (int i = 0; i < n_tgt; i++) {
    PDM_g_num_t arg_min = -1;
    double min_dist = 1E12;
    for (int j = 0; j < n_src; j++) {
      double dist = (tgt_coords[3*i+0] - src_coords[3*j+0])*(tgt_coords[3*i+0] - src_coords[3*j+0])
                  + (tgt_coords[3*i+1] - src_coords[3*j+1])*(tgt_coords[3*i+1] - src_coords[3*j+1])
                  + (tgt_coords[3*i+2] - src_coords[3*j+2])*(tgt_coords[3*i+2] - src_coords[3*j+2]);
      if (dist < min_dist) {
        arg_min = src_gnum[j];
        min_dist = dist;
      }
    }
    expected_closest_gnum[i] = arg_min;
    expected_closest_dist[i] = min_dist;
  }


  PDM_closest_point_t* clsp = PDM_closest_points_create (pdm_comm,
                                                         1,
                                                         PDM_OWNERSHIP_USER);

  PDM_closest_points_n_part_cloud_set (clsp, 1, 1);

  PDM_closest_points_src_cloud_set (clsp,
                                    0,
                                    n_src,
                                    src_coords,
                                    src_gnum);

  PDM_closest_points_tgt_cloud_set (clsp,
                                    0,
                                    n_tgt,
                                    tgt_coords,
                                    tgt_gnum);


  PDM_closest_points_compute (clsp);



  PDM_g_num_t *closest_src_gnum = NULL;
  double      *closest_src_dist = NULL;

  PDM_closest_points_get (clsp,
                          0,
                          &closest_src_gnum,
                          &closest_src_dist);


  PDM_closest_points_free (clsp);

  for (int i = 0; i < n_tgt; i++) {
    assert (PDM_ABS(closest_src_dist[i] - expected_closest_dist[i]) < 1E-6);
    CHECK (closest_src_gnum[i] == expected_closest_gnum[i]);
  }

  PDM_free(expected_closest_gnum);
  PDM_free(expected_closest_dist);
  PDM_free(closest_src_gnum);
  PDM_free(closest_src_dist);
}
