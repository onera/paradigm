#include <memory>
#include <vector>
#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_gnum.h"
#include "pdm_logging.h"
#include "pdm_mpi.h"
#include "pdm_gnum_location.h"

MPI_TEST_CASE("[pdm_gnum_location] - 1p",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  int n_part_in  = 1;
  int n_part_out = 1;
  std::vector<PDM_g_num_t> gnum_in = {8, 1, 12, 8, 4};
  int n_elt_in = gnum_in.size();

  std::vector<PDM_g_num_t> gnum_out = {1, 2, 3, 4};
  int n_elt_out = gnum_out.size();

  PDM_gnum_location_t *gl = PDM_gnum_location_create(n_part_in,
                                                     n_part_out,
                                                     pdm_comm,
                                                     PDM_OWNERSHIP_KEEP);

  PDM_gnum_location_elements_set(gl, 0, n_elt_in, gnum_in.data());

  PDM_gnum_location_requested_elements_set(gl, 0, n_elt_out, gnum_out.data());

  PDM_gnum_location_compute(gl);

  int *location_idx = NULL;
  int *location     = NULL;
  PDM_gnum_location_get(gl, 0, &location_idx, &location);

  if(0 == 1) {
    PDM_log_trace_array_int       (location_idx, n_elt_out+1, "location_idx ::");
    PDM_log_trace_connectivity_int(location_idx, location, n_elt_out, "location ::");
  }

  int p0_expected_location_idx[5] = {0, 3, 3, 3, 6};
  int p0_expected_location    [6] = {0, 0, 1, 0, 0, 4};

  MPI_CHECK_EQ_C_ARRAY(0, location_idx, p0_expected_location_idx, 5);
  MPI_CHECK_EQ_C_ARRAY(0, location    , p0_expected_location    , 6);


  PDM_gnum_location_free(gl);

}


MPI_TEST_CASE("[pdm_gnum_location] - 2p",2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  int n_part_in  = 1;
  int n_part_out = 1;
  std::vector<std::vector<PDM_g_num_t>> gnum_in = {{8, 1}, {12, 8, 4}};
  int n_elt_in = gnum_in[i_rank].size();

  std::vector<std::vector<PDM_g_num_t>> gnum_out = {{1, 2, 3, 8}, {3, 8}};
  int n_elt_out = gnum_out[i_rank].size();

  PDM_gnum_location_t *gl = PDM_gnum_location_create(n_part_in,
                                                     n_part_out,
                                                     pdm_comm,
                                                     PDM_OWNERSHIP_KEEP);

  PDM_gnum_location_elements_set(gl, 0, n_elt_in, gnum_in[i_rank].data());

  PDM_gnum_location_requested_elements_set(gl, 0, n_elt_out, gnum_out[i_rank].data());

  PDM_gnum_location_compute(gl);

  int *location_idx = NULL;
  int *location     = NULL;
  PDM_gnum_location_get(gl, 0, &location_idx, &location);

  if(0 == 1) {
    PDM_log_trace_array_int       (location_idx, n_elt_out+1, "location_idx ::");
    PDM_log_trace_connectivity_int(location_idx, location, n_elt_out, "location ::");
    PDM_log_trace_array_int       (location, location_idx[n_elt_out], "location ::");
  }

  int p0_expected_location_idx[5] = {0, 3, 3, 3, 9};
  int p0_expected_location    [9] = {0, 0, 1, 0, 0, 0, 1, 0, 1};

  int p1_expected_location_idx[3] = {0, 0, 6};
  int p1_expected_location    [6] = {0, 0, 0, 1, 0, 1};

  MPI_CHECK_EQ_C_ARRAY(0, location_idx, p0_expected_location_idx, 5);
  MPI_CHECK_EQ_C_ARRAY(0, location    , p0_expected_location    , 9);

  MPI_CHECK_EQ_C_ARRAY(1, location_idx, p1_expected_location_idx, 3);
  MPI_CHECK_EQ_C_ARRAY(1, location    , p1_expected_location    , 6);

  PDM_gnum_location_free(gl);

}
