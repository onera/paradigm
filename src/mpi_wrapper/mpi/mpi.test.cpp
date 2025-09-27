#include <vector>
#include <numeric>
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_array.h"
#include "pdm_doctest.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include "pdm_mpi.h"
#include <functional>


MPI_TEST_CASE("[PDM_MPI_Topo_test]", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  SUBCASE("Standard communicator has MPI_UNDEFINED topology") {
    int status;
    PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
    int err_code = PDM_MPI_Topo_test(comm, &status);

    CHECK(err_code == MPI_SUCCESS);
    CHECK(status == PDM_MPI_COMM_UNDEFINED);
  }

  SUBCASE("Standard communicator has MPI_UNDEFINED topology") {
    int status;
    int err_code = PDM_MPI_Topo_test(pdm_comm, &status);

    CHECK(err_code == MPI_SUCCESS);
    CHECK(status == PDM_MPI_COMM_UNDEFINED);
  }

  SUBCASE("Distributed graph communicator has MPI_DIST_GRAPH topology") {
    int status;
    PDM_MPI_Comm dist_graph_comm;

    int sources[] = {i_rank};
    int destinations[] = {(i_rank + 1) % n_rank};

    int err_code_create = PDM_MPI_Dist_graph_create_adjacent(pdm_comm,
                                                             1, sources,
                                                             1, destinations,
                                                             0, &dist_graph_comm);

    REQUIRE(err_code_create == MPI_SUCCESS);

    // On teste la topologie du nouveau communicateur
    int err_code_test = PDM_MPI_Topo_test(dist_graph_comm, &status);

    CHECK(err_code_test == MPI_SUCCESS);
    CHECK(status == PDM_MPI_DIST_GRAPH);

    // Nettoyage
    PDM_MPI_Comm_free(&dist_graph_comm);
  }
}

MPI_TEST_CASE("[PDM_MPI_Ialltoallv_p2p]", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  std::vector<std::vector<int>> send_buf = {{1, 2, 3}, {-1, -2, -3}};
  std::vector<std::vector<int>> send_n   = {{1, 2   }, {3, 0   }};
  std::vector<std::vector<int>> send_idx = {{0, 1   }, {0, 3   }};

  std::vector<int> recv_n(n_rank);
  PDM_MPI_Alltoall(send_n[i_rank].data(), 1, PDM_MPI_INT,
                   recv_n        .data(), 1, PDM_MPI_INT,
                   pdm_comm);

  std::vector<int> recv_idx(n_rank+1, 0);
  for(int i = 0; i < n_rank; ++i) {
    recv_idx[i+1] = recv_idx[i] + recv_n[i];
  }

  std::vector<int> recv_buf(recv_idx[n_rank], -10000);
  int n_request = 0;
  PDM_MPI_Request *requests;
  PDM_MPI_Ialltoallv_p2p(send_buf[i_rank].data(),
                         send_n  [i_rank].data(),
                         send_idx[i_rank].data(),
                         PDM_MPI_INT,
                         0,
                         NULL,
                         recv_buf.data(),
                         recv_n  .data(),
                         recv_idx.data(),
                         PDM_MPI_INT,
                         0,
                         NULL,
                         10,
                         pdm_comm,
                         &n_request,
                         &requests);

  for(int i = 0; i < n_request; ++i) {
    PDM_MPI_Wait(&requests[i]);
  }

  static int recv_buf_expected_p0[4] = {1, -1, -2, -3};
  static int recv_buf_expected_p1[2] = {2, 3};

  MPI_CHECK(0, n_request == 4); // 2 Send / 2 Recv
  MPI_CHECK(1, n_request == 2); // 1 Send / 1 Recv

  MPI_CHECK_EQ_C_ARRAY(0, recv_buf.data(), recv_buf_expected_p0, 4);
  MPI_CHECK_EQ_C_ARRAY(1, recv_buf.data(), recv_buf_expected_p1, 2);

  if(0 == 1) {
    PDM_log_trace_array_int(recv_buf.data(), recv_idx[n_rank], "recv_buf :");
  }

  PDM_free(requests);

  for(int i = 0; i < static_cast<int>(recv_buf.size()); ++i) {
    recv_buf[i] = -10000;
  }

  std::vector<int> n_active_send = {2, 1};
  std::vector<int> n_active_recv = {2, 1};
  std::vector<std::vector<int>> active_rank_send = {{0, 1}, {0}};
  std::vector<std::vector<int>> active_rank_recv = {{0, 1}, {0}};

  if(0 == 1) {
    PDM_log_trace_array_int(recv_n.data(), n_rank, "recv_n :");
  }

  // Same but with shortcut
  PDM_MPI_Ialltoallv_select_p2p(send_buf[i_rank].data(),
                                send_n  [i_rank].data(),
                                send_idx[i_rank].data(),
                                PDM_MPI_INT,
                                n_active_send   [i_rank],
                                active_rank_send[i_rank].data(),
                                recv_buf.data(),
                                recv_n  .data(),
                                recv_idx.data(),
                                PDM_MPI_INT,
                                n_active_recv   [i_rank],
                                active_rank_recv[i_rank].data(),
                                10,
                                pdm_comm,
                                &n_request,
                                &requests);

  for(int i = 0; i < n_request; ++i) {
    PDM_MPI_Wait(&requests[i]);
  }

  MPI_CHECK(0, n_request == 4); // 2 Send / 2 Recv
  MPI_CHECK(1, n_request == 2); // 1 Send / 1 Recv

  MPI_CHECK_EQ_C_ARRAY(0, recv_buf.data(), recv_buf_expected_p0, 4);
  MPI_CHECK_EQ_C_ARRAY(1, recv_buf.data(), recv_buf_expected_p1, 2);

  if(0 == 1) {
    PDM_log_trace_array_int(recv_buf.data(), recv_idx[n_rank], "recv_buf :");
  }

  PDM_free(requests);

}
