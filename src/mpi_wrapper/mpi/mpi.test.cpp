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



MPI_TEST_CASE("[PDM_MPI_Alltoallv_p2p]", 2) {

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
  PDM_MPI_Alltoallv_p2p(send_buf[i_rank].data(),
                        send_n  [i_rank].data(),
                        send_idx[i_rank].data(),
                        PDM_MPI_INT,
                        recv_buf.data(),
                        recv_n  .data(),
                        recv_idx.data(),
                        PDM_MPI_INT,
                        pdm_comm);

  static int recv_buf_expected_p0[4] = {1, -1, -2, -3};
  static int recv_buf_expected_p1[2] = {2, 3};

  MPI_CHECK_EQ_C_ARRAY(0, recv_buf.data(), recv_buf_expected_p0, 4);
  MPI_CHECK_EQ_C_ARRAY(1, recv_buf.data(), recv_buf_expected_p1, 2);

  if(0 == 1) {
    PDM_log_trace_array_int(recv_buf.data(), recv_idx[n_rank], "recv_buf :");
  }

}


MPI_TEST_CASE("[PDM_MPI_Alltoallv_l]", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  std::vector<std::vector<int>>    send_buf = {{1, 2, 3}, {-1, -2, -3}};
  std::vector<std::vector<int>>    send_n   = {{1, 2   }, {3, 0   }};
  std::vector<std::vector<size_t>> send_idx = {{0, 1   }, {0, 3   }};

  std::vector<int> recv_n(n_rank);
  PDM_MPI_Alltoall(send_n[i_rank].data(), 1, PDM_MPI_INT,
                   recv_n        .data(), 1, PDM_MPI_INT,
                   pdm_comm);

  std::vector<size_t> recv_idx(n_rank+1, 0);
  for(int i = 0; i < n_rank; ++i) {
    recv_idx[i+1] = recv_idx[i] + recv_n[i];
  }

  std::vector<int> recv_buf(recv_idx[n_rank], -10000);
  PDM_MPI_Alltoallv_l(send_buf[i_rank].data(),
                      send_n  [i_rank].data(),
                      send_idx[i_rank].data(),
                      PDM_MPI_INT,
                      recv_buf.data(),
                      recv_n  .data(),
                      recv_idx.data(),
                      PDM_MPI_INT,
                      pdm_comm);

  static int recv_buf_expected_p0[4] = {1, -1, -2, -3};
  static int recv_buf_expected_p1[2] = {2, 3};

  MPI_CHECK_EQ_C_ARRAY(0, recv_buf.data(), recv_buf_expected_p0, 4);
  MPI_CHECK_EQ_C_ARRAY(1, recv_buf.data(), recv_buf_expected_p1, 2);

  if(0 == 1) {
    PDM_log_trace_array_int(recv_buf.data(), recv_idx[n_rank], "recv_buf :");
  }

}



MPI_TEST_CASE("[PDM_MPI_Alltoallv_p2p_l]", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  std::vector<std::vector<int>>    send_buf = {{1, 2, 3}, {-1, -2, -3}};
  std::vector<std::vector<int>>    send_n   = {{1, 2   }, {3, 0   }};
  std::vector<std::vector<size_t>> send_idx = {{0, 1   }, {0, 3   }};

  std::vector<int> recv_n(n_rank);
  PDM_MPI_Alltoall(send_n[i_rank].data(), 1, PDM_MPI_INT,
                   recv_n        .data(), 1, PDM_MPI_INT,
                   pdm_comm);

  std::vector<size_t> recv_idx(n_rank+1, 0);
  for(int i = 0; i < n_rank; ++i) {
    recv_idx[i+1] = recv_idx[i] + recv_n[i];
  }

  std::vector<int> recv_buf(recv_idx[n_rank], -10000);
  PDM_MPI_Alltoallv_p2p_l(send_buf[i_rank].data(),
                          send_n  [i_rank].data(),
                          send_idx[i_rank].data(),
                          PDM_MPI_INT,
                          recv_buf.data(),
                          recv_n  .data(),
                          recv_idx.data(),
                          PDM_MPI_INT,
                          pdm_comm);

  static int recv_buf_expected_p0[4] = {1, -1, -2, -3};
  static int recv_buf_expected_p1[2] = {2, 3};

  MPI_CHECK_EQ_C_ARRAY(0, recv_buf.data(), recv_buf_expected_p0, 4);
  MPI_CHECK_EQ_C_ARRAY(1, recv_buf.data(), recv_buf_expected_p1, 2);

  if(0 == 1) {
    PDM_log_trace_array_int(recv_buf.data(), recv_idx[n_rank], "recv_buf :");
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


MPI_TEST_CASE("[PDM_MPI_Send_init/PDM_MPI_Recv_init]", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  const int MSG_SIZE = 3;
  const int TAG      = 99;

  // --- 1. Définition des tampons ---
  std::vector<int> send_data_p0 = {10, 11, 12};
  std::vector<int> send_data_p1 = {20, 21, 22};

  std::vector<int> recv_data(MSG_SIZE, -1);

  PDM_MPI_Request send_request = PDM_MPI_REQUEST_NULL;
  PDM_MPI_Request recv_request = PDM_MPI_REQUEST_NULL;

  if (i_rank == 0) {
    // P0 send to P1
    PDM_MPI_Send_init(send_data_p0.data(), MSG_SIZE, PDM_MPI_INT, 1, TAG, pdm_comm, &send_request);

    // P0 receive from P1
    PDM_MPI_Recv_init(recv_data.data(), MSG_SIZE, PDM_MPI_INT, 1, TAG, pdm_comm, &recv_request);

  } else { // i_rank == 1
    // P1 send to P0
    PDM_MPI_Recv_init(recv_data.data(), MSG_SIZE, PDM_MPI_INT, 0, TAG, pdm_comm, &recv_request);

    // P1 receive from P0
    PDM_MPI_Send_init(send_data_p1.data(), MSG_SIZE, PDM_MPI_INT, 0, TAG, pdm_comm, &send_request);
  }

  // Vérification que les requêtes sont bien initialisées (non nulles)
  MPI_CHECK(0, send_request != PDM_MPI_REQUEST_NULL);
  MPI_CHECK(1, send_request != PDM_MPI_REQUEST_NULL);
  MPI_CHECK(0, recv_request != PDM_MPI_REQUEST_NULL);
  MPI_CHECK(1, recv_request != PDM_MPI_REQUEST_NULL);

  // --- 3. Start communication ---
  PDM_MPI_Start(&recv_request);
  PDM_MPI_Start(&send_request);

  // --- 4. Waiting completion of send and recv ---
  PDM_MPI_Wait(&send_request);
  PDM_MPI_Wait(&recv_request);

  // --- 5. Check results ---
  std::vector<int> expected_recv_p0 = {20, 21, 22}; // P0 reçoit de P1
  std::vector<int> expected_recv_p1 = {10, 11, 12}; // P1 reçoit de P0

  MPI_CHECK_EQ_C_ARRAY(0, recv_data.data(), expected_recv_p0.data(), MSG_SIZE);
  MPI_CHECK_EQ_C_ARRAY(1, recv_data.data(), expected_recv_p1.data(), MSG_SIZE);

  // --- 6. Free persistent request ---
  // PDM_MPI_Request_free need to be called to free properly persistent request
  // After, the request must be equal to PDM_MPI_REQUEST_NULL
  PDM_MPI_Request_free(&send_request);
  PDM_MPI_Request_free(&recv_request);

  MPI_CHECK(0, send_request == PDM_MPI_REQUEST_NULL);
  MPI_CHECK(1, recv_request == PDM_MPI_REQUEST_NULL);

}


MPI_TEST_CASE("[PDM_MPI_Alltoallv_p2p_init]", 2) {

    PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
    int i_rank;
    int n_rank;
    PDM_MPI_Comm_rank(pdm_comm, &i_rank);
    PDM_MPI_Comm_size(pdm_comm, &n_rank);

    // Vérification de la taille du communicateur
    const int MSG_SIZE = 3;
    const int TAG      = 99;

    // --- 1. Define data and buffers ---
    // Data send by P0 (to P1) : {10, 11, 12}
    // Data send by P1 (to P0) : {20, 21, 22}

    // The initial send buffer should contain all the data !
    // For P0, we send 3 ints to P1, P0 send nothing to P0.
    // For P1, we send 3 ints to P0, P1 send nothing to P1.
    std::vector<int> send_data(MSG_SIZE, 0);

    // Recv data by P0 (from P1) : {20, 21, 22}
    // Recv data by P1 (from P0) : {10, 11, 12}
    std::vector<int> recv_data(MSG_SIZE, -1);

    // --- 2. Define count et displacement for alltoallv ---
    // P0 send 3 data to P1, receive 3 from P1
    // P1 send 3 data to P0, receive 3 from P0
    std::vector<int> sendcounts(n_rank, 0);
    std::vector<int> recvcounts(n_rank, 0);
    std::vector<int> sdispls(n_rank, 0);
    std::vector<int> rdispls(n_rank, 0);

    if (i_rank == 0) {
      sendcounts[1] = MSG_SIZE;
      recvcounts[1] = MSG_SIZE;
      send_data = {10, 11, 12};
    } else { // i_rank == 1
      sendcounts[0] = MSG_SIZE;
      recvcounts[0] = MSG_SIZE;
      send_data = {20, 21, 22};
    }

    // Define displacement - Trivial
    sdispls[0] = 0; sdispls[1] = 0;
    rdispls[0] = 0; rdispls[1] = 0;

    // --- 3. Init persistent request
    int n_requests_expected = 2; // 1 Send + 1 Recv by proces = 2 requests
    int n_requests_actual   = 0;
    PDM_MPI_Request *requests = NULL;

    PDM_MPI_Alltoallv_p2p_init(send_data.data(),
                               sendcounts.data(),
                               sdispls.data(),
                               PDM_MPI_INT,
                               recv_data.data(),
                               recvcounts.data(),
                               rdispls.data(),
                               PDM_MPI_INT,
                               TAG,
                               pdm_comm,
                               &n_requests_actual,
                               &requests);

    CHECK(n_requests_actual == n_requests_expected);

    // --- 6. Results validation ---
    std::vector<int> expected_recv_p0 = {20, 21, 22};
    std::vector<int> expected_recv_p1 = {10, 11, 12};

    for(int iter = 0; iter < 5; ++iter) {

      PDM_MPI_Startall(n_requests_actual, requests);

      PDM_MPI_Waitall(n_requests_actual, requests);

      MPI_CHECK_EQ_C_ARRAY(0, recv_data.data(), expected_recv_p0.data(), MSG_SIZE);
      MPI_CHECK_EQ_C_ARRAY(1, recv_data.data(), expected_recv_p1.data(), MSG_SIZE);

      for(int i = 0; i < static_cast<int>(send_data.size()); ++i) {
        send_data[i] += 1;
      }

      for(int i = 0; i < static_cast<int>(expected_recv_p0.size()); ++i) {
        expected_recv_p0[i] += 1;
      }

      for(int i = 0; i < static_cast<int>(expected_recv_p1.size()); ++i) {
        expected_recv_p1[i] += 1;
      }

    }

    for (int i = 0; i < n_requests_actual; ++i) {
      PDM_MPI_Request_free(&requests[i]);
    }

    PDM_free(requests);
}


MPI_TEST_CASE("[PDM_MPI_Partofactiverank]", 4) {

    PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
    int i_rank;
    int n_rank;
    PDM_MPI_Comm_rank(pdm_comm, &i_rank);
    PDM_MPI_Comm_size(pdm_comm, &n_rank);

    // --- SCÉNARIO 1 : Half of ranks is active (2 sur 4) ---
    // Rank 0 et 1 is active. Rank 2 et 3 sont inactive.
    std::vector<int> sendcounts1(n_rank, 0);
    std::vector<int> recvcounts1(n_rank, 0);
    double part_active_rank1 = 0.0;

    // Rangs 0 et 1 actifs :
    if (i_rank == 0) {
      sendcounts1[1] = 5;
    } else if (i_rank == 1) {
      recvcounts1[0] = 5;
    }

    // Appel de la fonction
    PDM_MPI_Partofactiverank(sendcounts1.data(),
                             recvcounts1.data(),
                             pdm_comm,
                             &part_active_rank1);

    CHECK(part_active_rank1 == doctest::Approx(0.25).epsilon(0.01));

    // --- SCÉNARIO 2 : All ranks is actives (4 sur 4) ---
    std::vector<int> sendcounts2(n_rank, 0);
    std::vector<int> recvcounts2(n_rank, 0);
    double part_active_rank2 = 0.0;

    // Rangs 0, 1, 2, 3 actifs :
    if (i_rank == 0) {
      sendcounts2[1] = 1;
    } else if (i_rank == 1) {
      recvcounts2[0] = 1;
    } else if (i_rank == 2) {
      sendcounts2[3] = 1;
    } else if (i_rank == 3) {
      recvcounts2[2] = 1;
    }

    // Appel de la fonction
    PDM_MPI_Partofactiverank(sendcounts2.data(),
                             recvcounts2.data(),
                             pdm_comm,
                             &part_active_rank2);

    CHECK(part_active_rank2 == doctest::Approx(0.25).epsilon(0.01));

    // --- SCÉNARIO 3 : One rank is active actif (1 sur 4) ---
    std::vector<int> sendcounts3(n_rank, 0);
    std::vector<int> recvcounts3(n_rank, 0);
    double part_active_rank3 = 0.0;

    // Rang 0 actif :
    if (i_rank == 0) {
      sendcounts3[1] = 1;
    }

    PDM_MPI_Partofactiverank(sendcounts3.data(),
                             recvcounts3.data(),
                             pdm_comm,
                             &part_active_rank3);
    CHECK(part_active_rank3 == doctest::Approx(0.25).epsilon(0.01));
}



MPI_TEST_CASE("[PDM_MPI_Isends/PDM_MPI_Irecvs]", 2) {

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

  std::vector<int> n_active_send = {2, 1};
  std::vector<int> n_active_recv = {2, 1};
  std::vector<std::vector<int>> active_rank_send = {{0, 1}, {0}};
  std::vector<std::vector<int>> active_rank_recv = {{0, 1}, {0}};

  std::vector<int> recv_buf(recv_idx[n_rank], -10000);
  PDM_MPI_Request *requests_send;
  PDM_MPI_Request *requests_recv;

  // Same but with shortcut
  PDM_MPI_Irecvs(recv_buf.data(),
                 recv_n  .data(),
                 recv_idx.data(),
                 PDM_MPI_INT,
                 n_active_recv   [i_rank],
                 active_rank_recv[i_rank].data(),
                 10,
                 pdm_comm,
                 &requests_recv);

  PDM_MPI_Isends(send_buf[i_rank].data(),
                 send_n  [i_rank].data(),
                 send_idx[i_rank].data(),
                 PDM_MPI_INT,
                 n_active_send   [i_rank],
                 active_rank_send[i_rank].data(),
                 10,
                 pdm_comm,
                 &requests_send);

  int recv_buf_expected_p0[4] = {1, -1, -2, -3};
  int recv_buf_expected_p1[2] = {2, 3};

  PDM_MPI_Waitall(n_active_recv[i_rank], requests_recv);
  PDM_MPI_Waitall(n_active_send[i_rank], requests_send);

  MPI_CHECK_EQ_C_ARRAY(0, recv_buf.data(), recv_buf_expected_p0, 4);
  MPI_CHECK_EQ_C_ARRAY(1, recv_buf.data(), recv_buf_expected_p1, 2);

  if(0 == 1) {
    PDM_log_trace_array_int(recv_buf.data(), recv_idx[n_rank], "recv_buf :");
  }

  PDM_free(requests_send);
  PDM_free(requests_recv);

}


MPI_TEST_CASE("[PDM_MPI_Sends_init/PDM_MPI_Recvs_init]", 2) {

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

  std::vector<int> n_active_send = {2, 1};
  std::vector<int> n_active_recv = {2, 1};
  std::vector<std::vector<int>> active_rank_send = {{0, 1}, {0}};
  std::vector<std::vector<int>> active_rank_recv = {{0, 1}, {0}};

  std::vector<int> recv_buf(recv_idx[n_rank], -10000);
  PDM_MPI_Request *requests_send;
  PDM_MPI_Request *requests_recv;

  // Same but with shortcut
  PDM_MPI_Sends_init(send_buf[i_rank].data(),
                     send_n  [i_rank].data(),
                     send_idx[i_rank].data(),
                     PDM_MPI_INT,
                     n_active_send   [i_rank],
                     active_rank_send[i_rank].data(),
                     10,
                     pdm_comm,
                     &requests_send);

  PDM_MPI_Recvs_init(recv_buf.data(),
                     recv_n  .data(),
                     recv_idx.data(),
                     PDM_MPI_INT,
                     n_active_recv   [i_rank],
                     active_rank_recv[i_rank].data(),
                     10,
                     pdm_comm,
                     &requests_recv);

  int recv_buf_expected_p0[4] = {1, -1, -2, -3};
  int recv_buf_expected_p1[2] = {2, 3};
  for(int i_iter = 0; i_iter < 5; ++i_iter) {

    PDM_MPI_Startall(n_active_recv[i_rank], requests_recv);
    PDM_MPI_Startall(n_active_send[i_rank], requests_send);

    PDM_MPI_Waitall(n_active_recv[i_rank], requests_recv);
    PDM_MPI_Waitall(n_active_send[i_rank], requests_send);

    MPI_CHECK_EQ_C_ARRAY(0, recv_buf.data(), recv_buf_expected_p0, 4);
    MPI_CHECK_EQ_C_ARRAY(1, recv_buf.data(), recv_buf_expected_p1, 2);


    // Fake buffer changement for all iteration
    for(int k = 0; k < recv_idx[n_rank]; ++k) {
      recv_buf[k] = -1;
    }

    for(int k = 0; k < static_cast<int>(send_buf[i_rank].size()); ++k) {
      send_buf[i_rank][k] += 1;
    }

    for(int i = 0; i < 4; ++i) {
      recv_buf_expected_p0[i] += 1;
    }

    for(int i = 0; i < 2; ++i) {
      recv_buf_expected_p1[i] += 1;
    }

  }

  if(0 == 1) {
    PDM_log_trace_array_int(recv_buf.data(), recv_idx[n_rank], "recv_buf :");
  }

  PDM_free(requests_send);
  PDM_free(requests_recv);

}
