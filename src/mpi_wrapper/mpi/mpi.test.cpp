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
  // P0 envoie [10, 11, 12] et P1 reçoit
  // P1 envoie [20, 21, 22] et P0 reçoit
  std::vector<int> send_data_p0 = {10, 11, 12};
  std::vector<int> send_data_p1 = {20, 21, 22};

  // Tampons de réception (taille maximum nécessaire)
  std::vector<int> recv_data(MSG_SIZE, -1);

  PDM_MPI_Request send_request = PDM_MPI_REQUEST_NULL;
  PDM_MPI_Request recv_request = PDM_MPI_REQUEST_NULL;

  // --- 2. Initialisation des communications persistantes (Send_init/Recv_init) ---
  if (i_rank == 0) {
    // P0: Envoie à P1, Reçoit de P1

    // P0 envoie à P1
    PDM_MPI_Send_init(send_data_p0.data(), MSG_SIZE, PDM_MPI_INT, 1, TAG, pdm_comm, &send_request);

    // P0 reçoit de P1
    PDM_MPI_Recv_init(recv_data.data(), MSG_SIZE, PDM_MPI_INT, 1, TAG, pdm_comm, &recv_request);

  } else { // i_rank == 1
    // P1: Reçoit de P0, Envoie à P0

    // P1 reçoit de P0
    PDM_MPI_Recv_init(recv_data.data(), MSG_SIZE, PDM_MPI_INT, 0, TAG, pdm_comm, &recv_request);

    // P1 envoie à P0
    PDM_MPI_Send_init(send_data_p1.data(), MSG_SIZE, PDM_MPI_INT, 0, TAG, pdm_comm, &send_request);
  }

  // Vérification que les requêtes sont bien initialisées (non nulles)
  MPI_CHECK(0, send_request != PDM_MPI_REQUEST_NULL);
  MPI_CHECK(1, send_request != PDM_MPI_REQUEST_NULL);
  MPI_CHECK(0, recv_request != PDM_MPI_REQUEST_NULL);
  MPI_CHECK(1, recv_request != PDM_MPI_REQUEST_NULL);

  // --- 3. Démarrage des communications (Start) ---

  // Démarrage de la communication de réception (essentiel avant le Start du Send)
  PDM_MPI_Start(&recv_request);

  // Démarrage de la communication d'envoi
  PDM_MPI_Start(&send_request);

  // --- 4. Attente des communications (Wait) ---

  // Attente de l'envoi et de la réception
  PDM_MPI_Wait(&send_request);
  PDM_MPI_Wait(&recv_request);

  // --- 5. Vérification des résultats ---

  std::vector<int> expected_recv_p0 = {20, 21, 22}; // P0 reçoit de P1
  std::vector<int> expected_recv_p1 = {10, 11, 12}; // P1 reçoit de P0

  MPI_CHECK_EQ_C_ARRAY(0, recv_data.data(), expected_recv_p0.data(), MSG_SIZE);
  MPI_CHECK_EQ_C_ARRAY(1, recv_data.data(), expected_recv_p1.data(), MSG_SIZE);

  // --- 6. Libération des requêtes persistantes (Très important) ---

  // PDM_MPI_Request_free doit être appelée pour nettoyer la requête persistante
  // Une fois libérée, la requête doit être PDM_MPI_REQUEST_NULL
  PDM_MPI_Request_free(&send_request);
  PDM_MPI_Request_free(&recv_request);

  MPI_CHECK(0, send_request == PDM_MPI_REQUEST_NULL);
  MPI_CHECK(1, recv_request == PDM_MPI_REQUEST_NULL);

}

MPI_TEST_CASE("[PDM_MPI_Partofactiverank]", 4) {

    PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
    int i_rank;
    int n_rank;
    PDM_MPI_Comm_rank(pdm_comm, &i_rank);
    PDM_MPI_Comm_size(pdm_comm, &n_rank);

    // --- SCÉNARIO 1 : Moitié des rangs actifs (2 sur 4) ---
    // Les rangs 0 et 1 sont actifs. Les rangs 2 et 3 sont inactifs.
    // Taux d'actifs attendu : 2 / 4 = 0.5
    std::vector<int> sendcounts1(n_rank, 0); // Les envois
    std::vector<int> recvcounts1(n_rank, 0); // Les réceptions
    double part_active_rank1 = 0.0;

    // Rangs 0 et 1 actifs :
    if (i_rank == 0) {
      sendcounts1[1] = 5; // P0 envoie à P1
    } else if (i_rank == 1) {
      recvcounts1[0] = 5; // P1 reçoit de P0
    }

    // Appel de la fonction
    PDM_MPI_Partofactiverank(sendcounts1.data(),
                             recvcounts1.data(),
                             pdm_comm,
                             &part_active_rank1);

    CHECK(part_active_rank1 == doctest::Approx(0.25).epsilon(0.01));

    // --- SCÉNARIO 2 : Tous les rangs actifs (4 sur 4) ---
    std::vector<int> sendcounts2(n_rank, 0);
    std::vector<int> recvcounts2(n_rank, 0);
    double part_active_rank2 = 0.0;

    // Rangs 0, 1, 2, 3 actifs :
    if (i_rank == 0) {
        sendcounts2[1] = 1; // Envoi
    } else if (i_rank == 1) {
        recvcounts2[0] = 1; // Réception
    } else if (i_rank == 2) {
        sendcounts2[3] = 1; // Envoi
    } else if (i_rank == 3) {
        recvcounts2[2] = 1; // Réception
    }

    // Appel de la fonction
    PDM_MPI_Partofactiverank(sendcounts2.data(),
                             recvcounts2.data(),
                             pdm_comm,
                             &part_active_rank2);

    CHECK(part_active_rank2 == doctest::Approx(0.25).epsilon(0.01));

    // --- SCÉNARIO 3 : Un seul rang actif (1 sur 4) ---
    // Seul le rang 0 est actif (il a un envoi > 0). Les autres sont inactifs.
    // Taux d'actifs attendu : 1 / 4 = 0.25

    std::vector<int> sendcounts3(n_rank, 0);
    std::vector<int> recvcounts3(n_rank, 0);
    double part_active_rank3 = 0.0;

    // Rang 0 actif :
    if (i_rank == 0) {
      sendcounts3[1] = 1;
    }

    // Appel de la fonction
    PDM_MPI_Partofactiverank(sendcounts3.data(),
                             recvcounts3.data(),
                             pdm_comm,
                             &part_active_rank3);
    CHECK(part_active_rank3 == doctest::Approx(0.25).epsilon(0.01));
}
