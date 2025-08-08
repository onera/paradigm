#include <vector>
#include <numeric>
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_array.h"
#include "pdm_doctest.h"
#include "pdm_exchange_helper.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include <functional>


MPI_TEST_CASE("[PDM_exchange_helper] - Create ", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  PDM_exchange_helper_t* exch_helper = PDM_exchange_helper_create(pdm_comm, 1);

  std::vector<std::vector<int>> send_buffer = {{1 ,  2,  3,  4, 5},
                                               {10, 20, 30, 40   }};

  std::vector<std::vector<int>> send_idx = {{0, 3, 5},
                                            {0, 1, 4}};

  std::vector<int> send_n  (n_rank  );
  std::vector<int> recv_n  (n_rank  );
  std::vector<int> recv_idx(n_rank+1);

  for(int i = 0; i < n_rank; ++i ) {
    send_n[i] = send_idx[i_rank][i+1] - send_idx[i_rank][i];
  }

  PDM_MPI_Alltoall(send_n.data(), 1, PDM_MPI_INT,
                   recv_n.data(), 1, PDM_MPI_INT, pdm_comm);

  recv_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i ) {
    recv_idx[i+1] = recv_idx[i] + recv_n[i];
  }
  std::vector<int> recv_buffer(recv_idx[n_rank]);

  // Persistent
  int request_id_0 = PDM_exchange_helper_exch_init(exch_helper,
                                                   PDM_MPI_COMM_KIND_P2P,
                                                   sizeof(int),
                                                   1,
                                                   send_idx   [i_rank].data(),
                                                   send_n             .data(),
                                                   send_buffer[i_rank].data(),
                                                   recv_idx           .data(),
                                                   recv_n             .data(),
                                                   recv_buffer        .data());

  CHECK(request_id_0 == 0);

  int request_id_1 = PDM_exchange_helper_exch_init(exch_helper,
                                                   PDM_MPI_COMM_KIND_P2P,
                                                   sizeof(int),
                                                   1,
                                                   send_idx   [i_rank].data(),
                                                   send_n             .data(),
                                                   send_buffer[i_rank].data(),
                                                   recv_idx           .data(),
                                                   recv_n             .data(),
                                                   recv_buffer        .data());

  CHECK(request_id_1 == 1);

  PDM_exchange_helper_exch_free(exch_helper, request_id_0);

  int request_id_2 = PDM_exchange_helper_exch_init(exch_helper,
                                                   PDM_MPI_COMM_KIND_P2P,
                                                   sizeof(int),
                                                   1,
                                                   send_idx   [i_rank].data(),
                                                   send_n             .data(),
                                                   send_buffer[i_rank].data(),
                                                   recv_idx           .data(),
                                                   recv_n             .data(),
                                                   recv_buffer        .data());

  CHECK(request_id_2 == 0);

  PDM_exchange_helper_exch_free(exch_helper, request_id_1);
  PDM_exchange_helper_exch_free(exch_helper, request_id_2);

  PDM_exchange_helper_free(exch_helper);
}


MPI_TEST_CASE("[PDM_exchange_helper] - Exch ", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  PDM_exchange_helper_t* exch_helper = PDM_exchange_helper_create(pdm_comm, 10);

  std::vector<std::vector<int>> send_buffer = {{1 ,  2,  3,  4, 5},
                                               {10, 20, 30, 40   }};

  std::vector<std::vector<int>> send_idx = {{0, 3, 5},
                                            {0, 1, 4}};

  std::vector<int> send_n  (n_rank  );
  std::vector<int> recv_n  (n_rank  );
  std::vector<int> recv_idx(n_rank+1);

  for(int i = 0; i < n_rank; ++i ) {
    send_n[i] = send_idx[i_rank][i+1] - send_idx[i_rank][i];
  }

  PDM_MPI_Alltoall(send_n.data(), 1, PDM_MPI_INT,
                   recv_n.data(), 1, PDM_MPI_INT, pdm_comm);

  recv_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i ) {
    recv_idx[i+1] = recv_idx[i] + recv_n[i];
  }
  std::vector<int> recv_buffer(recv_idx[n_rank]);

  std::vector<PDM_mpi_comm_kind_t> lexch_type = {PDM_MPI_COMM_KIND_P2P,
                                                 PDM_MPI_COMM_KIND_COLLECTIVE,
                                                 PDM_MPI_COMM_KIND_WIN_RMA};
  int n_type_exch = lexch_type.size();


  static int recv_buffer_expected_p0[4] = {1, 2, 3, 10};
  static int recv_buffer_expected_p1[5] = {4, 5, 20, 30, 40};

  for(int i_type_exch = 0; i_type_exch < n_type_exch-1; ++i_type_exch) { // RMA + bloquant -> Pas géré

    // Synchronous
    PDM_exchange_helper_exch(exch_helper,
                             lexch_type[i_type_exch],
                             1,
                             sizeof(int),
                             send_idx   [i_rank].data(),
                             send_n             .data(),
                             send_buffer[i_rank].data(),
                             recv_idx           .data(),
                             recv_n             .data(),
                             recv_buffer        .data());


    MPI_CHECK_EQ_C_ARRAY(0, recv_buffer, recv_buffer_expected_p0, recv_idx.back());
    MPI_CHECK_EQ_C_ARRAY(1, recv_buffer, recv_buffer_expected_p1, recv_idx.back());
  }

  static int recv_buffer_expected2_p0[4] = {11, 12, 13, 20};
  static int recv_buffer_expected2_p1[5] = {14, 15, 30, 40, 50};

  int n_try = 4;
  for(int i_try = 0; i_try < n_try; ++i_try) {
    for(int i_type_exch = 0; i_type_exch < n_type_exch; ++i_type_exch) {

#ifndef HAVE_MPI_COLLECTIVE_INIT_FUNC
      if(lexch_type[i_type_exch] == PDM_MPI_COMM_KIND_COLLECTIVE) {
        continue;
      }
#endif

      // Persistent
      int request_id = PDM_exchange_helper_exch_init(exch_helper,
                                                     lexch_type[i_type_exch],
                                                     sizeof(int),
                                                     1,
                                                     send_idx   [i_rank].data(),
                                                     send_n             .data(),
                                                     send_buffer[i_rank].data(),
                                                     recv_idx           .data(),
                                                     recv_n             .data(),
                                                     recv_buffer        .data());

      PDM_exchange_helper_exch_start(exch_helper, request_id);

      PDM_exchange_helper_exch_wait(exch_helper, request_id);

      if(0 == 1) {
        PDM_log_trace_array_int(recv_buffer.data(), recv_idx[n_rank], "recv_buffer ::");
      }

      MPI_CHECK_EQ_C_ARRAY(0, recv_buffer, recv_buffer_expected_p0, recv_idx.back());
      MPI_CHECK_EQ_C_ARRAY(1, recv_buffer, recv_buffer_expected_p1, recv_idx.back());

      // Changement des buffers !
      for(int i = 0; i < static_cast<int>(send_buffer[i_rank].size()); ++i) {
        send_buffer[i_rank][i] += 10;
      }

      PDM_exchange_helper_exch_start(exch_helper, request_id);

      PDM_exchange_helper_exch_wait(exch_helper, request_id);

      if(0 == 1) {
        PDM_log_trace_array_int(recv_buffer.data(), recv_idx[n_rank], "recv_buffer ::");
      }

      MPI_CHECK_EQ_C_ARRAY(0, recv_buffer, recv_buffer_expected2_p0, recv_idx.back());
      MPI_CHECK_EQ_C_ARRAY(1, recv_buffer, recv_buffer_expected2_p1, recv_idx.back());

      PDM_exchange_helper_exch_free(exch_helper, request_id);

      // Changement des buffers !
      for(int i = 0; i < static_cast<int>(send_buffer[i_rank].size()); ++i) {
        send_buffer[i_rank][i] -= 10;
      }
    }
  }


  for(int i_try = 0; i_try < n_try; ++i_try) {
    for(int i_type_exch = 0; i_type_exch < n_type_exch; ++i_type_exch) {

      // Persistent
      int request_id = PDM_exchange_helper_iexch(exch_helper,
                                                 lexch_type[i_type_exch],
                                                 sizeof(int),
                                                 1,
                                                 send_idx   [i_rank].data(),
                                                 send_n             .data(),
                                                 send_buffer[i_rank].data(),
                                                 recv_idx           .data(),
                                                 recv_n             .data(),
                                                 recv_buffer        .data());

      PDM_exchange_helper_exch_wait(exch_helper, request_id);

      if(0 == 1) {
        PDM_log_trace_array_int(recv_buffer.data(), recv_idx[n_rank], "recv_buffer ::");
      }

      MPI_CHECK_EQ_C_ARRAY(0, recv_buffer, recv_buffer_expected_p0, recv_idx.back());
      MPI_CHECK_EQ_C_ARRAY(1, recv_buffer, recv_buffer_expected_p1, recv_idx.back());

    }
  }

  /*
   * One-way
   */
  std::vector<int> send_active_rank = {0, 1};
  std::vector<int> recv_active_rank = {0, 1};

  int n_send_active_rank = send_active_rank.size();
  int n_recv_active_rank = recv_active_rank.size();

  std::fill(begin(recv_buffer), end(recv_buffer), -1000);

  int req_recv = PDM_exchange_helper_exch_one_way_init(exch_helper,
                                                       PDM_EXCHANGE_DIRECTION_RECV,
                                                       sizeof(int),
                                                       1,
                                                       n_recv_active_rank,
                                                       recv_active_rank.data(),
                                                       recv_idx        .data(),
                                                       recv_n          .data(),
                                                       11,
                                                       recv_buffer.data());

  int req_send = PDM_exchange_helper_exch_one_way_init(exch_helper,
                                                       PDM_EXCHANGE_DIRECTION_SEND,
                                                       sizeof(int),
                                                       1,
                                                       n_send_active_rank,
                                                       send_active_rank   .data(),
                                                       send_idx   [i_rank].data(),
                                                       send_n             .data(),
                                                       11,
                                                       send_buffer[i_rank].data());

  PDM_exchange_helper_exch_start(exch_helper, req_send);
  PDM_exchange_helper_exch_start(exch_helper, req_recv);

  PDM_exchange_helper_exch_wait(exch_helper, req_send);
  PDM_exchange_helper_exch_wait(exch_helper, req_recv);

  PDM_exchange_helper_exch_free(exch_helper, req_send);
  PDM_exchange_helper_exch_free(exch_helper, req_recv);

  MPI_CHECK_EQ_C_ARRAY(0, recv_buffer, recv_buffer_expected_p0, recv_idx.back());
  MPI_CHECK_EQ_C_ARRAY(1, recv_buffer, recv_buffer_expected_p1, recv_idx.back());


  PDM_exchange_helper_free(exch_helper);
}
